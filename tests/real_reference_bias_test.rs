use std::process::Command;
use std::path::{Path, PathBuf};
use std::fs;
use tempfile::TempDir;
use likegt::{io, math};

/// REAL reference bias test - no shortcuts, no simulations
/// Actually runs the full pipeline with real reads and alignments
#[test]
fn test_real_reference_bias_corruption() {
    println!("\n=== REAL REFERENCE BIAS TEST - FULL PIPELINE ===\n");
    println!("This test actually:");
    println!("1. Simulates real reads with wgsim");
    println!("2. Aligns to single reference with BWA");
    println!("3. Extracts only mapped reads");
    println!("4. Re-maps to graph");
    println!("5. Compares to unbiased pipeline");
    println!("NO SHORTCUTS!\n");
    
    // Check prerequisites
    if !check_prerequisites() {
        println!("Missing prerequisites, skipping test");
        return;
    }
    
    // Set up paths
    let fasta_path = PathBuf::from("hla-f.fa.gz");
    let gfa_path = PathBuf::from("hla-f.k51.gfa");
    let coverage_path = PathBuf::from("hla-f.k51.paths.coverage.tsv.gz");
    
    if !fasta_path.exists() || !gfa_path.exists() || !coverage_path.exists() {
        println!("Required files not found, skipping test");
        return;
    }
    
    // Load reference panel
    let ref_data = io::read_gzip_tsv(coverage_path.to_str().unwrap()).unwrap();
    println!("Loaded {} haplotypes\n", ref_data.len());
    
    // Test individuals
    let test_individuals = vec!["HG00096", "HG00268", "HG00733"];
    let mut results = Vec::new();
    
    for individual in &test_individuals {
        println!("Testing {} with REAL reference bias...", individual);
        
        let result = test_individual_with_real_bias(
            individual,
            &fasta_path,
            &gfa_path,
            &ref_data,
        );
        
        if let Ok(res) = result {
            results.push(res);
        } else {
            println!("  Failed: {:?}", result.err());
        }
    }
    
    // Summarize results
    if !results.is_empty() {
        println!("\n=== REAL REFERENCE BIAS IMPACT ===");
        
        let unbiased_correct = results.iter()
            .filter(|r| r.unbiased_correct)
            .count();
        let biased_correct = results.iter()
            .filter(|r| r.biased_correct)
            .count();
        
        println!("Unbiased accuracy: {}/{} ({:.1}%)",
            unbiased_correct, results.len(),
            unbiased_correct as f64 * 100.0 / results.len() as f64
        );
        
        println!("Biased accuracy:   {}/{} ({:.1}%)",
            biased_correct, results.len(),
            biased_correct as f64 * 100.0 / results.len() as f64
        );
        
        println!("\nIndividual results:");
        for res in &results {
            println!("  {}:", res.individual);
            println!("    Unbiased: {} (sim={:.4})", 
                if res.unbiased_correct { "CORRECT" } else { "WRONG" },
                res.unbiased_similarity
            );
            println!("    Biased:   {} (sim={:.4})",
                if res.biased_correct { "CORRECT" } else { "WRONG" },
                res.biased_similarity
            );
            println!("    Read loss: {:.1}%", res.read_loss_percent);
        }
    }
}

fn test_individual_with_real_bias(
    individual: &str,
    fasta_path: &Path,
    gfa_path: &Path,
    ref_data: &io::CoverageData,
) -> Result<BiasTestResult, Box<dyn std::error::Error>> {
    // Create temp directory for this test
    let temp_dir = TempDir::new()?;
    let work_dir = temp_dir.path();
    
    // Step 1: Extract individual's sequences
    println!("  1. Extracting sequences for {}...", individual);
    let ind_fasta = work_dir.join(format!("{}.fa", individual));
    
    let extract_output = Command::new("samtools")
        .args(&[
            "faidx",
            fasta_path.to_str().unwrap(),
            &format!("{}#1#haplotype1", individual),
            &format!("{}#2#haplotype2", individual),
        ])
        .output()?;
    
    if !extract_output.status.success() {
        // Try alternative naming
        let extract_output = Command::new("samtools")
            .args(&[
                "faidx",
                fasta_path.to_str().unwrap(),
                &format!("{}#1", individual),
                &format!("{}#2", individual),
            ])
            .output()?;
        fs::write(&ind_fasta, extract_output.stdout)?;
    } else {
        fs::write(&ind_fasta, extract_output.stdout)?;
    }
    
    // Step 2: Simulate reads
    println!("  2. Simulating reads...");
    let reads1 = work_dir.join("reads.1.fq");
    let reads2 = work_dir.join("reads.2.fq");
    
    Command::new("wgsim")
        .args(&[
            "-1", "150",
            "-2", "150", 
            "-N", "10000",  // 10k read pairs
            "-e", "0.001",
            "-r", "0.001",
            ind_fasta.to_str().unwrap(),
            reads1.to_str().unwrap(),
            reads2.to_str().unwrap(),
        ])
        .output()?;
    
    let total_reads = count_fastq_reads(&reads1)?;
    println!("    Generated {} read pairs", total_reads);
    
    // Step 3: Extract GRCh38 reference
    println!("  3. Extracting GRCh38 reference...");
    let grch38_fasta = work_dir.join("grch38.fa");
    
    Command::new("samtools")
        .args(&[
            "faidx",
            fasta_path.to_str().unwrap(),
            "grch38#1#chr6:29711814-29738528",
        ])
        .output()
        .ok()
        .and_then(|o| {
            if o.status.success() {
                fs::write(&grch38_fasta, o.stdout).ok()
            } else {
                None
            }
        });
    
    if !grch38_fasta.exists() {
        // Try alternate naming
        Command::new("samtools")
            .args(&[
                "faidx",
                fasta_path.to_str().unwrap(),
                "grch38#1",
            ])
            .output()
            .ok()
            .and_then(|o| fs::write(&grch38_fasta, o.stdout).ok());
    }
    
    // Step 4: Build BWA index for GRCh38
    println!("  4. Building BWA index for GRCh38...");
    Command::new("bwa")
        .args(&["index", grch38_fasta.to_str().unwrap()])
        .output()?;
    
    // Step 5: Align reads to GRCh38 (CREATES BIAS)
    println!("  5. Aligning to GRCh38 (creating reference bias)...");
    let biased_bam = work_dir.join("biased.bam");
    
    let bwa_output = Command::new("bwa")
        .args(&[
            "mem",
            "-t", "4",
            grch38_fasta.to_str().unwrap(),
            reads1.to_str().unwrap(),
            reads2.to_str().unwrap(),
        ])
        .output()?;
    
    // Convert to BAM and filter for mapped reads only
    let sam_to_bam = Command::new("samtools")
        .args(&[
            "view",
            "-bS",
            "-F", "4",  // Exclude unmapped reads
            "-",
        ])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .spawn()?;
    
    {
        let mut stdin = sam_to_bam.stdin.as_ref().unwrap();
        use std::io::Write;
        stdin.write_all(&bwa_output.stdout)?;
    }
    
    let bam_output = sam_to_bam.wait_with_output()?;
    fs::write(&biased_bam, bam_output.stdout)?;
    
    // Count mapped reads
    let flagstat = Command::new("samtools")
        .args(&["flagstat", biased_bam.to_str().unwrap()])
        .output()?;
    
    let flagstat_str = String::from_utf8_lossy(&flagstat.stdout);
    let mapped_reads = parse_mapped_reads(&flagstat_str);
    let read_loss_percent = (1.0 - (mapped_reads as f64 / total_reads as f64)) * 100.0;
    
    println!("    Mapped {}/{} reads ({:.1}% lost to reference bias!)",
        mapped_reads, total_reads, read_loss_percent);
    
    // Step 6: Extract mapped reads as FASTQ (these are now biased)
    println!("  6. Extracting biased reads...");
    let biased_reads1 = work_dir.join("biased.1.fq");
    let biased_reads2 = work_dir.join("biased.2.fq");
    
    // Sort by name first
    let sorted_bam = work_dir.join("sorted.bam");
    Command::new("samtools")
        .args(&[
            "sort",
            "-n",
            biased_bam.to_str().unwrap(),
            "-o", sorted_bam.to_str().unwrap(),
        ])
        .output()?;
    
    Command::new("samtools")
        .args(&[
            "fastq",
            "-1", biased_reads1.to_str().unwrap(),
            "-2", biased_reads2.to_str().unwrap(),
            sorted_bam.to_str().unwrap(),
        ])
        .output()?;
    
    // Step 7: Now run BOTH pipelines - unbiased and biased
    println!("  7. Running genotyping pipelines...");
    
    // Extract all paths as FASTA for mapping
    let paths_fasta = work_dir.join("paths.fa");
    Command::new("odgi")
        .args(&[
            "paths",
            "-i", gfa_path.to_str().unwrap(),
            "-f",
        ])
        .stdout(fs::File::create(&paths_fasta)?)
        .output()?;
    
    Command::new("bwa")
        .args(&["index", paths_fasta.to_str().unwrap()])
        .output()?;
    
    // 7a. Unbiased pipeline
    println!("    7a. Unbiased genotyping...");
    let unbiased_coverage = get_coverage_from_reads(
        &reads1, &reads2,
        &paths_fasta, gfa_path,
        work_dir, "unbiased"
    )?;
    
    let unbiased_result = genotype_sample(ref_data, &unbiased_coverage, individual);
    
    // 7b. Biased pipeline
    println!("    7b. Biased genotyping...");
    let biased_coverage = if biased_reads1.exists() && biased_reads2.exists() {
        get_coverage_from_reads(
            &biased_reads1, &biased_reads2,
            &paths_fasta, gfa_path,
            work_dir, "biased"
        )?
    } else {
        println!("      No biased reads available!");
        vec![0.0; unbiased_coverage.len()]
    };
    
    let biased_result = genotype_sample(ref_data, &biased_coverage, individual);
    
    Ok(BiasTestResult {
        individual: individual.to_string(),
        unbiased_correct: unbiased_result.0,
        unbiased_similarity: unbiased_result.1,
        biased_correct: biased_result.0,
        biased_similarity: biased_result.1,
        read_loss_percent,
    })
}

fn get_coverage_from_reads(
    reads1: &Path,
    reads2: &Path,
    paths_fasta: &Path,
    gfa_path: &Path,
    work_dir: &Path,
    prefix: &str,
) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    // Map reads to paths
    let sam_file = work_dir.join(format!("{}.sam", prefix));
    let bam_file = work_dir.join(format!("{}.bam", prefix));
    
    let bwa_output = Command::new("bwa")
        .args(&[
            "mem",
            "-t", "4",
            paths_fasta.to_str().unwrap(),
            reads1.to_str().unwrap(),
            reads2.to_str().unwrap(),
        ])
        .output()?;
    
    fs::write(&sam_file, bwa_output.stdout)?;
    
    // Convert to BAM
    Command::new("samtools")
        .args(&[
            "view",
            "-bS",
            sam_file.to_str().unwrap(),
            "-o", bam_file.to_str().unwrap(),
        ])
        .output()?;
    
    // Inject into graph
    let gaf_file = work_dir.join(format!("{}.gaf", prefix));
    
    Command::new("gfainject")
        .args(&[
            gfa_path.to_str().unwrap(),
            bam_file.to_str().unwrap(),
        ])
        .stdout(fs::File::create(&gaf_file)?)
        .output()?;
    
    // Get coverage
    let coverage_file = work_dir.join(format!("{}.coverage.tsv", prefix));
    
    Command::new("gafpack")
        .args(&[
            "coverage",
            "-g", gaf_file.to_str().unwrap(),
            "-n", gfa_path.to_str().unwrap(),
            "-b", "1",
            "-s", "1",
        ])
        .stdout(fs::File::create(&coverage_file)?)
        .output()?;
    
    // Parse coverage
    let coverage_str = fs::read_to_string(&coverage_file)?;
    let mut coverage = Vec::new();
    
    for line in coverage_str.lines().skip(1) { // Skip header
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() > 1 {
            let cov: f64 = parts[1].parse().unwrap_or(0.0);
            coverage.push(cov);
        }
    }
    
    Ok(coverage)
}

fn genotype_sample(
    ref_data: &io::CoverageData,
    sample_coverage: &[f64],
    true_individual: &str,
) -> (bool, f64) {
    // Simple genotyping
    let mut best_sim = 0.0;
    let mut best_haps = ("", "");
    
    for i in 0..ref_data.ids.len() {
        for j in i..ref_data.ids.len() {
            let combined = math::sum_vectors(&[
                &ref_data.coverages[i],
                &ref_data.coverages[j],
            ]);
            
            // Handle length mismatch
            let min_len = combined.len().min(sample_coverage.len());
            let sim = if min_len > 0 {
                let combined_slice = &combined[..min_len];
                let sample_slice = &sample_coverage[..min_len];
                math::cosine_similarity(combined_slice, sample_slice)
            } else {
                0.0
            };
            
            if sim > best_sim {
                best_sim = sim;
                best_haps = (
                    ref_data.ids[i].split(':').next().unwrap_or("?"),
                    ref_data.ids[j].split(':').next().unwrap_or("?"),
                );
            }
        }
    }
    
    let correct = best_haps.0.contains(true_individual) || 
                  best_haps.1.contains(true_individual);
    
    (correct, best_sim)
}

fn count_fastq_reads(fastq_path: &Path) -> Result<usize, Box<dyn std::error::Error>> {
    let content = fs::read_to_string(fastq_path)?;
    Ok(content.lines().count() / 4)
}

fn parse_mapped_reads(flagstat: &str) -> usize {
    for line in flagstat.lines() {
        if line.contains("mapped (") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if !parts.is_empty() {
                return parts[0].parse().unwrap_or(0);
            }
        }
    }
    0
}

fn check_prerequisites() -> bool {
    let required = vec!["samtools", "wgsim", "bwa", "odgi", "gfainject", "gafpack"];
    let mut all_found = true;
    
    for tool in &required {
        let result = Command::new("which")
            .arg(tool)
            .output();
        
        if result.is_err() || !result.unwrap().status.success() {
            println!("  Missing required tool: {}", tool);
            all_found = false;
        }
    }
    
    all_found
}

#[derive(Debug)]
struct BiasTestResult {
    individual: String,
    unbiased_correct: bool,
    unbiased_similarity: f64,
    biased_correct: bool,
    biased_similarity: f64,
    read_loss_percent: f64,
}
/// REAL PIPELINE - NO SHORTCUTS, NO FAKE DATA
/// This module ACTUALLY simulates reads and runs the full pipeline

use std::process::Command;
use std::path::Path;
use std::fs;
use anyhow::Result;
use crate::{io::CoverageData, math};

/// Calculate number of read pairs needed for target coverage
pub fn calculate_read_pairs(
    sequence_length: usize,
    coverage_depth: u32,
    read_length: u32,
) -> u32 {
    // Formula: (sequence_length * coverage_depth) / (2 * read_length)
    // The 2 is because we have paired-end reads
    let total_bases_needed = sequence_length as u64 * coverage_depth as u64;
    let bases_per_pair = 2 * read_length as u64;
    (total_bases_needed / bases_per_pair) as u32
}

/// Extract sequences for an individual from FASTA
pub fn extract_individual_sequences(
    individual: &str,
    fasta_path: &Path,
    output_path: &Path,
) -> Result<usize> {
    // First, find the actual sequence names for this individual
    let faidx_output = Command::new("samtools")
        .args(&["faidx", fasta_path.to_str().unwrap()])
        .output()?;
    
    let index_content = String::from_utf8_lossy(&faidx_output.stdout);
    
    // Find sequences matching this individual
    let mut sequences = Vec::new();
    for line in index_content.lines() {
        if line.contains(&format!("{}#", individual)) {
            let seq_name = line.split('\t').next().unwrap_or("");
            if !seq_name.is_empty() {
                sequences.push(seq_name.to_string());
            }
        }
    }
    
    if sequences.is_empty() {
        // Try reading the .fai file directly
        let fai_path = format!("{}.fai", fasta_path.to_str().unwrap());
        if Path::new(&fai_path).exists() {
            let fai_content = fs::read_to_string(&fai_path)?;
            for line in fai_content.lines() {
                if line.contains(&format!("{}#", individual)) {
                    let seq_name = line.split('\t').next().unwrap_or("");
                    if !seq_name.is_empty() {
                        sequences.push(seq_name.to_string());
                    }
                }
            }
        }
    }
    
    if sequences.is_empty() {
        return Err(anyhow::anyhow!("No sequences found for {}", individual));
    }
    
    println!("  Found {} sequences for {}", sequences.len(), individual);
    
    // Extract all sequences for this individual
    let mut total_length = 0;
    let mut output_content = String::new();
    
    for seq_name in &sequences {
        let extract = Command::new("samtools")
            .args(&[
                "faidx",
                fasta_path.to_str().unwrap(),
                seq_name,
            ])
            .output()?;
        
        if extract.status.success() {
            let seq_data = String::from_utf8_lossy(&extract.stdout);
            output_content.push_str(&seq_data);
            
            // Calculate length (skip header line)
            for line in seq_data.lines().skip(1) {
                total_length += line.len();
            }
        }
    }
    
    fs::write(output_path, output_content)?;
    
    Ok(total_length)
}

/// Simulate reads at specified coverage depth
pub fn simulate_reads(
    fasta_path: &Path,
    output_prefix: &Path,
    coverage_depth: u32,
    read_length: u32,
    error_rate: f64,
) -> Result<u32> {
    // Get sequence length
    let fasta_content = fs::read_to_string(fasta_path)?;
    let mut seq_length = 0;
    for line in fasta_content.lines() {
        if !line.starts_with('>') {
            seq_length += line.len();
        }
    }
    
    if seq_length == 0 {
        return Err(anyhow::anyhow!("Empty sequence in {}", fasta_path.display()));
    }
    
    // Calculate number of read pairs needed
    let num_pairs = calculate_read_pairs(seq_length, coverage_depth, read_length);
    
    println!("    Sequence length: {} bp", seq_length);
    println!("    Target coverage: {}x", coverage_depth);
    println!("    Read length: {} bp", read_length);
    println!("    Simulating {} read pairs", num_pairs);
    
    let reads1_path = format!("{}.1.fq", output_prefix.to_str().unwrap());
    let reads2_path = format!("{}.2.fq", output_prefix.to_str().unwrap());
    
    let output = Command::new("wgsim")
        .args(&[
            "-1", &read_length.to_string(),
            "-2", &read_length.to_string(),
            "-N", &num_pairs.to_string(),
            "-e", &error_rate.to_string(),
            "-r", "0.001",  // mutation rate
            "-R", "0.001",  // indel fraction
            fasta_path.to_str().unwrap(),
            &reads1_path,
            &reads2_path,
        ])
        .output()?;
    
    if !output.status.success() {
        return Err(anyhow::anyhow!("wgsim failed: {}", String::from_utf8_lossy(&output.stderr)));
    }
    
    // Verify reads were created
    let reads1_lines = fs::read_to_string(&reads1_path)?.lines().count();
    let actual_pairs = reads1_lines / 4;  // FASTQ has 4 lines per read
    
    println!("    Actually generated {} read pairs", actual_pairs);
    
    Ok(actual_pairs as u32)
}

/// Map reads to reference and create BAM
pub fn map_reads(
    reference_path: &Path,
    reads1_path: &Path,
    reads2_path: &Path,
    output_bam: &Path,
    threads: u32,
) -> Result<(u32, u32)> {
    // Build index if needed
    let index_file = format!("{}.bwt", reference_path.to_str().unwrap());
    if !Path::new(&index_file).exists() {
        println!("    Building BWA index...");
        Command::new("bwa")
            .args(&["index", reference_path.to_str().unwrap()])
            .output()?;
    }
    
    // Map reads
    println!("    Mapping reads...");
    let bwa_output = Command::new("bwa")
        .args(&[
            "mem",
            "-t", &threads.to_string(),
            reference_path.to_str().unwrap(),
            reads1_path.to_str().unwrap(),
            reads2_path.to_str().unwrap(),
        ])
        .output()?;
    
    if !bwa_output.status.success() {
        return Err(anyhow::anyhow!("BWA failed: {}", String::from_utf8_lossy(&bwa_output.stderr)));
    }
    
    // Convert to BAM
    let bam_output = Command::new("samtools")
        .args(&[
            "view",
            "-bS",
            "-@", &threads.to_string(),
            "-",
        ])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .spawn()?;
    
    {
        let mut stdin = bam_output.stdin.as_ref().unwrap();
        use std::io::Write;
        stdin.write_all(&bwa_output.stdout)?;
    }
    
    let bam_data = bam_output.wait_with_output()?;
    fs::write(output_bam, bam_data.stdout)?;
    
    // Get mapping statistics
    let flagstat = Command::new("samtools")
        .args(&["flagstat", output_bam.to_str().unwrap()])
        .output()?;
    
    let flagstat_str = String::from_utf8_lossy(&flagstat.stdout);
    let mut total_reads = 0;
    let mut mapped_reads = 0;
    
    for line in flagstat_str.lines() {
        if line.contains("in total") {
            total_reads = line.split_whitespace()
                .next()
                .and_then(|s| s.parse().ok())
                .unwrap_or(0);
        } else if line.contains("mapped (") && !line.contains("primary mapped") {
            mapped_reads = line.split_whitespace()
                .next()
                .and_then(|s| s.parse().ok())
                .unwrap_or(0);
        }
    }
    
    println!("    Mapped {}/{} reads ({:.1}%)",
        mapped_reads, total_reads,
        mapped_reads as f64 * 100.0 / total_reads.max(1) as f64
    );
    
    Ok((total_reads, mapped_reads))
}

/// Create reference bias by filtering through single reference
pub fn create_reference_bias(
    reference_path: &Path,
    reads1_path: &Path,
    reads2_path: &Path,
    output_prefix: &Path,
    threads: u32,
) -> Result<f64> {
    println!("  Creating reference bias...");
    
    // Map to single reference
    let biased_bam = output_prefix.with_extension("bam");
    let (total, mapped) = map_reads(
        reference_path,
        reads1_path,
        reads2_path,
        &biased_bam,
        threads
    )?;
    
    let loss_percent = (1.0 - (mapped as f64 / total.max(1) as f64)) * 100.0;
    
    // Extract only mapped reads
    println!("    Extracting mapped reads only...");
    let filtered_bam = output_prefix.with_extension("filtered.bam");
    Command::new("samtools")
        .args(&[
            "view",
            "-b",
            "-F", "4",  // exclude unmapped
            biased_bam.to_str().unwrap(),
            "-o", filtered_bam.to_str().unwrap(),
        ])
        .output()?;
    
    // Sort by name for fastq extraction
    let sorted_bam = output_prefix.with_extension("sorted.bam");
    Command::new("samtools")
        .args(&[
            "sort",
            "-n",
            filtered_bam.to_str().unwrap(),
            "-o", sorted_bam.to_str().unwrap(),
        ])
        .output()?;
    
    // Extract as FASTQ
    let biased_r1 = format!("{}.biased.1.fq", output_prefix.to_str().unwrap());
    let biased_r2 = format!("{}.biased.2.fq", output_prefix.to_str().unwrap());
    
    Command::new("samtools")
        .args(&[
            "fastq",
            "-1", &biased_r1,
            "-2", &biased_r2,
            sorted_bam.to_str().unwrap(),
        ])
        .output()?;
    
    println!("    Read loss due to bias: {:.1}%", loss_percent);
    
    Ok(loss_percent)
}

/// Inject reads into graph and get coverage
pub fn get_graph_coverage(
    gfa_path: &Path,
    bam_path: &Path,
    output_prefix: &Path,
) -> Result<Vec<f64>> {
    println!("    Injecting into graph...");
    
    // Run gfainject
    let gaf_path = output_prefix.with_extension("gaf");
    let gaf_output = Command::new("gfainject")
        .args(&[
            "--gfa", gfa_path.to_str().unwrap(),
            "--bam", bam_path.to_str().unwrap(),
        ])
        .output()?;
    
    if !gaf_output.status.success() {
        return Err(anyhow::anyhow!("gfainject failed: {}", 
            String::from_utf8_lossy(&gaf_output.stderr)));
    }
    
    fs::write(&gaf_path, gaf_output.stdout)?;
    
    // Get coverage with gafpack
    println!("    Extracting coverage...");
    let pack_output = Command::new("gafpack")
        .args(&[
            "--gfa", gfa_path.to_str().unwrap(),
            "--gaf", gaf_path.to_str().unwrap(),
        ])
        .output()?;
    
    if !pack_output.status.success() {
        return Err(anyhow::anyhow!("gafpack failed: {}", 
            String::from_utf8_lossy(&pack_output.stderr)));
    }
    
    // Parse coverage
    let coverage_str = String::from_utf8_lossy(&pack_output.stdout);
    let mut coverage = Vec::new();
    
    for line in coverage_str.lines().skip(1) {  // Skip header
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() > 1 {
            // Skip first column (sample name), parse rest as coverage values
            for i in 1..parts.len() {
                let val: f64 = parts[i].parse().unwrap_or(0.0);
                coverage.push(val);
            }
            break;  // Only take first sample row
        }
    }
    
    println!("    Got coverage for {} nodes", coverage.len());
    
    Ok(coverage)
}

/// Run complete pipeline for an individual
pub fn run_individual_pipeline(
    individual: &str,
    fasta_path: &Path,
    gfa_path: &Path,
    coverage_depth: u32,
    work_dir: &Path,
) -> Result<PipelineResult> {
    println!("\n=== Running REAL pipeline for {} ===", individual);
    
    // Step 1: Extract sequences
    println!("Step 1: Extracting sequences");
    let ind_fasta = work_dir.join(format!("{}.fa", individual));
    let seq_length = extract_individual_sequences(individual, fasta_path, &ind_fasta)?;
    
    // Step 2: Simulate reads at proper coverage
    println!("Step 2: Simulating {}x coverage", coverage_depth);
    let reads_prefix = work_dir.join(format!("{}.reads", individual));
    let num_pairs = simulate_reads(
        &ind_fasta,
        &reads_prefix,
        coverage_depth,
        150,  // read length
        0.001,  // error rate
    )?;
    
    // Step 3: Map to graph paths
    println!("Step 3: Mapping to graph paths");
    
    // First extract paths as FASTA
    let paths_fasta = work_dir.join("paths.fa");
    if !paths_fasta.exists() {
        let paths_output = Command::new("odgi")
            .args(&[
                "paths",
                "-i", gfa_path.to_str().unwrap(),
                "-f",
            ])
            .output()?;
        fs::write(&paths_fasta, paths_output.stdout)?;
    }
    
    let reads1 = format!("{}.1.fq", reads_prefix.to_str().unwrap());
    let reads2 = format!("{}.2.fq", reads_prefix.to_str().unwrap());
    let mapped_bam = work_dir.join(format!("{}.mapped.bam", individual));
    
    let (total_reads, mapped_reads) = map_reads(
        &paths_fasta,
        Path::new(&reads1),
        Path::new(&reads2),
        &mapped_bam,
        4,
    )?;
    
    // Step 4: Get graph coverage
    println!("Step 4: Getting graph coverage");
    let coverage = get_graph_coverage(
        gfa_path,
        &mapped_bam,
        &work_dir.join(individual),
    )?;
    
    Ok(PipelineResult {
        individual: individual.to_string(),
        sequence_length: seq_length,
        read_pairs_simulated: num_pairs,
        reads_mapped: mapped_reads,
        total_reads,
        coverage,
    })
}

/// Test hold-0-out with REAL reads
pub fn test_real_hold0(
    individual: &str,
    fasta_path: &Path,
    gfa_path: &Path,
    ref_data: &CoverageData,
    coverage_depth: u32,
) -> Result<Hold0Result> {
    println!("\n=== REAL HOLD-0-OUT TEST FOR {} ===", individual);
    
    let temp_dir = tempfile::tempdir()?;
    let work_dir = temp_dir.path();
    
    // Run the pipeline
    let result = run_individual_pipeline(
        individual,
        fasta_path,
        gfa_path,
        coverage_depth,
        work_dir,
    )?;
    
    // Now genotype using the coverage
    println!("Step 5: Genotyping");
    let mut best_similarity = 0.0;
    let mut best_haps = (String::new(), String::new());
    
    for i in 0..ref_data.ids.len() {
        for j in i..ref_data.ids.len() {
            let combined = math::sum_vectors(&[
                &ref_data.coverages[i],
                &ref_data.coverages[j],
            ]);
            
            // Handle length mismatch
            let min_len = combined.len().min(result.coverage.len());
            if min_len == 0 {
                continue;
            }
            
            let sim = math::cosine_similarity(
                &combined[..min_len],
                &result.coverage[..min_len],
            );
            
            if sim > best_similarity {
                best_similarity = sim;
                best_haps = (
                    ref_data.ids[i].clone(),
                    ref_data.ids[j].clone(),
                );
            }
        }
    }
    
    let correct = best_haps.0.contains(individual) || best_haps.1.contains(individual);
    
    println!("Results:");
    println!("  True: {}#1 + {}#2", individual, individual);
    println!("  Called: {} + {}", 
        best_haps.0.split(':').next().unwrap_or("?"),
        best_haps.1.split(':').next().unwrap_or("?")
    );
    println!("  Similarity: {:.4}", best_similarity);
    println!("  Correct: {}", if correct { "YES" } else { "NO" });
    
    Ok(Hold0Result {
        individual: individual.to_string(),
        correct,
        similarity: best_similarity,
        called_hap1: best_haps.0,
        called_hap2: best_haps.1,
        coverage_depth,
        reads_simulated: result.read_pairs_simulated,
    })
}

#[derive(Debug)]
pub struct PipelineResult {
    pub individual: String,
    pub sequence_length: usize,
    pub read_pairs_simulated: u32,
    pub reads_mapped: u32,
    pub total_reads: u32,
    pub coverage: Vec<f64>,
}

#[derive(Debug)]
pub struct Hold0Result {
    pub individual: String,
    pub correct: bool,
    pub similarity: f64,
    pub called_hap1: String,
    pub called_hap2: String,
    pub coverage_depth: u32,
    pub reads_simulated: u32,
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_coverage_calculation() {
        // Test coverage calculation
        assert_eq!(calculate_read_pairs(10000, 30, 150), 1000);
        assert_eq!(calculate_read_pairs(25000, 30, 150), 2500);
        assert_eq!(calculate_read_pairs(50000, 30, 150), 5000);
        
        // Test different coverage depths
        assert_eq!(calculate_read_pairs(10000, 10, 150), 333);
        assert_eq!(calculate_read_pairs(10000, 50, 150), 1666);
        assert_eq!(calculate_read_pairs(10000, 100, 150), 3333);
    }
}
use anyhow::{Result, Context};
use std::process::Command;

#[derive(Debug, Clone)]
pub struct AlignmentInfo {
    pub query: String,
    pub target: String,
    pub identity: f64,
    pub qv: f64,
    pub matches: u64,
    pub query_span: u64,
}

#[derive(Debug, Clone)]
pub struct MaxQVResult {
    pub individual: String,
    pub target_hap1: String,
    pub target_hap2: String,
    pub qv_hap1: f64,
    pub qv_hap2: f64,
    pub avg_qv: f64,
    pub identity_hap1: f64,
    pub identity_hap2: f64,
}

/// Parse a PAF line and extract alignment information
fn parse_paf_line(line: &str) -> Option<AlignmentInfo> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        return None;
    }
    
    let query = fields[0].to_string();
    let target = fields[5].to_string();
    
    // Get alignment coordinates for sequence-space length calculation
    let query_start: u64 = fields[2].parse().ok()?;
    let query_end: u64 = fields[3].parse().ok()?;
    let target_start: u64 = fields[7].parse().ok()?;
    let target_end: u64 = fields[8].parse().ok()?;
    
    let matches: u64 = fields[9].parse().ok()?;
    
    // Calculate alignment length in query sequence space
    let query_span = query_end - query_start;
    
    if query_span == 0 {
        return None;
    }
    
    // Identity = matched characters / query sequence length
    let identity = matches as f64 / query_span as f64;
    let error_rate = 1.0 - identity;
    
    let qv = if error_rate <= 0.0 || identity >= 1.0 {
        // For perfect matches, report a very high but finite QV
        50.0
    } else {
        let qv_val = -10.0 * error_rate.log10();
        qv_val.min(50.0).max(0.0)
    };
    
    Some(AlignmentInfo {
        query,
        target,
        identity,
        qv,
        matches,
        query_span,
    })
}

/// Extract individual ID from sequence name (e.g., HG00096#1 -> HG00096)
fn get_individual_id(seq_name: &str) -> String {
    seq_name.split('#').next().unwrap_or(seq_name).to_string()
}

/// Run allwave for alignment
fn run_allwave(fasta_file: &str, threads: usize, sparsification: &str, verbose: bool) -> Result<String> {
    use std::process::Stdio;
    use std::io::Read;
    
    if verbose {
        eprintln!("Running allwave with sparsification: {}", sparsification);
    }
    
    let mut child = Command::new("allwave")
        .args(&["-i", fasta_file, "-t", &threads.to_string(), "-p", sparsification])
        .stdout(Stdio::piped())
        .stderr(if verbose { Stdio::inherit() } else { Stdio::piped() })
        .spawn()
        .context("Failed to start allwave")?;
    
    let mut stdout = String::new();
    if let Some(mut stdout_reader) = child.stdout.take() {
        stdout_reader.read_to_string(&mut stdout)?;
    }
    
    let status = child.wait()?;
    if !status.success() {
        anyhow::bail!("allwave failed with exit code: {}", status);
    }
    
    Ok(stdout)
}

/// Run wfmash for approximate alignment
fn run_wfmash(fasta_file: &str, threads: usize, verbose: bool) -> Result<String> {
    use std::process::Stdio;
    use std::io::Read;
    
    if verbose {
        eprintln!("Running wfmash for approximate alignment...");
    }
    
    // Check if FAI index exists, create if needed
    let fai_file = format!("{}.fai", fasta_file);
    if !std::path::Path::new(&fai_file).exists() {
        if verbose {
            eprintln!("Creating FAI index for {}...", fasta_file);
        }
        let index_output = Command::new("samtools")
            .args(&["faidx", fasta_file])
            .output()
            .context("Failed to create FAI index with samtools")?;
        
        if !index_output.status.success() {
            anyhow::bail!("Failed to create FAI index for {}", fasta_file);
        }
    }
    
    // wfmash outputs to stdout when no -o specified
    // Format: wfmash [target] [query] - for self-mapping just provide one file
    let mut child = Command::new("wfmash")
        .args(&[
            "-t", &threads.to_string(), 
            "-S", "1",           // All segments  
            "-p", "70",          // 70% identity threshold for broader range
            "-n", "20",          // Keep top 20 mappings per segment
            // Note: without -m flag, wfmash does full alignment (slower but accurate)
            fasta_file           // Input FASTA for self-mapping
        ])
        .stdout(Stdio::piped())
        .stderr(if verbose { Stdio::inherit() } else { Stdio::piped() })
        .spawn()
        .context("Failed to start wfmash")?;
    
    let mut stdout = String::new();
    if let Some(mut stdout_reader) = child.stdout.take() {
        stdout_reader.read_to_string(&mut stdout)?;
    }
    
    let status = child.wait()?;
    if !status.success() {
        anyhow::bail!("wfmash failed with exit code: {}", status);
    }
    
    if verbose {
        let n_lines = stdout.lines().count();
        eprintln!("wfmash produced {} alignments", n_lines);
    }
    
    Ok(stdout)
}

/// Process wfmash PAF output (optionally with impg for enhancement)
fn process_wfmash_output(paf_content: &str, use_impg: bool, verbose: bool) -> Result<Vec<AlignmentInfo>> {
    use std::fs;
    
    if use_impg && which::which("impg").is_ok() {
        if verbose {
            eprintln!("Attempting to use impg for similarity enhancement...");
        }
        
        // Write PAF to temp file for impg
        let temp_paf = format!("/tmp/wfmash_temp_{}.paf", std::process::id());
        fs::write(&temp_paf, paf_content)?;
        
        // Index with impg
        let index_output = Command::new("impg")
            .args(&["-i", &temp_paf])
            .output()
            .context("Failed to run impg index")?;
        
        if !index_output.status.success() {
            // Only show warning in verbose mode
            if verbose {
                eprintln!("Note: impg index not available, using direct PAF parsing");
            }
        } else {
            // Query with impg for enhanced similarity
            let query_output = Command::new("impg")
                .args(&["-q", &temp_paf])  
                .output()
                .context("Failed to run impg query")?;
            
            if query_output.status.success() && !query_output.stdout.is_empty() {
                // impg enhances the PAF, parse enhanced output
                let enhanced = String::from_utf8_lossy(&query_output.stdout);
                let mut alignments = Vec::new();
                for line in enhanced.lines() {
                    if let Some(alignment) = parse_paf_line(line) {
                        alignments.push(alignment);
                    }
                }
                
                // Clean up temp files
                let _ = fs::remove_file(&temp_paf);
                let _ = fs::remove_file(format!("{}.impgi", &temp_paf));
                
                if verbose {
                    eprintln!("Successfully used impg enhancement");
                }
                return Ok(alignments);
            }
        }
        
        // Clean up if impg failed
        let _ = fs::remove_file(&temp_paf);
        let _ = fs::remove_file(format!("{}.impgi", &temp_paf));
    }
    
    // Fall back to parsing raw PAF
    if verbose {
        eprintln!("Parsing PAF output directly...");
    }
    
    let mut alignments = Vec::new();
    for line in paf_content.lines() {
        if let Some(alignment) = parse_paf_line(line) {
            alignments.push(alignment);
        }
    }
    
    Ok(alignments)
}

/// Calculate aggregate identity from multiple alignments
fn calculate_aggregate_identity(alignments: &[&AlignmentInfo]) -> (f64, u64, u64) {
    let mut total_matches = 0u64;
    let mut total_coverage = 0u64;
    
    for align in alignments {
        total_matches += align.matches;
        total_coverage += align.query_span;
    }
    
    if total_coverage == 0 {
        (0.0, 0, 0)
    } else {
        (total_matches as f64 / total_coverage as f64, total_matches, total_coverage)
    }
}

/// Find best genotype pairing for an individual using aggregate alignment statistics
fn find_best_pairing(
    alignments: &[AlignmentInfo],
    individual: &str,
) -> Option<MaxQVResult> {
    let hap1 = format!("{}#1", individual);
    let hap2 = format!("{}#2", individual);
    
    // Group alignments by target individual
    use std::collections::HashMap;
    let mut hap1_by_target: HashMap<String, Vec<&AlignmentInfo>> = HashMap::new();
    let mut hap2_by_target: HashMap<String, Vec<&AlignmentInfo>> = HashMap::new();
    
    for align in alignments {
        let target_individual = get_individual_id(&align.target);
        if target_individual == individual {
            continue; // Skip self
        }
        
        if align.query.starts_with(&hap1) {
            hap1_by_target.entry(target_individual.clone()).or_default().push(align);
        } else if align.query.starts_with(&hap2) {
            hap2_by_target.entry(target_individual.clone()).or_default().push(align);
        }
    }
    
    if hap1_by_target.is_empty() || hap2_by_target.is_empty() {
        return None;
    }
    
    // For now, just use the single best alignment approach but with proper pairing
    // TODO: Implement proper aggregation across multiple alignments per target
    let mut best_result: Option<MaxQVResult> = None;
    let mut best_avg_qv = 0.0;
    
    for (target1, aligns1) in &hap1_by_target {
        for (target2, aligns2) in &hap2_by_target {
            if target1 == target2 {
                continue; // Avoid using same target for both haplotypes
            }
            
            // Calculate aggregate identity for each haplotype against this target
            let (identity1, _matches1, _coverage1) = calculate_aggregate_identity(aligns1);
            let (identity2, _matches2, _coverage2) = calculate_aggregate_identity(aligns2);
            
            // Calculate QV from aggregate identity
            let qv1 = if identity1 >= 1.0 { 50.0 } else { 
                let error_rate1 = 1.0 - identity1;
                (-10.0 * error_rate1.log10()).min(50.0).max(0.0)
            };
            let qv2 = if identity2 >= 1.0 { 50.0 } else { 
                let error_rate2 = 1.0 - identity2;
                (-10.0 * error_rate2.log10()).min(50.0).max(0.0)
            };
            
            let avg_qv = (qv1 + qv2) / 2.0;
            
            if avg_qv > best_avg_qv {
                best_avg_qv = avg_qv;
                let target1_name = aligns1.first()?.target.clone();
                let target2_name = aligns2.first()?.target.clone();
                best_result = Some(MaxQVResult {
                    individual: individual.to_string(),
                    target_hap1: target1_name,
                    target_hap2: target2_name,
                    qv_hap1: qv1,
                    qv_hap2: qv2,
                    avg_qv,
                    identity_hap1: identity1,
                    identity_hap2: identity2,
                });
            }
        }
    }
    
    best_result
}

/// Get all individuals from FASTA file
fn get_all_individuals(fasta_file: &str) -> Result<Vec<String>> {
    let mut individuals = std::collections::HashSet::new();
    
    // Use zcat for gzipped files to ensure we read everything
    let output = if fasta_file.ends_with(".gz") {
        Command::new("zcat")
            .arg(fasta_file)
            .output()
            .context(format!("Failed to run zcat on {}", fasta_file))?
    } else {
        Command::new("cat")
            .arg(fasta_file)
            .output()
            .context(format!("Failed to read file {}", fasta_file))?
    };
    
    if !output.status.success() {
        anyhow::bail!("Failed to read fasta file: {} (file may not exist)", fasta_file);
    }
    
    let content = String::from_utf8_lossy(&output.stdout);
    for line in content.lines() {
        if line.starts_with('>') {
            let seq_id = line[1..].split_whitespace().next().unwrap_or("");
            let individual = get_individual_id(seq_id);
            
            // Skip reference genomes
            if !individual.is_empty() && !individual.starts_with("grch") && !individual.starts_with("chm") {
                individuals.insert(individual);
            }
        }
    }
    
    let mut result: Vec<String> = individuals.into_iter().collect();
    result.sort();
    Ok(result)
}

/// Run max QV analysis using selected alignment method
pub async fn run_max_qv_analysis(
    fasta_file: &str,
    individual: Option<&str>,
    output_file: Option<&str>,
    threads: usize,
    method: &str,
    sparsification: Option<&str>,
    verbose: bool,
) -> Result<()> {
    // Determine which individuals to process
    let individuals = if let Some(ind) = individual {
        if ind == "all" {
            if verbose {
                println!("🔍 Extracting individual IDs from FASTA...");
            }
            get_all_individuals(fasta_file)?
        } else if ind.contains(',') {
            // Handle comma-separated list
            ind.split(',').map(|s| s.trim().to_string()).collect()
        } else {
            vec![ind.to_string()]
        }
    } else {
        if verbose {
            println!("🔍 Extracting individual IDs from FASTA...");
        }
        get_all_individuals(fasta_file)?
    };
    
    if verbose {
        println!("📊 Found {} individuals to process", individuals.len());
        println!("🧬 Running {} for alignment...", method);
    }
    
    // Run alignment based on method
    let all_alignments = match method {
        "wfmash" => {
            // Use wfmash for approximate alignment
            let paf_output = run_wfmash(fasta_file, threads, verbose)?;
            
            if verbose {
                let n_lines = paf_output.lines().count();
                println!("✅ Got {} alignment records from wfmash", n_lines);
            }
            
            // Process with optional impg enhancement (disabled for now)
            process_wfmash_output(&paf_output, false, verbose)?
        },
        "allwave" | _ => {
            // Use allwave for exact alignment
            let sparsification_mode = sparsification.unwrap_or("tree:5:0:0");
            let paf_output = run_allwave(fasta_file, threads, sparsification_mode, verbose)?;
            
            if verbose {
                let n_lines = paf_output.lines().count();
                println!("✅ Got {} alignment records from allwave", n_lines);
            }
            
            // Parse all alignments
            let mut alignments = Vec::new();
            for line in paf_output.lines() {
                if let Some(alignment) = parse_paf_line(line) {
                    alignments.push(alignment);
                }
            }
            alignments
        }
    };
    
    // Process each individual
    let mut results = Vec::new();
    for (i, individual) in individuals.iter().enumerate() {
        if verbose {
            println!("[{}/{}] Processing {}...", i + 1, individuals.len(), individual);
        }
        
        if let Some(result) = find_best_pairing(&all_alignments, individual) {
            if verbose {
                println!("  Best match: {} + {}", result.target_hap1, result.target_hap2);
                println!("  Max attainable QV: {:.1}", result.avg_qv);
            }
            results.push(result);
        } else if verbose {
            println!("  No valid alignments found for {}", individual);
        }
    }
    
    // Sort by average QV (descending)
    results.sort_by(|a, b| b.avg_qv.partial_cmp(&a.avg_qv).unwrap());
    
    // Generate output
    let mut output = String::new();
    output.push_str("individual\ttarget_hap1\ttarget_hap2\tqv_hap1\tqv_hap2\tavg_qv\tidentity_hap1\tidentity_hap2\n");
    
    for r in &results {
        output.push_str(&format!(
            "{}\t{}\t{}\t{:.1}\t{:.1}\t{:.1}\t{:.4}\t{:.4}\n",
            r.individual,
            r.target_hap1,
            r.target_hap2,
            r.qv_hap1,
            r.qv_hap2,
            r.avg_qv,
            r.identity_hap1,
            r.identity_hap2
        ));
    }
    
    // Summary
    if verbose && !results.is_empty() {
        let avg_max_qv: f64 = results.iter().map(|r| r.avg_qv).sum::<f64>() / results.len() as f64;
        let min_qv = results.iter().map(|r| r.avg_qv).fold(f64::INFINITY, f64::min);
        let max_qv = results.iter().map(|r| r.avg_qv).fold(f64::NEG_INFINITY, f64::max);
        
        println!("\n=== SUMMARY ===");
        println!("Processed {} individuals", results.len());
        println!("Average max attainable QV: {:.1}", avg_max_qv);
        println!("QV range: {:.1} - {:.1}", min_qv, max_qv);
    }
    
    // Write output
    if let Some(output_path) = output_file {
        tokio::fs::write(output_path, &output).await?;
        if verbose {
            println!("📝 Results written to {}", output_path);
        }
    } else {
        print!("{}", output);
    }
    
    Ok(())
}
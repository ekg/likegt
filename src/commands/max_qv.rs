use anyhow::{Result, Context};
use std::process::Command;

#[derive(Debug, Clone)]
pub struct AlignmentInfo {
    pub query: String,
    pub target: String,
    pub identity: f64,
    pub qv: f64,
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
    let matches: u64 = fields[9].parse().ok()?;
    let align_len: u64 = fields[10].parse().ok()?;
    
    if align_len == 0 {
        return None;
    }
    
    let identity = matches as f64 / align_len as f64;
    let error_rate = 1.0 - identity;
    
    let qv = if error_rate <= 0.0 {
        60.0
    } else {
        let qv_val = -10.0 * error_rate.log10();
        qv_val.min(60.0).max(0.0)
    };
    
    Some(AlignmentInfo {
        query,
        target,
        identity,
        qv,
    })
}

/// Extract individual ID from sequence name (e.g., HG00096#1 -> HG00096)
fn get_individual_id(seq_name: &str) -> String {
    seq_name.split('#').next().unwrap_or(seq_name).to_string()
}

/// Run allwave for all-vs-all alignment
fn run_allwave(fasta_file: &str, threads: usize, verbose: bool) -> Result<String> {
    use std::process::Stdio;
    use std::io::Read;
    
    let mut child = Command::new("allwave")
        .args(&["-i", fasta_file, "-t", &threads.to_string()])
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

/// Find best genotype pairing for an individual
fn find_best_pairing(
    alignments: &[AlignmentInfo],
    individual: &str,
) -> Option<MaxQVResult> {
    let hap1 = format!("{}#1", individual);
    let hap2 = format!("{}#2", individual);
    
    // Get alignments for each haplotype (excluding self)
    let mut hap1_aligns: Vec<&AlignmentInfo> = alignments
        .iter()
        .filter(|a| {
            a.query.starts_with(&hap1) && 
            get_individual_id(&a.target) != individual
        })
        .collect();
    
    let mut hap2_aligns: Vec<&AlignmentInfo> = alignments
        .iter()
        .filter(|a| {
            a.query.starts_with(&hap2) && 
            get_individual_id(&a.target) != individual
        })
        .collect();
    
    if hap1_aligns.is_empty() || hap2_aligns.is_empty() {
        return None;
    }
    
    // Sort by QV (best first)
    hap1_aligns.sort_by(|a, b| b.qv.partial_cmp(&a.qv).unwrap());
    hap2_aligns.sort_by(|a, b| b.qv.partial_cmp(&a.qv).unwrap());
    
    // Find best pairing (avoiding same target twice)
    let mut best_result: Option<MaxQVResult> = None;
    let mut best_avg_qv = 0.0;
    
    let n_candidates = 10.min(hap1_aligns.len()).min(hap2_aligns.len());
    
    for i in 0..n_candidates {
        for j in 0..n_candidates {
            let align1 = hap1_aligns[i];
            let align2 = hap2_aligns[j];
            
            // Skip if same target used twice
            if align1.target == align2.target {
                continue;
            }
            
            let avg_qv = (align1.qv + align2.qv) / 2.0;
            
            if avg_qv > best_avg_qv {
                best_avg_qv = avg_qv;
                best_result = Some(MaxQVResult {
                    individual: individual.to_string(),
                    target_hap1: align1.target.clone(),
                    target_hap2: align2.target.clone(),
                    qv_hap1: align1.qv,
                    qv_hap2: align2.qv,
                    avg_qv,
                    identity_hap1: align1.identity,
                    identity_hap2: align2.identity,
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
            .context("Failed to run zcat")?
    } else {
        Command::new("cat")
            .arg(fasta_file)
            .output()
            .context("Failed to read file")?
    };
    
    if !output.status.success() {
        anyhow::bail!("Failed to read fasta file");
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

/// Run max QV analysis using allwave
pub async fn run_max_qv_analysis(
    fasta_file: &str,
    individual: Option<&str>,
    output_file: Option<&str>,
    threads: usize,
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
        println!("🧬 Running allwave for all-vs-all alignment...");
    }
    
    // Run allwave once for all alignments
    let paf_output = run_allwave(fasta_file, threads, verbose)?;
    
    if verbose {
        let n_lines = paf_output.lines().count();
        println!("✅ Got {} alignment records", n_lines);
    }
    
    // Parse all alignments
    let mut all_alignments = Vec::new();
    for line in paf_output.lines() {
        if let Some(alignment) = parse_paf_line(line) {
            all_alignments.push(alignment);
        }
    }
    
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
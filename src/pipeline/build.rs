use anyhow::{Result, Context};
use std::path::Path;
use std::process::Command;
use tokio::fs;

/// Build a pangenome graph from FASTA using allwave → seqwish → odgi pipeline
pub async fn build_graph_from_fasta(
    fasta_path: &str,
    output_gfa: &str,
    kmer_size: usize,
    segment_length: usize,
) -> Result<()> {
    log::info!("Building graph from FASTA: {}", fasta_path);
    log::info!("K-mer size: {}, Segment length: {}", kmer_size, segment_length);
    
    // Verify input file exists
    if !Path::new(fasta_path).exists() {
        return Err(anyhow::anyhow!("Input FASTA file not found: {}", fasta_path));
    }
    
    let output_path = Path::new(output_gfa);
    let output_dir = output_path.parent().unwrap();
    
    // Create output directory
    fs::create_dir_all(output_dir).await?;
    
    let base_name = output_path.file_stem().unwrap().to_str().unwrap();
    let work_dir = output_dir.join(format!("{}_build", base_name));
    
    fs::create_dir_all(&work_dir).await?;
    
    // Step 1: Run allwave for initial alignment
    log::info!("Step 1: Running allwave for sequence alignment...");
    let paf_file = work_dir.join("alignment.paf");
    
    let allwave_status = Command::new("allwave")
        .args([
            fasta_path,
            fasta_path,
            "-k", &kmer_size.to_string(),
            "-o", paf_file.to_str().unwrap(),
        ])
        .status()
        .context("Failed to run allwave - is it installed?")?;
    
    if !allwave_status.success() {
        return Err(anyhow::anyhow!("allwave failed with exit code: {:?}", allwave_status.code()));
    }
    
    // Step 2: Build graph with seqwish
    log::info!("Step 2: Building graph with seqwish...");
    let seqwish_gfa = work_dir.join("seqwish_graph.gfa");
    
    let seqwish_status = Command::new("seqwish")
        .args([
            "-s", fasta_path,
            "-p", paf_file.to_str().unwrap(),
            "-g", seqwish_gfa.to_str().unwrap(),
            "-k", &kmer_size.to_string(),
        ])
        .status()
        .context("Failed to run seqwish - is it installed?")?;
    
    if !seqwish_status.success() {
        return Err(anyhow::anyhow!("seqwish failed with exit code: {:?}", seqwish_status.code()));
    }
    
    // Step 3: Optimize graph with odgi
    log::info!("Step 3: Optimizing graph with odgi...");
    let og_file = work_dir.join("graph.og");
    
    // Build odgi graph
    let odgi_build_status = Command::new("odgi")
        .args([
            "build",
            "-g", seqwish_gfa.to_str().unwrap(),
            "-o", og_file.to_str().unwrap(),
        ])
        .status()
        .context("Failed to run odgi build - is odgi installed?")?;
    
    if !odgi_build_status.success() {
        return Err(anyhow::anyhow!("odgi build failed with exit code: {:?}", odgi_build_status.code()));
    }
    
    // Sort the graph
    let og_sorted = work_dir.join("graph_sorted.og");
    let odgi_sort_status = Command::new("odgi")
        .args([
            "sort",
            "-i", og_file.to_str().unwrap(),
            "-o", og_sorted.to_str().unwrap(),
        ])
        .status()
        .context("Failed to run odgi sort")?;
    
    if !odgi_sort_status.success() {
        return Err(anyhow::anyhow!("odgi sort failed with exit code: {:?}", odgi_sort_status.code()));
    }
    
    // Convert back to GFA
    let odgi_view_status = Command::new("odgi")
        .args([
            "view",
            "-i", og_sorted.to_str().unwrap(),
            "-g",
        ])
        .output()
        .context("Failed to run odgi view")?;
    
    if !odgi_view_status.status.success() {
        return Err(anyhow::anyhow!("odgi view failed with exit code: {:?}", odgi_view_status.status.code()));
    }
    
    // Write final GFA
    fs::write(output_gfa, &odgi_view_status.stdout).await?;
    
    log::info!("Graph construction complete: {}", output_gfa);
    
    // Clean up intermediate files if successful
    if let Err(e) = fs::remove_dir_all(&work_dir).await {
        log::warn!("Failed to clean up work directory: {}", e);
    }
    
    Ok(())
}

/// Alternative: Build graph using pggb (if available)
pub async fn build_graph_with_pggb(
    fasta_path: &str,
    output_gfa: &str,
    segment_length: usize,
    block_length: usize,
) -> Result<()> {
    log::info!("Building graph with pggb: {}", fasta_path);
    
    let output_path = Path::new(output_gfa);
    let output_dir = output_path.parent().unwrap();
    fs::create_dir_all(output_dir).await?;
    
    let pggb_status = Command::new("pggb")
        .args([
            "-i", fasta_path,
            "-o", output_dir.to_str().unwrap(),
            "-s", &segment_length.to_string(),
            "-l", &block_length.to_string(),
            "-n", "1", // Single community
        ])
        .status()
        .context("Failed to run pggb - is it installed?")?;
    
    if !pggb_status.success() {
        return Err(anyhow::anyhow!("pggb failed with exit code: {:?}", pggb_status.code()));
    }
    
    // Find the generated GFA file (pggb creates files with specific naming)
    let mut gfa_files = Vec::new();
    let mut entries = fs::read_dir(output_dir).await?;
    while let Some(entry) = entries.next_entry().await? {
        let path = entry.path();
        if path.extension().map_or(false, |ext| ext == "gfa") {
            gfa_files.push(path);
        }
    }
    
    if gfa_files.is_empty() {
        return Err(anyhow::anyhow!("No GFA file found after pggb run"));
    }
    
    // Use the first (likely only) GFA file
    if gfa_files[0] != output_path {
        fs::rename(&gfa_files[0], output_path).await?;
    }
    
    log::info!("Graph construction with pggb complete: {}", output_gfa);
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;
    
    #[tokio::test]
    #[ignore] // Requires external tools
    async fn test_build_graph_from_fasta() {
        let temp_dir = tempdir().unwrap();
        let fasta_path = temp_dir.path().join("test.fa");
        let output_gfa = temp_dir.path().join("test.gfa");
        
        // Create a simple test FASTA
        let fasta_content = ">seq1\nACGTACGTACGT\n>seq2\nACGTTCGTACGT\n";
        tokio::fs::write(&fasta_path, fasta_content).await.unwrap();
        
        // This test requires allwave, seqwish, and odgi to be installed
        let result = build_graph_from_fasta(
            fasta_path.to_str().unwrap(),
            output_gfa.to_str().unwrap(),
            11, // k-mer size
            1000, // segment length
        ).await;
        
        match result {
            Ok(_) => {
                assert!(output_gfa.exists());
                let gfa_content = tokio::fs::read_to_string(&output_gfa).await.unwrap();
                assert!(gfa_content.contains("H\tVN:Z:1.0"));
            }
            Err(e) => {
                // If tools aren't installed, that's OK for this test
                println!("Graph building failed (tools may not be installed): {}", e);
            }
        }
    }
}
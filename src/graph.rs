use std::process::Command;
use std::path::{Path, PathBuf};
use anyhow::{Result, Context};
use indicatif::{ProgressBar, ProgressStyle};

#[derive(Debug, Clone)]
pub struct GraphBuilder {
    input_fasta: PathBuf,
    output_dir: PathBuf,
    threads: usize,
    k_values: Vec<usize>,
}

impl GraphBuilder {
    pub fn new(input_fasta: impl AsRef<Path>, output_dir: impl AsRef<Path>) -> Self {
        Self {
            input_fasta: input_fasta.as_ref().to_path_buf(),
            output_dir: output_dir.as_ref().to_path_buf(),
            threads: 8,
            k_values: vec![51],  // Default k value
        }
    }
    
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }
    
    pub fn k_values(mut self, k_values: Vec<usize>) -> Self {
        self.k_values = k_values;
        self
    }
    
    pub fn build(&self) -> Result<GraphOutput> {
        // Create output directory
        std::fs::create_dir_all(&self.output_dir)?;
        
        let base_name = self.input_fasta
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| anyhow::anyhow!("Invalid input file name"))?
            .replace(".fa", "");
        
        // Step 1: Run allwave
        let paf_file = self.output_dir.join(format!("{}.paf", base_name));
        self.run_allwave(&paf_file)?;
        
        // Step 2: Build graphs for each k value
        let mut graphs = Vec::new();
        for &k in &self.k_values {
            let graph_info = self.build_graph_for_k(&base_name, &paf_file, k)?;
            graphs.push(graph_info);
        }
        
        Ok(GraphOutput {
            graphs,
            base_name: base_name.to_string(),
        })
    }
    
    fn run_allwave(&self, output_paf: &Path) -> Result<()> {
        log::info!("Running allwave alignment...");
        
        let output = Command::new("allwave")
            .args(&[
                "-i", self.input_fasta.to_str().unwrap(),
                "-t", &self.threads.to_string(),
                "-p", "none",
            ])
            .output()
            .context("Failed to run allwave")?;
        
        if !output.status.success() {
            return Err(anyhow::anyhow!(
                "allwave failed: {}",
                String::from_utf8_lossy(&output.stderr)
            ));
        }
        
        std::fs::write(output_paf, output.stdout)?;
        log::info!("Alignment complete: {}", output_paf.display());
        
        Ok(())
    }
    
    fn build_graph_for_k(&self, base_name: &str, paf_file: &Path, k: usize) -> Result<GraphInfo> {
        log::info!("Building graph with k={}...", k);
        
        let gfa_temp = self.output_dir.join(format!("{}.seqwish-k{}.gfa", base_name, k));
        let og_file = self.output_dir.join(format!("{}.k{}.og", base_name, k));
        let og_sorted = self.output_dir.join(format!("{}.k{}.sorted.og", base_name, k));
        let gfa_final = self.output_dir.join(format!("{}.k{}.gfa", base_name, k));
        let paths_file = self.output_dir.join(format!("{}.k{}.paths.tsv.gz", base_name, k));
        
        // Run seqwish
        log::info!("  Running seqwish...");
        let status = Command::new("seqwish")
            .args(&[
                "-s", self.input_fasta.to_str().unwrap(),
                "-g", gfa_temp.to_str().unwrap(),
                "-t", &self.threads.to_string(),
                "-p", paf_file.to_str().unwrap(),
                "-k", &k.to_string(),
                "-P",
            ])
            .status()
            .context("Failed to run seqwish")?;
        
        if !status.success() {
            return Err(anyhow::anyhow!("seqwish failed"));
        }
        
        // Build odgi file
        log::info!("  Building odgi graph...");
        let status = Command::new("odgi")
            .args(&[
                "build",
                "-g", gfa_temp.to_str().unwrap(),
                "-o", og_file.to_str().unwrap(),
            ])
            .status()
            .context("Failed to run odgi build")?;
        
        if !status.success() {
            return Err(anyhow::anyhow!("odgi build failed"));
        }
        
        // Sort the graph
        log::info!("  Sorting graph...");
        let status = Command::new("odgi")
            .args(&[
                "sort",
                "-i", og_file.to_str().unwrap(),
                "-p", "Ygs",
                "-o", og_sorted.to_str().unwrap(),
            ])
            .status()
            .context("Failed to run odgi sort")?;
        
        if !status.success() {
            return Err(anyhow::anyhow!("odgi sort failed"));
        }
        
        // Convert back to GFA
        log::info!("  Converting to GFA...");
        let output = Command::new("odgi")
            .args(&[
                "view",
                "-i", og_sorted.to_str().unwrap(),
                "-g",
            ])
            .output()
            .context("Failed to run odgi view")?;
        
        if !output.status.success() {
            return Err(anyhow::anyhow!("odgi view failed"));
        }
        
        std::fs::write(&gfa_final, output.stdout)?;
        
        // Extract paths coverage
        log::info!("  Extracting path coverage...");
        let output = Command::new("odgi")
            .args(&[
                "paths",
                "-i", og_sorted.to_str().unwrap(),
                "-L",
            ])
            .output()
            .context("Failed to run odgi paths")?;
        
        if !output.status.success() {
            return Err(anyhow::anyhow!("odgi paths failed"));
        }
        
        // Compress the paths output
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;
        
        let file = std::fs::File::create(&paths_file)?;
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(&output.stdout)?;
        encoder.finish()?;
        
        // Clean up temporary files
        let _ = std::fs::remove_file(&gfa_temp);
        let _ = std::fs::remove_file(&og_file);
        
        log::info!("Graph k={} complete", k);
        
        Ok(GraphInfo {
            k,
            gfa_path: gfa_final,
            og_path: og_sorted,
            paths_path: paths_file,
        })
    }
}

#[derive(Debug)]
pub struct GraphInfo {
    pub k: usize,
    pub gfa_path: PathBuf,
    pub og_path: PathBuf,
    pub paths_path: PathBuf,
}

#[derive(Debug)]
pub struct GraphOutput {
    pub graphs: Vec<GraphInfo>,
    pub base_name: String,
}

/// Extract sequences from graph (minus specified paths)
pub fn extract_sequences_without_paths(
    gfa_path: &Path,
    exclude_paths: &[String],
    output_fasta: &Path,
) -> Result<()> {
    log::info!("Extracting sequences from graph, excluding {} paths", exclude_paths.len());
    
    // First, get all path names from the graph
    let output = Command::new("odgi")
        .args(&[
            "paths",
            "-i", gfa_path.to_str().unwrap(),
            "-L",
        ])
        .output()
        .context("Failed to list paths")?;
    
    if !output.status.success() {
        return Err(anyhow::anyhow!("odgi paths failed"));
    }
    
    // Parse path names
    let all_paths: Vec<String> = String::from_utf8_lossy(&output.stdout)
        .lines()
        .map(|line| line.split('\t').next().unwrap_or("").to_string())
        .filter(|s| !s.is_empty())
        .collect();
    
    // Filter out excluded paths
    let keep_paths: Vec<String> = all_paths
        .into_iter()
        .filter(|p| !exclude_paths.contains(p))
        .collect();
    
    log::info!("Keeping {} paths", keep_paths.len());
    
    // Extract sequences for kept paths
    // This would need more complex processing with odgi
    // For now, we'll use a placeholder
    
    Ok(())
}
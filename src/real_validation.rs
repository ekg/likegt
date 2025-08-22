use std::process::Command;
use std::path::{Path, PathBuf};
use std::fs;
use std::io::{self, Write, BufRead, BufReader};
use anyhow::{Result, Context};

use crate::math::cosine_similarity;

/// Real validation using actual read simulation and mapping
pub struct RealValidator {
    graph_path: PathBuf,
    reference_fasta: PathBuf,
    coverage_depth: u32,
    read_length: u32,
    temp_dir: PathBuf,
}

impl RealValidator {
    pub fn new(
        graph_path: impl AsRef<Path>,
        reference_fasta: impl AsRef<Path>,
        coverage_depth: u32,
    ) -> Result<Self> {
        // Create temp directory for intermediate files
        let temp_dir = PathBuf::from(format!("temp_validation_{}", std::process::id()));
        fs::create_dir_all(&temp_dir)?;
        
        Ok(Self {
            graph_path: graph_path.as_ref().to_path_buf(),
            reference_fasta: reference_fasta.as_ref().to_path_buf(),
            coverage_depth,
            read_length: 150,
            temp_dir,
        })
    }
    
    /// Generate coverage for a sample using real read simulation
    pub fn generate_sample_coverage(&self, sample_name: &str) -> Result<Vec<f64>> {
        println!("Generating real coverage for {}...", sample_name);
        
        // Extract sample sequences
        let sample_fasta = self.temp_dir.join(format!("{}.fa", sample_name));
        self.extract_sample_sequences(sample_name, &sample_fasta)?;
        
        // Get sequence length
        let seq_length = self.get_sequence_length(&sample_fasta)?;
        
        // Calculate read pairs needed
        let total_bases = seq_length * self.coverage_depth as usize;
        let bases_per_pair = 2 * self.read_length as usize;
        let n_pairs = total_bases / bases_per_pair;
        
        println!("  Simulating {} read pairs for {}x coverage...", n_pairs, self.coverage_depth);
        
        // Simulate reads with wgsim
        let reads1 = self.temp_dir.join(format!("{}.1.fq", sample_name));
        let reads2 = self.temp_dir.join(format!("{}.2.fq", sample_name));
        
        let output = Command::new("wgsim")
            .args(&[
                "-1", &self.read_length.to_string(),
                "-2", &self.read_length.to_string(),
                "-N", &n_pairs.to_string(),
                "-e", "0.01",
                "-r", "0.001",
                sample_fasta.to_str().unwrap(),
                reads1.to_str().unwrap(),
                reads2.to_str().unwrap(),
            ])
            .output()
            .context("Failed to run wgsim")?;
        
        if !output.status.success() {
            anyhow::bail!("wgsim failed: {}", String::from_utf8_lossy(&output.stderr));
        }
        
        // Map reads with BWA
        println!("  Mapping reads...");
        let bam_file = self.temp_dir.join(format!("{}.bam", sample_name));
        
        // First ensure BWA index exists
        self.ensure_bwa_index()?;
        
        // Run BWA mem
        let bwa_output = Command::new("bwa")
            .args(&[
                "mem",
                "-t", "4",
                self.temp_dir.join("paths.fa").to_str().unwrap(),
                reads1.to_str().unwrap(),
                reads2.to_str().unwrap(),
            ])
            .output()
            .context("Failed to run bwa mem")?;
        
        if !bwa_output.status.success() {
            anyhow::bail!("bwa mem failed: {}", String::from_utf8_lossy(&bwa_output.stderr));
        }
        
        // Convert to BAM
        fs::write(&bam_file, bwa_output.stdout)?;
        
        // Convert BAM to GAF using gfainject
        println!("  Converting to GAF...");
        let gaf_file = self.temp_dir.join(format!("{}.gaf", sample_name));
        
        let gfa_path = self.graph_path.with_extension("gfa");
        if !gfa_path.exists() {
            anyhow::bail!("GFA file not found: {:?}", gfa_path);
        }
        
        let gfainject_output = Command::new("gfainject")
            .args(&[
                "--gfa", gfa_path.to_str().unwrap(),
                "--bam", bam_file.to_str().unwrap(),
            ])
            .output()
            .context("Failed to run gfainject")?;
        
        if !gfainject_output.status.success() {
            anyhow::bail!("gfainject failed: {}", String::from_utf8_lossy(&gfainject_output.stderr));
        }
        
        fs::write(&gaf_file, gfainject_output.stdout)?;
        
        // Get coverage with gafpack
        println!("  Computing coverage...");
        let gafpack_output = Command::new("gafpack")
            .args(&[
                "--gfa", gfa_path.to_str().unwrap(),
                "--gaf", gaf_file.to_str().unwrap(),
            ])
            .output()
            .context("Failed to run gafpack")?;
        
        if !gafpack_output.status.success() {
            anyhow::bail!("gafpack failed: {}", String::from_utf8_lossy(&gafpack_output.stderr));
        }
        
        // Parse coverage output
        let coverage = self.parse_coverage_output(&gafpack_output.stdout)?;
        
        // Clean up temp files
        let _ = fs::remove_file(reads1);
        let _ = fs::remove_file(reads2);
        let _ = fs::remove_file(bam_file);
        let _ = fs::remove_file(gaf_file);
        
        Ok(coverage)
    }
    
    /// Extract sequences for a sample from reference FASTA
    fn extract_sample_sequences(&self, sample_name: &str, output: &Path) -> Result<()> {
        // Get all sequence names for this sample
        let faidx_output = Command::new("samtools")
            .args(&["faidx", self.reference_fasta.to_str().unwrap()])
            .output()
            .context("Failed to run samtools faidx")?;
        
        let faidx_str = String::from_utf8_lossy(&faidx_output.stdout);
        let mut seq_names = Vec::new();
        
        for line in faidx_str.lines() {
            if line.starts_with(sample_name) {
                if let Some(name) = line.split('\t').next() {
                    seq_names.push(name.to_string());
                }
            }
        }
        
        if seq_names.is_empty() {
            anyhow::bail!("No sequences found for sample {}", sample_name);
        }
        
        // Extract sequences
        let mut args = vec!["faidx", self.reference_fasta.to_str().unwrap()];
        for name in &seq_names {
            args.push(name);
        }
        
        let extract_output = Command::new("samtools")
            .args(&args)
            .output()
            .context("Failed to extract sequences")?;
        
        if !extract_output.status.success() {
            anyhow::bail!("samtools faidx failed: {}", String::from_utf8_lossy(&extract_output.stderr));
        }
        
        fs::write(output, extract_output.stdout)?;
        Ok(())
    }
    
    /// Get total sequence length from FASTA file
    fn get_sequence_length(&self, fasta_path: &Path) -> Result<usize> {
        let content = fs::read_to_string(fasta_path)?;
        let mut total_length = 0;
        
        for line in content.lines() {
            if !line.starts_with('>') {
                total_length += line.trim().len();
            }
        }
        
        Ok(total_length)
    }
    
    /// Ensure BWA index exists for paths
    fn ensure_bwa_index(&self) -> Result<()> {
        let paths_fa = self.temp_dir.join("paths.fa");
        
        if !paths_fa.exists() {
            // Extract paths from graph
            let output = Command::new("odgi")
                .args(&[
                    "paths",
                    "-i", self.graph_path.to_str().unwrap(),
                    "-f",
                ])
                .output()
                .context("Failed to extract paths")?;
            
            if !output.status.success() {
                anyhow::bail!("odgi paths failed: {}", String::from_utf8_lossy(&output.stderr));
            }
            
            fs::write(&paths_fa, output.stdout)?;
            
            // Build BWA index
            let index_output = Command::new("bwa")
                .args(&["index", paths_fa.to_str().unwrap()])
                .output()
                .context("Failed to build BWA index")?;
            
            if !index_output.status.success() {
                anyhow::bail!("bwa index failed: {}", String::from_utf8_lossy(&index_output.stderr));
            }
        }
        
        Ok(())
    }
    
    /// Parse coverage output from gafpack
    fn parse_coverage_output(&self, output: &[u8]) -> Result<Vec<f64>> {
        let output_str = String::from_utf8_lossy(output);
        let lines: Vec<&str> = output_str.lines().collect();
        
        if lines.is_empty() {
            anyhow::bail!("Empty coverage output");
        }
        
        // Skip header line if present
        let data_line = if lines[0].starts_with("#sample") {
            if lines.len() < 2 {
                anyhow::bail!("No data in coverage output");
            }
            lines[1]
        } else {
            lines[0]
        };
        
        // Parse coverage values (skip sample name column)
        let mut coverage = Vec::new();
        let parts: Vec<&str> = data_line.split('\t').collect();
        
        for (i, part) in parts.iter().enumerate() {
            if i == 0 && !part.starts_with("node") {
                // Skip sample name column
                continue;
            }
            
            let value: f64 = part.parse()
                .with_context(|| format!("Failed to parse coverage value: {}", part))?;
            coverage.push(value);
        }
        
        Ok(coverage)
    }
    
    /// Run hold-0-out validation
    pub fn validate_hold0out(&self, samples: &[&str]) -> Result<f64> {
        println!("Running real hold-0-out validation...");
        
        let mut correct = 0;
        let mut total = 0;
        
        for &sample in samples {
            println!("\nTesting {}...", sample);
            
            // Generate coverage for this sample
            let test_coverage = self.generate_sample_coverage(sample)?;
            
            // Compare against all samples
            let mut best_match = String::new();
            let mut best_score = -1.0;
            
            for &ref_sample in samples {
                let ref_coverage = self.generate_sample_coverage(ref_sample)?;
                let score = cosine_similarity(&test_coverage, &ref_coverage);
                
                println!("  vs {}: {:.6}", ref_sample, score);
                
                if score > best_score {
                    best_score = score;
                    best_match = ref_sample.to_string();
                }
            }
            
            if best_match == sample {
                correct += 1;
                println!("  ✓ Correctly matched to itself");
            } else {
                println!("  ✗ Incorrectly matched to {}", best_match);
            }
            
            total += 1;
        }
        
        let accuracy = correct as f64 / total as f64;
        println!("\nHold-0-out accuracy: {}/{} = {:.2}%", correct, total, accuracy * 100.0);
        
        Ok(accuracy)
    }
}

impl Drop for RealValidator {
    fn drop(&mut self) {
        // Clean up temp directory
        let _ = fs::remove_dir_all(&self.temp_dir);
    }
}
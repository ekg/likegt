use anyhow::Result;
use std::path::Path;
use tokio::fs;

use crate::pipeline::build::build_graph_from_fasta;
use crate::pipeline::validate::{ValidationResults, ValidationConfig};

/// Run comprehensive validation tests on a FASTA file
/// 
/// This function:
/// 1. Builds a pangenome graph from the FASTA
/// 2. Runs hold-0-out validation (test individuals included in reference)
/// 3. Runs hold-2-out validation (test individuals excluded from reference) 
/// 4. Runs reference bias test (impact of reference choice)
/// 5. Generates a comprehensive report
pub async fn run_validation(
    fasta_path: &str,
    output_dir: &str,
    threads: usize,
    kmer_size: usize,
) -> Result<()> {
    log::info!("Starting validation for: {}", fasta_path);
    log::info!("Output directory: {}", output_dir);
    log::info!("Threads: {}, K-mer size: {}", threads, kmer_size);
    
    // Create output directory
    fs::create_dir_all(output_dir).await?;
    
    // Step 1: Build graph from FASTA
    let graph_path = Path::new(output_dir).join("graph.gfa");
    log::info!("Building graph from FASTA...");
    
    build_graph_from_fasta(
        fasta_path,
        graph_path.to_str().unwrap(),
        kmer_size,
        10000, // segment_length
    ).await?;
    
    // Step 2: Run validation tests
    let config = ValidationConfig {
        graph_path: graph_path.clone(),
        fasta_path: fasta_path.into(),
        output_dir: output_dir.into(),
        threads,
        kmer_size,
    };
    
    log::info!("Running hold-0-out validation...");
    let hold0_results = run_hold0_validation(&config).await?;
    
    log::info!("Running hold-2-out validation...");
    let hold2_results = run_hold2_validation(&config).await?;
    
    log::info!("Running reference bias test...");
    let bias_results = run_reference_bias_test(&config).await?;
    
    // Step 3: Generate comprehensive report
    let report = generate_validation_report(&hold0_results, &hold2_results, &bias_results)?;
    
    let report_path = Path::new(output_dir).join("validation_report.txt");
    fs::write(&report_path, &report).await?;
    
    log::info!("Validation complete. Report saved to: {}", report_path.display());
    println!("{}", report);
    
    Ok(())
}

async fn run_hold0_validation(config: &ValidationConfig) -> Result<ValidationResults> {
    // Hold-0-out: test individuals are included in the reference graph
    // Should achieve perfect or near-perfect accuracy
    
    log::info!("Hold-0-out: Testing with full reference graph");
    
    // Use the full graph for both reference and sample coverage
    let results = crate::pipeline::validate::run_hold0_test(
        &config.graph_path,
        &config.fasta_path,
        config.threads,
    ).await?;
    
    Ok(results)
}

async fn run_hold2_validation(config: &ValidationConfig) -> Result<ValidationResults> {
    // Hold-2-out: test individuals are excluded from the reference graph
    // More realistic scenario, accuracy should still be high for good graphs
    
    log::info!("Hold-2-out: Testing with reduced reference graph");
    
    let results = crate::pipeline::validate::run_hold2_test(
        &config.graph_path,
        &config.fasta_path,
        config.threads,
    ).await?;
    
    Ok(results)
}

async fn run_reference_bias_test(config: &ValidationConfig) -> Result<ValidationResults> {
    // Reference bias: test how choice of reference affects genotyping
    // Helps identify if the algorithm is biased toward specific haplotypes
    
    log::info!("Reference bias: Testing impact of reference selection");
    
    let results = crate::pipeline::validate::run_bias_test(
        &config.graph_path,
        &config.fasta_path,
        config.threads,
    ).await?;
    
    Ok(results)
}

fn generate_validation_report(
    hold0: &ValidationResults,
    hold2: &ValidationResults, 
    bias: &ValidationResults,
) -> Result<String> {
    let mut report = String::new();
    
    report.push_str("# Pangenome Graph Genotyping Validation Report\n\n");
    
    // Hold-0-out results
    report.push_str("## Hold-0-out Validation (Test individuals in reference)\n");
    report.push_str(&format!("- Total tests: {}\n", hold0.total_tests));
    report.push_str(&format!("- Correct genotypes: {}\n", hold0.correct_genotypes));
    report.push_str(&format!("- Accuracy: {:.2}%\n", hold0.accuracy * 100.0));
    report.push_str(&format!("- Average similarity: {:.4}\n", hold0.average_similarity));
    report.push_str(&format!("- Standard deviation: {:.4}\n\n", hold0.similarity_stddev));
    
    // Hold-2-out results  
    report.push_str("## Hold-2-out Validation (Test individuals excluded from reference)\n");
    report.push_str(&format!("- Total tests: {}\n", hold2.total_tests));
    report.push_str(&format!("- Correct genotypes: {}\n", hold2.correct_genotypes));
    report.push_str(&format!("- Accuracy: {:.2}%\n", hold2.accuracy * 100.0));
    report.push_str(&format!("- Average similarity: {:.4}\n", hold2.average_similarity));
    report.push_str(&format!("- Standard deviation: {:.4}\n\n", hold2.similarity_stddev));
    
    // Reference bias results
    report.push_str("## Reference Bias Test\n");
    report.push_str(&format!("- Total tests: {}\n", bias.total_tests));
    report.push_str(&format!("- Consistent results: {}\n", bias.correct_genotypes));
    report.push_str(&format!("- Consistency: {:.2}%\n", bias.accuracy * 100.0));
    report.push_str(&format!("- Bias variance: {:.4}\n\n", bias.similarity_stddev));
    
    // Overall assessment
    report.push_str("## Overall Assessment\n");
    
    if hold0.accuracy >= 0.95 && hold2.accuracy >= 0.80 && bias.accuracy >= 0.90 {
        report.push_str("✅ **EXCELLENT**: This graph is highly suitable for genotyping\n");
    } else if hold0.accuracy >= 0.90 && hold2.accuracy >= 0.70 && bias.accuracy >= 0.80 {
        report.push_str("⚠️  **GOOD**: This graph is suitable for genotyping with some limitations\n");
    } else if hold0.accuracy >= 0.80 || hold2.accuracy >= 0.60 {
        report.push_str("⚠️  **FAIR**: This graph may be suitable for genotyping but with caution\n");
    } else {
        report.push_str("❌ **POOR**: This graph is not suitable for accurate genotyping\n");
    }
    
    report.push_str("\n");
    report.push_str("### Recommendations\n");
    
    if hold0.accuracy < 0.95 {
        report.push_str("- Consider using a different k-mer size or graph construction parameters\n");
    }
    
    if hold2.accuracy < hold0.accuracy - 0.20 {
        report.push_str("- Large accuracy drop in hold-2-out suggests graph may be too sparse\n");
    }
    
    if bias.accuracy < 0.85 {
        report.push_str("- High reference bias detected - results may depend heavily on reference choice\n");
    }
    
    Ok(report)
}
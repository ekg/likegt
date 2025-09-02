use anyhow::Result;
use crate::{io, max_qv};

/// Run max QV analysis on a coverage matrix
pub async fn run_max_qv_analysis(
    coverage_file: &str,
    output_file: Option<&str>,
    ploidy: usize,
    verbose: bool,
) -> Result<()> {
    if verbose {
        println!("🔬 Computing maximum attainable QV for all samples...");
    }
    
    // Load coverage data
    let ref_data = io::read_gzip_tsv(coverage_file)?;
    
    if verbose {
        println!("📊 Loaded {} haplotypes from {}", ref_data.ids.len(), coverage_file);
    }
    
    // Compute max QV for all individuals
    let results = max_qv::compute_all_max_qvs(&ref_data, ploidy)?;
    
    if verbose {
        println!("✅ Computed max QV for {} individuals", results.len());
    }
    
    // Generate report
    let report = max_qv::generate_max_qv_report(&results);
    
    // Output results
    if let Some(output) = output_file {
        tokio::fs::write(output, &report).await?;
        println!("📝 Results written to {}", output);
    } else {
        println!("{}", report);
    }
    
    Ok(())
}

/// Compute max QV for a specific set of held-out individuals
pub async fn compute_max_qv_for_holdout(
    coverage_file: &str,
    held_out_individuals: &[String],
    ploidy: usize,
    verbose: bool,
) -> Result<Vec<max_qv::MaxQVResult>> {
    if verbose {
        println!("🔬 Computing max attainable QV for held-out individuals: {:?}", held_out_individuals);
    }
    
    // Load coverage data
    let ref_data = io::read_gzip_tsv(coverage_file)?;
    
    // Compute max QV for held-out individuals
    let results = max_qv::compute_max_qv_hold_n_out(&ref_data, held_out_individuals, ploidy)?;
    
    if verbose {
        for result in &results {
            println!("   {} -> Max QV: {:.1}, Best rank: {}, Best similarity: {:.4}",
                result.sample_id,
                result.max_attainable_qv,
                result.best_rank_achievable,
                result.best_non_self_similarity
            );
        }
    }
    
    Ok(results)
}
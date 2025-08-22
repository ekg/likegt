use anyhow::Result;
use itertools::Itertools;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Mutex;
use indicatif::{ProgressBar, ProgressStyle};

use crate::io::{read_gzip_tsv, write_genotype_result, write_sorted_combinations};
use crate::coverage::GenotypeData;
use crate::math::{sum_vectors, cosine_similarity};

#[derive(Debug, Clone)]
struct GenotypeResult {
    haplotypes: Vec<String>,
    similarity: f64,
}

pub fn run_genotyping(
    paths_file: &str,
    gaf_file: &str,
    output_dir: &str,
    sample_id: &str,
    ploidy: usize,
    threads: usize,
) -> Result<()> {
    log::info!("Starting genotyping for sample: {}", sample_id);
    log::info!("Ploidy: {}, Threads: {}", ploidy, threads);
    
    // Set up thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    
    // Read coverage data
    log::info!("Reading reference haplotype coverage from: {}", paths_file);
    let reference_data = read_gzip_tsv(paths_file)?;
    
    log::info!("Reading sample coverage from: {}", gaf_file);
    let sample_data = read_gzip_tsv(gaf_file)?;
    
    // Create genotype data
    let genotype_data = GenotypeData::new(
        reference_data,
        sample_data,
        sample_id.to_string(),
    )?;
    
    log::info!(
        "Loaded {} reference haplotypes with {} nodes",
        genotype_data.num_haplotypes(),
        genotype_data.num_nodes()
    );
    
    // Generate and evaluate combinations
    let results = evaluate_combinations(&genotype_data, ploidy)?;
    
    // Sort by similarity (descending)
    let mut sorted_results = results;
    sorted_results.sort_by(|a, b| {
        b.similarity.partial_cmp(&a.similarity).unwrap()
    });
    
    // Write results
    if !sorted_results.is_empty() {
        // Best genotype
        let best = &sorted_results[0];
        write_genotype_result(
            output_dir,
            sample_id,
            &best.haplotypes,
            best.similarity,
            ploidy,
        )?;
        
        // All combinations
        let all_combos: Vec<(Vec<String>, f64)> = sorted_results
            .iter()
            .map(|r| (r.haplotypes.clone(), r.similarity))
            .collect();
        
        write_sorted_combinations(
            output_dir,
            sample_id,
            &all_combos,
            ploidy,
        )?;
        
        log::info!(
            "Best genotype: {} (similarity: {:.6})",
            best.haplotypes.join(" + "),
            best.similarity
        );
    }
    
    Ok(())
}

fn evaluate_combinations(
    genotype_data: &GenotypeData,
    ploidy: usize,
) -> Result<Vec<GenotypeResult>> {
    let n = genotype_data.num_haplotypes();
    
    if n == 0 {
        return Err(anyhow::anyhow!("No haplotypes available for genotyping"));
    }
    
    // Calculate total number of combinations
    let total_combinations = if ploidy == 2 {
        n * (n + 1) / 2  // n choose 2 with replacement
    } else {
        // General formula for combinations with replacement
        calculate_combinations_with_replacement(n, ploidy)
    };
    
    log::info!("Evaluating {} combinations", total_combinations);
    
    // Progress bar
    let pb = ProgressBar::new(total_combinations as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")?
            .progress_chars("#>-")
    );
    
    let results = Mutex::new(Vec::new());
    let _seen_homozygous: Mutex<HashMap<usize, bool>> = Mutex::new(HashMap::new());
    
    // Generate all combinations with replacement
    let indices: Vec<usize> = (0..n).collect();
    
    // Use itertools to generate combinations with replacement
    let combinations: Vec<Vec<usize>> = indices
        .iter()
        .cloned()
        .combinations_with_replacement(ploidy)
        .collect();
    
    combinations.par_iter().for_each(|combo| {
        // Prepare coverage vectors for this combination
        let coverage_refs: Vec<&[f64]> = combo
            .iter()
            .map(|&idx| genotype_data.reference_haplotypes.coverages[idx].as_slice())
            .collect();
        
        // Sum the coverage vectors
        let combined_coverage = sum_vectors(&coverage_refs);
        
        // Calculate cosine similarity
        let similarity = cosine_similarity(
            &combined_coverage,
            &genotype_data.sample_coverage,
        );
        
        // Get haplotype names
        let haplotype_names: Vec<String> = combo
            .iter()
            .map(|&idx| genotype_data.reference_haplotypes.ids[idx].clone())
            .collect();
        
        // Store result
        results.lock().unwrap().push(GenotypeResult {
            haplotypes: haplotype_names,
            similarity,
        });
        
        pb.inc(1);
    });
    
    pb.finish_with_message("Evaluation complete");
    
    let final_results = results.into_inner().unwrap();
    log::info!("Evaluated {} unique combinations", final_results.len());
    
    Ok(final_results)
}

fn calculate_combinations_with_replacement(n: usize, k: usize) -> usize {
    // Formula: C(n+k-1, k) = (n+k-1)! / (k! * (n-1)!)
    let mut result = 1usize;
    for i in 0..k {
        result = result * (n + i) / (i + 1);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_combinations_count() {
        assert_eq!(calculate_combinations_with_replacement(3, 2), 6);
        assert_eq!(calculate_combinations_with_replacement(4, 2), 10);
        assert_eq!(calculate_combinations_with_replacement(5, 3), 35);
    }
}
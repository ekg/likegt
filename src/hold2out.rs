use anyhow::{Result, Context};
use std::path::Path;
use rayon::prelude::*;
use indicatif::{ProgressBar, ProgressStyle};
use crate::{io::CoverageData, math, sequence_qv};

#[derive(Debug, Clone)]
pub struct Hold2OutResult {
    pub individual: String,
    pub true_hap1: String,
    pub true_hap2: String,
    pub called_hap1: String,
    pub called_hap2: String,
    pub similarity: f64,
    pub rank: usize,
    pub correct: bool,
    pub sequence_qv: Option<f64>,
}

/// Perform hold-2-out validation for an individual
/// This removes the individual's haplotypes and tries to recover them
pub fn hold2out_individual(
    ref_data: &CoverageData,
    individual: &str,
) -> Result<Hold2OutResult> {
    // Find this individual's haplotypes
    let hap1_pattern = format!("{}#1", individual);
    let hap2_pattern = format!("{}#2", individual);
    
    let hap1_idx = ref_data.ids.iter()
        .position(|id| id.contains(&hap1_pattern))
        .context(format!("Haplotype 1 not found for {}", individual))?;
    
    let hap2_idx = ref_data.ids.iter()
        .position(|id| id.contains(&hap2_pattern))
        .context(format!("Haplotype 2 not found for {}", individual))?;
    
    // Create sample coverage by summing the two haplotypes
    let sample_coverage = math::sum_vectors(&[
        &ref_data.coverages[hap1_idx],
        &ref_data.coverages[hap2_idx],
    ]);
    
    // Create reference dataset WITHOUT this individual
    let mut held_out_ids = Vec::new();
    let mut held_out_coverages = Vec::new();
    
    for (i, id) in ref_data.ids.iter().enumerate() {
        if i != hap1_idx && i != hap2_idx {
            held_out_ids.push(id.clone());
            held_out_coverages.push(ref_data.coverages[i].clone());
        }
    }
    
    // Find best genotype from remaining haplotypes
    let mut results = Vec::new();
    
    for i in 0..held_out_ids.len() {
        for j in i..held_out_ids.len() {
            let combined = math::sum_vectors(&[
                &held_out_coverages[i],
                &held_out_coverages[j],
            ]);
            
            let similarity = math::cosine_similarity(&combined, &sample_coverage);
            
            results.push((
                vec![held_out_ids[i].clone(), held_out_ids[j].clone()],
                similarity,
            ));
        }
    }
    
    // Sort by similarity
    results.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    
    // Find rank of true genotype (should not be findable in hold-2-out!)
    let true_genotype = vec![ref_data.ids[hap1_idx].clone(), ref_data.ids[hap2_idx].clone()];
    let rank = results.iter()
        .position(|(haps, _)| {
            (haps[0] == true_genotype[0] && haps[1] == true_genotype[1]) ||
            (haps[0] == true_genotype[1] && haps[1] == true_genotype[0])
        })
        .map(|r| r + 1)
        .unwrap_or(0); // 0 means not found (expected!)
    
    // Best result should NOT be the true genotype in hold-2-out
    let best = &results[0];
    let correct = rank == 0; // In hold-2-out, "correct" means we didn't find the held-out individual
    
    Ok(Hold2OutResult {
        individual: individual.to_string(),
        true_hap1: ref_data.ids[hap1_idx].clone(),
        true_hap2: ref_data.ids[hap2_idx].clone(),
        called_hap1: best.0[0].clone(),
        called_hap2: best.0[1].clone(),
        similarity: best.1,
        rank,
        correct,
        sequence_qv: None,
    })
}

/// Run hold-2-out validation on all individuals
pub fn validate_all_hold2out(
    ref_data: &CoverageData,
) -> Result<Vec<Hold2OutResult>> {
    // Extract unique individuals
    let mut individuals = std::collections::HashSet::new();
    for id in &ref_data.ids {
        if let Some(ind) = id.split('#').next() {
            individuals.insert(ind.to_string());
        }
    }
    
    let mut all_individuals: Vec<_> = individuals.into_iter().collect();
    all_individuals.sort();
    
    println!("Running hold-2-out validation for {} individuals", all_individuals.len());
    
    // Create progress bar
    let pb = ProgressBar::new(all_individuals.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40} {pos}/{len} {msg}")
            .unwrap()
    );
    
    // Process all individuals in parallel
    let results: Vec<Result<Hold2OutResult>> = all_individuals
        .par_iter()
        .map(|individual| {
            pb.inc(1);
            pb.set_message(format!("Processing {}", individual));
            hold2out_individual(ref_data, individual)
        })
        .collect();
    
    pb.finish_with_message("Hold-2-out validation complete");
    
    // Collect successful results
    let mut successful_results = Vec::new();
    let mut failed_count = 0;
    
    for result in results {
        match result {
            Ok(r) => successful_results.push(r),
            Err(e) => {
                eprintln!("Failed: {}", e);
                failed_count += 1;
            }
        }
    }
    
    if failed_count > 0 {
        println!("Warning: {} individuals failed validation", failed_count);
    }
    
    Ok(successful_results)
}

/// Compare hold-0-out vs hold-2-out accuracy
pub fn compare_validation_methods(
    ref_data: &CoverageData,
) -> Result<()> {
    println!("\n=== COMPARING HOLD-0-OUT VS HOLD-2-OUT ===\n");
    
    // Run hold-0-out validation
    println!("Running hold-0-out validation...");
    let hold0_results = crate::validation::validate_all_individuals(ref_data)?;
    let hold0_accuracy = hold0_results.iter()
        .filter(|r| r.correct)
        .count() as f64 / hold0_results.len() as f64;
    
    // Run hold-2-out validation
    println!("Running hold-2-out validation...");
    let hold2_results = validate_all_hold2out(ref_data)?;
    
    // For hold-2-out, we need to check if the called genotype is "close"
    // to the true genotype based on similarity
    let mut hold2_close_calls = 0;
    let mut similarity_threshold = 0.95; // Consider "close" if similarity > 0.95
    
    for result in &hold2_results {
        if result.similarity > similarity_threshold {
            hold2_close_calls += 1;
        }
    }
    
    let hold2_close_rate = hold2_close_calls as f64 / hold2_results.len() as f64;
    
    // Print comparison
    println!("\n=== RESULTS COMPARISON ===");
    println!("Hold-0-out accuracy: {:.1}% ({}/{})", 
        hold0_accuracy * 100.0,
        hold0_results.iter().filter(|r| r.correct).count(),
        hold0_results.len()
    );
    
    println!("Hold-2-out close calls: {:.1}% ({}/{})",
        hold2_close_rate * 100.0,
        hold2_close_calls,
        hold2_results.len()
    );
    
    // Find cases where hold-2-out finds very different genotypes
    println!("\n=== HOLD-2-OUT DIVERGENT CASES ===");
    println!("Cases where hold-2-out similarity < 0.95:");
    
    let mut divergent_cases: Vec<_> = hold2_results.iter()
        .filter(|r| r.similarity < similarity_threshold)
        .collect();
    
    divergent_cases.sort_by(|a, b| a.similarity.partial_cmp(&b.similarity).unwrap());
    
    for (i, case) in divergent_cases.iter().take(10).enumerate() {
        println!("{}. {} - similarity: {:.4}", 
            i + 1,
            case.individual,
            case.similarity
        );
        println!("   True:   {} + {}", 
            case.true_hap1.split(':').next().unwrap_or("?"),
            case.true_hap2.split(':').next().unwrap_or("?")
        );
        println!("   Called: {} + {}", 
            case.called_hap1.split(':').next().unwrap_or("?"),
            case.called_hap2.split(':').next().unwrap_or("?")
        );
    }
    
    // Summary statistics
    println!("\n=== SUMMARY STATISTICS ===");
    
    let hold2_similarities: Vec<f64> = hold2_results.iter()
        .map(|r| r.similarity)
        .collect();
    
    let mean_similarity = hold2_similarities.iter().sum::<f64>() / hold2_similarities.len() as f64;
    let min_similarity = hold2_similarities.iter()
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap_or(&0.0);
    let max_similarity = hold2_similarities.iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap_or(&1.0);
    
    println!("Hold-2-out similarity distribution:");
    println!("  Mean:   {:.4}", mean_similarity);
    println!("  Min:    {:.4}", min_similarity);
    println!("  Max:    {:.4}", max_similarity);
    
    // Histogram of similarities
    println!("\nSimilarity histogram:");
    let mut bins = vec![0; 10];
    for sim in &hold2_similarities {
        let bin = ((sim * 10.0).floor() as usize).min(9);
        bins[bin] += 1;
    }
    
    for (i, count) in bins.iter().enumerate() {
        let start = i as f64 / 10.0;
        let end = (i + 1) as f64 / 10.0;
        let bar = "█".repeat((*count as f64 / 5.0).ceil() as usize);
        println!("  {:.1}-{:.1}: {:3} {}", start, end, count, bar);
    }
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_hold2out_removes_individual() {
        // Create test data
        let mut ref_data = CoverageData::new();
        
        // Add test haplotypes
        ref_data.ids.push("HG001#1".to_string());
        ref_data.coverages.push(vec![1.0, 0.0, 1.0, 0.0]);
        
        ref_data.ids.push("HG001#2".to_string());
        ref_data.coverages.push(vec![0.0, 1.0, 0.0, 1.0]);
        
        ref_data.ids.push("HG002#1".to_string());
        ref_data.coverages.push(vec![1.0, 1.0, 0.0, 0.0]);
        
        ref_data.ids.push("HG002#2".to_string());
        ref_data.coverages.push(vec![0.0, 0.0, 1.0, 1.0]);
        
        // Run hold-2-out for HG001
        let result = hold2out_individual(&ref_data, "HG001").unwrap();
        
        // Verify the called genotype is NOT HG001
        assert!(!result.called_hap1.contains("HG001"));
        assert!(!result.called_hap2.contains("HG001"));
        
        // Should call HG002 as the best match
        assert!(result.called_hap1.contains("HG002") || result.called_hap2.contains("HG002"));
        
        // Rank should be 0 (not found in held-out set)
        assert_eq!(result.rank, 0);
    }
}
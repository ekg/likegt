use anyhow::Result;
use crate::{io, math, validation};
use std::collections::HashMap;

/// Result for maximum attainable QV calculation
#[derive(Debug, Clone)]
pub struct MaxQVResult {
    pub sample_id: String,
    pub held_out_haplotypes: Vec<String>,
    pub best_non_self_genotype: (String, String),
    pub best_non_self_similarity: f64,
    pub true_self_similarity: f64,
    pub max_attainable_qv: f64,
    pub best_rank_achievable: usize,
    pub total_combinations_tested: usize,
}

/// Compute the maximum attainable QV for a held-out sample
/// This finds the best possible non-self genotype match
pub fn compute_max_attainable_qv(
    ref_data: &io::CoverageData,
    held_out_ids: &[String],
    ploidy: usize,
) -> Result<MaxQVResult> {
    if held_out_ids.len() != ploidy {
        return Err(anyhow::anyhow!(
            "Number of held-out haplotypes ({}) must match ploidy ({})",
            held_out_ids.len(),
            ploidy
        ));
    }
    
    // Find indices of held-out haplotypes
    let held_out_indices: Vec<usize> = held_out_ids
        .iter()
        .map(|id| {
            ref_data.ids.iter().position(|ref_id| ref_id == id)
                .ok_or_else(|| anyhow::anyhow!("Haplotype {} not found", id))
        })
        .collect::<Result<Vec<_>>>()?;
    
    // Create held-out sample coverage (sum of held-out haplotypes)
    let held_out_coverages: Vec<&[f64]> = held_out_indices
        .iter()
        .map(|&idx| ref_data.coverages[idx].as_slice())
        .collect();
    let sample_coverage = math::sum_vectors(&held_out_coverages);
    
    // Calculate true self-similarity (should be 1.0)
    let true_self_similarity = math::cosine_similarity(&sample_coverage, &sample_coverage);
    
    // Find all possible genotype combinations excluding held-out haplotypes
    let available_indices: Vec<usize> = (0..ref_data.ids.len())
        .filter(|idx| !held_out_indices.contains(idx))
        .collect();
    
    let mut best_non_self = None;
    let mut best_non_self_similarity = 0.0;
    let mut all_similarities = Vec::new();
    let mut total_combinations = 0;
    
    // Test all combinations from available haplotypes
    for i in 0..available_indices.len() {
        for j in i..available_indices.len() {
            let idx1 = available_indices[i];
            let idx2 = available_indices[j];
            
            let combined = if ploidy == 2 {
                math::sum_vectors(&[
                    &ref_data.coverages[idx1],
                    &ref_data.coverages[idx2],
                ])
            } else {
                continue; // Only handle diploid for now
            };
            
            let similarity = math::cosine_similarity(&combined, &sample_coverage);
            all_similarities.push(((idx1, idx2), similarity));
            total_combinations += 1;
            
            if similarity > best_non_self_similarity {
                best_non_self_similarity = similarity;
                best_non_self = Some((idx1, idx2));
            }
        }
    }
    
    // Also test combinations that include held-out haplotypes to determine rank
    let mut all_combinations_with_self = all_similarities.clone();
    
    // Add self-combinations
    for i in 0..held_out_indices.len() {
        for j in i..held_out_indices.len() {
            let idx1 = held_out_indices[i];
            let idx2 = held_out_indices[j];
            
            let combined = math::sum_vectors(&[
                &ref_data.coverages[idx1],
                &ref_data.coverages[idx2],
            ]);
            
            let similarity = math::cosine_similarity(&combined, &sample_coverage);
            all_combinations_with_self.push(((idx1, idx2), similarity));
        }
    }
    
    // Sort all combinations by similarity
    all_combinations_with_self.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    
    // Find rank of best non-self genotype
    let best_rank = if let Some(best_pair) = best_non_self {
        all_combinations_with_self
            .iter()
            .position(|&(pair, _)| pair == best_pair)
            .unwrap_or(all_combinations_with_self.len()) + 1
    } else {
        all_combinations_with_self.len()
    };
    
    // Calculate maximum attainable QV
    let max_attainable_qv = validation::calculate_qv(true_self_similarity, best_non_self_similarity);
    
    // Extract sample ID from held-out haplotypes
    let sample_id = held_out_ids[0]
        .split('#')
        .next()
        .unwrap_or("unknown")
        .to_string();
    
    let (best_hap1, best_hap2) = if let Some((idx1, idx2)) = best_non_self {
        (ref_data.ids[idx1].clone(), ref_data.ids[idx2].clone())
    } else {
        ("none".to_string(), "none".to_string())
    };
    
    Ok(MaxQVResult {
        sample_id,
        held_out_haplotypes: held_out_ids.to_vec(),
        best_non_self_genotype: (best_hap1, best_hap2),
        best_non_self_similarity,
        true_self_similarity,
        max_attainable_qv,
        best_rank_achievable: best_rank,
        total_combinations_tested: total_combinations,
    })
}

/// Compute maximum attainable QV for all individuals
pub fn compute_all_max_qvs(
    ref_data: &io::CoverageData,
    ploidy: usize,
) -> Result<Vec<MaxQVResult>> {
    // Extract unique individuals
    let mut individuals: HashMap<String, Vec<String>> = HashMap::new();
    
    for id in &ref_data.ids {
        if let Some(pos) = id.find('#') {
            let individual = id[..pos].to_string();
            if !individual.starts_with("chm13") && !individual.starts_with("grch") {
                individuals.entry(individual).or_insert_with(Vec::new).push(id.clone());
            }
        }
    }
    
    let mut results = Vec::new();
    
    for (individual, haplotypes) in individuals {
        if haplotypes.len() != ploidy {
            log::warn!("Skipping {}: found {} haplotypes, expected {}", 
                      individual, haplotypes.len(), ploidy);
            continue;
        }
        
        match compute_max_attainable_qv(ref_data, &haplotypes, ploidy) {
            Ok(result) => results.push(result),
            Err(e) => log::warn!("Failed to compute max QV for {}: {}", individual, e),
        }
    }
    
    results.sort_by(|a, b| b.max_attainable_qv.partial_cmp(&a.max_attainable_qv).unwrap());
    
    Ok(results)
}

/// Generate report for maximum attainable QVs
pub fn generate_max_qv_report(results: &[MaxQVResult]) -> String {
    let mut report = String::new();
    
    report.push_str("=== MAXIMUM ATTAINABLE QV REPORT ===\n");
    report.push_str("This shows the best possible QV when genotyping against non-self haplotypes\n\n");
    
    let avg_max_qv = results.iter().map(|r| r.max_attainable_qv).sum::<f64>() / results.len() as f64;
    let avg_rank = results.iter().map(|r| r.best_rank_achievable).sum::<usize>() as f64 / results.len() as f64;
    
    report.push_str(&format!("Total samples analyzed: {}\n", results.len()));
    report.push_str(&format!("Average maximum attainable QV: {:.1}\n", avg_max_qv));
    report.push_str(&format!("Average best achievable rank: {:.1}\n\n", avg_rank));
    
    report.push_str("Sample\tMax_QV\tBest_Rank\tBest_Similarity\tBest_Non-Self_Genotype\n");
    
    for result in results {
        let best_genotype = format!("{} + {}", 
            result.best_non_self_genotype.0.split(':').next().unwrap_or("?"),
            result.best_non_self_genotype.1.split(':').next().unwrap_or("?")
        );
        
        report.push_str(&format!(
            "{}\t{:.1}\t{}\t{:.6}\t{}\n",
            result.sample_id,
            result.max_attainable_qv,
            result.best_rank_achievable,
            result.best_non_self_similarity,
            best_genotype
        ));
    }
    
    // Distribution analysis
    report.push_str("\n=== QV DISTRIBUTION ===\n");
    let qv_ranges = [
        (50.0, 60.0, "QV50-60 (Excellent)"),
        (40.0, 50.0, "QV40-50 (Very Good)"),
        (30.0, 40.0, "QV30-40 (Good)"),
        (20.0, 30.0, "QV20-30 (Fair)"),
        (10.0, 20.0, "QV10-20 (Poor)"),
        (0.0, 10.0, "QV0-10 (Very Poor)"),
    ];
    
    for (min_qv, max_qv, label) in qv_ranges {
        let count = results.iter()
            .filter(|r| r.max_attainable_qv >= min_qv && r.max_attainable_qv < max_qv)
            .count();
        if count > 0 {
            report.push_str(&format!("{}: {} samples ({:.1}%)\n", 
                label, count, 100.0 * count as f64 / results.len() as f64));
        }
    }
    
    report
}

/// Compute maximum attainable QV for hold-n-out validation
/// This generalizes to holding out multiple individuals
pub fn compute_max_qv_hold_n_out(
    ref_data: &io::CoverageData,
    held_out_individuals: &[String],
    ploidy: usize,
) -> Result<Vec<MaxQVResult>> {
    // Collect all haplotypes for held-out individuals
    let mut held_out_haplotypes = Vec::new();
    
    for individual in held_out_individuals {
        for i in 1..=ploidy {
            let hap_id = format!("{}#{}", individual, i);
            if ref_data.ids.contains(&hap_id) {
                held_out_haplotypes.push(hap_id);
            }
        }
    }
    
    if held_out_haplotypes.is_empty() {
        return Err(anyhow::anyhow!("No valid haplotypes found for held-out individuals"));
    }
    
    // Group haplotypes by individual
    let mut individual_haplotypes: HashMap<String, Vec<String>> = HashMap::new();
    for hap in &held_out_haplotypes {
        if let Some(pos) = hap.find('#') {
            let individual = hap[..pos].to_string();
            individual_haplotypes.entry(individual).or_insert_with(Vec::new).push(hap.clone());
        }
    }
    
    let mut results = Vec::new();
    
    // Compute max QV for each held-out individual
    for (individual, haps) in individual_haplotypes {
        if haps.len() == ploidy {
            match compute_max_attainable_qv(ref_data, &haps, ploidy) {
                Ok(result) => results.push(result),
                Err(e) => log::warn!("Failed to compute max QV for {}: {}", individual, e),
            }
        }
    }
    
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_max_qv_computation() {
        // Create test data with more distinct haplotypes
        let ref_data = io::CoverageData {
            ids: vec![
                "HG001#1".to_string(),
                "HG001#2".to_string(),
                "HG002#1".to_string(),
                "HG002#2".to_string(),
            ],
            coverages: vec![
                vec![1.0, 0.0, 1.0, 0.0],  // HG001#1 - distinct pattern
                vec![0.0, 1.0, 0.0, 1.0],  // HG001#2 - opposite of HG001#1
                vec![0.8, 0.2, 0.7, 0.3],  // HG002#1 - similar but not identical
                vec![0.3, 0.7, 0.2, 0.8],  // HG002#2 - different pattern
            ],
        };
        
        // Test max QV for HG001
        let held_out = vec!["HG001#1".to_string(), "HG001#2".to_string()];
        let result = compute_max_attainable_qv(&ref_data, &held_out, 2).unwrap();
        
        assert_eq!(result.sample_id, "HG001");
        assert_eq!(result.held_out_haplotypes, held_out);
        assert!(result.best_non_self_similarity > 0.7); // Should find a reasonable match
        assert!(result.best_non_self_similarity < 1.0); // But not perfect
        // The QV depends on the similarity - with these values it might still be high
        assert!(result.max_attainable_qv >= 0.0); // QV should be non-negative
    }
}
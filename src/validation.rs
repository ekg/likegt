use crate::{io, math};
use anyhow::Result;
// use std::collections::HashMap;

/// Quality Value calculation for genotype calls
/// QV = -10 * log10(error_probability)
pub fn calculate_qv(correct_similarity: f64, called_similarity: f64) -> f64 {
    if called_similarity >= correct_similarity {
        // Perfect or better call
        60.0 // Cap at QV 60 (1 in 1,000,000 error rate)
    } else {
        // Calculate error based on similarity difference
        let error_prob = 1.0 - (called_similarity / correct_similarity).min(1.0);
        if error_prob <= 0.0 {
            60.0
        } else {
            (-10.0 * error_prob.log10()).max(0.0).min(60.0)
        }
    }
}

/// Hold-out validation result
#[derive(Debug, Clone)]
pub struct ValidationResult {
    pub sample_id: String,
    pub true_hap1: String,
    pub true_hap2: String,
    pub called_hap1: String,
    pub called_hap2: String,
    pub true_similarity: f64,
    pub called_similarity: f64,
    pub rank: usize, // Rank of true genotype in results
    pub qv: f64,
    pub correct: bool,
}

impl ValidationResult {
    pub fn is_perfect(&self) -> bool {
        self.correct && (self.called_similarity - self.true_similarity).abs() < 1e-10
    }
}

/// Perform hold-out validation for a specific individual
pub fn hold_out_validation(
    ref_data: &io::CoverageData,
    individual_id: &str,
    ploidy: usize,
) -> Result<ValidationResult> {
    // Find the individual's haplotypes
    let hap1_pattern = format!("{}#1", individual_id);
    let hap2_pattern = format!("{}#2", individual_id);
    
    let hap1_idx = ref_data.ids.iter().position(|id| id.contains(&hap1_pattern));
    let hap2_idx = ref_data.ids.iter().position(|id| id.contains(&hap2_pattern));
    
    if hap1_idx.is_none() || hap2_idx.is_none() {
        return Err(anyhow::anyhow!("Individual {} not found in dataset", individual_id));
    }
    
    let hap1_idx = hap1_idx.unwrap();
    let hap2_idx = hap2_idx.unwrap();
    
    // Create sample coverage as sum of the two haplotypes
    let sample_coverage = math::sum_vectors(&[
        &ref_data.coverages[hap1_idx],
        &ref_data.coverages[hap2_idx],
    ]);
    
    // Calculate true similarity (should be 1.0 for hold-0-out)
    let true_similarity = math::cosine_similarity(&sample_coverage, &sample_coverage);
    
    // Now find the best genotype from all combinations
    let mut results = Vec::new();
    
    for i in 0..ref_data.ids.len() {
        for j in i..ref_data.ids.len() {
            let combined = if ploidy == 2 {
                math::sum_vectors(&[
                    &ref_data.coverages[i],
                    &ref_data.coverages[j],
                ])
            } else {
                continue; // Only handle diploid for now
            };
            
            let similarity = math::cosine_similarity(&combined, &sample_coverage);
            results.push(((i, j), similarity));
        }
    }
    
    // Sort by similarity
    results.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    
    // Find rank of true genotype
    let true_pair = if hap1_idx <= hap2_idx {
        (hap1_idx, hap2_idx)
    } else {
        (hap2_idx, hap1_idx)
    };
    
    let rank = results.iter().position(|&(pair, _)| pair == true_pair).unwrap_or(results.len()) + 1;
    
    // Get the called (best) genotype
    let best = results[0];
    let called_hap1 = &ref_data.ids[best.0.0];
    let called_hap2 = &ref_data.ids[best.0.1];
    
    // Check if correct
    let correct = best.0 == true_pair;
    
    // Calculate QV
    let qv = calculate_qv(true_similarity, best.1);
    
    Ok(ValidationResult {
        sample_id: individual_id.to_string(),
        true_hap1: ref_data.ids[hap1_idx].clone(),
        true_hap2: ref_data.ids[hap2_idx].clone(),
        called_hap1: called_hap1.clone(),
        called_hap2: called_hap2.clone(),
        true_similarity,
        called_similarity: best.1,
        rank,
        qv,
        correct,
    })
}

/// Run systematic validation on all individuals
pub fn validate_all_individuals(
    ref_data: &io::CoverageData,
) -> Result<Vec<ValidationResult>> {
    // Extract unique individual IDs
    let mut individuals = std::collections::HashSet::new();
    for id in &ref_data.ids {
        if let Some(pos) = id.find('#') {
            let individual = &id[..pos];
            if !individual.starts_with("chm13") && !individual.starts_with("grch") {
                individuals.insert(individual.to_string());
            }
        }
    }
    
    let mut all_results = Vec::new();
    
    for individual in individuals {
        match hold_out_validation(ref_data, &individual, 2) {
            Ok(result) => all_results.push(result),
            Err(e) => log::warn!("Failed to validate {}: {}", individual, e),
        }
    }
    
    all_results.sort_by(|a, b| b.qv.partial_cmp(&a.qv).unwrap());
    
    Ok(all_results)
}

/// Generate accuracy report
pub fn generate_accuracy_report(results: &[ValidationResult]) -> String {
    let total = results.len();
    let correct = results.iter().filter(|r| r.correct).count();
    let perfect = results.iter().filter(|r| r.is_perfect()).count();
    
    let avg_qv = results.iter().map(|r| r.qv).sum::<f64>() / total as f64;
    let avg_rank = results.iter().map(|r| r.rank).sum::<usize>() as f64 / total as f64;
    
    let mut report = String::new();
    report.push_str(&format!("=== GENOTYPING ACCURACY REPORT ===\n"));
    report.push_str(&format!("Total individuals tested: {}\n", total));
    report.push_str(&format!("Correct genotypes: {} ({:.1}%)\n", correct, 100.0 * correct as f64 / total as f64));
    report.push_str(&format!("Perfect genotypes: {} ({:.1}%)\n", perfect, 100.0 * perfect as f64 / total as f64));
    report.push_str(&format!("Average QV: {:.1}\n", avg_qv));
    report.push_str(&format!("Average rank of true genotype: {:.1}\n", avg_rank));
    
    report.push_str("\n=== INDIVIDUAL RESULTS ===\n");
    report.push_str("Sample\tCorrect\tQV\tRank\tSimilarity\n");
    
    for result in results {
        report.push_str(&format!(
            "{}\t{}\t{:.1}\t{}\t{:.6}\n",
            result.sample_id,
            if result.correct { "YES" } else { "NO" },
            result.qv,
            result.rank,
            result.called_similarity
        ));
    }
    
    // Find problematic cases
    let failures: Vec<_> = results.iter().filter(|r| !r.correct).collect();
    if !failures.is_empty() {
        report.push_str("\n=== FAILED GENOTYPES ===\n");
        for fail in failures {
            report.push_str(&format!(
                "{}: Expected {} + {}, Got {} + {} (rank {})\n",
                fail.sample_id,
                fail.true_hap1.split(':').next().unwrap_or(&fail.true_hap1),
                fail.true_hap2.split(':').next().unwrap_or(&fail.true_hap2),
                fail.called_hap1.split(':').next().unwrap_or(&fail.called_hap1),
                fail.called_hap2.split(':').next().unwrap_or(&fail.called_hap2),
                fail.rank
            ));
        }
    }
    
    report
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_qv_calculation() {
        // Perfect call
        assert_eq!(calculate_qv(1.0, 1.0), 60.0);
        
        // Very good call
        let qv = calculate_qv(1.0, 0.99);
        assert!(qv > 19.0 && qv < 21.0); // ~QV20 for 1% error
        
        // Poor call
        let qv = calculate_qv(1.0, 0.9);
        assert!(qv > 9.0 && qv < 11.0); // ~QV10 for 10% error
        
        // Bad call
        let qv = calculate_qv(1.0, 0.5);
        assert!(qv > 2.0 && qv < 4.0); // ~QV3 for 50% error
    }
}
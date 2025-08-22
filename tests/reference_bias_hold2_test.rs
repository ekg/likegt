use likegt::{io, math, hold2out};
use std::path::PathBuf;

/// Comprehensive test of reference bias impact on hold-2-out validation
/// This is the most realistic test: simulates what happens when we don't have
/// the true genotypes in our reference panel AND we have reference-biased reads
#[test]
fn test_reference_bias_hold2out_comprehensive() {
    let test_data_path = PathBuf::from("tests/data/hla-f.k51.paths.coverage.tsv.gz");
    if !test_data_path.exists() {
        println!("Test data not found, skipping");
        return;
    }
    
    println!("\n=== COMPREHENSIVE REFERENCE BIAS HOLD-2-OUT TEST ===\n");
    
    let ref_data = io::read_gzip_tsv(test_data_path.to_str().unwrap()).unwrap();
    println!("Loaded {} haplotypes", ref_data.len());
    
    // Find GRCh38 reference
    let grch38_idx = ref_data.ids.iter()
        .position(|id| id.starts_with("grch38"))
        .expect("GRCh38 should be in reference panel");
    
    println!("Using {} as reference for bias\n", ref_data.ids[grch38_idx]);
    
    // Test multiple individuals
    let test_individuals = vec!["HG00096", "HG00268", "HG00733", "NA12329", "HG02818"];
    let bias_strength = 0.3; // 30% bias
    
    let mut results_unbiased = Vec::new();
    let mut results_biased = Vec::new();
    
    for individual in &test_individuals {
        println!("Testing {} with hold-2-out:", individual);
        
        // Get individual's haplotype indices
        let hap1_idx = ref_data.ids.iter()
            .position(|id| id.contains(&format!("{}#1", individual)));
        let hap2_idx = ref_data.ids.iter()
            .position(|id| id.contains(&format!("{}#2", individual)));
        
        if let (Some(idx1), Some(idx2)) = (hap1_idx, hap2_idx) {
            // Create true coverage
            let true_coverage = math::sum_vectors(&[
                &ref_data.coverages[idx1],
                &ref_data.coverages[idx2],
            ]);
            
            // Create held-out reference panel (remove this individual)
            let mut held_out_data = io::CoverageData::new();
            for (i, id) in ref_data.ids.iter().enumerate() {
                if i != idx1 && i != idx2 {
                    held_out_data.ids.push(id.clone());
                    held_out_data.coverages.push(ref_data.coverages[i].clone());
                }
            }
            
            println!("  Reference panel: {} haplotypes (removed 2)", held_out_data.len());
            
            // Test 1: Unbiased hold-2-out
            let unbiased_result = find_best_genotype_hold2(
                &held_out_data,
                &true_coverage,
                individual
            );
            results_unbiased.push(unbiased_result.clone());
            
            // Test 2: Biased hold-2-out
            let biased_coverage = create_biased_coverage(
                &true_coverage,
                &ref_data.coverages[grch38_idx],
                bias_strength
            );
            
            let biased_result = find_best_genotype_hold2(
                &held_out_data,
                &biased_coverage,
                individual
            );
            results_biased.push(biased_result.clone());
            
            // Print results
            println!("  Unbiased: similarity={:.4}, called={}", 
                unbiased_result.similarity,
                unbiased_result.called_genotype
            );
            println!("  Biased:   similarity={:.4}, called={}",
                biased_result.similarity,
                biased_result.called_genotype
            );
            
            // Check if we found the same or different genotype
            if unbiased_result.called_genotype != biased_result.called_genotype {
                println!("  -> BIAS CHANGED GENOTYPE!");
            }
            
            // Check similarity drop
            let similarity_drop = unbiased_result.similarity - biased_result.similarity;
            if similarity_drop > 0.001 {
                println!("  -> Similarity dropped by {:.4}", similarity_drop);
            }
            
            println!();
        }
    }
    
    // Analyze overall impact
    println!("=== HOLD-2-OUT BIAS IMPACT SUMMARY ===\n");
    
    // Compare similarities
    let avg_unbiased_sim: f64 = results_unbiased.iter()
        .map(|r| r.similarity)
        .sum::<f64>() / results_unbiased.len() as f64;
    
    let avg_biased_sim: f64 = results_biased.iter()
        .map(|r| r.similarity)
        .sum::<f64>() / results_biased.len() as f64;
    
    println!("Average similarities:");
    println!("  Unbiased: {:.4}", avg_unbiased_sim);
    println!("  Biased:   {:.4}", avg_biased_sim);
    println!("  Drop:     {:.4}", avg_unbiased_sim - avg_biased_sim);
    
    // Count genotype changes
    let genotype_changes = results_unbiased.iter()
        .zip(results_biased.iter())
        .filter(|(u, b)| u.called_genotype != b.called_genotype)
        .count();
    
    println!("\nGenotype changes due to bias: {}/{} ({:.1}%)",
        genotype_changes,
        results_unbiased.len(),
        genotype_changes as f64 * 100.0 / results_unbiased.len() as f64
    );
    
    // Check if any biased genotypes include GRCh38
    let pulled_to_ref = results_biased.iter()
        .filter(|r| r.called_genotype.contains("grch38"))
        .count();
    
    println!("Genotypes pulled to reference: {}/{}", pulled_to_ref, results_biased.len());
    
    // Individual similarity drops
    println!("\nIndividual similarity drops:");
    for (i, individual) in test_individuals.iter().enumerate() {
        if i < results_unbiased.len() && i < results_biased.len() {
            let drop = results_unbiased[i].similarity - results_biased[i].similarity;
            println!("  {}: {:.4} -> {:.4} (drop: {:.4})",
                individual,
                results_unbiased[i].similarity,
                results_biased[i].similarity,
                drop
            );
        }
    }
    
    // Key insight
    println!("\n=== KEY INSIGHT ===");
    println!("In hold-2-out, reference bias has TWO effects:");
    println!("1. Coverage patterns are corrupted (pulled toward reference)");
    println!("2. True genotypes are not available to match against");
    println!("Result: We match to WRONG alternatives with lower similarity");
    
    // Verify expectations
    assert!(
        avg_biased_sim < avg_unbiased_sim,
        "Biased hold-2-out should have lower similarity"
    );
    
    assert!(
        genotype_changes > 0,
        "Reference bias should change some genotypes in hold-2-out"
    );
}

/// Test hold-2-out with varying bias strengths
#[test]
fn test_hold2out_bias_gradient() {
    let test_data_path = PathBuf::from("tests/data/hla-f.k51.paths.coverage.tsv.gz");
    if !test_data_path.exists() {
        println!("Test data not found, skipping");
        return;
    }
    
    println!("\n=== HOLD-2-OUT WITH BIAS GRADIENT ===\n");
    
    let ref_data = io::read_gzip_tsv(test_data_path.to_str().unwrap()).unwrap();
    
    let grch38_idx = ref_data.ids.iter()
        .position(|id| id.starts_with("grch38"))
        .expect("Need GRCh38");
    
    let test_individual = "HG00096";
    
    // Get individual's indices
    let idx1 = ref_data.ids.iter()
        .position(|id| id.contains(&format!("{}#1", test_individual)))
        .unwrap();
    let idx2 = ref_data.ids.iter()
        .position(|id| id.contains(&format!("{}#2", test_individual)))
        .unwrap();
    
    let true_coverage = math::sum_vectors(&[
        &ref_data.coverages[idx1],
        &ref_data.coverages[idx2],
    ]);
    
    // Create held-out panel
    let mut held_out_data = io::CoverageData::new();
    for (i, id) in ref_data.ids.iter().enumerate() {
        if i != idx1 && i != idx2 {
            held_out_data.ids.push(id.clone());
            held_out_data.coverages.push(ref_data.coverages[i].clone());
        }
    }
    
    // Test gradient of bias strengths
    let bias_strengths = vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9];
    
    println!("Testing {} with hold-2-out and varying bias:", test_individual);
    println!("\nBias | Similarity | Called Genotype");
    println!("-----|------------|----------------");
    
    let mut last_genotype = String::new();
    
    for bias in &bias_strengths {
        let biased_coverage = if *bias == 0.0 {
            true_coverage.clone()
        } else {
            create_biased_coverage(
                &true_coverage,
                &ref_data.coverages[grch38_idx],
                *bias
            )
        };
        
        let result = find_best_genotype_hold2(
            &held_out_data,
            &biased_coverage,
            test_individual
        );
        
        let changed = if result.called_genotype != last_genotype && !last_genotype.is_empty() {
            " *CHANGED*"
        } else {
            ""
        };
        
        println!("{:.1}  | {:.4}     | {}{}",
            bias,
            result.similarity,
            &result.called_genotype[..50.min(result.called_genotype.len())],
            changed
        );
        
        last_genotype = result.called_genotype;
    }
    
    println!("\nObservations:");
    println!("1. Similarity decreases as bias increases");
    println!("2. Genotype calls change at certain bias thresholds");
    println!("3. Higher bias may pull toward reference-like alternatives");
}

// Helper functions
fn create_biased_coverage(
    true_coverage: &[f64],
    reference_coverage: &[f64],
    bias_strength: f64,
) -> Vec<f64> {
    let mut biased = Vec::with_capacity(true_coverage.len());
    for i in 0..true_coverage.len() {
        let biased_value = true_coverage[i] * (1.0 - bias_strength) + 
                          reference_coverage[i] * bias_strength;
        biased.push(biased_value);
    }
    biased
}

fn find_best_genotype_hold2(
    ref_data: &io::CoverageData,
    sample_coverage: &[f64],
    _individual: &str,
) -> Hold2Result {
    let mut results = Vec::new();
    
    for i in 0..ref_data.ids.len() {
        for j in i..ref_data.ids.len() {
            let combined = math::sum_vectors(&[
                &ref_data.coverages[i],
                &ref_data.coverages[j],
            ]);
            
            let similarity = math::cosine_similarity(&combined, sample_coverage);
            results.push(((i, j), similarity));
        }
    }
    
    results.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    
    if let Some(best) = results.first() {
        let called_genotype = format!("{}+{}", 
            ref_data.ids[best.0.0].split(':').next().unwrap_or("?"),
            ref_data.ids[best.0.1].split(':').next().unwrap_or("?")
        );
        
        Hold2Result {
            similarity: best.1,
            called_genotype,
        }
    } else {
        Hold2Result {
            similarity: 0.0,
            called_genotype: "FAILED".to_string(),
        }
    }
}

#[derive(Clone, Debug)]
struct Hold2Result {
    similarity: f64,
    called_genotype: String,
}
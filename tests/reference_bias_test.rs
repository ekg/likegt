use likegt::{io, math, validation, hold2out};
use std::path::PathBuf;

/// This test demonstrates how reference bias corrupts genotyping accuracy
/// 
/// The key insight: when reads are first aligned to a single reference (e.g., GRCh38),
/// they lose information about variations that differ from that reference.
/// This creates systematic bias toward reference-like genotypes.
#[test]
fn test_reference_bias_corrupts_accuracy() {
    // Skip if test data doesn't exist
    let test_data_path = PathBuf::from("tests/data/hla-f.k51.paths.coverage.tsv.gz");
    if !test_data_path.exists() {
        println!("Test data not found, skipping reference bias test");
        return;
    }
    
    println!("\n=== REFERENCE BIAS CORRUPTION TEST ===\n");
    
    // Load the reference panel
    let ref_data = io::read_gzip_tsv(test_data_path.to_str().unwrap()).unwrap();
    println!("Loaded {} haplotypes", ref_data.len());
    
    // Find GRCh38 reference in the panel
    let grch38_idx = ref_data.ids.iter()
        .position(|id| id.starts_with("grch38"))
        .expect("GRCh38 should be in reference panel");
    
    println!("Using {} as reference for bias", ref_data.ids[grch38_idx]);
    
    // Simulate reference bias by creating corrupted coverage
    // that is pulled toward the reference
    let mut biased_results = Vec::new();
    let mut unbiased_results = Vec::new();
    
    // Test a subset of individuals
    let test_individuals = ["HG00096", "HG00268", "HG00733", "NA12329", "HG02818"];
    
    for individual in &test_individuals {
        // Find this individual's haplotypes
        let hap1_pattern = format!("{}#1", individual);
        let hap2_pattern = format!("{}#2", individual);
        
        let hap1_idx = ref_data.ids.iter()
            .position(|id| id.contains(&hap1_pattern));
        let hap2_idx = ref_data.ids.iter()
            .position(|id| id.contains(&hap2_pattern));
        
        if let (Some(idx1), Some(idx2)) = (hap1_idx, hap2_idx) {
            // Create true coverage (unbiased)
            let true_coverage = math::sum_vectors(&[
                &ref_data.coverages[idx1],
                &ref_data.coverages[idx2],
            ]);
            
            // Create reference-biased coverage
            // This simulates what happens when reads are forced through GRCh38 first:
            // The coverage becomes a weighted average pulled toward the reference
            let bias_strength = 0.3; // 30% pull toward reference
            let biased_coverage = create_biased_coverage(
                &true_coverage,
                &ref_data.coverages[grch38_idx],
                bias_strength
            );
            
            // Test unbiased genotyping
            let unbiased_result = genotype_sample(&ref_data, &true_coverage, individual);
            unbiased_results.push(unbiased_result.clone());
            
            // Test biased genotyping
            let biased_result = genotype_sample(&ref_data, &biased_coverage, individual);
            biased_results.push(biased_result.clone());
            
            println!("\n{} results:", individual);
            println!("  Unbiased: {} (similarity: {:.4})", 
                if unbiased_result.correct { "CORRECT" } else { "WRONG" },
                unbiased_result.similarity
            );
            println!("  Biased:   {} (similarity: {:.4}, called: {})",
                if biased_result.correct { "CORRECT" } else { "WRONG" },
                biased_result.similarity,
                biased_result.called_genotype
            );
            
            // Check if bias pulled toward reference
            if biased_result.called_genotype.contains("grch38") && !unbiased_result.called_genotype.contains("grch38") {
                println!("  -> Reference bias pulled genotype toward GRCh38!");
            }
        }
    }
    
    // Calculate accuracy
    let unbiased_accuracy = unbiased_results.iter()
        .filter(|r| r.correct)
        .count() as f64 / unbiased_results.len() as f64;
    
    let biased_accuracy = biased_results.iter()
        .filter(|r| r.correct)
        .count() as f64 / biased_results.len() as f64;
    
    println!("\n=== IMPACT OF REFERENCE BIAS ===");
    println!("Unbiased accuracy: {:.1}%", unbiased_accuracy * 100.0);
    println!("Biased accuracy:   {:.1}%", biased_accuracy * 100.0);
    println!("Accuracy drop:     {:.1}%", (unbiased_accuracy - biased_accuracy) * 100.0);
    
    // Count how often bias pulls toward reference
    let pulled_to_ref = biased_results.iter()
        .filter(|r| r.called_genotype.contains("grch38"))
        .count();
    println!("Genotypes pulled to GRCh38: {}/{}", pulled_to_ref, biased_results.len());
    
    // We expect bias to reduce accuracy
    assert!(
        biased_accuracy <= unbiased_accuracy,
        "Reference bias should not improve accuracy"
    );
}

/// Create coverage that is biased toward a reference
fn create_biased_coverage(
    true_coverage: &[f64],
    reference_coverage: &[f64],
    bias_strength: f64,
) -> Vec<f64> {
    let mut biased = Vec::with_capacity(true_coverage.len());
    
    for i in 0..true_coverage.len() {
        // Weighted average: pull coverage toward reference
        let biased_value = true_coverage[i] * (1.0 - bias_strength) + 
                          reference_coverage[i] * bias_strength;
        biased.push(biased_value);
    }
    
    biased
}

/// Simple genotyping function for testing
fn genotype_sample(
    ref_data: &io::CoverageData,
    sample_coverage: &[f64],
    true_individual: &str,
) -> GenotypeResult {
    let mut results = Vec::new();
    
    // Find all pairwise combinations
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
        let called_hap1 = &ref_data.ids[best.0.0];
        let called_hap2 = &ref_data.ids[best.0.1];
        
        let correct = called_hap1.contains(true_individual) || called_hap2.contains(true_individual);
        
        let called_genotype = format!("{}+{}", 
            called_hap1.split(':').next().unwrap_or("?"),
            called_hap2.split(':').next().unwrap_or("?")
        );
        
        GenotypeResult {
            correct,
            similarity: best.1,
            called_genotype,
        }
    } else {
        GenotypeResult {
            correct: false,
            similarity: 0.0,
            called_genotype: "FAILED".to_string(),
        }
    }
}

#[derive(Clone, Debug)]
struct GenotypeResult {
    correct: bool,
    similarity: f64,
    called_genotype: String,
}

#[test]
fn test_reference_bias_hold2out() {
    let test_data_path = PathBuf::from("tests/data/hla-f.k51.paths.coverage.tsv.gz");
    if !test_data_path.exists() {
        println!("Test data not found, skipping");
        return;
    }
    
    println!("\n=== REFERENCE BIAS WITH HOLD-2-OUT ===\n");
    
    let ref_data = io::read_gzip_tsv(test_data_path.to_str().unwrap()).unwrap();
    
    // Find GRCh38 reference
    let grch38_idx = ref_data.ids.iter()
        .position(|id| id.starts_with("grch38"));
    
    if grch38_idx.is_none() {
        println!("GRCh38 not found in reference panel");
        return;
    }
    
    let grch38_idx = grch38_idx.unwrap();
    let bias_strength = 0.3;
    
    // Test individual that we'll hold out
    let test_individual = "HG00096";
    
    // Regular hold-2-out (unbiased)
    let unbiased_result = hold2out::hold2out_individual(&ref_data, test_individual).unwrap();
    
    // Now simulate biased hold-2-out
    // First, get the individual's true coverage
    let hap1_idx = ref_data.ids.iter()
        .position(|id| id.contains(&format!("{}#1", test_individual)))
        .unwrap();
    let hap2_idx = ref_data.ids.iter()
        .position(|id| id.contains(&format!("{}#2", test_individual)))
        .unwrap();
    
    let true_coverage = math::sum_vectors(&[
        &ref_data.coverages[hap1_idx],
        &ref_data.coverages[hap2_idx],
    ]);
    
    // Apply reference bias
    let biased_coverage = create_biased_coverage(
        &true_coverage,
        &ref_data.coverages[grch38_idx],
        bias_strength
    );
    
    // Create held-out reference panel (without the test individual)
    let mut held_out_data = io::CoverageData::new();
    for (i, id) in ref_data.ids.iter().enumerate() {
        if i != hap1_idx && i != hap2_idx {
            held_out_data.ids.push(id.clone());
            held_out_data.coverages.push(ref_data.coverages[i].clone());
        }
    }
    
    // Find best genotype with biased coverage
    let biased_result = genotype_sample(&held_out_data, &biased_coverage, test_individual);
    
    println!("Hold-2-out for {}:", test_individual);
    println!("  Unbiased similarity: {:.4}", unbiased_result.similarity);
    println!("  Biased similarity:   {:.4}", biased_result.similarity);
    println!("  Biased called:       {}", biased_result.called_genotype);
    
    // With bias, similarity should generally be lower
    // because we're matching against corrupted coverage
    assert!(
        biased_result.similarity <= unbiased_result.similarity + 0.01, // Small tolerance
        "Biased coverage should not improve similarity"
    );
}

#[test]
fn test_varying_bias_strength() {
    let test_data_path = PathBuf::from("tests/data/hla-f.k51.paths.coverage.tsv.gz");
    if !test_data_path.exists() {
        println!("Test data not found, skipping");
        return;
    }
    
    println!("\n=== EFFECT OF VARYING BIAS STRENGTH ===\n");
    
    let ref_data = io::read_gzip_tsv(test_data_path.to_str().unwrap()).unwrap();
    
    // Find GRCh38 and a test individual
    let grch38_idx = ref_data.ids.iter()
        .position(|id| id.starts_with("grch38"))
        .expect("Need GRCh38 in panel");
    
    let test_individual = "HG00096";
    let hap1_idx = ref_data.ids.iter()
        .position(|id| id.contains(&format!("{}#1", test_individual)))
        .unwrap();
    let hap2_idx = ref_data.ids.iter()
        .position(|id| id.contains(&format!("{}#2", test_individual)))
        .unwrap();
    
    let true_coverage = math::sum_vectors(&[
        &ref_data.coverages[hap1_idx],
        &ref_data.coverages[hap2_idx],
    ]);
    
    // Test different bias strengths
    let bias_strengths = vec![0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9];
    
    println!("Testing {} with varying reference bias:", test_individual);
    println!("Bias | Correct | Similarity | Called Genotype");
    println!("-----|---------|------------|----------------");
    
    for bias in &bias_strengths {
        let biased_coverage = create_biased_coverage(
            &true_coverage,
            &ref_data.coverages[grch38_idx],
            *bias
        );
        
        let result = genotype_sample(&ref_data, &biased_coverage, test_individual);
        
        println!("{:.1}  | {:7} | {:.4}     | {}",
            bias,
            if result.correct { "YES" } else { "NO" },
            result.similarity,
            result.called_genotype
        );
    }
    
    println!("\nAs bias increases:");
    println!("1. Accuracy typically decreases");
    println!("2. Genotypes shift toward the reference");
    println!("3. Similarity scores may decrease due to distorted coverage");
}
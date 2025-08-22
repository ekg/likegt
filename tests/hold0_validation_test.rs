use likegt::{io, math, coverage};

/// This test validates that hold-0-out genotyping works correctly:
/// When we create a sample by summing coverage from two haplotypes in the graph,
/// we should always recover those exact haplotypes as the best genotype
#[test]
fn test_hold0_always_correct() {
    // Create a realistic pangenome graph scenario
    let mut ref_data = io::CoverageData::new();
    
    // Simulate 20 haplotypes with realistic coverage patterns
    // Each haplotype covers different nodes in the graph
    let num_haplotypes = 20;
    let num_nodes = 200;
    
    for i in 0..num_haplotypes {
        ref_data.ids.push(format!("haplotype_{:03}", i));
        
        let mut coverage = vec![0.0; num_nodes];
        
        // Each haplotype has a unique pattern of covered nodes
        // Simulate that each haplotype covers ~30% of nodes
        for j in 0..num_nodes {
            // Use a deterministic pattern based on haplotype ID
            let hash = (i * 37 + j * 13) % 100;
            if hash < 30 {
                // This node is covered by this haplotype
                coverage[j] = 1.0;
            }
        }
        
        ref_data.coverages.push(coverage);
    }
    
    // Test multiple scenarios
    let test_pairs = vec![
        (0, 1),   // Different haplotypes
        (5, 10),  // Different haplotypes
        (3, 3),   // Homozygous
        (15, 19), // Different haplotypes
        (7, 7),   // Homozygous
    ];
    
    for (true_hap1, true_hap2) in test_pairs {
        println!("Testing pair: haplotype_{:03} + haplotype_{:03}", true_hap1, true_hap2);
        
        // Create perfect sample as sum of two haplotypes
        let sample_coverage = math::sum_vectors(&[
            &ref_data.coverages[true_hap1],
            &ref_data.coverages[true_hap2],
        ]);
        
        // Find the best genotype by testing all combinations
        let mut results = Vec::new();
        
        for i in 0..num_haplotypes {
            for j in i..num_haplotypes {
                let combined = math::sum_vectors(&[
                    &ref_data.coverages[i],
                    &ref_data.coverages[j],
                ]);
                
                let similarity = math::cosine_similarity(&combined, &sample_coverage);
                
                results.push(((i, j), similarity));
            }
        }
        
        // Sort by similarity (descending)
        results.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        
        // The best result should be our true genotype
        let best = results[0];
        
        assert_eq!(
            best.0, 
            (true_hap1.min(true_hap2), true_hap1.max(true_hap2)),
            "Failed to recover correct genotype for pair ({}, {}). Got ({}, {}) instead with similarity {}",
            true_hap1, true_hap2, best.0.0, best.0.1, best.1
        );
        
        // With perfect data, we should get perfect similarity (1.0)
        assert!(
            (best.1 - 1.0).abs() < 1e-10,
            "Expected perfect similarity (1.0) for true genotype, got {}",
            best.1
        );
        
        // Verify that no other genotype has perfect similarity
        if results.len() > 1 {
            assert!(
                results[1].1 < 1.0 - 1e-10,
                "Multiple genotypes with perfect similarity detected"
            );
        }
    }
}

/// Test that adding noise to coverage still recovers correct genotype
/// as long as the noise is not too large
#[test]
fn test_hold0_with_small_noise() {
    let mut ref_data = io::CoverageData::new();
    
    // Create distinct haplotypes with more unique patterns
    for i in 0..10 {
        ref_data.ids.push(format!("hap{}", i));
        
        let mut coverage = vec![0.0; 50];
        // Each haplotype has a more distinct coverage pattern
        for j in 0..50 {
            // Use different patterns for each haplotype
            let pattern = match i {
                0 => j % 5 == 0,
                1 => j % 7 == 0,
                2 => j % 3 == 0,
                3 => j % 11 == 0,
                4 => j % 13 == 0,
                5 => (j + 1) % 5 == 0,
                6 => (j + 2) % 7 == 0,
                7 => j % 4 == 0,
                8 => (j + 3) % 9 == 0,
                9 => j % 6 == 0,
                _ => false,
            };
            
            if pattern {
                coverage[j] = 10.0;
            }
        }
        ref_data.coverages.push(coverage);
    }
    
    // True genotype
    let true_hap1 = 2;
    let true_hap2 = 7;
    
    // Create sample with small noise
    let mut sample_coverage = math::sum_vectors(&[
        &ref_data.coverages[true_hap1],
        &ref_data.coverages[true_hap2],
    ]);
    
    // Add small random noise (up to 5% of signal)
    for i in 0..sample_coverage.len() {
        let noise = ((i * 17) % 10) as f64 / 100.0 - 0.05; // -0.05 to +0.05
        sample_coverage[i] = (sample_coverage[i] + sample_coverage[i] * noise).max(0.0);
    }
    
    // Find best genotype
    let mut best_similarity = 0.0;
    let mut best_pair = (0, 0);
    
    for i in 0..10 {
        for j in i..10 {
            let combined = math::sum_vectors(&[
                &ref_data.coverages[i],
                &ref_data.coverages[j],
            ]);
            let similarity = math::cosine_similarity(&combined, &sample_coverage);
            
            if similarity > best_similarity {
                best_similarity = similarity;
                best_pair = (i, j);
            }
        }
    }
    
    // Should still recover the correct genotype
    assert_eq!(
        best_pair, 
        (true_hap1, true_hap2),
        "Failed to recover correct genotype with small noise"
    );
    
    // Similarity should be very high but not perfect due to noise
    assert!(
        best_similarity > 0.99,
        "Similarity too low with small noise: {}",
        best_similarity
    );
    assert!(
        best_similarity < 1.0,
        "Similarity should not be perfect with noise: {}",
        best_similarity
    );
}
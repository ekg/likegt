use likegt::{io, math, coverage};
use std::path::PathBuf;

#[test]
fn test_integration_small_example() {
    // Create test reference haplotypes
    let mut ref_data = io::CoverageData::new();
    
    // Add 5 haplotypes with distinct patterns
    for i in 0..5 {
        ref_data.ids.push(format!("hap{}", i));
        let mut coverage = vec![0.0; 10];
        // Create unique pattern for each haplotype
        for j in 0..10 {
            if (i + j) % 3 == 0 {
                coverage[j] = 10.0 + i as f64;
            }
        }
        ref_data.coverages.push(coverage);
    }
    
    // Create a sample that's the sum of hap1 and hap3
    let mut sample_data = io::CoverageData::new();
    sample_data.ids.push("test_sample".to_string());
    
    let sample_cov = math::sum_vectors(&[
        &ref_data.coverages[1],
        &ref_data.coverages[3],
    ]);
    sample_data.coverages.push(sample_cov.clone());
    
    // Create genotype data
    let genotype_data = coverage::GenotypeData::new(
        ref_data.clone(),
        sample_data,
        "test_sample".to_string(),
    ).unwrap();
    
    // Calculate similarities for all combinations
    let mut results = Vec::new();
    for i in 0..5 {
        for j in i..5 {
            let combined = math::sum_vectors(&[
                &ref_data.coverages[i],
                &ref_data.coverages[j],
            ]);
            let similarity = math::cosine_similarity(&combined, &sample_cov);
            results.push(((i, j), similarity));
        }
    }
    
    // Sort by similarity
    results.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    
    // The best should be hap1 + hap3
    assert_eq!(results[0].0, (1, 3));
    assert!((results[0].1 - 1.0).abs() < 1e-10);
}

#[test]
#[cfg(unix)]
fn test_read_test_data_if_exists() {
    // This test only runs if test data exists
    let test_data_path = PathBuf::from("tests/data/hla-f.k51.paths.coverage.tsv.gz");
    
    if test_data_path.exists() {
        let data = io::read_gzip_tsv(test_data_path.to_str().unwrap()).unwrap();
        
        // Should have read the HLA data
        assert!(data.len() > 0);
        assert!(data.ids.len() == data.coverages.len());
        
        // Each coverage vector should have the same length
        if !data.coverages.is_empty() {
            let expected_len = data.coverages[0].len();
            for cov in &data.coverages {
                assert_eq!(cov.len(), expected_len);
            }
        }
    }
}

#[test]
fn test_perfect_genotyping_scenario() {
    // Test that when we have perfect coverage data,
    // we always recover the correct genotype
    
    let mut ref_data = io::CoverageData::new();
    
    // Create 10 random haplotypes
    for i in 0..10 {
        ref_data.ids.push(format!("haplotype_{:03}", i));
        
        // Create random but reproducible coverage pattern
        let mut coverage = Vec::new();
        for j in 0..100 {
            // Use a simple hash-like function for reproducibility
            let value = ((i * 31 + j * 17) % 100) as f64 / 10.0;
            coverage.push(value);
        }
        ref_data.coverages.push(coverage);
    }
    
    // Test multiple true genotypes
    let test_cases = vec![
        (0, 1),
        (2, 5),
        (3, 3), // Homozygous
        (7, 9),
    ];
    
    for (hap1_idx, hap2_idx) in test_cases {
        // Create perfect sample
        let sample_cov = math::sum_vectors(&[
            &ref_data.coverages[hap1_idx],
            &ref_data.coverages[hap2_idx],
        ]);
        
        // Find best genotype
        let mut best_similarity = 0.0;
        let mut best_pair = (0, 0);
        
        for i in 0..10 {
            for j in i..10 {
                let combined = math::sum_vectors(&[
                    &ref_data.coverages[i],
                    &ref_data.coverages[j],
                ]);
                let similarity = math::cosine_similarity(&combined, &sample_cov);
                
                if similarity > best_similarity {
                    best_similarity = similarity;
                    best_pair = (i, j);
                }
            }
        }
        
        // Should recover the true genotype with perfect similarity
        assert_eq!(best_pair, (hap1_idx, hap2_idx));
        assert!(
            (best_similarity - 1.0).abs() < 1e-10,
            "Expected perfect similarity for pair ({}, {}), got {}",
            hap1_idx, hap2_idx, best_similarity
        );
    }
}
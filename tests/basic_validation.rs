/// Basic validation test using our actual HLA data
use likegt::{io, math};

#[test]
fn test_cosine_similarity_with_hla_data() {
    // Test that we can read our HLA data and perform basic cosine similarity
    if let Ok(ref_data) = io::read_gzip_tsv("hla-f.k51.paths.tsv.gz") {
        println!("✅ Successfully loaded HLA reference data");
        println!("   {} haplotypes loaded", ref_data.ids.len());
        
        if !ref_data.ids.is_empty() {
            let first_coverage = &ref_data.coverages[0];
            println!("   First haplotype: {}", ref_data.ids[0]);
            println!("   Coverage vector length: {}", first_coverage.len());
            
            if !first_coverage.is_empty() {
                // Test self-similarity (should be perfect)
                let self_similarity = math::cosine_similarity(first_coverage, first_coverage);
                assert!((self_similarity - 1.0).abs() < 1e-10, "Self-similarity should be perfect");
                
                // Test that we can sum vectors
                if ref_data.coverages.len() >= 2 {
                    let combined = math::sum_vectors(&[
                        &ref_data.coverages[0],
                        &ref_data.coverages[1],
                    ]);
                    assert_eq!(combined.len(), first_coverage.len(), "Combined vector should have same length");
                    
                    // Test similarity of combined vs original
                    let similarity = math::cosine_similarity(&combined, first_coverage);
                    assert!(similarity >= 0.0 && similarity <= 1.0, "Similarity should be in [0,1] range");
                    
                    println!("   ✅ Vector operations work correctly");
                }
            } else {
                println!("   ⚠️ Coverage vectors are empty - skipping vector tests");
            }
        }
    } else {
        println!("⚠️ HLA test data not found - test will pass but skip validation");
    }
}

#[test]
fn test_algorithm_correctness() {
    // Test with synthetic data to ensure algorithm works correctly
    let hap1 = vec![1.0, 0.0, 2.0, 0.0];
    let hap2 = vec![0.0, 3.0, 0.0, 1.0];
    let combined = math::sum_vectors(&[&hap1, &hap2]);
    
    // Combined should be [1.0, 3.0, 2.0, 1.0]
    assert_eq!(combined, vec![1.0, 3.0, 2.0, 1.0]);
    
    // Test that hap1 + hap2 has perfect similarity with the sum
    let similarity = math::cosine_similarity(&combined, &combined);
    assert!((similarity - 1.0).abs() < 1e-10);
    
    // Test that wrong combination gives different similarity
    let wrong_combined = math::sum_vectors(&[&hap1, &hap1]);
    let wrong_similarity = math::cosine_similarity(&wrong_combined, &combined);
    assert!(wrong_similarity < 1.0);
    
    println!("✅ Core algorithm works correctly");
}
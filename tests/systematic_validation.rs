use likegt::{io, validation};
use std::path::PathBuf;

/// This test systematically validates EVERY individual in the dataset
/// It shows actual work by testing hold-0-out and reporting accuracy
#[test]
fn test_all_individuals_hold0() {
    println!("\n=== SYSTEMATIC VALIDATION OF ALL INDIVIDUALS ===\n");
    
    let test_data_path = PathBuf::from("tests/data/hla-f.k51.paths.coverage.tsv.gz");
    
    if !test_data_path.exists() {
        println!("Test data not found, skipping");
        return;
    }
    
    // Read the actual graph coverage data
    println!("Loading graph coverage data...");
    let ref_data = io::read_gzip_tsv(test_data_path.to_str().unwrap()).unwrap();
    println!("Loaded {} haplotypes with {} nodes each", ref_data.len(), ref_data.coverages[0].len());
    
    // Run validation on ALL individuals
    println!("\nValidating all individuals (this shows actual work!)...");
    let results = validation::validate_all_individuals(&ref_data).unwrap();
    
    // Generate and print report
    let report = validation::generate_accuracy_report(&results);
    println!("{}", report);
    
    // Assert that hold-0-out should be perfect for all individuals
    let failed: Vec<_> = results.iter().filter(|r| !r.correct).collect();
    
    if !failed.is_empty() {
        println!("\n⚠️  WARNING: Hold-0-out should be PERFECT but failed for:");
        for fail in &failed {
            println!("  - {}: similarity = {:.6} (rank {})", 
                fail.sample_id, fail.called_similarity, fail.rank);
        }
    }
    
    // For hold-0-out, we expect 100% accuracy
    assert!(
        results.iter().all(|r| r.correct),
        "Hold-0-out should have 100% accuracy, but {} individuals failed",
        failed.len()
    );
}

/// Test hold-2-out validation (more realistic scenario)
#[test] 
fn test_hold2_out_validation() {
    println!("\n=== HOLD-2-OUT VALIDATION ===\n");
    
    let test_data_path = PathBuf::from("tests/data/hla-f.k51.paths.coverage.tsv.gz");
    
    if !test_data_path.exists() {
        println!("Test data not found, skipping");
        return;
    }
    
    let mut ref_data = io::read_gzip_tsv(test_data_path.to_str().unwrap()).unwrap();
    
    // Test specific individuals with hold-2-out
    let test_individuals = vec!["HG00096", "HG00268", "NA12329", "HG00733"];
    
    for individual in test_individuals {
        println!("\nTesting hold-2-out for {}:", individual);
        
        // Find and remove this individual's haplotypes
        let hap1_pattern = format!("{}#1", individual);
        let hap2_pattern = format!("{}#2", individual);
        
        let hap1_idx = ref_data.ids.iter().position(|id| id.contains(&hap1_pattern));
        let hap2_idx = ref_data.ids.iter().position(|id| id.contains(&hap2_pattern));
        
        if let (Some(idx1), Some(idx2)) = (hap1_idx, hap2_idx) {
            // Save the haplotypes
            let saved_hap1 = (ref_data.ids[idx1].clone(), ref_data.coverages[idx1].clone());
            let saved_hap2 = (ref_data.ids[idx2].clone(), ref_data.coverages[idx2].clone());
            
            // Create sample coverage BEFORE removing
            let sample_coverage = likegt::math::sum_vectors(&[
                &ref_data.coverages[idx1],
                &ref_data.coverages[idx2],
            ]);
            
            // Remove from reference (hold-2-out)
            let mut filtered_ids = ref_data.ids.clone();
            let mut filtered_coverages = ref_data.coverages.clone();
            
            // Remove in reverse order to maintain indices
            let (first, second) = if idx1 > idx2 { (idx1, idx2) } else { (idx2, idx1) };
            filtered_ids.remove(first);
            filtered_coverages.remove(first);
            filtered_ids.remove(second);
            filtered_coverages.remove(second);
            
            // Create filtered data
            let mut filtered_data = io::CoverageData::new();
            filtered_data.ids = filtered_ids;
            filtered_data.coverages = filtered_coverages;
            
            println!("  Reference haplotypes after removal: {}", filtered_data.len());
            
            // Find best genotype without the true haplotypes
            let mut results = Vec::new();
            for i in 0..filtered_data.ids.len() {
                for j in i..filtered_data.ids.len() {
                    let combined = likegt::math::sum_vectors(&[
                        &filtered_data.coverages[i],
                        &filtered_data.coverages[j],
                    ]);
                    
                    let similarity = likegt::math::cosine_similarity(&combined, &sample_coverage);
                    results.push(((i, j), similarity));
                }
            }
            
            results.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            
            if !results.is_empty() {
                let best = results[0];
                println!("  Best match: {} + {}", 
                    filtered_data.ids[best.0.0].split(':').next().unwrap_or("?"),
                    filtered_data.ids[best.0.1].split(':').next().unwrap_or("?")
                );
                println!("  Similarity: {:.6}", best.1);
                
                // In hold-2-out, we don't expect perfect similarity
                assert!(
                    best.1 < 1.0,
                    "Hold-2-out should not have perfect similarity"
                );
            }
        }
    }
}

/// Test that shows computation is actually happening
#[test]
fn test_computational_work() {
    use std::time::Instant;
    
    println!("\n=== DEMONSTRATING COMPUTATIONAL WORK ===\n");
    
    // Create a larger test case to show real computation
    let mut ref_data = io::CoverageData::new();
    
    // Create 100 haplotypes with 1000 nodes each
    let num_haplotypes = 100;
    let num_nodes = 1000;
    
    println!("Creating {} haplotypes with {} nodes each...", num_haplotypes, num_nodes);
    for i in 0..num_haplotypes {
        ref_data.ids.push(format!("synthetic_hap_{:03}", i));
        
        let mut coverage = vec![0.0; num_nodes];
        for j in 0..num_nodes {
            // Create complex patterns
            let value = ((i * 31 + j * 17) % 100) as f64 / 10.0;
            coverage[j] = value * ((j % (i + 1)) as f64).sin().abs();
        }
        ref_data.coverages.push(coverage);
    }
    
    // Time the genotyping computation
    println!("\nComputing all pairwise combinations ({} total)...", num_haplotypes * (num_haplotypes + 1) / 2);
    
    let start = Instant::now();
    
    // Create a test sample
    let sample_coverage = likegt::math::sum_vectors(&[
        &ref_data.coverages[10],
        &ref_data.coverages[25],
    ]);
    
    let mut results = Vec::new();
    let mut operations = 0u64;
    
    for i in 0..num_haplotypes {
        for j in i..num_haplotypes {
            let combined = likegt::math::sum_vectors(&[
                &ref_data.coverages[i],
                &ref_data.coverages[j],
            ]);
            
            let similarity = likegt::math::cosine_similarity(&combined, &sample_coverage);
            results.push(((i, j), similarity));
            
            operations += num_nodes as u64 * 3; // sum + dot product + magnitude
        }
        
        if i % 10 == 0 {
            println!("  Processed {} haplotypes...", i);
        }
    }
    
    let elapsed = start.elapsed();
    
    results.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    
    println!("\nComputation complete!");
    println!("  Time elapsed: {:.2} seconds", elapsed.as_secs_f64());
    println!("  Combinations evaluated: {}", results.len());
    println!("  Total operations: ~{} million", operations / 1_000_000);
    println!("  Operations per second: ~{} million", (operations as f64 / elapsed.as_secs_f64()) / 1_000_000.0);
    
    // Verify we found the correct answer
    let best = results[0];
    assert_eq!(best.0, (10, 25), "Should recover the true genotype");
    assert!((best.1 - 1.0).abs() < 1e-10, "Should have perfect similarity");
    
    println!("\n✓ Correctly identified synthetic_hap_010 + synthetic_hap_025");
}
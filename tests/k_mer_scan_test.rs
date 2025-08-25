use likegt::{io, validation};
use std::path::PathBuf;

#[test]
fn test_k_mer_size_accuracy_scan() {
    println!("\n=== K-MER SIZE ACCURACY SCAN ===\n");
    
    let k_values = vec![
        (51, "tests/data/hla-f.k51.paths.coverage.tsv.gz"),
        (101, "tests/data/k101/hla-f.k101.paths.coverage.tsv.gz"), 
        (179, "tests/data/k179/hla-f.k179.paths.coverage.tsv.gz"),
    ];
    
    let mut results = Vec::new();
    
    for (k, path) in k_values {
        let test_data_path = PathBuf::from(path);
        
        if !test_data_path.exists() {
            println!("⚠️  K={}: Test data not found at {}", k, path);
            continue;
        }
        
        println!("🔍 Running hold-0-out validation for K={}...", k);
        
        match io::read_gzip_tsv(test_data_path.to_str().unwrap()) {
            Ok(ref_data) => {
                println!("   Loaded {} haplotypes with {} nodes", 
                    ref_data.ids.len(), ref_data.coverages[0].len());
                
                match validation::validate_all_individuals(&ref_data) {
                    Ok(validation_results) => {
                        let total = validation_results.len();
                        let correct = validation_results.iter().filter(|r| r.correct).count();
                        let accuracy = 100.0 * correct as f64 / total as f64;
                        let avg_qv = validation_results.iter().map(|r| r.qv).sum::<f64>() / total as f64;
                        let avg_rank = validation_results.iter().map(|r| r.rank).sum::<usize>() as f64 / total as f64;
                        
                        println!("   ✅ K={}: {}/{} correct ({:.1}%), Avg QV: {:.1}, Avg rank: {:.1}", 
                            k, correct, total, accuracy, avg_qv, avg_rank);
                        
                        // Find failed cases
                        let failures: Vec<_> = validation_results.iter().filter(|r| !r.correct).collect();
                        if !failures.is_empty() {
                            println!("   ❌ Failed cases: {}", 
                                failures.iter()
                                    .map(|f| f.sample_id.as_str())
                                    .collect::<Vec<_>>()
                                    .join(", ")
                            );
                        }
                        
                        results.push((k, accuracy, avg_qv, correct, total, failures.len()));
                    }
                    Err(e) => {
                        println!("   ❌ K={}: Validation failed: {}", k, e);
                    }
                }
            }
            Err(e) => {
                println!("   ❌ K={}: Failed to load data: {}", k, e);
            }
        }
        println!();
    }
    
    // Summary comparison
    println!("=== K-MER SIZE COMPARISON SUMMARY ===");
    println!("K-mer\tAccuracy\tAvg QV\tCorrect\tTotal\tFailed");
    println!("-----\t--------\t------\t-------\t-----\t------");
    
    for (k, accuracy, avg_qv, correct, total, failed) in &results {
        println!("{}\t{:.1}%\t\t{:.1}\t{}\t{}\t{}", 
            k, accuracy, avg_qv, correct, total, failed);
    }
    
    // Analysis
    if results.len() > 1 {
        let best_k = results.iter().max_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        if let Some((best_k_val, best_acc, _, _, _, _)) = best_k {
            println!("\n🏆 Best performing K-mer size: {} ({:.1}% accuracy)", best_k_val, best_acc);
        }
        
        println!("\n📊 Trends:");
        for i in 1..results.len() {
            let (k_prev, acc_prev, _, _, _, fail_prev) = results[i-1];
            let (k_curr, acc_curr, _, _, _, fail_curr) = results[i];
            let acc_change = acc_curr - acc_prev;
            let fail_change = fail_curr as i32 - fail_prev as i32;
            
            if acc_change > 0.1 {
                println!("   ⬆️ K {} → {}: Accuracy improved by {:.1}% ({} fewer failures)", 
                    k_prev, k_curr, acc_change, -fail_change);
            } else if acc_change < -0.1 {
                println!("   ⬇️ K {} → {}: Accuracy decreased by {:.1}% ({} more failures)",
                    k_prev, k_curr, -acc_change, fail_change);
            } else {
                println!("   ➡️ K {} → {}: Similar accuracy ({:.1}% vs {:.1}%)",
                    k_prev, k_curr, acc_prev, acc_curr);
            }
        }
    }
}

#[test]
fn test_k311_if_available() {
    println!("\n=== TESTING K=311 (if available) ===\n");
    
    // Look for K=311 data in various possible locations
    let possible_paths = vec![
        "tests/data/k311/hla-f.k311.paths.coverage.tsv.gz",
        "tests/data/hla-f.k311.paths.coverage.tsv.gz",
        "hla-f.k311.paths.coverage.tsv.gz",
    ];
    
    let mut found = false;
    for path in possible_paths {
        let test_data_path = PathBuf::from(path);
        if test_data_path.exists() {
            println!("📁 Found K=311 data at: {}", path);
            
            match io::read_gzip_tsv(path) {
                Ok(ref_data) => {
                    println!("   Loaded {} haplotypes with {} nodes", 
                        ref_data.ids.len(), ref_data.coverages[0].len());
                    
                    match validation::validate_all_individuals(&ref_data) {
                        Ok(validation_results) => {
                            let total = validation_results.len();
                            let correct = validation_results.iter().filter(|r| r.correct).count();
                            let accuracy = 100.0 * correct as f64 / total as f64;
                            let avg_qv = validation_results.iter().map(|r| r.qv).sum::<f64>() / total as f64;
                            
                            println!("   ✅ K=311: {}/{} correct ({:.1}%), Avg QV: {:.1}", 
                                correct, total, accuracy, avg_qv);
                            
                            // Compare with K=51 baseline
                            println!("   📈 Compared to K=51 (89.2% accuracy):");
                            if accuracy > 89.2 {
                                println!("      🎯 IMPROVED by {:.1} percentage points!", accuracy - 89.2);
                            } else if accuracy < 89.2 {
                                println!("      📉 Decreased by {:.1} percentage points", 89.2 - accuracy);
                            } else {
                                println!("      ➡️  Same accuracy as K=51");
                            }
                        }
                        Err(e) => {
                            println!("   ❌ Validation failed: {}", e);
                        }
                    }
                }
                Err(e) => {
                    println!("   ❌ Failed to load data: {}", e);
                }
            }
            found = true;
            break;
        }
    }
    
    if !found {
        println!("ℹ️  K=311 data not found. You can generate it with:");
        println!("   seqwish -s hla-f.fa.gz -k 311 -p hla-f.k311");
        println!("   odgi build -g hla-f.k311.gfa -o hla-f.k311.og");
        println!("   odgi paths -i hla-f.k311.og -H | cut -f 1,4- | gzip > hla-f.k311.paths.coverage.tsv.gz");
    }
}
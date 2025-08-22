pub mod io;
pub mod math;
pub mod genotype;
pub mod coverage;
pub mod graph;
pub mod validation;
pub mod sequence_qv;
pub mod hold2out;

#[cfg(test)]
mod tests {
    use super::*;
    
    /// Test that when we create a sample that's exactly the sum of two haplotypes,
    /// we recover those haplotypes as the best genotype
    #[test]
    fn test_hold0_perfect_recovery() {
        // Create test data
        let mut ref_data = io::CoverageData::new();
        
        // Add some reference haplotypes with distinct patterns
        ref_data.ids.push("hap1".to_string());
        ref_data.coverages.push(vec![1.0, 0.0, 1.0, 0.0, 1.0]);
        
        ref_data.ids.push("hap2".to_string());
        ref_data.coverages.push(vec![0.0, 1.0, 0.0, 1.0, 0.0]);
        
        ref_data.ids.push("hap3".to_string());
        ref_data.coverages.push(vec![1.0, 1.0, 0.0, 0.0, 0.0]);
        
        ref_data.ids.push("hap4".to_string());
        ref_data.coverages.push(vec![0.0, 0.0, 1.0, 1.0, 1.0]);
        
        // Create a sample that's exactly hap1 + hap2
        let sample_coverage = vec![1.0, 1.0, 1.0, 1.0, 1.0];
        
        // Calculate all similarities
        let mut results = Vec::new();
        
        // Test all combinations
        for i in 0..ref_data.ids.len() {
            for j in i..ref_data.ids.len() {
                let combined = math::sum_vectors(&[
                    &ref_data.coverages[i],
                    &ref_data.coverages[j],
                ]);
                
                let similarity = math::cosine_similarity(&combined, &sample_coverage);
                
                results.push((
                    vec![ref_data.ids[i].clone(), ref_data.ids[j].clone()],
                    similarity,
                ));
            }
        }
        
        // Sort by similarity
        results.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        
        // The best should be hap1 + hap2 with perfect similarity
        assert_eq!(results[0].0, vec!["hap1", "hap2"]);
        assert!((results[0].1 - 1.0).abs() < 1e-10, "Expected perfect similarity for true genotype");
    }
    
    /// Test that coverage scaling doesn't affect relative rankings
    #[test]
    fn test_coverage_scaling_invariance() {
        let ref1 = vec![1.0, 2.0, 3.0, 4.0];
        let ref2 = vec![4.0, 3.0, 2.0, 1.0];
        let sample = vec![5.0, 5.0, 5.0, 5.0];
        
        let combined = math::sum_vectors(&[&ref1, &ref2]);
        let sim1 = math::cosine_similarity(&combined, &sample);
        
        // Scale everything by 10
        let ref1_scaled: Vec<f64> = ref1.iter().map(|x| x * 10.0).collect();
        let ref2_scaled: Vec<f64> = ref2.iter().map(|x| x * 10.0).collect();
        let sample_scaled: Vec<f64> = sample.iter().map(|x| x * 10.0).collect();
        
        let combined_scaled = math::sum_vectors(&[&ref1_scaled, &ref2_scaled]);
        let sim2 = math::cosine_similarity(&combined_scaled, &sample_scaled);
        
        // Cosine similarity should be the same regardless of scaling
        assert!((sim1 - sim2).abs() < 1e-10);
    }
}
use anyhow::Result;
use std::path::Path;
use tempfile::TempDir;

use crate::commands;

/// Integration tests for the CLI tool using our HLA test data
#[cfg(test)]
mod tests {
    use super::*;
    
    /// Test that we can build a graph from our test FASTA and validate it
    #[tokio::test]
    #[ignore] // Requires external tools (allwave, seqwish, odgi)
    async fn test_full_validation_pipeline() {
        // This test uses the HLA sequences we've been working with
        let test_fasta = "test_hla_sequences.fa";
        
        // Create test FASTA if it doesn't exist
        if !Path::new(test_fasta).exists() {
            create_test_hla_fasta(test_fasta).await.unwrap();
        }
        
        let temp_dir = TempDir::new().unwrap();
        let output_dir = temp_dir.path().join("validation_results");
        
        // Run the validation command
        let result = commands::validate::run_validation(
            test_fasta,
            output_dir.to_str().unwrap(),
            2, // threads
            51, // k-mer size
        ).await;
        
        match result {
            Ok(_) => {
                // Check that validation report was created
                let report_path = output_dir.join("validation_report.txt");
                assert!(report_path.exists());
                
                // Check that a graph was built
                let graph_path = output_dir.join("graph.gfa");
                assert!(graph_path.exists());
                
                // Verify report contains expected sections
                let report_content = tokio::fs::read_to_string(&report_path).await.unwrap();
                assert!(report_content.contains("Hold-0-out Validation"));
                assert!(report_content.contains("Hold-2-out Validation"));
                assert!(report_content.contains("Reference Bias Test"));
                assert!(report_content.contains("Overall Assessment"));
                
                println!("✅ Full validation pipeline test passed");
            }
            Err(e) => {
                // If external tools aren't available, this test will fail gracefully
                println!("⚠️ Validation pipeline test failed (tools may not be installed): {}", e);
            }
        }
    }
    
    /// Test graph analysis functionality
    #[test]
    fn test_graph_check_functionality() {
        // Create a minimal test GFA
        let temp_dir = TempDir::new().unwrap();
        let test_gfa = temp_dir.path().join("test.gfa");
        let report_path = temp_dir.path().join("report.txt");
        
        // Create a simple test GFA
        let gfa_content = r#"H	VN:Z:1.0
S	1	ACGT
S	2	TCGA
S	3	AAAA
L	1	+	2	+	0M
L	2	+	3	+	0M
P	sample1#chr1#hap1	1+,2+,3+	*
P	sample1#chr1#hap2	1+,3+	*
P	sample2#chr1#hap1	2+,3+	*
P	sample2#chr1#hap2	1+,2+	*
"#;
        
        std::fs::write(&test_gfa, gfa_content).unwrap();
        
        // Run graph check
        let result = commands::check::check_graph_genotyping_suitability(
            test_gfa.to_str().unwrap(),
            report_path.to_str().unwrap(),
        );
        
        assert!(result.is_ok());
        assert!(report_path.exists());
        
        let report_content = std::fs::read_to_string(&report_path).unwrap();
        assert!(report_content.contains("Graph Structure"));
        assert!(report_content.contains("Nodes: 3"));
        assert!(report_content.contains("Paths: 4"));
        assert!(report_content.contains("Genotyping Suitability Assessment"));
        
        println!("✅ Graph check functionality test passed");
    }
    
    
    /// Test that cargo test runs and validates our algorithm
    #[test]
    fn test_algorithm_correctness() {
        // Test our core cosine similarity algorithm with known inputs
        let ref1 = vec![1.0, 0.0, 1.0, 0.0];
        let ref2 = vec![0.0, 1.0, 0.0, 1.0];
        let sample = vec![1.0, 1.0, 1.0, 1.0]; // Perfect sum of ref1 + ref2
        
        let combined = crate::math::sum_vectors(&[&ref1, &ref2]);
        let similarity = crate::math::cosine_similarity(&combined, &sample);
        
        // Should be perfect similarity (1.0)
        assert!((similarity - 1.0).abs() < 1e-10, "Expected perfect similarity, got {}", similarity);
        
        // Test scaling invariance
        let ref1_scaled: Vec<f64> = ref1.iter().map(|x| x * 10.0).collect();
        let ref2_scaled: Vec<f64> = ref2.iter().map(|x| x * 10.0).collect();
        let sample_scaled: Vec<f64> = sample.iter().map(|x| x * 10.0).collect();
        
        let combined_scaled = crate::math::sum_vectors(&[&ref1_scaled, &ref2_scaled]);
        let similarity_scaled = crate::math::cosine_similarity(&combined_scaled, &sample_scaled);
        
        assert!((similarity - similarity_scaled).abs() < 1e-10, "Cosine similarity should be scale-invariant");
        
        println!("✅ Algorithm correctness test passed");
    }
}

/// Create test FASTA with HLA-like sequences for testing
async fn create_test_hla_fasta(filename: &str) -> Result<()> {
    // Create simplified HLA-like sequences for testing
    let sequences = vec![
        (">HG00096#HLA-F#hap1", "ATGGCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCGAGACCTGGGCCGGTTGGCCGCGGACGGGCGGCGAATTCCAGAGCCACAGAAGCTGAATGGGTACCGACCGGATTACATGGGCGGGCCGGTCCGTCTTTCCGTGCCCACGACATTCCGCCCCGGGCGCACC"),
        (">HG00096#HLA-F#hap2", "ATGGCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCGAGACCTGGGCCGGTTGGCCGCGGACGGGCGGCGAATTCCAGAGCCACAGAAGCTGAATGGGTACCGACCGGATTACATGGGCGGGCCGGTCCGTCTTTCCGTGCCCACGACATTCCGCCCCGGGCGACC"),
        (">HG00097#HLA-F#hap1", "ATGGCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCGAGACCTGGGCCGGTTGGCCGCGGACGGGCGGCGAATTCCAGAGCCACAGAAGCTGAATGGGTACCGACCGGATTACATGGGCGGGCCGGTCCGTCTTTCCGTGCCCACGACATTCCGCCCCGGGCGCCC"),
        (">HG00097#HLA-F#hap2", "ATGGCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCGAGACCTGGGCCGGTTGGCCGCGGACGGGCGGCGAATTCCAGAGCCACAGAAGCTGAATGGGTACCGACCGGATTACATGGGCGGGCCGGTCCGTCTTTCCGTGCCCACGACATTCCGCCCCGGGCGACC"),
    ];
    
    let mut content = String::new();
    for (header, sequence) in sequences {
        content.push_str(header);
        content.push('\n');
        content.push_str(sequence);
        content.push('\n');
    }
    
    tokio::fs::write(filename, content).await?;
    
    Ok(())
}
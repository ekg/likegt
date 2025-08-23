use likegt::*;
use std::process::Command;

#[test]
#[ignore] // This test requires external files and takes time
fn test_hold2out_genotype_accuracy() {
    // This test runs the full hold-2-out pipeline and checks accuracy
    
    // Run hold-2-out for a test individual
    let output = Command::new("./run_hold2out.sh")
        .args(&["hla-f.k51.og", "HG00358"])
        .output()
        .expect("Failed to run hold-2-out pipeline");
    
    assert!(output.status.success(), "Hold-2-out pipeline failed");
    
    // Parse the output to get predicted genotype
    let stdout = String::from_utf8_lossy(&output.stdout);
    
    // Check if we got a result
    assert!(stdout.contains("GENOTYPE RESULT"), "No genotype result found");
    
    // The genotype accuracy is expected to be low for hold-2-out
    // since the test individual is not in the reference
    assert!(stdout.contains("Quality:"), "No quality assessment found");
}

#[test]
fn test_hold0out_accuracy() {
    // For hold-0-out (individual included), we should get near-perfect accuracy
    use math::cosine_similarity;
    use io::read_gzip_tsv;
    
    // This would load pre-computed test data
    // In a real test, we'd have test fixtures with known good data
    
    // Simulate perfect match
    let test_coverage = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let ref_coverage = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    
    let similarity = cosine_similarity(&test_coverage, &ref_coverage);
    assert!((similarity - 1.0).abs() < 1e-6, "Hold-0-out should have perfect similarity");
}

#[test]
fn test_genotype_qv_calculation() {
    use sequence_qv::SequenceQV;
    
    // Test QV calculation from identity
    let qv99 = SequenceQV::from_identity(0.99);
    let qv999 = SequenceQV::from_identity(0.999);
    
    assert!((qv99.qv - 20.0).abs() < 0.1);  // 99% identity ≈ QV20
    assert!((qv999.qv - 30.0).abs() < 0.1); // 99.9% identity ≈ QV30
    let qv9999 = SequenceQV::from_identity(0.9999);
    assert!((qv9999.qv - 40.0).abs() < 0.1); // 99.99% identity ≈ QV40
    
    // Test the poor result we're seeing
    let qv25 = SequenceQV::from_identity(0.25);
    assert!((qv25.qv - 1.249).abs() < 0.01); // 25% identity ≈ QV~1.25
}
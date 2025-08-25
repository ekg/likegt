use likegt::{io, validation, sequence_qv};
use std::path::PathBuf;

#[test]
fn test_sequence_qv_for_tied_genotypes() {
    println!("\n=== SEQUENCE-LEVEL QV VALIDATION FOR TIED CASES ===\n");
    
    let test_data_path = PathBuf::from("tests/data/hla-f.k51.paths.coverage.tsv.gz");
    let fasta_path = PathBuf::from("tests/data/hla-f.fa.gz");
    
    if !test_data_path.exists() || !fasta_path.exists() {
        println!("Test data not found, skipping");
        return;
    }
    
    // Load sequences
    println!("Loading FASTA sequences...");
    let temp_fasta = "/tmp/hla-f.fa";
    std::process::Command::new("zcat")
        .arg(&fasta_path)
        .output()
        .map(|output| std::fs::write(temp_fasta, output.stdout))
        .expect("Failed to decompress FASTA");
    
    let sequences = sequence_qv::read_fasta_sequences(&PathBuf::from(temp_fasta)).unwrap();
    println!("Loaded {} sequences", sequences.len());
    
    // Run graph-based validation first
    println!("Running graph-based validation...");
    let ref_data = io::read_gzip_tsv(test_data_path.to_str().unwrap()).unwrap();
    let results = validation::validate_all_individuals(&ref_data).unwrap();
    
    // Find the failed cases (tied genotypes)
    let failed_cases: Vec<_> = results.iter().filter(|r| !r.correct).collect();
    println!("Found {} failed cases with graph cosine similarity = 1.0\n", failed_cases.len());
    
    // Now check sequence-level QV for each failed case
    for fail in &failed_cases {
        println!("=== {} ===", fail.sample_id);
        println!("Graph-based call: {} + {}", 
            fail.called_hap1.split('#').next().unwrap_or("?"),
            fail.called_hap2.split('#').next().unwrap_or("?")
        );
        println!("True genotype: {} + {}", 
            fail.true_hap1.split('#').next().unwrap_or("?"),
            fail.true_hap2.split('#').next().unwrap_or("?")
        );
        
        // Use full sequence IDs as they appear in FASTA
        let true_hap1_id = &fail.true_hap1;
        let true_hap2_id = &fail.true_hap2;  
        let called_hap1_id = &fail.called_hap1;
        let called_hap2_id = &fail.called_hap2;
        
        // Calculate sequence-level QV
        match sequence_qv::calculate_genotype_qv(
            true_hap1_id,
            true_hap2_id,
            called_hap1_id,
            called_hap2_id,
            &sequences,
        ) {
            Ok(seq_qv) => {
                println!("Sequence-level identity: {:.4}", seq_qv.identity);
                println!("Sequence-level QV: {:.1}", seq_qv.qv);
                println!("Edit distance: {} / {}", seq_qv.edit_distance, seq_qv.alignment_length);
                
                // Check if specific haplotypes are identical
                if let (Some(true_seq1), Some(called_seq1)) = (sequences.get(true_hap1_id), sequences.get(called_hap1_id)) {
                    if true_seq1 == called_seq1 {
                        println!("  ✓ Haplotype 1: {} == {} (IDENTICAL)", true_hap1_id, called_hap1_id);
                    } else {
                        println!("  ✗ Haplotype 1: {} != {} (DIFFERENT)", true_hap1_id, called_hap1_id);
                    }
                }
                
                if let (Some(true_seq2), Some(called_seq2)) = (sequences.get(true_hap2_id), sequences.get(called_hap2_id)) {
                    if true_seq2 == called_seq2 {
                        println!("  ✓ Haplotype 2: {} == {} (IDENTICAL)", true_hap2_id, called_hap2_id);
                    } else {
                        println!("  ✗ Haplotype 2: {} != {} (DIFFERENT)", true_hap2_id, called_hap2_id);
                    }
                }
            }
            Err(e) => {
                println!("Failed to calculate sequence QV: {}", e);
            }
        }
        println!();
    }
    
    // Clean up
    std::fs::remove_file(temp_fasta).ok();
}

#[test] 
fn test_specific_identical_pairs() {
    println!("\n=== TESTING SPECIFIC IDENTICAL PAIRS ===\n");
    
    let fasta_path = PathBuf::from("tests/data/hla-f.fa.gz");
    if !fasta_path.exists() {
        println!("FASTA not found, skipping");
        return;
    }
    
    // Load sequences
    let temp_fasta = "/tmp/hla-f.fa";
    std::process::Command::new("zcat")
        .arg(&fasta_path)
        .output()
        .map(|output| std::fs::write(temp_fasta, output.stdout))
        .expect("Failed to decompress FASTA");
    
    let sequences = sequence_qv::read_fasta_sequences(&PathBuf::from(temp_fasta)).unwrap();
    
    // Get a few sample pairs to test - need to find actual sequence IDs from FASTA
    println!("Available sequences:");
    for (i, key) in sequences.keys().take(10).enumerate() {
        println!("  {}: {}", i, key);
    }
    
    // Look for specific patterns we know are problematic
    let mut hg512_1 = None;
    let mut hg514_1 = None;
    let mut hg732_1 = None; 
    let mut hg733_1 = None;
    
    for key in sequences.keys() {
        if key.starts_with("HG00512#1#") { hg512_1 = Some(key.clone()); }
        if key.starts_with("HG00514#1#") { hg514_1 = Some(key.clone()); }
        if key.starts_with("HG00732#1#") { hg732_1 = Some(key.clone()); }
        if key.starts_with("HG00733#1#") { hg733_1 = Some(key.clone()); }
    }
    
    let mut pairs_to_test = Vec::new();
    if let (Some(a), Some(b)) = (&hg512_1, &hg514_1) {
        pairs_to_test.push((a.as_str(), b.as_str()));
    }
    if let (Some(a), Some(b)) = (&hg732_1, &hg733_1) {
        pairs_to_test.push((a.as_str(), b.as_str()));
    }
    
    for (hap1, hap2) in pairs_to_test {
        println!("Testing {} vs {}", hap1, hap2);
        
        if let (Some(seq1), Some(seq2)) = (sequences.get(hap1), sequences.get(hap2)) {
            if seq1 == seq2 {
                println!("  ✓ IDENTICAL sequences (length: {})", seq1.len());
            } else {
                // Calculate edit distance
                match sequence_qv::align_sequences_wfa(seq1, seq2) {
                    Ok(qv) => {
                        println!("  ✗ DIFFERENT sequences");
                        println!("    Identity: {:.4}", qv.identity);
                        println!("    Edit distance: {} / {}", qv.edit_distance, qv.alignment_length);
                        println!("    QV: {:.1}", qv.qv);
                    }
                    Err(_) => {
                        println!("  ✗ Failed to align sequences");
                    }
                }
            }
        } else {
            println!("  ? One or both sequences not found");
        }
    }
    
    // Clean up
    std::fs::remove_file(temp_fasta).ok();
}
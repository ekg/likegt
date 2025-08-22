use likegt::{io, real_pipeline};
use std::path::Path;

fn main() {
    println!("=== TESTING REAL PIPELINE WITH ACTUAL 30x COVERAGE ===\n");
    
    let fasta_path = Path::new("hla-f.fa.gz");
    let gfa_path = Path::new("hla-f.k51.gfa");
    let coverage_path = Path::new("hla-f.k51.paths.coverage.tsv.gz");
    
    // Load reference coverage
    let ref_data = io::read_gzip_tsv(coverage_path.to_str().unwrap())
        .expect("Failed to load reference data");
    
    println!("Loaded {} reference haplotypes\n", ref_data.len());
    
    // Test HG00096 with REAL 30x coverage
    let result = real_pipeline::test_real_hold0(
        "HG00096",
        fasta_path,
        gfa_path,
        &ref_data,
        30,  // ACTUAL 30x coverage
    );
    
    match result {
        Ok(res) => {
            println!("\n=== FINAL RESULTS ===");
            println!("Individual: {}", res.individual);
            println!("Coverage depth: {}x", res.coverage_depth);
            println!("Reads simulated: {}", res.reads_simulated);
            println!("Correct: {}", res.correct);
            println!("Similarity: {:.4}", res.similarity);
        }
        Err(e) => {
            println!("ERROR: {}", e);
        }
    }
}
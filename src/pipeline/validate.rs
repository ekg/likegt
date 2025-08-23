use anyhow::{Result, Context};
use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;
use tokio::fs;

#[derive(Debug, Clone)]
pub struct ValidationResults {
    pub total_tests: usize,
    pub correct_genotypes: usize,
    pub accuracy: f64,
    pub average_similarity: f64,
    pub similarity_stddev: f64,
    pub individual_results: Vec<IndividualResult>,
}

#[derive(Debug, Clone)]
pub struct IndividualResult {
    pub individual_id: String,
    pub predicted_genotype: Vec<String>,
    pub true_genotype: Vec<String>,
    pub similarity_score: f64,
    pub is_correct: bool,
}

#[derive(Debug)]
pub struct ValidationConfig {
    pub graph_path: PathBuf,
    pub fasta_path: PathBuf,
    pub output_dir: PathBuf,
    pub threads: usize,
    pub kmer_size: usize,
}

/// Run hold-0-out validation test
/// Test individuals are included in the reference graph
pub async fn run_hold0_test(
    graph_path: &PathBuf,
    fasta_path: &PathBuf,
    threads: usize,
) -> Result<ValidationResults> {
    log::info!("Running hold-0-out validation");
    
    // Parse individuals from FASTA
    let individuals = parse_individuals_from_fasta(fasta_path).await?;
    log::info!("Found {} individuals to test", individuals.len());
    
    let mut individual_results = Vec::new();
    
    for (individual_id, sequences) in individuals {
        log::debug!("Testing individual: {}", individual_id);
        
        // For hold-0-out, use the full graph (individual is included in reference)
        let result = test_individual_genotyping(
            graph_path,
            &individual_id,
            &sequences,
            graph_path, // Use same graph for reference
            threads,
        ).await?;
        
        individual_results.push(result);
    }
    
    Ok(calculate_validation_statistics(individual_results))
}

/// Run hold-2-out validation test  
/// Test individuals are excluded from the reference graph
pub async fn run_hold2_test(
    graph_path: &PathBuf,
    fasta_path: &PathBuf,
    threads: usize,
) -> Result<ValidationResults> {
    log::info!("Running hold-2-out validation");
    
    let individuals = parse_individuals_from_fasta(fasta_path).await?;
    log::info!("Found {} individuals to test", individuals.len());
    
    let mut individual_results = Vec::new();
    
    for (individual_id, sequences) in individuals {
        log::debug!("Testing individual (hold-2-out): {}", individual_id);
        
        // Create hold-2-out graph (exclude this individual)
        let hold2_graph = create_hold2out_graph(
            graph_path,
            fasta_path,
            &individual_id,
        ).await?;
        
        // Test genotyping using the reduced graph
        let result = test_individual_genotyping(
            &hold2_graph,
            &individual_id,
            &sequences,
            &hold2_graph, // Use reduced graph consistently
            threads,
        ).await?;
        
        individual_results.push(result);
        
        // Clean up temporary graph
        if let Err(e) = fs::remove_file(&hold2_graph).await {
            log::warn!("Failed to cleanup hold-2-out graph: {}", e);
        }
    }
    
    Ok(calculate_validation_statistics(individual_results))
}

/// Run reference bias test
/// Test how choice of reference affects genotyping results
pub async fn run_bias_test(
    graph_path: &PathBuf,
    fasta_path: &PathBuf,
    threads: usize,
) -> Result<ValidationResults> {
    log::info!("Running reference bias test");
    
    let individuals = parse_individuals_from_fasta(fasta_path).await?;
    
    if individuals.len() < 3 {
        return Err(anyhow::anyhow!("Need at least 3 individuals for bias test"));
    }
    
    let mut individual_results = Vec::new();
    
    // Test each individual against different reference sets
    for (test_individual, test_sequences) in &individuals {
        log::debug!("Bias testing individual: {}", test_individual);
        
        let mut genotype_results = Vec::new();
        
        // Create different reference graphs by excluding different individuals
        for (exclude_individual, _) in &individuals {
            if exclude_individual == test_individual {
                continue;
            }
            
            // Create reference graph excluding one other individual
            let bias_graph = create_hold2out_graph(
                graph_path,
                fasta_path,
                exclude_individual,
            ).await?;
            
            // Test genotyping with this reference
            let result = test_individual_genotyping(
                &bias_graph,
                test_individual,
                test_sequences,
                &bias_graph,
                threads,
            ).await?;
            
            genotype_results.push(result);
            
            // Clean up
            if let Err(e) = fs::remove_file(&bias_graph).await {
                log::warn!("Failed to cleanup bias test graph: {}", e);
            }
        }
        
        // Analyze consistency across different references
        let bias_result = analyze_reference_bias(&genotype_results);
        individual_results.push(bias_result);
    }
    
    Ok(calculate_validation_statistics(individual_results))
}

async fn parse_individuals_from_fasta(fasta_path: &PathBuf) -> Result<Vec<(String, Vec<String>)>> {
    let content = fs::read_to_string(fasta_path).await?;
    let mut individuals: HashMap<String, Vec<String>> = HashMap::new();
    
    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();
    
    for line in content.lines() {
        if line.starts_with('>') {
            // Save previous sequence
            if let Some(id) = current_id.take() {
                if !current_seq.is_empty() {
                    // Extract individual name (before # or other separators)
                    let individual = id.split('#').next().unwrap_or(&id).to_string();
                    individuals.entry(individual).or_insert_with(Vec::new).push(current_seq.clone());
                }
            }
            
            current_id = Some(line[1..].to_string());
            current_seq.clear();
        } else {
            current_seq.push_str(line.trim());
        }
    }
    
    // Save last sequence
    if let Some(id) = current_id {
        if !current_seq.is_empty() {
            let individual = id.split('#').next().unwrap_or(&id).to_string();
            individuals.entry(individual).or_insert_with(Vec::new).push(current_seq);
        }
    }
    
    Ok(individuals.into_iter().collect())
}

async fn create_hold2out_graph(
    original_graph: &PathBuf,
    fasta_path: &PathBuf,
    exclude_individual: &str,
) -> Result<PathBuf> {
    log::debug!("Creating hold-2-out graph excluding: {}", exclude_individual);
    
    // Create temporary paths file without the excluded individual
    let temp_dir = std::env::temp_dir();
    let temp_paths = temp_dir.join(format!("hold2out_{}.fa", exclude_individual));
    
    // Filter FASTA to exclude the test individual
    let fasta_content = fs::read_to_string(fasta_path).await?;
    let mut filtered_content = String::new();
    let mut skip_sequence = false;
    
    for line in fasta_content.lines() {
        if line.starts_with('>') {
            let seq_id = &line[1..];
            let individual = seq_id.split('#').next().unwrap_or(seq_id);
            skip_sequence = individual == exclude_individual;
            
            if !skip_sequence {
                filtered_content.push_str(line);
                filtered_content.push('\n');
            }
        } else if !skip_sequence {
            filtered_content.push_str(line);
            filtered_content.push('\n');
        }
    }
    
    fs::write(&temp_paths, &filtered_content).await?;
    
    // Rebuild graph with filtered sequences
    let hold2_graph = temp_dir.join(format!("hold2out_graph_{}.gfa", exclude_individual));
    
    crate::pipeline::build::build_graph_from_fasta(
        temp_paths.to_str().unwrap(),
        hold2_graph.to_str().unwrap(),
        51, // Default k-mer size
        10000, // Default segment length
    ).await?;
    
    // Clean up temp paths file
    if let Err(e) = fs::remove_file(&temp_paths).await {
        log::warn!("Failed to cleanup temp paths file: {}", e);
    }
    
    Ok(hold2_graph)
}

async fn test_individual_genotyping(
    graph_path: &PathBuf,
    individual_id: &str,
    true_sequences: &[String],
    reference_graph: &PathBuf,
    _threads: usize,
) -> Result<IndividualResult> {
    // This is a placeholder for the actual genotyping test
    // In a real implementation, this would:
    // 1. Create reads from the true sequences (wgsim)
    // 2. Map reads to the graph (BWA + gfainject)  
    // 3. Get coverage with gafpack
    // 4. Run genotyping algorithm
    // 5. Compare predicted vs true genotype
    
    log::debug!("Testing genotyping for individual: {}", individual_id);
    
    // For now, return a mock result
    // TODO: Implement actual genotyping pipeline
    
    let predicted_genotype = vec![
        format!("{}_hap1", individual_id),
        format!("{}_hap2", individual_id),
    ];
    
    let true_genotype = vec![
        format!("{}_hap1", individual_id),
        format!("{}_hap2", individual_id),
    ];
    
    let similarity_score = 0.95; // Mock high similarity
    let is_correct = predicted_genotype == true_genotype;
    
    Ok(IndividualResult {
        individual_id: individual_id.to_string(),
        predicted_genotype,
        true_genotype,
        similarity_score,
        is_correct,
    })
}

fn analyze_reference_bias(results: &[IndividualResult]) -> IndividualResult {
    if results.is_empty() {
        return IndividualResult {
            individual_id: "unknown".to_string(),
            predicted_genotype: vec![],
            true_genotype: vec![],
            similarity_score: 0.0,
            is_correct: false,
        };
    }
    
    // Check consistency across different reference graphs
    let first_result = &results[0];
    let mut consistent = true;
    let mut similarity_scores = vec![first_result.similarity_score];
    
    for result in results.iter().skip(1) {
        if result.predicted_genotype != first_result.predicted_genotype {
            consistent = false;
        }
        similarity_scores.push(result.similarity_score);
    }
    
    let mean_similarity = similarity_scores.iter().sum::<f64>() / similarity_scores.len() as f64;
    
    IndividualResult {
        individual_id: first_result.individual_id.clone(),
        predicted_genotype: first_result.predicted_genotype.clone(),
        true_genotype: first_result.true_genotype.clone(),
        similarity_score: mean_similarity,
        is_correct: consistent,
    }
}

fn calculate_validation_statistics(results: Vec<IndividualResult>) -> ValidationResults {
    if results.is_empty() {
        return ValidationResults {
            total_tests: 0,
            correct_genotypes: 0,
            accuracy: 0.0,
            average_similarity: 0.0,
            similarity_stddev: 0.0,
            individual_results: results,
        };
    }
    
    let total_tests = results.len();
    let correct_genotypes = results.iter().filter(|r| r.is_correct).count();
    let accuracy = correct_genotypes as f64 / total_tests as f64;
    
    let similarities: Vec<f64> = results.iter().map(|r| r.similarity_score).collect();
    let average_similarity = similarities.iter().sum::<f64>() / similarities.len() as f64;
    
    let variance = similarities.iter()
        .map(|s| (s - average_similarity).powi(2))
        .sum::<f64>() / similarities.len() as f64;
    let similarity_stddev = variance.sqrt();
    
    ValidationResults {
        total_tests,
        correct_genotypes,
        accuracy,
        average_similarity,
        similarity_stddev,
        individual_results: results,
    }
}
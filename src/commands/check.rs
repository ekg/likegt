use anyhow::Result;
use std::collections::HashMap;
use std::fs;

/// Check if a graph is suitable for genotyping
/// 
/// This analyzes graph properties to determine genotyping suitability:
/// - Node count and complexity
/// - Path count and diversity  
/// - Connectivity and structure
/// - Sequence length distribution
pub fn check_graph_genotyping_suitability(gfa_path: &str, output_path: &str) -> Result<()> {
    log::info!("Analyzing graph: {}", gfa_path);
    
    let gfa_content = fs::read_to_string(gfa_path)?;
    let analysis = analyze_gfa_structure(&gfa_content)?;
    
    let report = generate_suitability_report(&analysis);
    
    fs::write(output_path, &report)?;
    
    log::info!("Analysis complete. Report saved to: {}", output_path);
    println!("{}", report);
    
    Ok(())
}

#[derive(Debug)]
struct GraphAnalysis {
    node_count: usize,
    edge_count: usize,
    path_count: usize,
    total_sequence_length: usize,
    node_length_stats: NodeLengthStats,
    path_length_stats: PathLengthStats,
    complexity_metrics: ComplexityMetrics,
}

#[derive(Debug)]
struct NodeLengthStats {
    min: usize,
    max: usize,
    mean: f64,
    median: usize,
}

#[derive(Debug)]
struct PathLengthStats {
    min: usize,
    max: usize,
    mean: f64,
    paths_per_sample: HashMap<String, usize>,
}

#[derive(Debug)]
struct ComplexityMetrics {
    branching_factor: f64,
    path_diversity: f64,
    structural_complexity: f64,
}

fn analyze_gfa_structure(gfa_content: &str) -> Result<GraphAnalysis> {
    let mut node_count = 0;
    let mut edge_count = 0;
    let mut path_count = 0;
    let mut total_sequence_length = 0;
    let mut node_lengths = Vec::new();
    let mut path_lengths = Vec::new();
    let mut paths_per_sample = HashMap::new();
    
    for line in gfa_content.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() {
            continue;
        }
        
        match fields[0] {
            "S" => {
                // Segment (node)
                node_count += 1;
                if fields.len() > 2 {
                    let seq_len = fields[2].len();
                    total_sequence_length += seq_len;
                    node_lengths.push(seq_len);
                }
            }
            "L" => {
                // Link (edge)
                edge_count += 1;
            }
            "P" => {
                // Path
                path_count += 1;
                if fields.len() > 2 {
                    let path_nodes = fields[2].split(',').count();
                    path_lengths.push(path_nodes);
                    
                    // Extract sample name (assuming format like "sample#chr#hap")
                    let path_name = fields[1];
                    if let Some(sample_name) = path_name.split('#').next() {
                        *paths_per_sample.entry(sample_name.to_string()).or_insert(0) += 1;
                    }
                }
            }
            _ => {}
        }
    }
    
    // Calculate statistics
    let node_length_stats = calculate_node_stats(&node_lengths);
    let path_length_stats = calculate_path_stats(&path_lengths, paths_per_sample);
    let complexity_metrics = calculate_complexity_metrics(
        node_count, 
        edge_count, 
        path_count,
        &node_lengths,
    );
    
    Ok(GraphAnalysis {
        node_count,
        edge_count,
        path_count,
        total_sequence_length,
        node_length_stats,
        path_length_stats,
        complexity_metrics,
    })
}

fn calculate_node_stats(lengths: &[usize]) -> NodeLengthStats {
    if lengths.is_empty() {
        return NodeLengthStats { min: 0, max: 0, mean: 0.0, median: 0 };
    }
    
    let mut sorted_lengths = lengths.to_vec();
    sorted_lengths.sort_unstable();
    
    NodeLengthStats {
        min: sorted_lengths[0],
        max: sorted_lengths[sorted_lengths.len() - 1],
        mean: lengths.iter().sum::<usize>() as f64 / lengths.len() as f64,
        median: sorted_lengths[sorted_lengths.len() / 2],
    }
}

fn calculate_path_stats(lengths: &[usize], paths_per_sample: HashMap<String, usize>) -> PathLengthStats {
    if lengths.is_empty() {
        return PathLengthStats { 
            min: 0, max: 0, mean: 0.0, 
            paths_per_sample 
        };
    }
    
    let mut sorted_lengths = lengths.to_vec();
    sorted_lengths.sort_unstable();
    
    PathLengthStats {
        min: sorted_lengths[0],
        max: sorted_lengths[sorted_lengths.len() - 1],
        mean: lengths.iter().sum::<usize>() as f64 / lengths.len() as f64,
        paths_per_sample,
    }
}

fn calculate_complexity_metrics(
    node_count: usize,
    edge_count: usize, 
    path_count: usize,
    node_lengths: &[usize],
) -> ComplexityMetrics {
    // Branching factor: average number of edges per node
    let branching_factor = if node_count > 0 {
        edge_count as f64 / node_count as f64
    } else {
        0.0
    };
    
    // Path diversity: how many paths relative to nodes
    let path_diversity = if node_count > 0 {
        path_count as f64 / node_count as f64
    } else {
        0.0
    };
    
    // Structural complexity: coefficient of variation in node lengths
    let structural_complexity = if !node_lengths.is_empty() {
        let mean = node_lengths.iter().sum::<usize>() as f64 / node_lengths.len() as f64;
        let variance = node_lengths.iter()
            .map(|&x| (x as f64 - mean).powi(2))
            .sum::<f64>() / node_lengths.len() as f64;
        let stddev = variance.sqrt();
        
        if mean > 0.0 {
            stddev / mean
        } else {
            0.0
        }
    } else {
        0.0
    };
    
    ComplexityMetrics {
        branching_factor,
        path_diversity,
        structural_complexity,
    }
}

fn generate_suitability_report(analysis: &GraphAnalysis) -> String {
    let mut report = String::new();
    
    report.push_str("# Graph Genotyping Suitability Analysis\n\n");
    
    // Basic statistics
    report.push_str("## Graph Structure\n");
    report.push_str(&format!("- Nodes: {}\n", analysis.node_count));
    report.push_str(&format!("- Edges: {}\n", analysis.edge_count));
    report.push_str(&format!("- Paths: {}\n", analysis.path_count));
    report.push_str(&format!("- Total sequence length: {} bp\n\n", analysis.total_sequence_length));
    
    // Node statistics
    report.push_str("## Node Length Distribution\n");
    report.push_str(&format!("- Min: {} bp\n", analysis.node_length_stats.min));
    report.push_str(&format!("- Max: {} bp\n", analysis.node_length_stats.max));
    report.push_str(&format!("- Mean: {:.1} bp\n", analysis.node_length_stats.mean));
    report.push_str(&format!("- Median: {} bp\n\n", analysis.node_length_stats.median));
    
    // Path statistics
    report.push_str("## Path Statistics\n");
    report.push_str(&format!("- Unique samples: {}\n", analysis.path_length_stats.paths_per_sample.len()));
    report.push_str(&format!("- Average path length: {:.1} nodes\n", analysis.path_length_stats.mean));
    report.push_str("- Paths per sample:\n");
    
    for (sample, count) in &analysis.path_length_stats.paths_per_sample {
        report.push_str(&format!("  - {}: {} paths\n", sample, count));
    }
    report.push_str("\n");
    
    // Complexity metrics
    report.push_str("## Complexity Metrics\n");
    report.push_str(&format!("- Branching factor: {:.2}\n", analysis.complexity_metrics.branching_factor));
    report.push_str(&format!("- Path diversity: {:.2}\n", analysis.complexity_metrics.path_diversity));
    report.push_str(&format!("- Structural complexity: {:.2}\n\n", analysis.complexity_metrics.structural_complexity));
    
    // Suitability assessment
    report.push_str("## Genotyping Suitability Assessment\n");
    
    let mut score = 0;
    let mut recommendations = Vec::new();
    
    // Check minimum requirements
    if analysis.node_count < 100 {
        recommendations.push("Graph may be too simple (< 100 nodes)");
    } else if analysis.node_count > 100000 {
        recommendations.push("Graph may be too complex (> 100k nodes) - consider chunking");
    } else {
        score += 1;
    }
    
    if analysis.path_count < 4 {
        recommendations.push("Too few paths (< 4) for meaningful genotyping");
    } else if analysis.path_count >= 10 {
        score += 1;
    }
    
    if analysis.complexity_metrics.branching_factor < 1.5 {
        recommendations.push("Low branching factor suggests linear structure");
    } else if analysis.complexity_metrics.branching_factor <= 4.0 {
        score += 1;
    } else {
        recommendations.push("High branching factor may complicate genotyping");
    }
    
    if analysis.complexity_metrics.path_diversity >= 0.1 {
        score += 1;
    } else {
        recommendations.push("Low path diversity may limit genotyping resolution");
    }
    
    // Overall assessment
    match score {
        4 => report.push_str("✅ **EXCELLENT**: Graph appears highly suitable for genotyping\n"),
        3 => report.push_str("✅ **GOOD**: Graph appears suitable for genotyping\n"),
        2 => report.push_str("⚠️  **MODERATE**: Graph may work but with limitations\n"),
        _ => report.push_str("❌ **POOR**: Graph may not be suitable for genotyping\n"),
    }
    
    if !recommendations.is_empty() {
        report.push_str("\n### Recommendations:\n");
        for rec in recommendations {
            report.push_str(&format!("- {}\n", rec));
        }
    }
    
    report
}
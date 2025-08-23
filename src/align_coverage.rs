use anyhow::Result;
use std::collections::HashSet;

/// Aligns coverage data between odgi and gafpack outputs
/// 
/// odgi paths -H outputs only visited nodes (subset)
/// gafpack outputs all nodes in the graph
/// 
/// This function ensures both have the same nodes by either:
/// 1. Padding odgi output with zeros for missing nodes, or
/// 2. Filtering gafpack output to only include visited nodes
pub fn align_coverage_dimensions(
    odgi_headers: &[String],
    odgi_coverage: &[Vec<f64>],
    gafpack_headers: &[String],
    gafpack_coverage: &[f64],
) -> Result<(Vec<String>, Vec<Vec<f64>>, Vec<f64>)> {
    // Extract node IDs from headers
    let odgi_nodes: HashSet<String> = odgi_headers
        .iter()
        .filter(|h| h.starts_with("node."))
        .cloned()
        .collect();
    
    let gafpack_nodes: HashSet<String> = gafpack_headers
        .iter()
        .filter(|h| h.starts_with("node."))
        .cloned()
        .collect();
    
    log::info!(
        "Aligning coverage: odgi has {} nodes, gafpack has {} nodes",
        odgi_nodes.len(),
        gafpack_nodes.len()
    );
    
    // Strategy 1: If odgi is subset, use only those nodes
    if odgi_nodes.len() < gafpack_nodes.len() {
        log::info!("Using odgi nodes as reference (visited nodes only)");
        
        // Find indices of odgi nodes in gafpack
        // Note: we need to account for non-node columns like "path.name" or "sample"
        let mut gafpack_indices = Vec::new();
        let mut node_column_start = 0;
        
        // Find where node columns start in gafpack
        for (i, header) in gafpack_headers.iter().enumerate() {
            if header.starts_with("node.") {
                node_column_start = i;
                break;
            }
        }
        
        // For each odgi node column, find its position in gafpack
        for (_odgi_idx, odgi_header) in odgi_headers.iter().enumerate() {
            if odgi_header.starts_with("node.") {
                // Find this node in gafpack headers
                if let Some(gafpack_idx) = gafpack_headers.iter().position(|h| h == odgi_header) {
                    // The actual data index is relative to where nodes start
                    let data_idx = gafpack_idx - node_column_start;
                    if data_idx < gafpack_coverage.len() {
                        gafpack_indices.push(data_idx);
                    } else {
                        log::warn!("Index {} out of bounds for gafpack coverage", data_idx);
                        return Err(anyhow::anyhow!(
                            "Index out of bounds when aligning coverage"
                        ));
                    }
                } else {
                    // Node in odgi but not in gafpack - shouldn't happen
                    log::warn!("Node {} in odgi but not in gafpack", odgi_header);
                    return Err(anyhow::anyhow!(
                        "Node {} found in odgi but not in gafpack",
                        odgi_header
                    ));
                }
            }
        }
        
        // Filter gafpack coverage to match odgi nodes
        let filtered_gafpack: Vec<f64> = gafpack_indices
            .iter()
            .map(|&idx| gafpack_coverage[idx])
            .collect();
        
        // Return aligned data
        Ok((
            odgi_headers.to_vec(),
            odgi_coverage.to_vec(),
            filtered_gafpack,
        ))
    } else if gafpack_nodes.len() < odgi_nodes.len() {
        // This shouldn't happen in practice
        log::warn!("Gafpack has fewer nodes than odgi - unexpected");
        Err(anyhow::anyhow!(
            "Unexpected: gafpack has fewer nodes ({}) than odgi ({})",
            gafpack_nodes.len(),
            odgi_nodes.len()
        ))
    } else {
        // Same number of nodes
        log::info!("Same number of nodes in both files");
        Ok((
            odgi_headers.to_vec(),
            odgi_coverage.to_vec(),
            gafpack_coverage.to_vec(),
        ))
    }
}

/// Alternative approach: Use gafpack for both reference and sample
/// This ensures consistent dimensions
pub fn use_gafpack_for_both(
    reference_gaf: &str,
    sample_gaf: &str,
) -> Result<()> {
    log::info!("Using gafpack for both reference and sample coverage");
    log::info!("Reference GAF: {}", reference_gaf);
    log::info!("Sample GAF: {}", sample_gaf);
    
    // This would require creating GAF files for reference paths
    // Then using gafpack on both
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_align_coverage() {
        let odgi_headers = vec![
            "path.name".to_string(),
            "node.1".to_string(),
            "node.3".to_string(),
        ];
        
        let odgi_coverage = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
        ];
        
        let gafpack_headers = vec![
            "sample".to_string(),
            "node.1".to_string(),
            "node.2".to_string(),
            "node.3".to_string(),
        ];
        
        // Gafpack coverage for a single sample (3 nodes: node.1, node.2, node.3)
        let gafpack_coverage = vec![2.0, 0.0, 3.0];
        
        let (headers, _ref_cov, sample_cov) = align_coverage_dimensions(
            &odgi_headers,
            &odgi_coverage,
            &gafpack_headers,
            &gafpack_coverage,
        ).unwrap();
        
        // We should get back odgi headers (which has only node.1 and node.3)
        assert_eq!(headers.len(), 3); // path.name, node.1, node.3
        
        // Sample coverage should be filtered to match odgi nodes
        // Original gafpack: [2.0, 0.0, 3.0] for [node.1, node.2, node.3]
        // Filtered to match odgi: [2.0, 3.0] for [node.1, node.3]
        assert_eq!(sample_cov.len(), 2);
        assert_eq!(sample_cov, vec![2.0, 3.0]);
    }
}
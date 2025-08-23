use crate::io::CoverageData;
use anyhow::Result;
use std::collections::HashSet;

/// Coverage vectors for genotyping
#[derive(Debug)]
pub struct GenotypeData {
    pub reference_haplotypes: CoverageData,
    pub sample_coverage: Vec<f64>,
    pub sample_id: String,
}

impl GenotypeData {
    pub fn new(
        reference_haplotypes: CoverageData,
        sample_data: CoverageData,
        sample_id: String,
    ) -> Result<Self> {
        // Ensure sample data has exactly one coverage vector
        if sample_data.coverages.is_empty() {
            return Err(anyhow::anyhow!("No sample coverage data found"));
        }
        
        if sample_data.coverages.len() > 1 {
            log::warn!("Multiple sample coverage vectors found, using first one");
        }
        
        let sample_coverage = sample_data.coverages[0].clone();
        
        // Handle dimension mismatch between reference and sample
        // Reference from odgi paths -H only includes visited nodes
        // Sample from gafpack includes all nodes in the graph
        let ref_len = reference_haplotypes.coverages[0].len();
        let sample_len = sample_coverage.len();
        
        let aligned_sample_coverage = if ref_len != sample_len {
            log::warn!(
                "Dimension mismatch: reference has {} nodes, sample has {} nodes",
                ref_len,
                sample_len
            );
            
            if ref_len < sample_len {
                // Truncate sample to reference size (assume first nodes match)
                log::info!("Truncating sample coverage to match reference dimensions");
                sample_coverage[..ref_len].to_vec()
            } else {
                // Pad sample with zeros
                log::info!("Padding sample coverage to match reference dimensions");
                let mut padded = sample_coverage.clone();
                padded.resize(ref_len, 0.0);
                padded
            }
        } else {
            sample_coverage
        };
        
        // Verify all reference vectors have same length
        for (i, cov) in reference_haplotypes.coverages.iter().enumerate() {
            if cov.len() != ref_len {
                return Err(anyhow::anyhow!(
                    "Reference coverage vectors have inconsistent lengths: {} has {} nodes, expected {}",
                    reference_haplotypes.ids[i],
                    cov.len(),
                    ref_len
                ));
            }
        }
        
        Ok(Self {
            reference_haplotypes,
            sample_coverage: aligned_sample_coverage,
            sample_id,
        })
    }
    
    pub fn filter_blacklist(&mut self, blacklist: &HashSet<String>) {
        if blacklist.is_empty() {
            return;
        }
        
        let mut keep_indices = Vec::new();
        
        for (i, id) in self.reference_haplotypes.ids.iter().enumerate() {
            let should_exclude = blacklist.iter().any(|pattern| id.contains(pattern));
            if !should_exclude {
                keep_indices.push(i);
            }
        }
        
        let filtered_ids: Vec<String> = keep_indices
            .iter()
            .map(|&i| self.reference_haplotypes.ids[i].clone())
            .collect();
            
        let filtered_coverages: Vec<Vec<f64>> = keep_indices
            .iter()
            .map(|&i| self.reference_haplotypes.coverages[i].clone())
            .collect();
        
        let removed = self.reference_haplotypes.ids.len() - filtered_ids.len();
        if removed > 0 {
            log::info!("Filtered out {} haplotypes using blacklist", removed);
        }
        
        self.reference_haplotypes.ids = filtered_ids;
        self.reference_haplotypes.coverages = filtered_coverages;
    }
    
    pub fn num_haplotypes(&self) -> usize {
        self.reference_haplotypes.len()
    }
    
    pub fn num_nodes(&self) -> usize {
        self.sample_coverage.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::CoverageData;
    
    #[test]
    fn test_genotype_data_creation() {
        let mut ref_data = CoverageData::new();
        ref_data.ids.push("hap1".to_string());
        ref_data.coverages.push(vec![1.0, 2.0, 3.0]);
        ref_data.ids.push("hap2".to_string());
        ref_data.coverages.push(vec![4.0, 5.0, 6.0]);
        
        let mut sample_data = CoverageData::new();
        sample_data.ids.push("sample".to_string());
        sample_data.coverages.push(vec![5.0, 7.0, 9.0]);
        
        let genotype_data = GenotypeData::new(
            ref_data,
            sample_data,
            "test_sample".to_string()
        ).unwrap();
        
        assert_eq!(genotype_data.num_haplotypes(), 2);
        assert_eq!(genotype_data.num_nodes(), 3);
        assert_eq!(genotype_data.sample_coverage, vec![5.0, 7.0, 9.0]);
    }
    
    #[test]
    fn test_genotype_data_mismatched_lengths() {
        let mut ref_data = CoverageData::new();
        ref_data.ids.push("hap1".to_string());
        ref_data.coverages.push(vec![1.0, 2.0, 3.0]);
        
        let mut sample_data = CoverageData::new();
        sample_data.ids.push("sample".to_string());
        sample_data.coverages.push(vec![5.0, 7.0]); // Wrong length
        
        let result = GenotypeData::new(
            ref_data,
            sample_data,
            "test_sample".to_string()
        );
        
        // Now we handle dimension mismatch by truncating/padding
        assert!(result.is_ok());
        let data = result.unwrap();
        // Sample should be padded to match reference length
        assert_eq!(data.sample_coverage.len(), 3);
    }
    
    #[test]
    fn test_filter_blacklist() {
        let mut ref_data = CoverageData::new();
        ref_data.ids.push("keep1".to_string());
        ref_data.coverages.push(vec![1.0, 2.0]);
        ref_data.ids.push("exclude_this".to_string());
        ref_data.coverages.push(vec![3.0, 4.0]);
        ref_data.ids.push("keep2".to_string());
        ref_data.coverages.push(vec![5.0, 6.0]);
        
        let mut sample_data = CoverageData::new();
        sample_data.ids.push("sample".to_string());
        sample_data.coverages.push(vec![7.0, 8.0]);
        
        let mut genotype_data = GenotypeData::new(
            ref_data,
            sample_data,
            "test".to_string()
        ).unwrap();
        
        let mut blacklist = HashSet::new();
        blacklist.insert("exclude".to_string());
        
        genotype_data.filter_blacklist(&blacklist);
        
        assert_eq!(genotype_data.num_haplotypes(), 2);
        assert_eq!(genotype_data.reference_haplotypes.ids, vec!["keep1", "keep2"]);
    }
    
    #[test]
    fn test_empty_blacklist() {
        let mut ref_data = CoverageData::new();
        ref_data.ids.push("hap1".to_string());
        ref_data.coverages.push(vec![1.0, 2.0]);
        ref_data.ids.push("hap2".to_string());
        ref_data.coverages.push(vec![3.0, 4.0]);
        
        let mut sample_data = CoverageData::new();
        sample_data.ids.push("sample".to_string());
        sample_data.coverages.push(vec![5.0, 6.0]);
        
        let mut genotype_data = GenotypeData::new(
            ref_data,
            sample_data,
            "test".to_string()
        ).unwrap();
        
        let blacklist = HashSet::new();
        genotype_data.filter_blacklist(&blacklist);
        
        assert_eq!(genotype_data.num_haplotypes(), 2);
    }
}
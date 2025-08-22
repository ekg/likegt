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
        
        // Verify all coverage vectors have same length
        let expected_len = sample_coverage.len();
        for (i, cov) in reference_haplotypes.coverages.iter().enumerate() {
            if cov.len() != expected_len {
                return Err(anyhow::anyhow!(
                    "Coverage vector length mismatch: {} has {} nodes, expected {}",
                    reference_haplotypes.ids[i],
                    cov.len(),
                    expected_len
                ));
            }
        }
        
        Ok(Self {
            reference_haplotypes,
            sample_coverage,
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
        
        assert!(result.is_err());
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
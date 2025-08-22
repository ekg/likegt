#[cfg(test)]
mod tests {
    use super::super::*;
    use std::io::Write;
    use tempfile::tempdir;
    
    #[test]
    fn test_coverage_data_creation() {
        let mut data = CoverageData::new();
        assert_eq!(data.len(), 0);
        assert!(data.is_empty());
        
        data.ids.push("test".to_string());
        data.coverages.push(vec![1.0, 2.0, 3.0]);
        
        assert_eq!(data.len(), 1);
        assert!(!data.is_empty());
    }
    
    #[test]
    fn test_read_write_gzip_tsv() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test.tsv.gz");
        
        // Create test data
        use flate2::write::GzEncoder;
        use flate2::Compression;
        
        let file = std::fs::File::create(&file_path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        
        // Write header and data
        writeln!(encoder, "path\tnode1\tnode2\tnode3").unwrap();
        writeln!(encoder, "hap1\t1.0\t2.0\t3.0").unwrap();
        writeln!(encoder, "hap2\t4.0\t5.0\t6.0").unwrap();
        encoder.finish().unwrap();
        
        // Read it back
        let data = read_gzip_tsv(file_path.to_str().unwrap()).unwrap();
        
        assert_eq!(data.len(), 2);
        assert_eq!(data.ids[0], "hap1");
        assert_eq!(data.ids[1], "hap2");
        assert_eq!(data.coverages[0], vec![1.0, 2.0, 3.0]);
        assert_eq!(data.coverages[1], vec![4.0, 5.0, 6.0]);
    }
    
    #[test]
    fn test_write_genotype_result() {
        let dir = tempdir().unwrap();
        let output_dir = dir.path().to_str().unwrap();
        
        write_genotype_result(
            output_dir,
            "sample1",
            &vec!["hap1".to_string(), "hap2".to_string()],
            0.987654,
            2,
        ).unwrap();
        
        let output_file = dir.path().join("sample1.genotype.tsv");
        assert!(output_file.exists());
        
        let contents = std::fs::read_to_string(output_file).unwrap();
        assert!(contents.contains("sample1"));
        assert!(contents.contains("hap1"));
        assert!(contents.contains("hap2"));
        assert!(contents.contains("0.987654"));
    }
}
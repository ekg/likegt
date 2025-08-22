use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use flate2::read::GzDecoder;
use csv::ReaderBuilder;
use anyhow::{Result, Context};

#[derive(Debug, Clone)]
pub struct CoverageData {
    pub ids: Vec<String>,
    pub coverages: Vec<Vec<f64>>,
}

impl CoverageData {
    pub fn new() -> Self {
        Self {
            ids: Vec::new(),
            coverages: Vec::new(),
        }
    }
    
    pub fn len(&self) -> usize {
        self.ids.len()
    }
    
    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }
}

/// Read a gzip-compressed TSV file containing coverage data
/// Format: first column is ID, remaining columns are coverage values
pub fn read_gzip_tsv(path: &str) -> Result<CoverageData> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open file: {}", path))?;
    
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);
    
    let mut csv_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(reader);
    
    let mut data = CoverageData::new();
    
    for result in csv_reader.records() {
        let record = result.with_context(|| "Failed to read CSV record")?;
        
        if record.is_empty() {
            continue;
        }
        
        // First field is the ID
        let id = record.get(0)
            .ok_or_else(|| anyhow::anyhow!("Missing ID field"))?
            .to_string();
        
        // Remaining fields are coverage values
        let mut coverage = Vec::new();
        for i in 1..record.len() {
            let value = record.get(i)
                .ok_or_else(|| anyhow::anyhow!("Missing coverage field at index {}", i))?
                .parse::<f64>()
                .with_context(|| format!("Failed to parse coverage value: {}", record.get(i).unwrap()))?;
            coverage.push(value);
        }
        
        data.ids.push(id);
        data.coverages.push(coverage);
    }
    
    log::info!("Read {} coverage vectors from {}", data.len(), path);
    
    Ok(data)
}

/// Write results to TSV files
pub fn write_genotype_result(
    output_dir: &str,
    sample_id: &str,
    haplotypes: &[String],
    similarity: f64,
    ploidy: usize,
) -> Result<()> {
    use std::fs;
    use std::io::Write;
    
    // Create output directory if it doesn't exist
    fs::create_dir_all(output_dir)?;
    
    let output_path = Path::new(output_dir).join(format!("{}.genotype.tsv", sample_id));
    let mut file = fs::File::create(&output_path)?;
    
    // Write header
    let mut header = vec!["sample.id".to_string()];
    for i in 1..=ploidy {
        header.push(format!("haplotype.{}", i));
    }
    header.push("cosine.similarity".to_string());
    writeln!(file, "{}", header.join("\t"))?;
    
    // Write data
    let mut row = vec![sample_id.to_string()];
    for hap in haplotypes {
        row.push(hap.clone());
    }
    row.push(format!("{:.16}", similarity));
    writeln!(file, "{}", row.join("\t"))?;
    
    log::info!("Wrote genotype results to {}", output_path.display());
    
    Ok(())
}

/// Write all combinations sorted by similarity
pub fn write_sorted_combinations(
    output_dir: &str,
    sample_id: &str,
    combinations: &[(Vec<String>, f64)],
    ploidy: usize,
) -> Result<()> {
    use std::fs;
    use std::io::Write;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    
    let output_path = Path::new(output_dir).join(format!("{}.sorted_combos.tsv.gz", sample_id));
    let file = fs::File::create(&output_path)?;
    let mut encoder = GzEncoder::new(file, Compression::default());
    
    // Write header
    let mut header = Vec::new();
    for i in 1..=ploidy {
        header.push(format!("haplotype.{}", i));
    }
    header.push("cosine.similarity".to_string());
    writeln!(encoder, "{}", header.join("\t"))?;
    
    // Write combinations
    for (haplotypes, similarity) in combinations {
        let mut row = Vec::new();
        for hap in haplotypes {
            row.push(hap.clone());
        }
        row.push(format!("{:.16}", similarity));
        writeln!(encoder, "{}", row.join("\t"))?;
    }
    
    encoder.finish()?;
    log::info!("Wrote {} combinations to {}", combinations.len(), output_path.display());
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_coverage_data_creation() {
        let mut data = CoverageData::new();
        assert_eq!(data.len(), 0);
        
        data.ids.push("test".to_string());
        data.coverages.push(vec![1.0, 2.0, 3.0]);
        
        assert_eq!(data.len(), 1);
    }
}
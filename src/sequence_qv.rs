use anyhow::Result;
use lib_wfa2::affine_wavefront::{
    AffineWavefronts, AlignmentScope, AlignmentSpan, AlignmentStatus, HeuristicStrategy, MemoryMode,
};
use std::path::Path;
use bio::io::fasta;

/// Calculate QV based on actual sequence alignment using biWFA
pub struct SequenceQV {
    pub edit_distance: usize,
    pub alignment_length: usize,
    pub identity: f64,
    pub qv: f64,
}

impl SequenceQV {
    /// Calculate QV from identity percentage
    pub fn from_identity(identity: f64) -> Self {
        // QV = -10 * log10(1 - identity)
        // Cap at QV 60 for perfect matches
        let error_rate = (1.0 - identity).max(0.000001);
        let qv = (-10.0 * error_rate.log10()).min(60.0).max(0.0);
        
        Self {
            edit_distance: 0,
            alignment_length: 0,
            identity,
            qv,
        }
    }
    
    /// Calculate QV from edit distance and alignment length
    pub fn from_edit_distance(edit_distance: usize, alignment_length: usize) -> Self {
        let identity = if alignment_length > 0 {
            1.0 - (edit_distance as f64 / alignment_length as f64)
        } else {
            1.0
        };
        
        let error_rate = (1.0 - identity).max(0.000001);
        let qv = (-10.0 * error_rate.log10()).min(60.0).max(0.0);
        
        Self {
            edit_distance,
            alignment_length,
            identity,
            qv,
        }
    }
}

/// Align two sequences using biWFA
pub fn align_sequences_wfa(seq1: &str, seq2: &str) -> Result<SequenceQV> {
    // Configure WFA parameters (similar to allwave)
    let mut aligner = AffineWavefronts::with_penalties_and_memory_mode(
        0,  // match score (0 for affine)
        4,  // mismatch penalty
        6,  // gap opening penalty
        2,  // gap extension penalty
        MemoryMode::Medium,
    );
    
    // Configure alignment settings
    aligner.set_alignment_scope(AlignmentScope::Alignment);
    aligner.set_alignment_span(AlignmentSpan::End2End);
    aligner.set_heuristic(&HeuristicStrategy::None);
    
    // Perform alignment
    let status = aligner.align(seq1.as_bytes(), seq2.as_bytes());
    
    match status {
        AlignmentStatus::Completed => {
            let cigar = aligner.cigar();
            
            // Parse CIGAR to calculate identity
            let mut matches = 0;
            let mut mismatches = 0;
            let mut insertions = 0;
            let mut deletions = 0;
            
            for &op in cigar {
                match op {
                    b'M' => matches += 1,      // Match in WFA2
                    b'X' => mismatches += 1,   // Mismatch
                    b'I' => deletions += 1,    // WFA2 'I' means standard 'D'
                    b'D' => insertions += 1,   // WFA2 'D' means standard 'I'
                    _ => {}
                }
            }
            
            let alignment_length = matches + mismatches + insertions + deletions;
            let edit_distance = mismatches + insertions + deletions;
            let _identity = matches as f64 / alignment_length.max(1) as f64;
            
            Ok(SequenceQV::from_edit_distance(edit_distance, alignment_length))
        }
        _ => {
            // Alignment failed, return 0 identity
            Ok(SequenceQV::from_identity(0.0))
        }
    }
}

/// Read a FASTA file and return sequences as a HashMap
pub fn read_fasta_sequences(path: &Path) -> Result<std::collections::HashMap<String, String>> {
    let reader = fasta::Reader::from_file(path)?;
    let mut sequences = std::collections::HashMap::new();
    
    for record in reader.records() {
        let record = record?;
        let id = record.id().to_string();
        let seq = String::from_utf8(record.seq().to_vec())?;
        sequences.insert(id, seq);
    }
    
    Ok(sequences)
}

/// Calculate sequence-based QV for a genotype call
pub fn calculate_genotype_qv(
    true_hap1: &str,
    true_hap2: &str,
    called_hap1: &str,
    called_hap2: &str,
    sequences: &std::collections::HashMap<String, String>,
) -> Result<SequenceQV> {
    // Get sequences
    let true_seq1 = sequences.get(true_hap1)
        .ok_or_else(|| anyhow::anyhow!("Sequence not found: {}", true_hap1))?;
    let true_seq2 = sequences.get(true_hap2)
        .ok_or_else(|| anyhow::anyhow!("Sequence not found: {}", true_hap2))?;
    let called_seq1 = sequences.get(called_hap1)
        .ok_or_else(|| anyhow::anyhow!("Sequence not found: {}", called_hap1))?;
    let called_seq2 = sequences.get(called_hap2)
        .ok_or_else(|| anyhow::anyhow!("Sequence not found: {}", called_hap2))?;
    
    // Try both possible pairings and choose the best
    // Pairing 1: true1-called1, true2-called2
    let qv1_1 = align_sequences_wfa(true_seq1, called_seq1)?;
    let qv1_2 = align_sequences_wfa(true_seq2, called_seq2)?;
    let total_identity1 = (qv1_1.identity + qv1_2.identity) / 2.0;
    
    // Pairing 2: true1-called2, true2-called1
    let qv2_1 = align_sequences_wfa(true_seq1, called_seq2)?;
    let qv2_2 = align_sequences_wfa(true_seq2, called_seq1)?;
    let total_identity2 = (qv2_1.identity + qv2_2.identity) / 2.0;
    
    // Return the best pairing
    if total_identity1 >= total_identity2 {
        let total_edit = qv1_1.edit_distance + qv1_2.edit_distance;
        let total_len = qv1_1.alignment_length + qv1_2.alignment_length;
        Ok(SequenceQV::from_edit_distance(total_edit, total_len))
    } else {
        let total_edit = qv2_1.edit_distance + qv2_2.edit_distance;
        let total_len = qv2_1.alignment_length + qv2_2.alignment_length;
        Ok(SequenceQV::from_edit_distance(total_edit, total_len))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_qv_calculation() {
        // Perfect match
        let qv = SequenceQV::from_identity(1.0);
        assert_eq!(qv.qv, 60.0);
        
        // 99% identity = QV ~20
        let qv = SequenceQV::from_identity(0.99);
        assert!(qv.qv > 19.0 && qv.qv < 21.0);
        
        // 90% identity = QV ~10
        let qv = SequenceQV::from_identity(0.90);
        assert!(qv.qv > 9.0 && qv.qv < 11.0);
        
        // 50% identity = QV ~3
        let qv = SequenceQV::from_identity(0.50);
        assert!(qv.qv > 2.0 && qv.qv < 4.0);
    }
    
    #[test]
    fn test_edit_distance_qv() {
        // No edits
        let qv = SequenceQV::from_edit_distance(0, 1000);
        assert_eq!(qv.identity, 1.0);
        assert_eq!(qv.qv, 60.0);
        
        // 1 edit in 100 bases
        let qv = SequenceQV::from_edit_distance(1, 100);
        assert_eq!(qv.identity, 0.99);
        
        // 10 edits in 100 bases
        let qv = SequenceQV::from_edit_distance(10, 100);
        assert_eq!(qv.identity, 0.90);
    }
    
    #[test]
    fn test_wfa_alignment() {
        // Test perfect match
        let seq1 = "ACGTACGTACGT";
        let seq2 = "ACGTACGTACGT";
        let qv = align_sequences_wfa(seq1, seq2).unwrap();
        assert_eq!(qv.identity, 1.0);
        assert_eq!(qv.edit_distance, 0);
        
        // Test with mismatches
        let seq1 = "ACGTACGTACGT";
        let seq2 = "ACGTACCTACGT";  // One mismatch
        let qv = align_sequences_wfa(seq1, seq2).unwrap();
        assert!(qv.identity > 0.9);
        assert_eq!(qv.edit_distance, 1);
    }
}
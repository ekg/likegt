#!/usr/bin/env python3
"""
Check the accuracy of genotype assignment and calculate QV.
Compare predicted genotype sequences to actual test sequences.
"""

import subprocess
import sys
from pathlib import Path

def get_sequence_from_fasta(fasta_file, seq_name):
    """Extract a sequence from a FASTA file"""
    cmd = f"samtools faidx {fasta_file} {seq_name}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        # Try with zcat if it's compressed
        cmd = f"zcat {fasta_file} | awk -v name='>{seq_name}' 'BEGIN{{p=0}} /^>/{{if($0==name){{p=1}}else{{p=0}}}} p' | grep -v '^>'"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    # Return just the sequence, no header
    lines = result.stdout.strip().split('\n')
    return ''.join([line for line in lines if not line.startswith('>')])

def calculate_identity(seq1, seq2):
    """Calculate sequence identity between two sequences"""
    if len(seq1) != len(seq2):
        # Align if different lengths - for now just truncate to shorter
        min_len = min(len(seq1), len(seq2))
        seq1 = seq1[:min_len]
        seq2 = seq2[:min_len]
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    identity = matches / len(seq1)
    return identity

def identity_to_qv(identity):
    """Convert sequence identity to QV score"""
    import math
    if identity >= 1.0:
        return 60.0  # Cap at QV60
    elif identity <= 0:
        return 0.0
    else:
        # QV = -10 * log10(1 - identity)
        return -10 * math.log10(1 - identity)

def main():
    if len(sys.argv) < 2:
        print("Usage: check_genotype_accuracy.py <test_individual> [predicted_hap1] [predicted_hap2]")
        print("Example: check_genotype_accuracy.py HG00171 HG02769#1#haplotype1-0000020:29626886-29650800 HG03248#2#haplotype2-0000148:29762683-29786601")
        sys.exit(1)
    
    test_ind = sys.argv[1]
    
    # If predictions provided, use them; otherwise read from last run
    if len(sys.argv) >= 4:
        pred_hap1 = sys.argv[2]
        pred_hap2 = sys.argv[3]
    else:
        # Try to parse from last genotyping output
        print("No predictions provided, using test data...")
        pred_hap1 = "HG02769#1#haplotype1-0000020:29626886-29650800"
        pred_hap2 = "HG03248#2#haplotype2-0000148:29762683-29786601"
    
    ref_fasta = "hla-f.fa.gz"
    
    print(f"Test individual: {test_ind}")
    print(f"Predicted haplotype 1: {pred_hap1}")
    print(f"Predicted haplotype 2: {pred_hap2}")
    print()
    
    # Get actual test sequences
    print("Extracting actual test sequences...")
    test_hap1_name = f"{test_ind}#1#haplotype1"
    test_hap2_name = f"{test_ind}#2#haplotype2"
    
    # Find exact sequence names
    cmd = f"zcat {ref_fasta} | grep '^>{test_ind}#' | sed 's/>// '"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    test_seq_names = result.stdout.strip().split('\n')
    
    if len(test_seq_names) < 2:
        print(f"Error: Could not find sequences for {test_ind}")
        sys.exit(1)
    
    print(f"Found test sequences: {test_seq_names}")
    
    # Extract sequences
    test_seq1 = get_sequence_from_fasta(ref_fasta, test_seq_names[0])
    test_seq2 = get_sequence_from_fasta(ref_fasta, test_seq_names[1])
    pred_seq1 = get_sequence_from_fasta(ref_fasta, pred_hap1)
    pred_seq2 = get_sequence_from_fasta(ref_fasta, pred_hap2)
    
    print(f"\nSequence lengths:")
    print(f"  Test hap1: {len(test_seq1)} bp")
    print(f"  Test hap2: {len(test_seq2)} bp")
    print(f"  Pred hap1: {len(pred_seq1)} bp")
    print(f"  Pred hap2: {len(pred_seq2)} bp")
    
    # Calculate best matching (test could be in different order)
    # Compare all combinations
    id_11 = calculate_identity(test_seq1, pred_seq1)
    id_12 = calculate_identity(test_seq1, pred_seq2)
    id_21 = calculate_identity(test_seq2, pred_seq1)
    id_22 = calculate_identity(test_seq2, pred_seq2)
    
    # Find best assignment
    if id_11 + id_22 >= id_12 + id_21:
        # Order 1: test1->pred1, test2->pred2
        best_id1 = id_11
        best_id2 = id_22
        assignment = "straight"
    else:
        # Order 2: test1->pred2, test2->pred1
        best_id1 = id_12
        best_id2 = id_21
        assignment = "crossed"
    
    avg_identity = (best_id1 + best_id2) / 2
    
    print(f"\n=== SEQUENCE IDENTITY ===")
    print(f"Assignment: {assignment}")
    print(f"Haplotype 1 identity: {best_id1:.4f}")
    print(f"Haplotype 2 identity: {best_id2:.4f}")
    print(f"Average identity: {avg_identity:.4f}")
    
    # Calculate QV
    qv1 = identity_to_qv(best_id1)
    qv2 = identity_to_qv(best_id2)
    avg_qv = identity_to_qv(avg_identity)
    
    print(f"\n=== QUALITY VALUES (QV) ===")
    print(f"Haplotype 1 QV: {qv1:.1f}")
    print(f"Haplotype 2 QV: {qv2:.1f}")
    print(f"Average QV: {avg_qv:.1f}")
    
    # Interpretation
    print(f"\n=== INTERPRETATION ===")
    if avg_identity > 0.99:
        print("✓ EXCELLENT: Near-perfect genotype assignment (>99% identity)")
    elif avg_identity > 0.95:
        print("✓ GOOD: High-quality genotype assignment (>95% identity)")
    elif avg_identity > 0.90:
        print("⚠ FAIR: Reasonable genotype assignment (>90% identity)")
    else:
        print("✗ POOR: Low-quality genotype assignment (<90% identity)")
        print("  The predicted genotype is significantly different from the actual sequence")

if __name__ == "__main__":
    main()
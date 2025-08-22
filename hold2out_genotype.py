#!/usr/bin/env python3
"""
Proper hold-2-out genotyping following cosigt approach.
Compares test sample coverage against reference combinations.
"""

import numpy as np
import sys

def cosine_similarity(a, b):
    """Calculate cosine similarity between two vectors"""
    dot = np.dot(a, b)
    mag_a = np.linalg.norm(a)
    mag_b = np.linalg.norm(b)
    return dot / (mag_a * mag_b) if mag_a * mag_b > 0 else 0.0

def load_test_coverage(filename):
    """Load coverage from gafpack output"""
    with open(filename, 'r') as f:
        header = f.readline()
        line = f.readline().strip().split('\t')
        # Column 0 is sample name, rest are node coverages
        return np.array([float(x) for x in line[1:]])

def load_reference_coverage(filename):
    """Load reference coverage matrix from odgi paths -H output"""
    ref_matrix = []
    ref_names = []
    
    with open(filename, 'r') as f:
        header = f.readline().strip().split('\t')
        
        # Detect format - if column 2 is "path.length", we need to skip it
        if len(header) > 2 and header[1] == "path.length":
            start_col = 3  # Skip name, length, step.count
        else:
            start_col = 1  # Skip only name
            
        for line in f:
            parts = line.strip().split('\t')
            ref_names.append(parts[0])
            # Extract node coverage values
            coverage = [float(x) for x in parts[start_col:]]
            ref_matrix.append(coverage)
    
    return np.array(ref_matrix), ref_names

def find_best_genotype(test_coverage, ref_matrix, ref_names):
    """Find best matching pair of haplotypes"""
    n_refs = len(ref_names)
    best_score = -1
    best_pair = None
    
    # Test all pairs
    for i in range(n_refs):
        for j in range(i, n_refs):
            # Combine two haplotypes
            combined = ref_matrix[i] + ref_matrix[j]
            score = cosine_similarity(test_coverage, combined)
            
            if score > best_score:
                best_score = score
                best_pair = (ref_names[i], ref_names[j])
    
    return best_pair, best_score

def main():
    if len(sys.argv) != 3:
        print("Usage: hold2out_genotype.py <test_coverage.tsv> <reference_coverage.tsv>")
        sys.exit(1)
    
    test_file = sys.argv[1]
    ref_file = sys.argv[2]
    
    print("Loading test coverage...")
    test_coverage = load_test_coverage(test_file)
    print(f"  Test coverage: {len(test_coverage)} nodes")
    print(f"  Non-zero nodes: {np.count_nonzero(test_coverage)}")
    
    print("\nLoading reference coverage...")
    ref_matrix, ref_names = load_reference_coverage(ref_file)
    print(f"  Reference: {len(ref_names)} paths x {ref_matrix.shape[1]} nodes")
    
    # Ensure dimensions match
    if len(test_coverage) != ref_matrix.shape[1]:
        # Try truncating to smaller dimension
        min_dim = min(len(test_coverage), ref_matrix.shape[1])
        print(f"  Warning: Dimension mismatch, truncating to {min_dim} nodes")
        test_coverage = test_coverage[:min_dim]
        ref_matrix = ref_matrix[:, :min_dim]
    
    print("\nFinding best genotype...")
    best_pair, best_score = find_best_genotype(test_coverage, ref_matrix, ref_names)
    
    print(f"\n=== GENOTYPE RESULT ===")
    print(f"Haplotype 1: {best_pair[0]}")
    print(f"Haplotype 2: {best_pair[1]}")
    print(f"Similarity score: {best_score:.4f}")
    
    # Parse individual names
    ind1 = best_pair[0].split('#')[0] if '#' in best_pair[0] else best_pair[0]
    ind2 = best_pair[1].split('#')[0] if '#' in best_pair[1] else best_pair[1]
    
    if ind1 == ind2:
        print(f"\nGenotype: {ind1} (homozygous)")
    else:
        print(f"\nGenotype: {ind1}/{ind2} (heterozygous)")
    
    # Quality assessment
    if best_score > 0.8:
        print("Quality: HIGH (>0.8)")
    elif best_score > 0.5:
        print("Quality: MEDIUM (0.5-0.8)")
    else:
        print("Quality: LOW (<0.5)")

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
Test hold-0-out genotyping with REAL read coverage.
This should achieve near-perfect accuracy when a sample is included in the reference.
"""

import numpy as np

def cosine_similarity(a, b):
    """Calculate cosine similarity between two vectors"""
    dot = np.dot(a, b)
    mag_a = np.linalg.norm(a)
    mag_b = np.linalg.norm(b)
    if mag_a * mag_b > 0:
        return dot / (mag_a * mag_b)
    return 0.0

def load_coverage_matrix(filename):
    """Load coverage matrix from TSV file"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Parse header
    header = lines[0].strip().split('\t')
    n_nodes = len(header) - 1  # Subtract 'sample' column
    
    # Parse data
    samples = []
    coverage_matrix = []
    
    for line in lines[1:]:
        parts = line.strip().split('\t')
        sample_name = parts[0]
        values = [float(x) for x in parts[1:]]
        samples.append(sample_name)
        coverage_matrix.append(values)
    
    return np.array(coverage_matrix), samples

# Load the reference matrix
print("Loading reference coverage matrix...")
ref_matrix, ref_samples = load_coverage_matrix("test_reference_coverage/test_reference_matrix.tsv")
print(f"Reference matrix shape: {ref_matrix.shape}")
print(f"Samples: {ref_samples}")

# For each sample, combine its haplotypes and test
individuals = {}
for i, sample in enumerate(ref_samples):
    ind_name = sample.split('_')[0]
    if ind_name not in individuals:
        individuals[ind_name] = []
    individuals[ind_name].append(i)

print("\n=== HOLD-0-OUT TEST (Sample in Reference) ===")
for ind_name, indices in individuals.items():
    print(f"\nTesting {ind_name}...")
    
    # Combine haplotypes (sum coverage)
    test_coverage = np.sum(ref_matrix[indices], axis=0)
    
    # Test against all reference combinations
    best_match = None
    best_score = -1
    
    for ref_ind, ref_indices in individuals.items():
        ref_coverage = np.sum(ref_matrix[ref_indices], axis=0)
        score = cosine_similarity(test_coverage, ref_coverage)
        
        print(f"  vs {ref_ind}: {score:.6f}")
        
        if score > best_score:
            best_score = score
            best_match = ref_ind
    
    correct = best_match == ind_name
    print(f"  Best match: {best_match} (score: {best_score:.6f}) {'✓' if correct else '✗'}")

# Also test with just individual haplotypes
print("\n=== INDIVIDUAL HAPLOTYPE TEST ===")
for i, sample in enumerate(ref_samples):
    print(f"\nTesting {sample}...")
    test_coverage = ref_matrix[i]
    
    best_match = None
    best_score = -1
    
    for j, ref_sample in enumerate(ref_samples):
        ref_coverage = ref_matrix[j]
        score = cosine_similarity(test_coverage, ref_coverage)
        
        if j < 3 or abs(j - i) <= 1:  # Show nearby samples
            print(f"  vs {ref_sample}: {score:.6f}")
        
        if score > best_score:
            best_score = score
            best_match = ref_sample
    
    correct = best_match == sample
    print(f"  Best match: {best_match} (score: {best_score:.6f}) {'✓' if correct else '✗'}")
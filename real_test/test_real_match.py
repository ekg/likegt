#!/usr/bin/env python3
import numpy as np

def cosine_similarity(a, b):
    """Calculate cosine similarity between two vectors"""
    dot = np.dot(a, b)
    mag_a = np.linalg.norm(a)
    mag_b = np.linalg.norm(b)
    if mag_a * mag_b > 0:
        return dot / (mag_a * mag_b)
    return 0.0

def load_coverage(filename):
    """Load coverage from TSV file"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Skip header line
    header = lines[0].strip().split('\t')
    n_nodes = len(header) - 1  # Subtract 'sample' column
    
    # Parse data lines (should be just one for gafpack output)
    coverage = []
    for line in lines[1:]:
        parts = line.strip().split('\t')
        # Skip the sample name, get the coverage values
        values = [float(x) for x in parts[1:]]
        coverage.extend(values)
    
    return np.array(coverage)

# Load the coverage files
print("Loading HG00096 coverage...")
coverage1 = load_coverage("HG00096.coverage.tsv")
print(f"Coverage vector shape: {coverage1.shape}")

print("\nLoading HG00268 coverage...")
coverage2 = load_coverage("HG00268.coverage.tsv")
print(f"Coverage vector shape: {coverage2.shape}")

# Show some statistics
print(f"\nHG00096 coverage stats:")
print(f"  Min: {coverage1.min()}, Max: {coverage1.max()}, Mean: {coverage1.mean():.2f}")
print(f"  Non-zero nodes: {np.count_nonzero(coverage1)}/{len(coverage1)}")

print(f"\nHG00268 coverage stats:")
print(f"  Min: {coverage2.min()}, Max: {coverage2.max()}, Mean: {coverage2.mean():.2f}")
print(f"  Non-zero nodes: {np.count_nonzero(coverage2)}/{len(coverage2)}")

# Calculate similarities
sim_self = cosine_similarity(coverage1, coverage1)
sim_cross = cosine_similarity(coverage1, coverage2)

print(f"\n=== COSINE SIMILARITY (WITH ACTUAL COVERAGE) ===")
print(f"HG00096 vs HG00096 (self): {sim_self:.6f}")
print(f"HG00096 vs HG00268: {sim_cross:.6f}")

# Test with binary presence/absence (like the old reference)
binary1 = (coverage1 > 0).astype(float)
binary2 = (coverage2 > 0).astype(float)

sim_binary_self = cosine_similarity(binary1, binary1)
sim_binary_cross = cosine_similarity(binary1, binary2)

print(f"\n=== BINARY (PRESENCE/ABSENCE) SIMILARITY ===")
print(f"HG00096 vs HG00096 (self): {sim_binary_self:.6f}")
print(f"HG00096 vs HG00268: {sim_binary_cross:.6f}")

# Show overlap statistics
overlap = np.logical_and(coverage1 > 0, coverage2 > 0).sum()
only1 = np.logical_and(coverage1 > 0, coverage2 == 0).sum()
only2 = np.logical_and(coverage1 == 0, coverage2 > 0).sum()
neither = np.logical_and(coverage1 == 0, coverage2 == 0).sum()

print(f"\n=== NODE OVERLAP STATS ===")
print(f"Nodes covered by both: {overlap}")
print(f"Nodes only in HG00096: {only1}")
print(f"Nodes only in HG00268: {only2}")
print(f"Nodes in neither: {neither}")
print(f"Total nodes: {len(coverage1)}")
#!/usr/bin/env python3
"""
Demonstrate reference bias with a simple synthetic example.
Shows how pre-aligning to a single reference corrupts coverage.
"""

import numpy as np

def cosine_similarity(a, b):
    """Calculate cosine similarity"""
    dot = np.dot(a, b)
    mag_a = np.linalg.norm(a)
    mag_b = np.linalg.norm(b)
    if mag_a * mag_b > 0:
        return dot / (mag_a * mag_b)
    return 0.0

print("=== REFERENCE BIAS DEMONSTRATION ===")
print("Simulating how pre-alignment to single reference corrupts coverage\n")

# Simulate a pangenome with 100 nodes
n_nodes = 100

# Create coverage for sample HG00096 (covers nodes specific to its haplotype)
# Let's say it covers nodes 0-39 and 60-79 (60 nodes total)
hg00096_true = np.zeros(n_nodes)
hg00096_true[0:40] = 30  # 30x coverage on these nodes
hg00096_true[60:80] = 30  # 30x coverage on these nodes

# Create coverage for reference (e.g., GRCh38)
# Reference only has nodes 0-29 and 70-89 (50 nodes, partially overlapping)
grch38_nodes = np.zeros(n_nodes)
grch38_nodes[0:30] = 1
grch38_nodes[70:90] = 1

print("True HG00096 coverage:")
print(f"  Covers nodes: {np.where(hg00096_true > 0)[0].tolist()[:10]}... (60 nodes total)")
print(f"  Total coverage: {hg00096_true.sum():.0f}")

print("\nGRCh38 reference:")
print(f"  Has nodes: {np.where(grch38_nodes > 0)[0].tolist()[:10]}... (50 nodes total)")

# Simulate biased coverage: reads can only map to nodes present in reference
hg00096_biased = hg00096_true * grch38_nodes

print("\n=== AFTER PRE-ALIGNING TO GRCh38 ===")
print("Biased HG00096 coverage (only nodes in GRCh38):")
print(f"  Covers nodes: {np.where(hg00096_biased > 0)[0].tolist()[:10]}... ({np.count_nonzero(hg00096_biased)} nodes)")
print(f"  Total coverage: {hg00096_biased.sum():.0f}")

# Calculate loss
lost_nodes = (hg00096_true > 0) & (hg00096_biased == 0)
coverage_loss = (hg00096_true.sum() - hg00096_biased.sum()) / hg00096_true.sum() * 100

print("\n=== IMPACT OF REFERENCE BIAS ===")
print(f"Nodes lost due to bias: {lost_nodes.sum()} / {np.count_nonzero(hg00096_true)}")
print(f"Coverage loss: {coverage_loss:.1f}%")

# Calculate similarity
similarity = cosine_similarity(hg00096_true, hg00096_biased)
print(f"Cosine similarity (true vs biased): {similarity:.4f}")

# Show genotyping impact
print("\n=== GENOTYPING IMPACT ===")

# Create another sample (HG00268) with different coverage
hg00268_true = np.zeros(n_nodes)
hg00268_true[20:50] = 30  # Different pattern
hg00268_true[70:90] = 30

# Calculate similarities for genotyping
print("\nTrue genotyping (unbiased):")
self_sim_true = cosine_similarity(hg00096_true, hg00096_true)
cross_sim_true = cosine_similarity(hg00096_true, hg00268_true)
print(f"  HG00096 vs itself: {self_sim_true:.4f}")
print(f"  HG00096 vs HG00268: {cross_sim_true:.4f}")
print(f"  Discrimination ratio: {self_sim_true / cross_sim_true:.2f}x")

print("\nBiased genotyping (after GRCh38 pre-alignment):")
# Bias HG00268 too
hg00268_biased = hg00268_true * grch38_nodes

self_sim_biased = cosine_similarity(hg00096_biased, hg00096_biased)
cross_sim_biased = cosine_similarity(hg00096_biased, hg00268_biased)
print(f"  HG00096 vs itself: {self_sim_biased:.4f}")
print(f"  HG00096 vs HG00268: {cross_sim_biased:.4f}")

if cross_sim_biased > 0:
    print(f"  Discrimination ratio: {self_sim_biased / cross_sim_biased:.2f}x")
else:
    print(f"  Discrimination ratio: undefined (no overlap!)")

print("\n=== CONCLUSION ===")
if coverage_loss > 20:
    print(f"⚠️  SEVERE BIAS: {coverage_loss:.0f}% of coverage lost!")
    print("   Pre-aligning to single reference destroys genotyping accuracy")
elif coverage_loss > 5:
    print(f"⚠️  MODERATE BIAS: {coverage_loss:.0f}% of coverage lost")
    print("   This will reduce genotyping accuracy")
else:
    print(f"✓  Minimal bias: {coverage_loss:.0f}% coverage lost")

print("\nThis demonstrates why we must map directly to the pangenome,")
print("not pre-align to a single reference first!")
#!/usr/bin/env python3
import gzip
import csv
import numpy as np
from collections import defaultdict

def investigate_coverage_similarity():
    """Investigate why some haplotypes have identical or near-identical coverage"""
    
    print("Loading coverage data...")
    coverages = {}
    
    with gzip.open('tests/data/hla-f.k51.paths.coverage.tsv.gz', 'rt') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        
        for row in reader:
            path_name = row[0]
            coverage = [int(x) for x in row[1:]]
            coverages[path_name] = np.array(coverage)
    
    print(f"Loaded {len(coverages)} haplotypes")
    
    # Find pairs with identical coverage
    print("\nFinding haplotypes with identical or near-identical coverage...")
    
    identical_pairs = []
    near_identical_pairs = []
    
    haplotypes = list(coverages.keys())
    
    for i in range(len(haplotypes)):
        for j in range(i+1, len(haplotypes)):
            hap1, hap2 = haplotypes[i], haplotypes[j]
            cov1, cov2 = coverages[hap1], coverages[hap2]
            
            # Check if identical
            if np.array_equal(cov1, cov2):
                identical_pairs.append((hap1, hap2))
            # Check if very similar (>99.9% similar)
            elif np.sum(cov1 == cov2) / len(cov1) > 0.999:
                near_identical_pairs.append((hap1, hap2, np.sum(cov1 == cov2) / len(cov1)))
    
    if identical_pairs:
        print(f"\nFound {len(identical_pairs)} pairs with IDENTICAL coverage:")
        for hap1, hap2 in identical_pairs[:10]:  # Show first 10
            print(f"  {hap1[:50]} == {hap2[:50]}")
    
    if near_identical_pairs:
        print(f"\nFound {len(near_identical_pairs)} pairs with >99.9% similar coverage:")
        for hap1, hap2, sim in near_identical_pairs[:10]:  # Show first 10
            print(f"  {hap1[:30]} ~ {hap2[:30]} ({sim:.4f})")
    
    # Analyze the failed cases from our test
    print("\n=== Analyzing Failed Genotype Cases ===")
    
    failed_cases = [
        ("HG01890#1", "HG01890#2", "HG01505#2", "HG01890#1"),
        ("HG00733#1", "HG00733#2", "HG00732#1", "HG00733#2"),
        ("HG03065#1", "HG03065#2", "HG03065#1", "HG03065#1"),
    ]
    
    for true_hap1, true_hap2, called_hap1, called_hap2 in failed_cases:
        # Find the full names
        true1 = [k for k in coverages.keys() if k.startswith(true_hap1)][0]
        true2 = [k for k in coverages.keys() if k.startswith(true_hap2)][0]
        called1 = [k for k in coverages.keys() if k.startswith(called_hap1)][0]
        called2 = [k for k in coverages.keys() if k.startswith(called_hap2)][0]
        
        print(f"\nCase: {true_hap1} + {true_hap2}")
        print(f"  Called: {called_hap1} + {called_hap2}")
        
        # Compare the coverage patterns
        true_sum = coverages[true1] + coverages[true2]
        called_sum = coverages[called1] + coverages[called2]
        
        # Calculate cosine similarity
        dot = np.dot(true_sum, called_sum)
        mag1 = np.linalg.norm(true_sum)
        mag2 = np.linalg.norm(called_sum)
        
        if mag1 > 0 and mag2 > 0:
            cosine_sim = dot / (mag1 * mag2)
            print(f"  Cosine similarity: {cosine_sim:.10f}")
        
        # Check if they're identical
        if np.array_equal(true_sum, called_sum):
            print(f"  ⚠️  Coverage sums are IDENTICAL!")
            
            # Check individual haplotypes
            if np.array_equal(coverages[true1], coverages[called1]):
                print(f"    {true_hap1} has identical coverage to {called_hap1}")
            if np.array_equal(coverages[true2], coverages[called2]):
                print(f"    {true_hap2} has identical coverage to {called_hap2}")
    
    # Calculate statistics
    print("\n=== Coverage Statistics ===")
    
    # How many nodes does each haplotype cover?
    coverage_counts = {}
    for name, cov in coverages.items():
        coverage_counts[name] = np.sum(cov > 0)
    
    counts = list(coverage_counts.values())
    print(f"Nodes covered per haplotype:")
    print(f"  Mean: {np.mean(counts):.1f}")
    print(f"  Min: {np.min(counts)}")
    print(f"  Max: {np.max(counts)}")
    print(f"  Std: {np.std(counts):.1f}")
    
    # Find haplotypes that cover exactly the same nodes
    coverage_patterns = defaultdict(list)
    for name, cov in coverages.items():
        pattern = tuple(cov > 0)  # Binary pattern of coverage
        coverage_patterns[pattern].append(name)
    
    identical_coverage_groups = [group for group in coverage_patterns.values() if len(group) > 1]
    
    if identical_coverage_groups:
        print(f"\n{len(identical_coverage_groups)} groups of haplotypes cover exactly the same nodes:")
        for group in identical_coverage_groups[:5]:  # Show first 5
            print(f"  Group with {len(group)} haplotypes:")
            for hap in group[:3]:  # Show first 3 in group
                print(f"    - {hap[:60]}")

if __name__ == "__main__":
    investigate_coverage_similarity()
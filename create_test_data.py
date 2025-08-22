#!/usr/bin/env python3
import gzip
import csv
import random
import os

# Create synthetic test data for genotyping
def create_test_data():
    """Create synthetic coverage data for testing"""
    
    # Number of nodes in the graph
    num_nodes = 1000
    # Number of haplotypes
    num_haplotypes = 10
    
    # Create reference haplotype coverage
    print("Creating reference haplotype coverage...")
    with gzip.open('test_paths.tsv.gz', 'wt') as f:
        writer = csv.writer(f, delimiter='\t')
        # Write header
        header = ['path'] + [f'node_{i}' for i in range(num_nodes)]
        writer.writerow(header)
        
        # Write haplotype coverage (random values)
        for hap_id in range(num_haplotypes):
            row = [f'haplotype_{hap_id:03d}']
            # Create distinct patterns for each haplotype
            random.seed(hap_id)  # Reproducible patterns
            coverage = []
            for i in range(num_nodes):
                # Create some pattern - some nodes have high coverage, others low
                if random.random() < 0.3:  # 30% of nodes have coverage
                    coverage.append(random.uniform(5, 50))
                else:
                    coverage.append(0.0)
            row.extend(coverage)
            writer.writerow(row)
    
    # Create sample coverage (mix of two haplotypes for diploid)
    print("Creating sample coverage...")
    hap1_idx = 2  # Use haplotype_002
    hap2_idx = 5  # Use haplotype_005
    
    with gzip.open('test_sample.tsv.gz', 'wt') as f:
        writer = csv.writer(f, delimiter='\t')
        # Write header
        header = ['sample'] + [f'node_{i}' for i in range(num_nodes)]
        writer.writerow(header)
        
        # Generate sample coverage as sum of two haplotypes with some noise
        row = ['test_sample']
        sample_coverage = []
        
        # Read back the reference coverages
        ref_coverages = []
        with gzip.open('test_paths.tsv.gz', 'rt') as ref_f:
            reader = csv.reader(ref_f, delimiter='\t')
            next(reader)  # Skip header
            for ref_row in reader:
                ref_coverages.append([float(x) for x in ref_row[1:]])
        
        # Create sample as mix of hap1 and hap2
        for i in range(num_nodes):
            cov = ref_coverages[hap1_idx][i] + ref_coverages[hap2_idx][i]
            # Add some noise
            if cov > 0:
                cov += random.gauss(0, 2)
                cov = max(0, cov)
            sample_coverage.append(cov)
        
        row.extend(sample_coverage)
        writer.writerow(row)
    
    print(f"Created test data:")
    print(f"  - test_paths.tsv.gz: {num_haplotypes} reference haplotypes")
    print(f"  - test_sample.tsv.gz: Sample (mix of haplotype_002 + haplotype_005)")
    print(f"  - Each with {num_nodes} nodes")

if __name__ == "__main__":
    create_test_data()
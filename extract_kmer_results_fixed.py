#!/usr/bin/env python3
import re
import sys

def extract_results(filename, k_value):
    """Extract valid TSV results from k-mer validation output files"""
    results = []
    header_found = False
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            # Skip debug and error lines
            if line.startswith('DEBUG:') or line.startswith('Error processing') or 'Blocking waiting' in line:
                continue
            
            # Look for header
            if line.startswith('sample\t'):
                header_found = True
                results.append(f"k_value\t{line}")
                continue
            
            # Extract actual sample results (starts with sample ID, not warnings/messages)
            if header_found and '\t' in line:
                # Check if it looks like a valid result line (has expected number of columns)
                cols = line.split('\t')
                if len(cols) >= 10:  # Expected number of columns
                    # Check if first column looks like a sample ID
                    if cols[0] and not cols[0].startswith('warning:') and not 'Finished' in cols[0]:
                        results.append(f"{k_value}\t{line}")
    
    return results

if __name__ == "__main__":
    k_values = [0, 25, 51, 101, 179, 311]
    all_results = []
    
    # Add header
    all_results.append("k_value\tsample\ttrue_hap1\ttrue_hap2\tcalled_hap1\tcalled_hap2\tsimilarity\tqv\talignment\tbias_loss\ttime")
    
    for k in k_values:
        if k == 51:
            # Use existing k51 results
            filename = "hla_f_final.tsv"
        elif k == 25:
            # Use original k25 results
            filename = "hla_f_k25_results.tsv"
        else:
            # Use fixed results
            filename = f"hla_f_k{k}_results_fixed.tsv"
        
        try:
            results = extract_results(filename, k)
            # Skip header from individual files
            if len(results) > 1:
                all_results.extend(results[1:])
            print(f"Extracted {len(results)-1 if len(results) > 1 else 0} results for k={k}")
        except FileNotFoundError:
            print(f"File not found: {filename}")
    
    # Write combined results
    with open("hla_f_all_kmers_complete.tsv", 'w') as f:
        f.write('\n'.join(all_results) + '\n')
    
    print(f"Combined {len(all_results)-1} total results")
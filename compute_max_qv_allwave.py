#!/usr/bin/env python3
"""
Compute maximum attainable QV using allwave alignments.
For each held-out individual, finds the best non-self alignment.
"""

import sys
import subprocess
import argparse
from collections import defaultdict
import math
import gzip

def parse_paf_line(line):
    """Parse a PAF format line and extract key alignment metrics."""
    fields = line.strip().split('\t')
    if len(fields) < 12:
        return None
    
    return {
        'query': fields[0],
        'query_len': int(fields[1]),
        'query_start': int(fields[2]),
        'query_end': int(fields[3]),
        'strand': fields[4],
        'target': fields[5],
        'target_len': int(fields[6]),
        'target_start': int(fields[7]),
        'target_end': int(fields[8]),
        'matches': int(fields[9]),
        'alignment_len': int(fields[10]),
        'mapping_quality': int(fields[11]),
        # Parse additional tags
        'tags': {tag.split(':')[0]: tag.split(':')[2] for tag in fields[12:] if ':' in tag}
    }

def compute_qv_from_alignment(alignment):
    """
    Compute QV from alignment statistics.
    QV = -10 * log10(error_rate)
    """
    if alignment['alignment_len'] == 0:
        return 0.0
    
    # Identity = matches / alignment_length
    identity = alignment['matches'] / alignment['alignment_len']
    
    # Error rate = 1 - identity
    error_rate = 1.0 - identity
    
    if error_rate <= 0:
        return 60.0  # Cap at QV 60 for perfect matches
    
    # QV = -10 * log10(error_rate)
    qv = -10.0 * math.log10(error_rate)
    
    # Cap between 0 and 60
    return min(60.0, max(0.0, qv))

def extract_individual_id(sequence_name):
    """Extract individual ID from sequence name (e.g., HG00096#1 -> HG00096)."""
    if '#' in sequence_name:
        return sequence_name.split('#')[0]
    return sequence_name

def run_allwave_all(fasta_file, threads=4):
    """
    Run allwave for all-vs-all alignment.
    Returns PAF output as string.
    """
    cmd = [
        'allwave',
        '-i', fasta_file,
        '-t', str(threads)
    ]
    
    print(f"Running: {' '.join(cmd)}", file=sys.stderr)
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running allwave: {e}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        return None

def find_best_genotype_pair(alignments, query_hap1, query_hap2):
    """
    Find the best pair of target haplotypes for a diploid genotype.
    Tests both possible pairings and returns the best.
    """
    # Group alignments by query
    hap1_alignments = [a for a in alignments if a['query'] == query_hap1]
    hap2_alignments = [a for a in alignments if a['query'] == query_hap2]
    
    if not hap1_alignments or not hap2_alignments:
        return None
    
    # Sort by QV (best first)
    hap1_alignments.sort(key=lambda x: x['qv'], reverse=True)
    hap2_alignments.sort(key=lambda x: x['qv'], reverse=True)
    
    # Find best pairing
    best_pairing = None
    best_avg_qv = 0
    
    # Try top N targets for each haplotype
    n_candidates = min(10, len(hap1_alignments), len(hap2_alignments))
    
    for i in range(n_candidates):
        for j in range(n_candidates):
            target1 = hap1_alignments[i]['target']
            target2 = hap2_alignments[j]['target']
            
            # Skip if same target used twice (for diploid)
            if target1 == target2:
                continue
            
            # Average QV for this pairing
            avg_qv = (hap1_alignments[i]['qv'] + hap2_alignments[j]['qv']) / 2.0
            
            if avg_qv > best_avg_qv:
                best_avg_qv = avg_qv
                best_pairing = {
                    'target_hap1': target1,
                    'target_hap2': target2,
                    'qv_hap1': hap1_alignments[i]['qv'],
                    'qv_hap2': hap2_alignments[j]['qv'],
                    'avg_qv': avg_qv,
                    'identity_hap1': hap1_alignments[i]['identity'],
                    'identity_hap2': hap2_alignments[j]['identity']
                }
    
    return best_pairing

def process_individual(fasta_file, individual, paf_output=None, threads=4):
    """
    Process a single individual to find max attainable QV.
    """
    # If no PAF provided, run allwave
    if paf_output is None:
        paf_output = run_allwave_all(fasta_file, threads)
        if not paf_output:
            return None
    
    # Parse alignments, filtering out self-alignments
    alignments = []
    for line in paf_output.strip().split('\n'):
        if not line:
            continue
        alignment = parse_paf_line(line)
        if alignment:
            # Skip self-alignments (same individual)
            query_individual = extract_individual_id(alignment['query'])
            target_individual = extract_individual_id(alignment['target'])
            
            # Only keep alignments where query is from our individual
            # but target is from a different individual
            if query_individual == individual and target_individual != individual:
                alignment['qv'] = compute_qv_from_alignment(alignment)
                alignment['identity'] = alignment['matches'] / alignment['alignment_len'] if alignment['alignment_len'] > 0 else 0
                alignments.append(alignment)
    
    # Find alignments for this individual's haplotypes
    hap1 = f"{individual}#1"
    hap2 = f"{individual}#2"
    
    # Find best genotype pairing
    best_genotype = find_best_genotype_pair(alignments, hap1, hap2)
    
    return best_genotype

def get_all_individuals(fasta_file):
    """Extract all unique individual IDs from FASTA file."""
    individuals = set()
    
    if fasta_file.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    with opener(fasta_file, mode) as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line[1:].strip().split()[0]
                individual = extract_individual_id(seq_id)
                if individual and not individual.startswith('grch') and not individual.startswith('chm'):
                    individuals.add(individual)
    
    return sorted(individuals)

def main():
    parser = argparse.ArgumentParser(description='Compute max attainable QV using allwave')
    parser.add_argument('-f', '--fasta', required=True, help='Input FASTA file')
    parser.add_argument('-i', '--individual', help='Individual to test (default: all)')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('-o', '--output', help='Output file (default: stdout)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Determine which individuals to process
    if args.individual and args.individual != 'all':
        individuals = [args.individual]
    else:
        if args.verbose:
            print("Extracting individual IDs from FASTA...", file=sys.stderr)
        individuals = get_all_individuals(args.fasta)
        if args.verbose:
            print(f"Found {len(individuals)} individuals to process", file=sys.stderr)
    
    # Run allwave once for all alignments
    if args.verbose:
        print(f"\nRunning allwave for all-vs-all alignment...", file=sys.stderr)
    paf_output = run_allwave_all(args.fasta, args.threads)
    
    if not paf_output:
        print("Error: allwave failed to produce output", file=sys.stderr)
        return 1
    
    if args.verbose:
        n_lines = len(paf_output.strip().split('\n'))
        print(f"Got {n_lines} alignment records", file=sys.stderr)
    
    # Process each individual
    results = []
    for i, individual in enumerate(individuals, 1):
        if args.verbose:
            print(f"\n[{i}/{len(individuals)}] Processing {individual}...", file=sys.stderr)
        
        result = process_individual(args.fasta, individual, paf_output=paf_output, threads=args.threads)
        
        if result:
            results.append({
                'individual': individual,
                **result
            })
            
            if args.verbose:
                print(f"  Best match: {result['target_hap1']} + {result['target_hap2']}", file=sys.stderr)
                print(f"  Max attainable QV: {result['avg_qv']:.1f}", file=sys.stderr)
        else:
            if args.verbose:
                print(f"  No valid alignments found for {individual}", file=sys.stderr)
    
    # Output results
    output = sys.stdout
    if args.output:
        output = open(args.output, 'w')
    
    # Header
    print("individual\ttarget_hap1\ttarget_hap2\tqv_hap1\tqv_hap2\tavg_qv\tidentity_hap1\tidentity_hap2", file=output)
    
    # Sort by average QV (descending)
    results.sort(key=lambda x: x['avg_qv'], reverse=True)
    
    # Results
    for r in results:
        print(f"{r['individual']}\t{r['target_hap1']}\t{r['target_hap2']}\t"
              f"{r['qv_hap1']:.1f}\t{r['qv_hap2']:.1f}\t{r['avg_qv']:.1f}\t"
              f"{r['identity_hap1']:.4f}\t{r['identity_hap2']:.4f}", file=output)
    
    if args.output:
        output.close()
    
    # Summary
    if args.verbose and results:
        avg_max_qv = sum(r['avg_qv'] for r in results) / len(results)
        print(f"\n=== SUMMARY ===", file=sys.stderr)
        print(f"Processed {len(results)} individuals", file=sys.stderr)
        print(f"Average max attainable QV: {avg_max_qv:.1f}", file=sys.stderr)
        print(f"QV range: {min(r['avg_qv'] for r in results):.1f} - {max(r['avg_qv'] for r in results):.1f}", file=sys.stderr)

if __name__ == '__main__':
    main()
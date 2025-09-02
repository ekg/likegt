#!/usr/bin/env python3
"""
Simplified version to compute max attainable QV using allwave.
"""

import sys
import subprocess
import argparse
import math
import gzip
from collections import defaultdict

def parse_paf_line(line):
    """Parse PAF line."""
    fields = line.strip().split('\t')
    if len(fields) < 12:
        return None
    return {
        'query': fields[0],
        'query_len': int(fields[1]),
        'target': fields[5],
        'target_len': int(fields[6]),
        'matches': int(fields[9]),
        'align_len': int(fields[10]),
    }

def compute_qv(identity):
    """Compute QV from identity."""
    error_rate = 1.0 - identity
    if error_rate <= 0:
        return 60.0
    qv = -10.0 * math.log10(error_rate)
    return min(60.0, max(0.0, qv))

def get_individual_id(seq_name):
    """Extract individual ID."""
    return seq_name.split('#')[0] if '#' in seq_name else seq_name

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', required=True)
    parser.add_argument('-i', '--individual', required=True)
    parser.add_argument('-t', '--threads', type=int, default=4)
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()
    
    # Run allwave
    if args.verbose:
        print(f"Running allwave...", file=sys.stderr)
    
    cmd = ['allwave', '-i', args.fasta, '-t', str(args.threads)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running allwave: {result.stderr}", file=sys.stderr)
        return 1
    
    # Parse alignments for this individual
    hap1 = f"{args.individual}#1"
    hap2 = f"{args.individual}#2"
    
    hap1_alignments = []
    hap2_alignments = []
    
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        
        aln = parse_paf_line(line)
        if not aln:
            continue
        
        # Check if this alignment involves our individual
        query_id = get_individual_id(aln['query'])
        target_id = get_individual_id(aln['target'])
        
        # Skip self-alignments
        if query_id == args.individual and target_id == args.individual:
            continue
        
        # Calculate identity and QV
        identity = aln['matches'] / aln['align_len'] if aln['align_len'] > 0 else 0
        qv = compute_qv(identity)
        
        # Store alignments for our haplotypes
        if aln['query'].startswith(hap1) and target_id != args.individual:
            hap1_alignments.append((aln['target'], identity, qv))
        elif aln['query'].startswith(hap2) and target_id != args.individual:
            hap2_alignments.append((aln['target'], identity, qv))
        # Also check if our haplotypes are targets (allwave might only output one direction)
        elif aln['target'].startswith(hap1) and query_id != args.individual:
            hap1_alignments.append((aln['query'], identity, qv))
        elif aln['target'].startswith(hap2) and query_id != args.individual:
            hap2_alignments.append((aln['query'], identity, qv))
    
    if args.verbose:
        print(f"Found {len(hap1_alignments)} alignments for {hap1}", file=sys.stderr)
        print(f"Found {len(hap2_alignments)} alignments for {hap2}", file=sys.stderr)
    
    if not hap1_alignments or not hap2_alignments:
        print(f"No valid alignments found for {args.individual}", file=sys.stderr)
        return 1
    
    # Sort by QV
    hap1_alignments.sort(key=lambda x: x[2], reverse=True)
    hap2_alignments.sort(key=lambda x: x[2], reverse=True)
    
    # Find best pairing (avoid using same target twice)
    best_qv = 0
    best_pair = None
    
    for t1, id1, qv1 in hap1_alignments[:10]:  # Check top 10
        for t2, id2, qv2 in hap2_alignments[:10]:
            if t1 != t2:  # Different targets
                avg_qv = (qv1 + qv2) / 2.0
                if avg_qv > best_qv:
                    best_qv = avg_qv
                    best_pair = (t1, t2, qv1, qv2, id1, id2)
    
    if best_pair:
        t1, t2, qv1, qv2, id1, id2 = best_pair
        print(f"\n{args.individual} max attainable QV: {best_qv:.1f}")
        print(f"  {hap1} -> {t1[:50]}... (QV: {qv1:.1f}, identity: {id1:.4f})")
        print(f"  {hap2} -> {t2[:50]}... (QV: {qv2:.1f}, identity: {id2:.4f})")
    else:
        print(f"Could not find valid pairing for {args.individual}")
        return 1
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
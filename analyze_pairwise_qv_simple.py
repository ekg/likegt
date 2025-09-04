#!/usr/bin/env python3
"""
Analyze pairwise QV distribution from allwave alignments.
Shows the full distribution of QV values from all pairwise alignments,
excluding self-alignments.
"""

import sys
import subprocess
import gzip
import math
from pathlib import Path
from collections import defaultdict, Counter
import statistics

def get_individual_id(haplotype_name):
    """Extract individual ID from haplotype name."""
    # Format: HG00096#1#JAHBCB010000008.1#0#3341634
    parts = haplotype_name.split('#')
    if len(parts) >= 2:
        return parts[0]
    return haplotype_name.split('.')[0]

def run_allwave(fasta_path, threads=4):
    """Run allwave and return PAF output."""
    cmd = ['allwave', '-i', fasta_path, '-t', str(threads), '-p', 'none']
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running allwave: {result.stderr}", file=sys.stderr)
        return None
    return result.stdout

def parse_paf_line(line):
    """Parse a PAF line and calculate QV."""
    fields = line.strip().split('\t')
    if len(fields) < 12:
        return None
    
    query = fields[0]
    query_len = int(fields[1])
    query_start = int(fields[2])
    query_end = int(fields[3])
    strand = fields[4]
    target = fields[5]
    target_len = int(fields[6])
    target_start = int(fields[7])
    target_end = int(fields[8])
    matches = int(fields[9])
    align_len = int(fields[10])
    mapq = int(fields[11])
    
    if align_len == 0:
        return None
    
    identity = matches / align_len
    error_rate = 1.0 - identity
    
    if error_rate <= 0.0:
        qv = 60.0
    else:
        qv = -10.0 * math.log10(error_rate)
        qv = min(60.0, max(0.0, qv))
    
    return {
        'query': query,
        'target': target,
        'matches': matches,
        'align_len': align_len,
        'identity': identity,
        'qv': qv,
        'query_individual': get_individual_id(query),
        'target_individual': get_individual_id(target)
    }

def categorize_qv(qv):
    """Categorize QV into bins."""
    if qv <= 17:
        return 'Very Low (≤17)'
    elif qv <= 23:
        return 'Low (17-23)'
    elif qv <= 33:
        return 'Mid (23-33)'
    else:
        return 'High (>33)'

def analyze_pairwise_qv(fasta_path, output_prefix, threads=4):
    """Analyze pairwise QV distribution from a FASTA file."""
    
    print(f"\n{'='*60}")
    print(f"Processing {Path(fasta_path).name}")
    print(f"{'='*60}\n")
    
    # Run allwave
    print("Running allwave...")
    paf_output = run_allwave(fasta_path, threads)
    if not paf_output:
        return None
    
    # Parse alignments
    alignments = []
    self_alignments = []
    
    for line in paf_output.strip().split('\n'):
        if not line:
            continue
        
        alignment = parse_paf_line(line)
        if not alignment:
            continue
        
        # Check if self-alignment
        if alignment['query_individual'] == alignment['target_individual']:
            self_alignments.append(alignment)
        else:
            alignments.append(alignment)
    
    print(f"Total alignments: {len(alignments) + len(self_alignments)}")
    print(f"Non-self alignments: {len(alignments)}")
    print(f"Self alignments: {len(self_alignments)} (excluded)")
    
    if not alignments:
        print("No valid non-self alignments found!")
        return None
    
    # Extract QV values
    qv_values = [a['qv'] for a in alignments]
    qv_values.sort()
    
    # Calculate statistics
    mean_qv = statistics.mean(qv_values)
    median_qv = statistics.median(qv_values)
    stdev_qv = statistics.stdev(qv_values) if len(qv_values) > 1 else 0
    min_qv = min(qv_values)
    max_qv = max(qv_values)
    q25 = qv_values[len(qv_values)//4]
    q75 = qv_values[3*len(qv_values)//4]
    
    print(f"\nPairwise QV Statistics (excluding self):")
    print(f"  Mean QV: {mean_qv:.2f}")
    print(f"  Median QV: {median_qv:.2f}")
    print(f"  Std Dev: {stdev_qv:.2f}")
    print(f"  Min QV: {min_qv:.2f}")
    print(f"  Max QV: {max_qv:.2f}")
    print(f"  Q25: {q25:.2f}")
    print(f"  Q75: {q75:.2f}")
    
    # QV categories
    category_counts = Counter()
    for qv in qv_values:
        category_counts[categorize_qv(qv)] += 1
    
    print(f"\nQV Category Distribution:")
    for cat in ['Very Low (≤17)', 'Low (17-23)', 'Mid (23-33)', 'High (>33)']:
        count = category_counts.get(cat, 0)
        pct = 100 * count / len(qv_values) if len(qv_values) > 0 else 0
        print(f"  {cat}: {count} ({pct:.1f}%)")
    
    # Create histogram
    print(f"\nQV Histogram (bins of 5):")
    bins = list(range(0, 65, 5))
    hist = defaultdict(int)
    
    for qv in qv_values:
        for i in range(len(bins)-1):
            if bins[i] <= qv < bins[i+1]:
                hist[f"{bins[i]}-{bins[i+1]}"] = hist.get(f"{bins[i]}-{bins[i+1]}", 0) + 1
                break
        if qv >= 60:
            hist["60-65"] = hist.get("60-65", 0) + 1
    
    max_count = max(hist.values()) if hist else 1
    for i in range(len(bins)-1):
        bin_label = f"{bins[i]:2d}-{bins[i+1]:2d}"
        count = hist.get(f"{bins[i]}-{bins[i+1]}", 0)
        bar_width = int(50 * count / max_count)
        bar = '#' * bar_width
        print(f"  {bin_label}: {bar} {count}")
    
    # Save detailed results
    output_file = f"{output_prefix}_pairwise_report.txt"
    with open(output_file, 'w') as f:
        f.write(f"Pairwise QV Analysis Report\n")
        f.write(f"{'='*60}\n\n")
        f.write(f"Input: {Path(fasta_path).name}\n")
        f.write(f"Output prefix: {output_prefix}\n\n")
        
        f.write(f"Alignment Summary:\n")
        f.write(f"  Total alignments: {len(alignments) + len(self_alignments)}\n")
        f.write(f"  Non-self alignments: {len(alignments)}\n")
        f.write(f"  Self alignments: {len(self_alignments)} (excluded)\n\n")
        
        f.write(f"QV Statistics (excluding self):\n")
        f.write(f"  Mean: {mean_qv:.2f}\n")
        f.write(f"  Median: {median_qv:.2f}\n")
        f.write(f"  Std Dev: {stdev_qv:.2f}\n")
        f.write(f"  Min: {min_qv:.2f}\n")
        f.write(f"  Max: {max_qv:.2f}\n")
        f.write(f"  Q25: {q25:.2f}\n")
        f.write(f"  Q75: {q75:.2f}\n\n")
        
        f.write(f"QV Category Distribution:\n")
        for cat in ['Very Low (≤17)', 'Low (17-23)', 'Mid (23-33)', 'High (>33)']:
            count = category_counts.get(cat, 0)
            pct = 100 * count / len(qv_values) if len(qv_values) > 0 else 0
            f.write(f"  {cat}: {count} ({pct:.1f}%)\n")
        
        f.write(f"\n{'='*60}\n")
        f.write(f"Self-alignment verification:\n")
        f.write(f"  Total self-alignments found and excluded: {len(self_alignments)}\n")
        if self_alignments:
            self_qvs = [a['qv'] for a in self_alignments[:5]]
            f.write(f"  Sample self-alignment QVs (first 5): {self_qvs}\n")
            f.write(f"  Mean self-alignment QV: {statistics.mean([a['qv'] for a in self_alignments]):.2f}\n")
        f.write(f"\n")
        
        f.write(f"Note: Self-alignments are automatically excluded from all statistics above.\n")
        f.write(f"This ensures that the analysis reflects true pairwise comparisons between\n")
        f.write(f"different individuals, not self-similarity.\n")
    
    print(f"\nReport saved to: {output_file}")
    
    # Save TSV with all alignments
    tsv_file = f"{output_prefix}_pairwise_alignments.tsv"
    with open(tsv_file, 'w') as f:
        f.write("query\ttarget\tquery_individual\ttarget_individual\tidentity\tqv\tmatches\talign_len\tis_self\n")
        for a in alignments:
            f.write(f"{a['query']}\t{a['target']}\t{a['query_individual']}\t{a['target_individual']}\t")
            f.write(f"{a['identity']:.4f}\t{a['qv']:.2f}\t{a['matches']}\t{a['align_len']}\tFalse\n")
        # Also write self-alignments for verification
        for a in self_alignments:
            f.write(f"{a['query']}\t{a['target']}\t{a['query_individual']}\t{a['target_individual']}\t")
            f.write(f"{a['identity']:.4f}\t{a['qv']:.2f}\t{a['matches']}\t{a['align_len']}\tTrue\n")
    
    print(f"Alignments saved to: {tsv_file}")
    
    return {
        'file': Path(fasta_path).name,
        'mean_qv': mean_qv,
        'median_qv': median_qv,
        'n_alignments': len(alignments),
        'n_self_excluded': len(self_alignments)
    }

def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_pairwise_qv_simple.py <fasta.gz> [output_prefix] [threads]")
        print("\nExample:")
        print("  python analyze_pairwise_qv_simple.py chr6.fasta.gz chr6_analysis 8")
        sys.exit(1)
    
    fasta_path = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else Path(fasta_path).stem
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    
    if not Path(fasta_path).exists():
        print(f"Error: {fasta_path} not found!")
        sys.exit(1)
    
    analyze_pairwise_qv(fasta_path, output_prefix, threads)

if __name__ == "__main__":
    main()
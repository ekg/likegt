#!/usr/bin/env python3
"""
Analyze pairwise QV distribution from allwave alignments.
Shows the full distribution of QV values from all pairwise alignments,
excluding self-alignments.
"""

import sys
import subprocess
import gzip
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import tempfile
import os

def parse_fasta(fasta_path):
    """Parse FASTA file and return sequences."""
    sequences = {}
    current_name = None
    current_seq = []
    
    if fasta_path.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    with opener(fasta_path, mode) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            sequences[current_name] = ''.join(current_seq)
    
    return sequences

def get_individual_id(haplotype_name):
    """Extract individual ID from haplotype name."""
    # Format: HG00096#1#JAHBCB010000008.1#0#3341634
    parts = haplotype_name.split('#')
    if len(parts) >= 2:
        return parts[0]
    return haplotype_name.split('.')[0]

def run_allwave(fasta_path, threads=4):
    """Run allwave and return PAF output."""
    cmd = ['allwave', '-t', str(threads), fasta_path]
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
        qv = -10.0 * np.log10(error_rate)
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

def analyze_pairwise_qv(fasta_path, output_prefix, threads=4):
    """Analyze pairwise QV distribution from a FASTA file."""
    
    print(f"Processing {fasta_path}...")
    
    # Run allwave
    print("Running allwave...")
    paf_output = run_allwave(fasta_path, threads)
    if not paf_output:
        return
    
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
        return
    
    # Create DataFrame
    df = pd.DataFrame(alignments)
    
    # Calculate statistics
    print("\nPairwise QV Statistics (excluding self):")
    print(f"  Mean QV: {df['qv'].mean():.2f}")
    print(f"  Median QV: {df['qv'].median():.2f}")
    print(f"  Std Dev: {df['qv'].std():.2f}")
    print(f"  Min QV: {df['qv'].min():.2f}")
    print(f"  Max QV: {df['qv'].max():.2f}")
    print(f"  Q25: {df['qv'].quantile(0.25):.2f}")
    print(f"  Q75: {df['qv'].quantile(0.75):.2f}")
    
    # QV categories
    df['qv_category'] = pd.cut(df['qv'], 
                                bins=[-np.inf, 17, 23, 33, np.inf],
                                labels=['Very Low (≤17)', 'Low (17-23)', 'Mid (23-33)', 'High (>33)'])
    
    print("\nQV Category Distribution:")
    for cat, count in df['qv_category'].value_counts().items():
        pct = 100 * count / len(df)
        print(f"  {cat}: {count} ({pct:.1f}%)")
    
    # Create visualizations
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'Pairwise QV Distribution: {Path(fasta_path).stem}', fontsize=14, fontweight='bold')
    
    # 1. Histogram of QV values
    ax = axes[0, 0]
    ax.hist(df['qv'], bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    ax.axvline(df['qv'].mean(), color='red', linestyle='--', label=f'Mean: {df["qv"].mean():.1f}')
    ax.axvline(df['qv'].median(), color='green', linestyle='--', label=f'Median: {df["qv"].median():.1f}')
    ax.set_xlabel('Quality Value (QV)')
    ax.set_ylabel('Count')
    ax.set_title('QV Distribution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. Violin plot by category
    ax = axes[0, 1]
    category_colors = {'Very Low (≤17)': '#d62728', 'Low (17-23)': '#ff7f0e', 
                      'Mid (23-33)': '#ffff00', 'High (>33)': '#2ca02c'}
    df_cat = df.groupby('qv_category')['qv'].apply(list).to_dict()
    positions = []
    labels = []
    data_to_plot = []
    colors = []
    for i, (cat, values) in enumerate(df_cat.items()):
        if len(values) > 0:
            data_to_plot.append(values)
            positions.append(i)
            labels.append(f"{cat}\n(n={len(values)})")
            colors.append(category_colors.get(cat, 'gray'))
    
    if data_to_plot:
        parts = ax.violinplot(data_to_plot, positions=positions, widths=0.7, showmeans=True)
        for pc, color in zip(parts['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
    
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('QV')
    ax.set_title('QV by Category')
    ax.grid(True, alpha=0.3, axis='y')
    
    # 3. Cumulative distribution
    ax = axes[1, 0]
    sorted_qv = np.sort(df['qv'])
    cumulative = np.arange(1, len(sorted_qv) + 1) / len(sorted_qv)
    ax.plot(sorted_qv, cumulative, linewidth=2, color='darkblue')
    ax.set_xlabel('Quality Value (QV)')
    ax.set_ylabel('Cumulative Probability')
    ax.set_title('Cumulative Distribution')
    ax.grid(True, alpha=0.3)
    ax.axhline(0.5, color='red', linestyle='--', alpha=0.5)
    ax.axvline(df['qv'].median(), color='red', linestyle='--', alpha=0.5)
    
    # 4. Box plot with outliers
    ax = axes[1, 1]
    bp = ax.boxplot([df['qv']], patch_artist=True)
    bp['boxes'][0].set_facecolor('lightblue')
    ax.set_ylabel('QV')
    ax.set_xticklabels(['All Pairwise\nAlignments'])
    ax.set_title(f'Box Plot (n={len(df)})')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add text with outlier info
    q1, q3 = df['qv'].quantile([0.25, 0.75])
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    outliers = df[(df['qv'] < lower_bound) | (df['qv'] > upper_bound)]
    ax.text(0.5, 0.95, f'Outliers: {len(outliers)} ({100*len(outliers)/len(df):.1f}%)',
            transform=ax.transAxes, ha='center', va='top', fontsize=10)
    
    plt.tight_layout()
    
    # Save plots
    png_file = f"{output_prefix}_pairwise_qv_dist.png"
    pdf_file = f"{output_prefix}_pairwise_qv_dist.pdf"
    
    plt.savefig(png_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.savefig(pdf_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved plots: {png_file}, {pdf_file}")
    
    # Save detailed statistics to CSV
    stats_file = f"{output_prefix}_pairwise_stats.csv"
    df.to_csv(stats_file, index=False)
    print(f"Saved statistics: {stats_file}")
    
    # Create summary report
    summary_file = f"{output_prefix}_pairwise_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"Pairwise QV Analysis Report\n")
        f.write(f"===========================\n\n")
        f.write(f"Input file: {fasta_path}\n")
        f.write(f"Total alignments: {len(alignments) + len(self_alignments)}\n")
        f.write(f"Non-self alignments: {len(alignments)}\n")
        f.write(f"Self alignments: {len(self_alignments)} (excluded)\n\n")
        
        f.write(f"QV Statistics (excluding self):\n")
        f.write(f"  Mean: {df['qv'].mean():.2f}\n")
        f.write(f"  Median: {df['qv'].median():.2f}\n")
        f.write(f"  Std Dev: {df['qv'].std():.2f}\n")
        f.write(f"  Min: {df['qv'].min():.2f}\n")
        f.write(f"  Max: {df['qv'].max():.2f}\n")
        f.write(f"  Q25: {df['qv'].quantile(0.25):.2f}\n")
        f.write(f"  Q75: {df['qv'].quantile(0.75):.2f}\n\n")
        
        f.write(f"QV Category Distribution:\n")
        for cat, count in df['qv_category'].value_counts().items():
            pct = 100 * count / len(df)
            f.write(f"  {cat}: {count} ({pct:.1f}%)\n")
    
    print(f"Saved summary: {summary_file}")
    
    return df

def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_pairwise_qv.py <fasta.gz> [output_prefix] [threads]")
        print("\nExample:")
        print("  python analyze_pairwise_qv.py chr6.fasta.gz chr6_analysis 8")
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
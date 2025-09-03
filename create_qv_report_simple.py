#!/usr/bin/env python3
"""
Create QV distribution report with stacked bar plots for each genomic region
(No pandas dependency version)
"""

import os
import csv
from pathlib import Path

def categorize_qv(qv):
    """Categorize QV into quality bins"""
    if qv <= 17:
        return 'Very Low'
    elif qv <= 23:
        return 'Low'
    elif qv <= 33:
        return 'Mid'
    else:
        return 'High'

def process_region(tsv_file):
    """Process a single region's TSV file"""
    region_name = Path(tsv_file).stem.replace('_max_qv', '')
    
    qv_values = []
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            qv_values.append(float(row['avg_qv']))
    
    # Count categories
    categories = {'Very Low': 0, 'Low': 0, 'Mid': 0, 'High': 0}
    for qv in qv_values:
        cat = categorize_qv(qv)
        categories[cat] += 1
    
    total = len(qv_values)
    percentages = {cat: 100 * count / total for cat, count in categories.items()}
    avg_qv = sum(qv_values) / len(qv_values) if qv_values else 0
    
    return region_name, percentages, avg_qv

def create_ascii_bar(percentage, width=50):
    """Create an ASCII bar for a percentage"""
    filled = int(percentage * width / 100)
    return '█' * filled + '░' * (width - filled)

def create_text_report(results):
    """Create a text-based visualization and summary"""
    report = []
    report.append("=" * 100)
    report.append("QV DISTRIBUTION REPORT - Maximum Attainable Quality Values")
    report.append("=" * 100)
    report.append("")
    
    # Sort by region name
    results.sort(key=lambda x: x[0])
    
    # Create stacked bar visualization in text
    report.append("Stacked Bar Chart (Quality Distribution by Region)")
    report.append("-" * 100)
    report.append("Region                              0%    25%    50%    75%    100%  Avg QV")
    report.append("                                    |------+------+------+------|")
    
    for region, percentages, avg_qv in results:
        # Create stacked representation
        very_low = percentages['Very Low']
        low = percentages['Low']
        mid = percentages['Mid']
        high = percentages['High']
        
        # Build the bar
        bar_width = 50
        vl_width = int(very_low * bar_width / 100)
        l_width = int(low * bar_width / 100)
        m_width = int(mid * bar_width / 100)
        h_width = bar_width - vl_width - l_width - m_width  # Ensure we fill exactly
        
        bar = '▓' * vl_width + '▒' * l_width + '░' * m_width + '█' * h_width
        
        # Format region name to fixed width
        region_display = region[:35].ljust(35)
        
        report.append(f"{region_display} {bar} {avg_qv:5.1f}")
    
    report.append("")
    report.append("Legend: ▓ = Very Low (≤17)  ▒ = Low (17-23)  ░ = Mid (23-33)  █ = High (>33)")
    report.append("")
    
    # Detailed breakdown
    report.append("=" * 100)
    report.append("DETAILED BREAKDOWN BY REGION")
    report.append("=" * 100)
    
    # Sort by average QV for detailed view
    results_sorted = sorted(results, key=lambda x: x[2], reverse=True)
    
    for region, percentages, avg_qv in results_sorted:
        report.append(f"\n{region}")
        report.append("-" * len(region))
        report.append(f"Average QV: {avg_qv:.1f}")
        report.append(f"  Very Low (≤17): {percentages['Very Low']:5.1f}% {create_ascii_bar(percentages['Very Low'], 20)}")
        report.append(f"  Low (17-23):    {percentages['Low']:5.1f}% {create_ascii_bar(percentages['Low'], 20)}")
        report.append(f"  Mid (23-33):    {percentages['Mid']:5.1f}% {create_ascii_bar(percentages['Mid'], 20)}")
        report.append(f"  High (>33):     {percentages['High']:5.1f}% {create_ascii_bar(percentages['High'], 20)}")
    
    # Overall statistics
    report.append("")
    report.append("=" * 100)
    report.append("OVERALL STATISTICS")
    report.append("=" * 100)
    
    all_avg_qvs = [r[2] for r in results]
    mean_qv = sum(all_avg_qvs) / len(all_avg_qvs)
    min_qv = min(all_avg_qvs)
    max_qv = max(all_avg_qvs)
    
    # Calculate std dev manually
    variance = sum((qv - mean_qv) ** 2 for qv in all_avg_qvs) / len(all_avg_qvs)
    std_dev = variance ** 0.5
    
    report.append(f"Number of regions analyzed: {len(results)}")
    report.append(f"Mean QV across all regions: {mean_qv:.1f}")
    report.append(f"Standard Deviation: {std_dev:.1f}")
    report.append(f"Range: {min_qv:.1f} - {max_qv:.1f}")
    report.append("")
    
    # Region rankings
    report.append("Regions ranked by average QV (best to worst):")
    for i, (region, _, avg_qv) in enumerate(results_sorted, 1):
        quality = "High" if avg_qv > 33 else "Mid" if avg_qv > 23 else "Low" if avg_qv > 17 else "Very Low"
        report.append(f"  {i:2}. {region:40} {avg_qv:5.1f} ({quality})")
    
    return "\n".join(report)

def main():
    # Directory containing TSV files
    results_dir = "max_qv_results"
    
    # Find all TSV files
    tsv_files = sorted(Path(results_dir).glob("*_max_qv.tsv"))
    
    if not tsv_files:
        print(f"No TSV files found in {results_dir}")
        return
    
    print(f"Processing {len(tsv_files)} regions...")
    
    # Process each file
    results = []
    for tsv_file in tsv_files:
        region_name, percentages, avg_qv = process_region(tsv_file)
        results.append((region_name, percentages, avg_qv))
        print(f"  Processed: {region_name} (avg QV: {avg_qv:.1f})")
    
    # Create text report
    report = create_text_report(results)
    
    # Save report
    report_file = "qv_distribution_report.txt"
    with open(report_file, 'w') as f:
        f.write(report)
    print(f"\nSaved report to {report_file}")
    
    # Also print to console
    print("\n" + report)

if __name__ == "__main__":
    main()
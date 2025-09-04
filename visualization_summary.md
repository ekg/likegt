# QV Analysis Visualization Summary

## Generated Visualizations

### 1. Max-QV Analysis (Best Attainable QV)
- **Files**: 
  - `qv_report_final_qv_distribution.pdf` / `.png`
  - `qv_report_final_qv_summary.csv`
- **Description**: Shows the distribution of maximum attainable QV values when each individual is held out and genotyped using the best matching non-self haplotypes
- **Key Finding**: Chr16 shows highest mean max-QV (50.2) with true all-vs-all alignments

### 2. Pairwise QV Distribution Analysis
- **Files**: 
  - `pairwise_chr16_pairwise_qv_analysis.pdf` / `.png` - 4-panel comprehensive analysis
  - `pairwise_chr16_pairwise_histogram.pdf` / `.png` - Detailed histogram
  - `pairwise_chr16_pairwise_summary.csv`
- **Description**: Full distribution of all pairwise alignment QV values (excluding self)
- **Key Stats for chr16**:
  - 17,162 non-self alignments (vs 1,277 with sparsification)
  - Mean QV: 30.2
  - 43.1% high quality (QV > 33)
  - 7.1% very low quality (QV ≤ 17)

## Important Fix Applied

Both the Rust code (`src/commands/max_qv.rs`) and Python analysis script now use `allwave -p none` to ensure true all-vs-all alignments. This gives:
- ~13x more alignments
- Proper detection of all self-alignments (130 for 65 individuals)
- More accurate QV estimates

## Visualization Features

### R Scripts Created:
1. **`create_qv_plot.R`**: Creates stacked bar charts showing QV category distribution across regions
2. **`create_pairwise_qv_plot.R`**: Creates comprehensive pairwise QV analysis with:
   - Violin plots by region
   - Density distributions
   - Cumulative distributions
   - Category breakdowns
   - Detailed histograms

Both scripts use ggplot2 and consistent color scheme:
- Green: High QV (>33)
- Yellow: Mid QV (23-33)
- Orange: Low QV (17-23)
- Red: Very Low QV (≤17)
#!/bin/bash

# Regenerate plots from existing TSV files
# Usage: ./regenerate_plots.sh <tsv_dir>

set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 <tsv_directory>"
    echo "Example: $0 likegt_results/"
    echo ""
    echo "This will regenerate plots from existing *_max_qv.tsv files"
    exit 1
fi

TSV_DIR="$1"

echo "========================================="
echo "Regenerating plots from TSV files in $TSV_DIR"
echo "========================================="
echo

# Check if TSV files exist
if ! ls "$TSV_DIR"/*_max_qv.tsv 1> /dev/null 2>&1; then
    echo "❌ No *_max_qv.tsv files found in $TSV_DIR"
    exit 1
fi

echo "Found $(ls -1 $TSV_DIR/*_max_qv.tsv | wc -l) TSV files"
echo

# Generate standard bar plot
echo "Creating standard bar plot..."
Rscript create_qv_plot.R "$TSV_DIR/combined_report" --dir "$TSV_DIR"
echo "✅ Saved: $TSV_DIR/combined_report_qv_distribution.png"

# Generate violin plot
echo "Creating violin plot..."
Rscript create_qv_violin_plot.R "$TSV_DIR/combined_violin" --dir "$TSV_DIR"
echo "✅ Saved: $TSV_DIR/combined_violin_plot.png"

echo
echo "✅ Complete! Plots regenerated in $TSV_DIR/"
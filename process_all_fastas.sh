#!/bin/bash

# Process all FASTA files in a directory and create plots
# Usage: ./process_all_fastas.sh <input_dir> [output_dir] [neighbors]

set -e

# Check arguments
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_dir> [output_dir] [neighbors]"
    echo ""
    echo "Arguments:"
    echo "  input_dir   Directory containing FASTA files"
    echo "  output_dir  Output directory (default: likegt_results)"
    echo "  neighbors   Number of nearest neighbors for tree sparsification (default: 10)"
    echo "              Use 'none' for exact all-vs-all (slower)"
    echo ""
    echo "Examples:"
    echo "  $0 loci/                    # Use 10 neighbors"
    echo "  $0 loci/ loci.n5/ 5         # Use 5 neighbors, output to loci.n5/"  
    echo "  $0 loci/ loci.n20/ 20       # Use 20 neighbors, output to loci.n20/"
    echo "  $0 loci/ exact/ none        # Exact all-vs-all (slow)"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="${2:-likegt_results}"  # Default to likegt_results if not specified
NEIGHBORS="${3:-10}"               # Default to 10 neighbors
THREADS=8

# Set sparsification strategy
if [ "$NEIGHBORS" = "none" ]; then
    SPARSIFICATION="none"
    echo "Using exact all-vs-all alignment (no sparsification)"
else
    SPARSIFICATION="tree:${NEIGHBORS}:0:0"
    echo "Using tree sparsification with $NEIGHBORS nearest neighbors"
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "========================================="
echo "Processing all FASTA files in $INPUT_DIR"
echo "========================================="
echo

# Step 1: Run max-qv on all FASTA files
echo "Step 1: Running max-qv analysis..."
for fasta in "$INPUT_DIR"/*.fasta.gz "$INPUT_DIR"/*.fasta "$INPUT_DIR"/*.fa.gz "$INPUT_DIR"/*.fa; do
    if [ -f "$fasta" ]; then
        basename=$(basename "$fasta" | sed 's/\.\(fasta\|fa\)\(\.gz\)\?$//')
        output_file="$OUTPUT_DIR/${basename}_max_qv.tsv"
        
        echo "Processing: $basename"
        cargo run --release -- max-qv \
            -f "$fasta" \
            -i all \
            -o "$output_file" \
            -t "$THREADS" \
            -p "$SPARSIFICATION" 2>/dev/null
        
        if [ -f "$output_file" ]; then
            echo "  ✅ Saved to $output_file"
        else
            echo "  ❌ Failed"
        fi
    fi
done

echo
echo "Step 2: Creating visualizations..."

# Step 2: Generate combined plot
if ls "$OUTPUT_DIR"/*_max_qv.tsv 1> /dev/null 2>&1; then
    Rscript create_qv_plot.R "$OUTPUT_DIR/combined_report" --dir "$OUTPUT_DIR"
    echo "✅ Combined plot saved to $OUTPUT_DIR/combined_report_qv_distribution.png"
else
    echo "❌ No TSV files found to plot"
    exit 1
fi

echo
echo "========================================="
echo "✅ Complete!"
echo "========================================="
echo
echo "Results in $OUTPUT_DIR/:"
echo "  - Individual TSV files: $(ls -1 $OUTPUT_DIR/*_max_qv.tsv 2>/dev/null | wc -l) files"
echo "  - Combined plot: $OUTPUT_DIR/combined_report_qv_distribution.png"
echo "  - Summary stats: $OUTPUT_DIR/combined_report_qv_summary.csv"
echo
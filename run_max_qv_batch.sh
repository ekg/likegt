#!/bin/bash

# Process all FASTA files in cositest/alleles_chiara_analysis
FASTA_DIR="../cositest/alleles_chiara_analysis"
OUTPUT_DIR="max_qv_results"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "Processing FASTA files from $FASTA_DIR..."
echo "Results will be saved to $OUTPUT_DIR/"
echo

# Process each FASTA file
for fasta in "$FASTA_DIR"/*.fasta.gz; do
    if [ -f "$fasta" ]; then
        # Get basename without path and extension
        basename=$(basename "$fasta" .fasta.gz)
        output_file="$OUTPUT_DIR/${basename}_max_qv.tsv"
        
        echo "Processing: $basename"
        echo "  Input: $fasta"
        echo "  Output: $output_file"
        
        # Run max-qv with all individuals, save to TSV
        cargo run --release -- max-qv \
            -f "$fasta" \
            -i all \
            -o "$output_file" \
            -t 4 \
            -v 2>&1 | grep -E "Found|Running|Got|SUMMARY|Average" | sed 's/^/  /'
        
        # Check if output was created
        if [ -f "$output_file" ]; then
            line_count=$(wc -l < "$output_file")
            echo "  ✅ Completed: $line_count lines written"
        else
            echo "  ❌ Failed to create output file"
        fi
        echo
    fi
done

echo "All files processed. Results in $OUTPUT_DIR/"
ls -lh "$OUTPUT_DIR"/*.tsv 2>/dev/null | tail -10
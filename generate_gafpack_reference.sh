#!/bin/bash
# Generate reference coverage matrix using gafpack for consistency

set -euo pipefail

echo "=== GENERATING REFERENCE WITH GAFPACK ==="
echo "This ensures reference and test use same node numbering"

GRAPH="hla-f.k51.gfa"
REF_FASTA="hla-f.fa.gz"
OUTPUT="reference_gafpack.tsv"

# Extract paths for BWA
odgi build -g $GRAPH -o temp.og
odgi paths -i temp.og -f > all_paths.fa
bwa index all_paths.fa 2>/dev/null

# Process each haplotype
echo "Processing haplotypes..."
first=true

zcat $REF_FASTA | grep "^>" | sed 's/>//' | while read path_name; do
    echo "  $path_name"
    
    # Extract sequence
    zcat $REF_FASTA | awk -v name=">$path_name" '
        /^>/ { if ($0 == name) p=1; else p=0 }
        p && !/^>/ { print }
    ' > temp_seq.fa
    echo ">$path_name" > temp.fa
    cat temp_seq.fa >> temp.fa
    
    # Simulate minimal reads (just enough to mark presence)
    # For path coverage, we just need 1x coverage
    seq_len=$(cat temp_seq.fa | tr -d '\n' | wc -c)
    n_pairs=$((seq_len / 300 + 1))  # Minimal coverage
    
    wgsim -1 150 -2 150 -N $n_pairs temp.fa temp.1.fq temp.2.fq 2>/dev/null
    
    # Map and convert
    bwa mem -t 4 all_paths.fa temp.1.fq temp.2.fq 2>/dev/null | \
        samtools view -b - > temp.bam
    gfainject --gfa $GRAPH --bam temp.bam > temp.gaf 2>/dev/null
    gafpack --gfa $GRAPH --gaf temp.gaf > temp_coverage.tsv
    
    if [ "$first" = true ]; then
        # First iteration - include header
        cat temp_coverage.tsv > $OUTPUT
        first=false
    else
        # Subsequent iterations - append data only
        tail -n +2 temp_coverage.tsv >> $OUTPUT
    fi
    
    # Clean up
    rm -f temp.fa temp_seq.fa temp.1.fq temp.2.fq temp.bam temp.gaf temp_coverage.tsv
done

# Clean up
rm -f temp.og all_paths.fa*

echo "Generated reference matrix: $OUTPUT"
wc -l $OUTPUT
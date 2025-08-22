#!/bin/bash
set -e

echo "=== Building graphs for all k values ==="
echo "This will test k = 0, 25, 51, 101, 179, 311"

# Check if PAF exists, if not create it
if [ ! -f "tests/data/hla-f.paf" ]; then
    echo "Creating PAF alignment..."
    allwave -i tests/data/hla-f.fa.gz -t 8 -p none > tests/data/hla-f.paf
else
    echo "Using existing PAF file"
fi

# Build graphs for each k value
for k in 0 25 51 101 179 311; do
    echo ""
    echo "=== Building graph with k=$k ==="
    
    OUTPUT_DIR="tests/data/k${k}"
    mkdir -p $OUTPUT_DIR
    
    # Check if already built
    if [ -f "$OUTPUT_DIR/hla-f.k${k}.paths.coverage.tsv.gz" ]; then
        echo "Graph for k=$k already exists, skipping..."
        continue
    fi
    
    echo "Running seqwish with k=$k..."
    seqwish -s tests/data/hla-f.fa.gz \
            -g $OUTPUT_DIR/hla-f.seqwish-k${k}.gfa \
            -t 8 \
            -p tests/data/hla-f.paf \
            -k $k \
            -P
    
    echo "Building odgi graph..."
    odgi build -g $OUTPUT_DIR/hla-f.seqwish-k${k}.gfa \
               -o $OUTPUT_DIR/hla-f.k${k}.og
    
    echo "Sorting graph..."
    odgi sort -i $OUTPUT_DIR/hla-f.k${k}.og \
              -p Ygs \
              -o $OUTPUT_DIR/hla-f.k${k}.sorted.og
    
    echo "Converting to GFA..."
    odgi view -i $OUTPUT_DIR/hla-f.k${k}.sorted.og -g > $OUTPUT_DIR/hla-f.k${k}.gfa
    
    echo "Extracting path coverage matrix..."
    odgi paths -i $OUTPUT_DIR/hla-f.k${k}.sorted.og -H | \
         cut -f 1,4- | gzip > $OUTPUT_DIR/hla-f.k${k}.paths.coverage.tsv.gz
    
    # Get statistics
    echo "Graph statistics for k=$k:"
    odgi stats -i $OUTPUT_DIR/hla-f.k${k}.sorted.og -S
    
    # Clean up intermediate files
    rm -f $OUTPUT_DIR/hla-f.seqwish-k${k}.gfa
    
    echo "Completed k=$k"
done

echo ""
echo "=== All graphs built successfully ==="
ls -lh tests/data/k*/hla-f.*.paths.coverage.tsv.gz
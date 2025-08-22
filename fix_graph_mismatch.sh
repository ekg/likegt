#!/bin/bash
# Fix the graph mismatch issue - use consistent graph throughout

set -euo pipefail

echo "=== FIXING GRAPH MISMATCH ==="

# Use the SAME graph for everything
GRAPH_OG="hla-f.k51.og"
GRAPH_GFA="hla-f.k51.gfa"

# 1. Generate reference coverage with ALL nodes (not just visited ones)
echo "1. Generating reference coverage for ALL nodes..."
odgi paths -i $GRAPH_OG -c > reference_all_nodes.tsv

echo "   Reference coverage shape:"
head -1 reference_all_nodes.tsv | wc -w

# 2. Test with a sample using the SAME graph
echo "2. Simulating test reads..."
TEST_IND="HG00096"

# Extract test sequences
zcat hla-f.fa.gz | awk -v ind="$TEST_IND" '
    /^>/ { if ($0 ~ "^>"ind"#") {print_it=1} else {print_it=0} }
    print_it {print}
' > ${TEST_IND}_test.fa

# Simulate reads
seq_len=$(grep -v "^>" ${TEST_IND}_test.fa | tr -d '\n' | wc -c)
n_pairs=$((seq_len * 30 / 300))
wgsim -1 150 -2 150 -N $n_pairs ${TEST_IND}_test.fa \
    ${TEST_IND}_test.1.fq ${TEST_IND}_test.2.fq 2>/dev/null

# 3. Map to the GRAPH directly (using vg or GraphAligner would be better, but we'll use gfainject)
echo "3. Mapping to graph..."

# Extract paths for BWA mapping
odgi paths -i $GRAPH_OG -f > graph_paths.fa
bwa index graph_paths.fa 2>/dev/null

bwa mem -t 4 graph_paths.fa ${TEST_IND}_test.1.fq ${TEST_IND}_test.2.fq 2>/dev/null | \
    samtools view -b - > ${TEST_IND}_test.bam

# Convert to GAF using the SAME graph
gfainject --gfa $GRAPH_GFA --bam ${TEST_IND}_test.bam > ${TEST_IND}_test.gaf 2>/dev/null

# Get coverage using the SAME graph
gafpack --gfa $GRAPH_GFA --gaf ${TEST_IND}_test.gaf > ${TEST_IND}_test_coverage.tsv

echo "4. Checking dimensions..."
ref_nodes=$(head -1 reference_all_nodes.tsv | wc -w)
test_nodes=$(head -1 ${TEST_IND}_test_coverage.tsv | wc -w)

echo "   Reference nodes: $ref_nodes"
echo "   Test nodes: $test_nodes"

if [ "$ref_nodes" -eq "$test_nodes" ]; then
    echo "✓ Dimensions match!"
else
    echo "✗ STILL MISMATCHED"
    echo "   This means odgi paths -c and gafpack are using different node sets"
fi

# Clean up
rm -f ${TEST_IND}_test.fa ${TEST_IND}_test.1.fq ${TEST_IND}_test.2.fq 
rm -f ${TEST_IND}_test.bam ${TEST_IND}_test.gaf
rm -f graph_paths.fa*
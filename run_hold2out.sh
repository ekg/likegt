#!/bin/bash
# Complete hold-2-out genotyping pipeline
# Following the cosigt approach

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: $0 <graph.og> <test_individual> [reference.fa.gz]"
    echo "Example: $0 hla-f.k51.og HG00096 hla-f.fa.gz"
    exit 1
fi

GRAPH="$1"
TEST_IND="$2"
REF_FASTA="${3:-hla-f.fa.gz}"
COVERAGE_DEPTH=30
READ_LENGTH=150

echo "=== HOLD-2-OUT GENOTYPING ==="
echo "Graph: $GRAPH"
echo "Test individual: $TEST_IND"
echo ""

# Step 1: Remove test individual from graph
echo "1. Creating hold-2-out graph..."
odgi paths -i "$GRAPH" -L | grep -v "^${TEST_IND}#" > kept_paths.txt
odgi paths -i "$GRAPH" -K kept_paths.txt -o hold2out_temp.og
odgi prune -i hold2out_temp.og -o hold2out.og
rm hold2out_temp.og

# Step 2: Extract reference coverage matrix
echo "2. Extracting reference coverage..."
odgi paths -i hold2out.og -H > hold2out_reference.tsv

# Step 3: Convert to GFA for gfainject
echo "3. Converting to GFA..."
odgi view -i hold2out.og -g > hold2out.gfa

# Step 4: Extract paths for BWA mapping
echo "4. Preparing BWA index..."
odgi paths -i hold2out.og -f > hold2out_paths.fa
bwa index hold2out_paths.fa 2>/dev/null

# Step 5: Simulate reads from test individual
echo "5. Simulating reads from $TEST_IND..."
# Get sequences for test individual
zcat "$REF_FASTA" | awk -v ind="$TEST_IND" '
    /^>/ { if ($0 ~ "^>"ind"#") {print_it=1} else {print_it=0} }
    print_it {print}
' > ${TEST_IND}.fa

seq_len=$(grep -v "^>" ${TEST_IND}.fa | tr -d '\n' | wc -c)
n_pairs=$((seq_len * COVERAGE_DEPTH / (2 * READ_LENGTH)))

wgsim -1 $READ_LENGTH -2 $READ_LENGTH -N $n_pairs -e 0.01 -r 0.001 \
    ${TEST_IND}.fa ${TEST_IND}.1.fq ${TEST_IND}.2.fq 2>/dev/null

# Step 6: Map reads
echo "6. Mapping reads..."
bwa mem -t 4 hold2out_paths.fa ${TEST_IND}.1.fq ${TEST_IND}.2.fq 2>/dev/null | \
    samtools view -b - > ${TEST_IND}.bam

# Step 7: Convert to GAF
echo "7. Converting to GAF..."
gfainject --gfa hold2out.gfa --bam ${TEST_IND}.bam > ${TEST_IND}.gaf 2>/dev/null

# Step 8: Get coverage
echo "8. Computing coverage..."
gafpack --gfa hold2out.gfa --gaf ${TEST_IND}.gaf > ${TEST_IND}_coverage.tsv

# Step 9: Run genotyping
echo "9. Genotyping..."
echo ""
python3 hold2out_genotype.py ${TEST_IND}_coverage.tsv hold2out_reference.tsv

# Cleanup
rm -f kept_paths.txt ${TEST_IND}.fa ${TEST_IND}.1.fq ${TEST_IND}.2.fq ${TEST_IND}.bam ${TEST_IND}.gaf
rm -f hold2out_paths.fa*

echo ""
echo "Hold-2-out genotyping complete!"
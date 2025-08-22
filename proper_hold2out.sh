#!/bin/bash
set -euo pipefail

# Proper hold-2-out: Remove test individual's paths and prune uncovered nodes

GRAPH="hla-f.k51.og"
TEST_INDIVIDUAL="HG00096"

echo "=== PROPER HOLD-2-OUT TEST ==="
echo "Testing $TEST_INDIVIDUAL with paths removed from graph"
echo ""

# Step 1: Get list of all paths except test individual
echo "1. Creating list of paths to keep (excluding $TEST_INDIVIDUAL)..."
odgi paths -i $GRAPH -L | grep -v "^${TEST_INDIVIDUAL}#" > kept_paths.txt
echo "   Keeping $(wc -l < kept_paths.txt) paths"

# Step 2: Extract subgraph with only kept paths
echo "2. Extracting subgraph without $TEST_INDIVIDUAL..."
odgi paths -i $GRAPH -K kept_paths.txt -o hold2out_temp.og

# Step 3: Prune nodes with 0 coverage (not visited by any remaining path)
echo "3. Pruning nodes with 0 coverage..."
odgi prune -i hold2out_temp.og -o hold2out_pruned.og

# Step 4: Get stats on the pruned graph
echo "4. Graph statistics:"
echo "   Original graph:"
odgi stats -i $GRAPH -S | head -2

echo "   Hold-2-out graph (pruned):"
odgi stats -i hold2out_pruned.og -S | head -2

# Step 5: Convert to GFA for use with gafpack
echo "5. Converting to GFA..."
odgi view -i hold2out_pruned.og -g > hold2out.gfa

# Step 6: Extract coverage matrix for the pruned graph
echo "6. Extracting coverage matrix for pruned graph..."
odgi paths -i hold2out_pruned.og -H > hold2out_coverage.tsv

# Check dimensions
n_nodes_gfa=$(grep "^S" hold2out.gfa | wc -l)
n_cols=$(($(head -1 hold2out_coverage.tsv | wc -w) - 1))

echo ""
echo "=== DIMENSION CHECK ==="
echo "Pruned GFA nodes: $n_nodes_gfa"
echo "Coverage matrix columns: $n_cols"

if [ "$n_nodes_gfa" -eq "$n_cols" ]; then
    echo "✓ Dimensions match!"
else
    echo "⚠️  DIMENSION MISMATCH!"
fi

# Step 7: Now simulate reads from test individual and map to pruned graph
echo ""
echo "7. Simulating reads from $TEST_INDIVIDUAL..."

# Extract test sequences
samtools faidx hla-f.fa.gz \
    $(zcat hla-f.fa.gz | grep "^>${TEST_INDIVIDUAL}#" | sed 's/>//' | head -2) \
    > ${TEST_INDIVIDUAL}.fa 2>/dev/null

# Get sequence length and calculate read pairs
seq_len=$(grep -v "^>" ${TEST_INDIVIDUAL}.fa | tr -d '\n' | wc -c)
n_pairs=$((seq_len * 30 / 300))  # 30x coverage, 150bp paired reads

echo "   Simulating $n_pairs read pairs..."
wgsim -1 150 -2 150 -N $n_pairs ${TEST_INDIVIDUAL}.fa \
    ${TEST_INDIVIDUAL}.1.fq ${TEST_INDIVIDUAL}.2.fq 2>/dev/null

# Extract paths from pruned graph for mapping
echo "8. Mapping reads to hold-2-out graph..."
odgi paths -i hold2out_pruned.og -f > hold2out_paths.fa
bwa index hold2out_paths.fa 2>/dev/null

bwa mem -t 4 hold2out_paths.fa ${TEST_INDIVIDUAL}.1.fq ${TEST_INDIVIDUAL}.2.fq 2>/dev/null | \
    samtools view -b - > ${TEST_INDIVIDUAL}_hold2out.bam

# Convert to GAF and get coverage
echo "9. Computing coverage on pruned graph..."
gfainject --gfa hold2out.gfa --bam ${TEST_INDIVIDUAL}_hold2out.bam > ${TEST_INDIVIDUAL}_hold2out.gaf 2>/dev/null
gafpack --gfa hold2out.gfa --gaf ${TEST_INDIVIDUAL}_hold2out.gaf > ${TEST_INDIVIDUAL}_hold2out_coverage.tsv

# Step 8: Compare dimensions and coverage
echo ""
echo "=== FINAL CHECK ==="
test_cols=$(($(head -1 ${TEST_INDIVIDUAL}_hold2out_coverage.tsv | wc -w) - 1))
echo "Test coverage columns: $test_cols"
echo "Reference coverage columns: $n_cols"

if [ "$test_cols" -eq "$n_cols" ]; then
    echo "✓ DIMENSIONS MATCH! Ready for genotyping"
    
    # Quick similarity check
    python3 << 'EOF'
import numpy as np

def cosine_similarity(a, b):
    dot = np.dot(a, b)
    mag_a = np.linalg.norm(a)
    mag_b = np.linalg.norm(b)
    return dot / (mag_a * mag_b) if mag_a * mag_b > 0 else 0.0

# Load reference (without test individual)
ref_matrix = []
with open('hold2out_coverage.tsv', 'r') as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if not parts[0].startswith('HG00096#'):  # Skip if somehow still there
            ref_matrix.append([float(x) for x in parts[1:]])

if len(ref_matrix) > 0:
    # Sum all reference paths (simplified - should do proper combinations)
    ref_sum = np.sum(ref_matrix, axis=0)
    
    # Load test coverage
    with open('HG00096_hold2out_coverage.tsv', 'r') as f:
        f.readline()  # Skip header
        test_parts = f.readline().strip().split('\t')
        test_coverage = np.array([float(x) for x in test_parts[1:]])
    
    print(f"\nReference shape: {ref_sum.shape}")
    print(f"Test shape: {test_coverage.shape}")
    
    if ref_sum.shape == test_coverage.shape:
        sim = cosine_similarity(ref_sum, test_coverage)
        print(f"Similarity (test vs all ref): {sim:.4f}")
EOF
else
    echo "✗ DIMENSION MISMATCH - Something went wrong!"
fi

echo ""
echo "Hold-2-out test complete. The test individual was properly removed"
echo "from the graph and nodes were pruned before testing."
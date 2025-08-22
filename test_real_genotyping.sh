#!/bin/bash
# Test REAL genotyping with actual read coverage

set -euo pipefail

echo "=== REAL GENOTYPING TEST ==="
echo "Using actual read coverage, not path presence"
echo ""

# First, generate reference coverage for a few individuals
INDIVIDUALS="HG00096 HG00268 HG00733"
COVERAGE=30

# We'll store coverage vectors
mkdir -p real_test
cd real_test

# Extract paths once
odgi paths -i ../hla-f.k51.gfa -f > paths.fa
bwa index paths.fa 2>/dev/null

# Generate reference coverage
echo "Generating reference coverage..."
for IND in $INDIVIDUALS; do
    echo "  $IND..."
    
    # Extract sequences
    samtools faidx ../hla-f.fa.gz \
        $(samtools faidx ../hla-f.fa.gz 2>/dev/null | grep "^${IND}#" | cut -f1) \
        > ${IND}.fa 2>/dev/null
    
    # Simulate reads
    SEQ_LEN=$(grep -v "^>" ${IND}.fa | tr -d '\n' | wc -c)
    READ_PAIRS=$((SEQ_LEN * COVERAGE / 300))
    
    wgsim -1 150 -2 150 -N $READ_PAIRS ${IND}.fa ${IND}.1.fq ${IND}.2.fq 2>/dev/null
    
    # Map and get coverage
    bwa mem -t 4 paths.fa ${IND}.1.fq ${IND}.2.fq 2>/dev/null | \
        samtools view -bS - > ${IND}.bam
    
    gfainject --gfa ../hla-f.k51.gfa --bam ${IND}.bam > ${IND}.gaf 2>/dev/null
    gafpack --gfa ../hla-f.k51.gfa --gaf ${IND}.gaf > ${IND}.coverage.tsv
done

echo ""
echo "Now testing hold-0-out (sample should match itself)..."

# For HG00096, combine its two haplotype coverages
echo "Test: HG00096"

# The coverage we got is already the combined diploid coverage
# Now find best match among all three
python3 << 'EOF'
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

# Read coverage vectors
individuals = ["HG00096", "HG00268", "HG00733"]
coverages = {}

for ind in individuals:
    with open(f"{ind}.coverage.tsv") as f:
        lines = f.readlines()
        if len(lines) >= 2:
            # Skip header, get coverage values
            values = lines[1].strip().split('\t')[1:]
            coverages[ind] = np.array([float(v) for v in values])

# Test hold-0-out for HG00096
test_cov = coverages["HG00096"]

# Compare to all (including itself)
for ind in individuals:
    ref_cov = coverages[ind]
    
    # Ensure same length
    min_len = min(len(test_cov), len(ref_cov))
    sim = cosine_similarity(
        test_cov[:min_len].reshape(1, -1),
        ref_cov[:min_len].reshape(1, -1)
    )[0, 0]
    
    print(f"  {ind}: similarity = {sim:.4f}")

best = max(individuals, key=lambda x: cosine_similarity(
    test_cov[:min(len(test_cov), len(coverages[x]))].reshape(1, -1),
    coverages[x][:min(len(test_cov), len(coverages[x]))].reshape(1, -1)
)[0, 0])

print(f"\nBest match: {best}")
print(f"Correct: {'YES' if best == 'HG00096' else 'NO'}")
EOF

echo ""
echo "This is how genotyping SHOULD work - with real read coverage!"
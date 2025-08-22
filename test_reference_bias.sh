#!/bin/bash
set -euo pipefail

# Test reference bias by comparing:
# 1. Direct mapping to pangenome (unbiased)
# 2. Pre-aligning to single reference then projecting (biased)

GRAPH="../hla-f.k51.gfa"
REF_FASTA="../hla-f.fa.gz"
OUTPUT_DIR="reference_bias_test"
COVERAGE_DEPTH=30
READ_LENGTH=150

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "=== REFERENCE BIAS TEST ==="
echo "Comparing unbiased vs biased mapping for HG00096"
echo ""

# Extract HG00096 sequences
echo "Extracting HG00096 sequences..."
samtools faidx "$REF_FASTA" \
    $(samtools faidx "$REF_FASTA" 2>/dev/null | grep "^HG00096#" | cut -f1) \
    > HG00096.fa 2>/dev/null

# Get sequence length
SEQ_LEN=$(grep -v "^>" HG00096.fa | tr -d '\n' | wc -c)
TOTAL_BASES=$((SEQ_LEN * COVERAGE_DEPTH))
BASES_PER_PAIR=$((2 * READ_LENGTH))
N_PAIRS=$((TOTAL_BASES / BASES_PER_PAIR))

echo "Simulating $N_PAIRS read pairs for ${COVERAGE_DEPTH}x coverage..."
wgsim -1 $READ_LENGTH -2 $READ_LENGTH -N $N_PAIRS -e 0.01 -r 0.001 \
    HG00096.fa HG00096.1.fq HG00096.2.fq 2>/dev/null

# === UNBIASED: Direct mapping to pangenome ===
echo ""
echo "1. UNBIASED: Mapping directly to pangenome..."

# Extract all paths for mapping
odgi paths -i $GRAPH -f > all_paths.fa 2>/dev/null || {
    # If GRAPH is .gfa not .og, use the paths directly
    odgi build -g $GRAPH -o temp.og
    odgi paths -i temp.og -f > all_paths.fa
    rm temp.og
}

bwa index all_paths.fa 2>/dev/null

# Map reads
bwa mem -t 4 all_paths.fa HG00096.1.fq HG00096.2.fq 2>/dev/null | \
    samtools view -b - > HG00096_unbiased.bam

# Convert to GAF and get coverage
gfainject --gfa $GRAPH --bam HG00096_unbiased.bam > HG00096_unbiased.gaf 2>/dev/null
gafpack --gfa $GRAPH --gaf HG00096_unbiased.gaf > HG00096_unbiased.coverage.tsv

# === BIASED: Pre-align to single reference ===
echo "2. BIASED: Pre-aligning to grch38 reference first..."

# Extract GRCh38 as the single reference
samtools faidx "$REF_FASTA" "grch38#1#chr6:29711814-29738528" > grch38.fa 2>/dev/null

# Build index for single reference
bwa index grch38.fa 2>/dev/null

# Map to single reference
bwa mem -t 4 grch38.fa HG00096.1.fq HG00096.2.fq 2>/dev/null | \
    samtools view -b - > HG00096_to_grch38.bam

# Now re-map the aligned reads to pangenome
echo "   Re-mapping grch38-aligned reads to pangenome..."
samtools fastq HG00096_to_grch38.bam 2>/dev/null | \
    awk '{if(NR%4==1) print "@"substr($0,2); else print}' > HG00096_from_grch38.fq

# Split into paired reads (approximate)
awk 'NR%8<4' HG00096_from_grch38.fq > HG00096_from_grch38.1.fq
awk 'NR%8>=4' HG00096_from_grch38.fq > HG00096_from_grch38.2.fq

# Map biased reads to pangenome
bwa mem -t 4 all_paths.fa HG00096_from_grch38.1.fq HG00096_from_grch38.2.fq 2>/dev/null | \
    samtools view -b - > HG00096_biased.bam

# Convert to GAF and get coverage
gfainject --gfa $GRAPH --bam HG00096_biased.bam > HG00096_biased.gaf 2>/dev/null
gafpack --gfa $GRAPH --gaf HG00096_biased.gaf > HG00096_biased.coverage.tsv

# === ANALYSIS ===
echo ""
echo "3. Analyzing coverage differences..."

python3 << 'EOF'
import numpy as np

def load_coverage(filename):
    """Load coverage from TSV file"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Skip header, parse values
    if len(lines) > 1:
        parts = lines[1].strip().split('\t')
        return np.array([float(x) for x in parts[1:]])
    return np.array([])

def cosine_similarity(a, b):
    """Calculate cosine similarity"""
    dot = np.dot(a, b)
    mag_a = np.linalg.norm(a)
    mag_b = np.linalg.norm(b)
    if mag_a * mag_b > 0:
        return dot / (mag_a * mag_b)
    return 0.0

# Load coverage
print("Loading coverage data...")
unbiased = load_coverage("HG00096_unbiased.coverage.tsv")
biased = load_coverage("HG00096_biased.coverage.tsv")

print(f"Unbiased coverage shape: {unbiased.shape}")
print(f"Biased coverage shape: {biased.shape}")

# Compare coverage
print("\n=== COVERAGE STATISTICS ===")
print(f"Unbiased:")
print(f"  Total coverage: {unbiased.sum():.0f}")
print(f"  Mean per node: {unbiased.mean():.1f}")
print(f"  Nodes with coverage: {np.count_nonzero(unbiased)}")

print(f"\nBiased (pre-aligned to grch38):")
print(f"  Total coverage: {biased.sum():.0f}")
print(f"  Mean per node: {biased.mean():.1f}")
print(f"  Nodes with coverage: {np.count_nonzero(biased)}")

# Calculate similarity
similarity = cosine_similarity(unbiased, biased)
print(f"\n=== SIMILARITY ===")
print(f"Cosine similarity between unbiased and biased: {similarity:.4f}")

# Show which nodes lost coverage
lost_coverage = (unbiased > 0) & (biased == 0)
gained_coverage = (unbiased == 0) & (biased > 0)

print(f"\n=== COVERAGE CHANGES ===")
print(f"Nodes that lost all coverage due to bias: {lost_coverage.sum()}")
print(f"Nodes that gained coverage due to bias: {gained_coverage.sum()}")

# Calculate coverage loss percentage
if unbiased.sum() > 0:
    coverage_loss = (1 - biased.sum() / unbiased.sum()) * 100
    print(f"Total coverage loss: {coverage_loss:.1f}%")

# Show impact on genotyping
print("\n=== GENOTYPING IMPACT ===")
if similarity < 0.95:
    print("⚠️  SEVERE BIAS DETECTED!")
    print(f"   The biased approach has similarity of only {similarity:.2f} with unbiased")
    print(f"   This will cause incorrect genotyping results!")
elif similarity < 0.99:
    print("⚠️  MODERATE BIAS DETECTED")
    print(f"   The biased approach has similarity of {similarity:.2f} with unbiased")
    print(f"   This may affect genotyping accuracy")
else:
    print("✓  Minimal bias detected")
    print(f"   Similarity: {similarity:.4f}")
EOF

echo ""
echo "=== CONCLUSION ==="
echo "Reference bias test complete. Check the statistics above to see how"
echo "pre-aligning to a single reference corrupts the coverage distribution."
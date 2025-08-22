#!/bin/bash
set -e

echo "=== Validating hold-0-out genotyping ==="
echo "Testing that we ALWAYS get the correct answer when simulating from graph paths"

# Select random paths to test
PATHS=(
    "HG00096#1#haplotype1-0000024:29677236-29701492"
    "HG00096#2#haplotype2-0000117:29686034-29710088"
    "HG00268#1#haplotype1-0000024:141273481-141297761"
    "NA12329#1#haplotype1-0000031:142215709-142239790"
)

# Test different combinations
echo "Testing path pair: ${PATHS[0]} + ${PATHS[1]}"

# Extract sequences from the graph paths themselves (not from original FASTA)
echo "Extracting sequences from graph paths..."
odgi paths -i hla-f.k51.sorted.og -E | grep -A1 "^>${PATHS[0]}" | tail -1 > path1.seq
odgi paths -i hla-f.k51.sorted.og -E | grep -A1 "^>${PATHS[1]}" | tail -1 > path2.seq

# Convert to proper FASTA
echo ">${PATHS[0]}" > path1.fa
cat path1.seq >> path1.fa
echo ">${PATHS[1]}" > path2.fa
cat path2.seq >> path2.fa

# Simulate reads with high coverage and no errors
echo "Simulating reads (high coverage, no errors)..."
wgsim -N 5000 -1 150 -2 150 -r 0 -R 0 -X 0 -e 0 path1.fa path1_reads.1.fq path1_reads.2.fq
wgsim -N 5000 -1 150 -2 150 -r 0 -R 0 -X 0 -e 0 path2.fa path2_reads.1.fq path2_reads.2.fq

# Combine reads
cat path1_reads.1.fq path2_reads.1.fq > combined.1.fq
cat path1_reads.2.fq path2_reads.2.fq > combined.2.fq

# Map to graph paths (not original sequences)
echo "Mapping reads to graph..."
bwa mem -t 8 hla-f.k51.paths.fa combined.1.fq combined.2.fq | samtools view -bS - > combined.bam

# Convert to GAF
echo "Converting to GAF..."
gfainject --gfa hla-f.k51.gfa --bam combined.bam > combined.gaf

# Compute coverage
echo "Computing coverage..."
gafpack --gfa hla-f.k51.gfa --gaf combined.gaf --len-scale --weight-queries | gzip > combined.gafpack.gz

# Run genotyping
echo "Running genotyping..."
RUST_LOG=info cargo run --release -- genotype \
    --paths hla-f.k51.paths.coverage.tsv.gz \
    --gaf combined.gafpack.gz \
    --output validate_results \
    --id validate_sample \
    --ploidy 2

echo ""
echo "=== RESULTS ==="
echo "Expected: ${PATHS[0]} + ${PATHS[1]}"
echo "Got:"
cat validate_results/validate_sample.genotype.tsv

echo ""
echo "Checking if correct genotype is in top 5:"
zcat validate_results/validate_sample.sorted_combos.tsv.gz | head -6

echo ""
echo "Searching for exact match:"
zcat validate_results/validate_sample.sorted_combos.tsv.gz | grep -E "${PATHS[0]}.*${PATHS[1]}|${PATHS[1]}.*${PATHS[0]}" | head -1 || echo "Not found in results!"
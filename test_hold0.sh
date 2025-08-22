#!/bin/bash
set -e

echo "Testing hold-0-out genotyping..."

# We'll use the first two sequences as our "sample"
# Extract two sequences to simulate reads from
odgi paths -i hla-f.k51.sorted.og -L | head -2 | cut -f1 > test_paths.txt

echo "Selected paths for simulation:"
cat test_paths.txt

# Extract those sequences
samtools faidx hla-f.fa.gz $(cat test_paths.txt | head -1) > test_seq1.fa
samtools faidx hla-f.fa.gz $(cat test_paths.txt | tail -1) > test_seq2.fa

# Simulate reads from these sequences
echo "Simulating reads..."
wgsim -N 1000 -1 150 -2 150 -r 0 -R 0 -X 0 test_seq1.fa test_seq1.1.fq test_seq1.2.fq
wgsim -N 1000 -1 150 -2 150 -r 0 -R 0 -X 0 test_seq2.fa test_seq2.1.fq test_seq2.2.fq

# Combine reads
cat test_seq1.1.fq test_seq2.1.fq > test_sample.1.fq
cat test_seq1.2.fq test_seq2.2.fq > test_sample.2.fq

# Map reads to graph sequences
echo "Extracting graph sequences..."
odgi paths -i hla-f.k51.sorted.og -f > hla-f.k51.paths.fa

echo "Mapping reads..."
bwa index hla-f.k51.paths.fa
bwa mem -t 8 hla-f.k51.paths.fa test_sample.1.fq test_sample.2.fq | samtools view -bS - > test_sample.bam

# Convert to GAF
echo "Converting to GAF..."
gfainject --gfa hla-f.k51.gfa --bam test_sample.bam > test_sample.gaf

# Get coverage with gafpack
echo "Computing coverage..."
gafpack --gfa hla-f.k51.gfa --gaf test_sample.gaf --len-scale --weight-queries | gzip > test_sample.gafpack.gz

# Run genotyping
echo "Running genotyping..."
RUST_LOG=info cargo run -- genotype \
    --paths hla-f.k51.paths.coverage.tsv.gz \
    --gaf test_sample.gafpack.gz \
    --output hold0_results \
    --id test_sample \
    --ploidy 2

echo "Expected paths:"
cat test_paths.txt

echo "Results:"
cat hold0_results/test_sample.genotype.tsv
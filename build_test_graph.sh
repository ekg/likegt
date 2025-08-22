#!/bin/bash
set -e

# Build a test graph from HLA-F sequences
echo "Building test graph from HLA-F sequences..."

# Step 1: All-vs-all alignment
echo "Running allwave..."
allwave -i hla-f.fa.gz -t 8 -p none > hla-f.paf

# Step 2: Build graph with k=51 for testing
k=51
echo "Building graph with k=$k..."
seqwish -s hla-f.fa.gz -g hla-f.seqwish-k$k.gfa -t 8 -p hla-f.paf -k $k -P

# Step 3: Sort the graph
echo "Sorting graph..."
odgi build -g hla-f.seqwish-k$k.gfa -o hla-f.k$k.og
odgi sort -i hla-f.k$k.og -p Ygs -o hla-f.k$k.sorted.og
odgi view -i hla-f.k$k.sorted.og -g > hla-f.k$k.gfa

# Step 4: Extract path coverage for reference haplotypes
echo "Extracting path coverage..."
odgi paths -i hla-f.k$k.sorted.og -L | gzip > hla-f.k$k.paths.tsv.gz

echo "Graph construction complete!"
echo "Graph: hla-f.k$k.gfa"
echo "Paths: hla-f.k$k.paths.tsv.gz"
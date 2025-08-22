#!/bin/bash
set -euo pipefail

# Generate real coverage matrix for all individuals in the reference graph
# This replaces the fake binary path presence matrix

GRAPH="../hla-f.gfa.k51.B50.og"
REF_FASTA="../hla-f.fa.gz"
OUTPUT_DIR="real_coverage_matrix"
COVERAGE_DEPTH=30
READ_LENGTH=150

mkdir -p "$OUTPUT_DIR"

# Extract all haplotypes from the reference
echo "Extracting haplotypes from reference..."
odgi paths -i "$GRAPH" -f > "$OUTPUT_DIR/all_paths.fa"

# Build BWA index
echo "Building BWA index..."
cd "$OUTPUT_DIR"
bwa index all_paths.fa

# Process each individual
echo "Processing individuals..."
zcat "$REF_FASTA" | grep "^>" | sed 's/>//' | while read -r path_name; do
    individual=$(echo "$path_name" | cut -d'#' -f1)
    haplotype=$(echo "$path_name" | cut -d'#' -f2)
    
    echo "Processing $individual haplotype $haplotype..."
    
    # Extract this haplotype's sequence
    samtools faidx ../hla-f.fa.gz "$path_name" > "${individual}_${haplotype}.fa"
    
    # Calculate read pairs needed for coverage
    seq_length=$(grep -v "^>" "${individual}_${haplotype}.fa" | tr -d '\n' | wc -c)
    total_bases=$((seq_length * COVERAGE_DEPTH))
    bases_per_pair=$((2 * READ_LENGTH))
    n_pairs=$((total_bases / bases_per_pair))
    
    echo "  Sequence length: $seq_length bp"
    echo "  Simulating $n_pairs read pairs for ${COVERAGE_DEPTH}x coverage..."
    
    # Simulate reads
    wgsim -1 $READ_LENGTH -2 $READ_LENGTH -N $n_pairs -e 0.01 -r 0.001 \
        "${individual}_${haplotype}.fa" \
        "${individual}_${haplotype}.1.fq" \
        "${individual}_${haplotype}.2.fq" 2>/dev/null
    
    # Map reads
    bwa mem -t 4 all_paths.fa \
        "${individual}_${haplotype}.1.fq" \
        "${individual}_${haplotype}.2.fq" 2>/dev/null | \
        samtools view -b - > "${individual}_${haplotype}.bam"
    
    # Convert to GAF
    samtools view "${individual}_${haplotype}.bam" | \
        awk '{print $1"\t"length($10)"\t0\t"length($10)"\t+\t"$3"\t"length($10)"\t0\t"length($10)"\t"length($10)"\t"length($10)"\t60"}' > \
        "${individual}_${haplotype}.gaf"
    
    # Get coverage
    gafpack coverage -g ../hla-f.gfa.k51.B50.og -a "${individual}_${haplotype}.gaf" \
        -o "${individual}_${haplotype}.coverage.tsv"
    
    # Clean up intermediate files
    rm "${individual}_${haplotype}.1.fq" "${individual}_${haplotype}.2.fq" 
    rm "${individual}_${haplotype}.bam" "${individual}_${haplotype}.gaf"
done

# Combine all coverage files into matrix
echo "Combining coverage files into matrix..."
python3 << 'EOF'
import os
import glob

# Get all coverage files
coverage_files = sorted(glob.glob("*.coverage.tsv"))
print(f"Found {len(coverage_files)} coverage files")

# Read header from first file
with open(coverage_files[0], 'r') as f:
    header = f.readline().strip()
    
# Write combined matrix
with open("reference_coverage_matrix.tsv", 'w') as out:
    # Write header
    out.write(header + '\n')
    
    # Write each sample's coverage
    for cf in coverage_files:
        sample_name = cf.replace('.coverage.tsv', '')
        with open(cf, 'r') as f:
            f.readline()  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                # Replace sample column with actual sample name
                parts[0] = sample_name
                out.write('\t'.join(parts) + '\n')

print("Created reference_coverage_matrix.tsv")
EOF

echo "Done! Reference coverage matrix saved to $OUTPUT_DIR/reference_coverage_matrix.tsv"
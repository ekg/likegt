#!/bin/bash
set -euo pipefail

# Test version - just process 2 individuals to verify the approach works

GRAPH="hla-f.k51.og"
REF_FASTA="hla-f.fa.gz"
OUTPUT_DIR="test_reference_coverage"
COVERAGE_DEPTH=30
READ_LENGTH=150

mkdir -p "$OUTPUT_DIR"

# Extract paths for testing
echo "Extracting paths..."
odgi paths -i "$GRAPH" -f > "$OUTPUT_DIR/all_paths.fa"

# Build BWA index
echo "Building BWA index..."
cd "$OUTPUT_DIR"
bwa index all_paths.fa 2>/dev/null

# Process just HG00096 and HG00268 for testing
for individual in HG00096 HG00268; do
    for hap_num in 1 2; do
        # Find the exact path name
        path_name=$(zcat "../$REF_FASTA" | grep "^>${individual}#${hap_num}#" | sed 's/>//' | head -1)
        
        if [ -z "$path_name" ]; then
            echo "Warning: No path found for ${individual} haplotype ${hap_num}"
            continue
        fi
        
        echo "Processing $path_name..."
        
        # Extract this haplotype's sequence
        samtools faidx "../$REF_FASTA" "$path_name" > "${individual}_hap${hap_num}.fa"
        
        # Calculate read pairs needed for coverage
        seq_length=$(grep -v "^>" "${individual}_hap${hap_num}.fa" | tr -d '\n' | wc -c)
        total_bases=$((seq_length * COVERAGE_DEPTH))
        bases_per_pair=$((2 * READ_LENGTH))
        n_pairs=$((total_bases / bases_per_pair))
        
        echo "  Sequence length: $seq_length bp"
        echo "  Simulating $n_pairs read pairs for ${COVERAGE_DEPTH}x coverage..."
        
        # Simulate reads
        wgsim -1 $READ_LENGTH -2 $READ_LENGTH -N $n_pairs -e 0.01 -r 0.001 \
            "${individual}_hap${hap_num}.fa" \
            "${individual}_hap${hap_num}.1.fq" \
            "${individual}_hap${hap_num}.2.fq" 2>/dev/null
        
        # Map reads to pangenome
        echo "  Mapping reads..."
        bwa mem -t 4 all_paths.fa \
            "${individual}_hap${hap_num}.1.fq" \
            "${individual}_hap${hap_num}.2.fq" 2>/dev/null | \
            samtools view -b - > "${individual}_hap${hap_num}.bam"
        
        # Convert to GAF
        echo "  Converting to GAF..."
        gfainject --gfa "../hla-f.k51.gfa" --bam "${individual}_hap${hap_num}.bam" > \
            "${individual}_hap${hap_num}.gaf" 2>/dev/null
        
        # Get coverage
        echo "  Computing coverage..."
        gafpack --gfa "../hla-f.k51.gfa" --gaf "${individual}_hap${hap_num}.gaf" > \
            "${individual}_hap${hap_num}.coverage.tsv"
        
        # Clean up intermediate files
        rm "${individual}_hap${hap_num}.1.fq" "${individual}_hap${hap_num}.2.fq" 
        rm "${individual}_hap${hap_num}.bam"
    done
done

# Combine coverage files  
echo "Creating combined matrix..."
python3 << 'EOF'
import os
import glob

# Get all coverage files
coverage_files = sorted(glob.glob("*_hap*.coverage.tsv"))
print(f"Found {len(coverage_files)} coverage files: {coverage_files}")

if len(coverage_files) == 0:
    print("ERROR: No coverage files found!")
    exit(1)

# Read header from first file
with open(coverage_files[0], 'r') as f:
    header = f.readline().strip()
    
# Write combined matrix
with open("test_reference_matrix.tsv", 'w') as out:
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

print("Created test_reference_matrix.tsv")

# Show some stats
import numpy as np

all_coverage = []
for cf in coverage_files:
    with open(cf, 'r') as f:
        f.readline()  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            values = [float(x) for x in parts[1:]]
            all_coverage.append(values)
            
all_coverage = np.array(all_coverage)
print(f"\nMatrix shape: {all_coverage.shape}")
print(f"Coverage range: {all_coverage.min():.0f} - {all_coverage.max():.0f}")
print(f"Mean coverage per node: {all_coverage.mean():.1f}")
print(f"Non-zero fraction: {(all_coverage > 0).mean():.3f}")
EOF

echo "Done! Test reference matrix saved to $OUTPUT_DIR/test_reference_matrix.tsv"
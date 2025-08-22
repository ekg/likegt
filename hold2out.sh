#!/bin/bash
# Hold-2-out validation script
# Remove an individual's haplotypes, simulate reads, map, and genotype

set -euo pipefail

# Arguments
INDIVIDUAL=${1:-HG00096}
K=${2:-51}
COVERAGE=${3:-30}
READ_LENGTH=${4:-150}
THREADS=${5:-8}

echo "=== Hold-2-out validation for $INDIVIDUAL ==="
echo "Parameters: k=$K, coverage=${COVERAGE}x, read_length=$READ_LENGTH"

# Input files
FASTA="hla-f.fa.gz"
GFA="hla-f.k${K}.gfa"
REFERENCE_PATHS="hla-f.k${K}.paths.coverage.tsv.gz"

# Check required files exist
if [[ ! -f "$FASTA" ]]; then
    echo "Error: $FASTA not found"
    exit 1
fi

if [[ ! -f "$GFA" ]]; then
    echo "Error: $GFA not found. Run: ./build_all_graphs.sh first"
    exit 1
fi

# Create working directory
WORKDIR="hold2out_${INDIVIDUAL}_k${K}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

echo "Working in $WORKDIR"

# Step 1: Extract the individual's sequences
echo "Step 1: Extracting sequences for $INDIVIDUAL"
samtools faidx "../$FASTA" "${INDIVIDUAL}#1#hla_f#0" > "${INDIVIDUAL}.hap1.fa" 2>/dev/null || \
    samtools faidx "../$FASTA" "${INDIVIDUAL}#1#HLA-F#0" > "${INDIVIDUAL}.hap1.fa" 2>/dev/null || \
    samtools faidx "../$FASTA" "${INDIVIDUAL}#1" > "${INDIVIDUAL}.hap1.fa" 2>/dev/null || \
    echo "Warning: Could not extract haplotype 1"

samtools faidx "../$FASTA" "${INDIVIDUAL}#2#hla_f#0" > "${INDIVIDUAL}.hap2.fa" 2>/dev/null || \
    samtools faidx "../$FASTA" "${INDIVIDUAL}#2#HLA-F#0" > "${INDIVIDUAL}.hap2.fa" 2>/dev/null || \
    samtools faidx "../$FASTA" "${INDIVIDUAL}#2" > "${INDIVIDUAL}.hap2.fa" 2>/dev/null || \
    echo "Warning: Could not extract haplotype 2"

# Combine haplotypes
cat "${INDIVIDUAL}.hap1.fa" "${INDIVIDUAL}.hap2.fa" > "${INDIVIDUAL}.diploid.fa"

# Check if we got sequences
if [[ ! -s "${INDIVIDUAL}.diploid.fa" ]]; then
    echo "Error: Failed to extract sequences for $INDIVIDUAL"
    echo "Available sequences:"
    samtools faidx "../$FASTA"
    exit 1
fi

echo "Extracted $(grep -c "^>" "${INDIVIDUAL}.diploid.fa") sequences"

# Step 2: Create graph without this individual
echo "Step 2: Creating graph without $INDIVIDUAL"
# Remove paths from GFA
grep -v "${INDIVIDUAL}#" "../$GFA" > "graph_without_${INDIVIDUAL}.gfa" || true

# Count remaining paths
REMAINING_PATHS=$(grep "^P" "graph_without_${INDIVIDUAL}.gfa" | wc -l || echo 0)
echo "Graph has $REMAINING_PATHS paths (removed 2)"

# Step 3: Simulate reads
echo "Step 3: Simulating reads at ${COVERAGE}x coverage"
# Calculate number of reads needed
SEQ_LENGTH=$(grep -v "^>" "${INDIVIDUAL}.diploid.fa" | tr -d '\n' | wc -c)
NUM_READS=$((SEQ_LENGTH * COVERAGE / READ_LENGTH / 2))  # /2 for paired-end

echo "Sequence length: $SEQ_LENGTH bp"
echo "Simulating $NUM_READS read pairs"

# Use wgsim to simulate reads
wgsim -1 $READ_LENGTH -2 $READ_LENGTH \
      -N $NUM_READS \
      -e 0.001 \
      -r 0.001 \
      -R 0.001 \
      "${INDIVIDUAL}.diploid.fa" \
      "${INDIVIDUAL}.reads.1.fq" \
      "${INDIVIDUAL}.reads.2.fq" \
      2> wgsim.log

echo "Generated $(wc -l < "${INDIVIDUAL}.reads.1.fq") / 4 = $(($(wc -l < "${INDIVIDUAL}.reads.1.fq") / 4)) reads"

# Step 4: Map reads to graph
echo "Step 4: Mapping reads to graph"

# First extract paths as FASTA for BWA index
odgi paths -i "graph_without_${INDIVIDUAL}.gfa" -f > paths_without_${INDIVIDUAL}.fa

# Build BWA index
bwa index paths_without_${INDIVIDUAL}.fa 2>/dev/null

# Map reads
bwa mem -t $THREADS paths_without_${INDIVIDUAL}.fa \
    "${INDIVIDUAL}.reads.1.fq" \
    "${INDIVIDUAL}.reads.2.fq" \
    2> bwa.log | \
    samtools view -bS - > "${INDIVIDUAL}.mapped.bam"

# Check mapping
samtools flagstat "${INDIVIDUAL}.mapped.bam"

# Step 5: Inject alignments into graph
echo "Step 5: Injecting alignments into graph"
gfainject "graph_without_${INDIVIDUAL}.gfa" \
          "${INDIVIDUAL}.mapped.bam" \
          > "${INDIVIDUAL}.gaf" 2> gfainject.log || {
    echo "gfainject failed. Trying with sorted BAM..."
    samtools sort "${INDIVIDUAL}.mapped.bam" > "${INDIVIDUAL}.sorted.bam"
    samtools index "${INDIVIDUAL}.sorted.bam"
    gfainject "graph_without_${INDIVIDUAL}.gfa" \
              "${INDIVIDUAL}.sorted.bam" \
              > "${INDIVIDUAL}.gaf" 2> gfainject.retry.log
}

# Step 6: Get coverage from GAF
echo "Step 6: Extracting coverage from GAF"
gafpack coverage -g "${INDIVIDUAL}.gaf" \
                 -n "graph_without_${INDIVIDUAL}.gfa" \
                 -b 1 -s 1 \
                 > "${INDIVIDUAL}.coverage.tsv" 2> gafpack.log

# Check coverage output
if [[ ! -s "${INDIVIDUAL}.coverage.tsv" ]]; then
    echo "Warning: Coverage file is empty"
    echo "GAF file size: $(wc -l < "${INDIVIDUAL}.gaf") lines"
    head "${INDIVIDUAL}.gaf"
fi

# Step 7: Run genotyping
echo "Step 7: Running genotyping"

# Need to prepare reference coverage without the individual
zcat "../$REFERENCE_PATHS" | \
    grep -v "${INDIVIDUAL}#" | \
    gzip > "reference_without_${INDIVIDUAL}.tsv.gz" || true

# Run likegt genotyping
cargo run --release -- genotype \
    --paths "reference_without_${INDIVIDUAL}.tsv.gz" \
    --gaf "${INDIVIDUAL}.coverage.tsv" \
    --output "${INDIVIDUAL}.genotype.json" \
    --id "$INDIVIDUAL" \
    --ploidy 2 \
    --threads $THREADS 2>&1 | tee genotype.log || {
    echo "Genotyping failed. Debug info:"
    echo "Reference paths: $(zcat "reference_without_${INDIVIDUAL}.tsv.gz" | wc -l) lines"
    echo "Coverage: $(wc -l < "${INDIVIDUAL}.coverage.tsv") lines"
}

# Step 8: Evaluate results
echo "Step 8: Evaluating results"
if [[ -f "${INDIVIDUAL}.genotype.json" ]]; then
    echo "Genotype results:"
    jq . "${INDIVIDUAL}.genotype.json"
    
    # Extract called genotype
    CALLED_HAP1=$(jq -r '.haplotypes[0]' "${INDIVIDUAL}.genotype.json")
    CALLED_HAP2=$(jq -r '.haplotypes[1]' "${INDIVIDUAL}.genotype.json")
    SIMILARITY=$(jq -r '.similarity' "${INDIVIDUAL}.genotype.json")
    
    echo ""
    echo "=== RESULTS ==="
    echo "True genotype:   ${INDIVIDUAL}#1, ${INDIVIDUAL}#2"
    echo "Called genotype: $CALLED_HAP1, $CALLED_HAP2"
    echo "Similarity:      $SIMILARITY"
    
    # Check if correct
    if [[ "$CALLED_HAP1" == *"${INDIVIDUAL}"* ]] || [[ "$CALLED_HAP2" == *"${INDIVIDUAL}"* ]]; then
        echo "ERROR: Called genotype contains the held-out individual!"
    else
        echo "Success: Held-out individual not in called genotype (expected for hold-2-out)"
    fi
else
    echo "No genotype results found"
fi

echo ""
echo "Hold-2-out validation complete for $INDIVIDUAL"
echo "Results saved in $WORKDIR"
#!/bin/bash
# Generate REAL reference coverage from simulated reads
# NOT path presence bullshit

set -euo pipefail

FASTA="hla-f.fa.gz"
GFA="hla-f.k51.gfa"
COVERAGE=30
THREADS=8
OUTDIR="real_reference_coverage"

mkdir -p $OUTDIR

echo "=== GENERATING REAL REFERENCE COVERAGE ==="
echo "This will take a while - simulating ${COVERAGE}x reads for ALL individuals"
echo ""

# Extract all paths as reference for mapping
echo "Extracting graph paths..."
odgi paths -i $GFA -f > $OUTDIR/paths.fa
bwa index $OUTDIR/paths.fa 2>/dev/null

# Get list of individuals
samtools faidx $FASTA
INDIVIDUALS=$(grep "^>" hla-f.fa.gz.fai | cut -f1 | cut -d'#' -f1 | sort -u)

# Initialize coverage matrix file
echo -e "individual\t$(odgi paths -i $GFA -L | head -1 | cut -f1)" > $OUTDIR/coverage_matrix.tsv

for IND in $INDIVIDUALS; do
    echo "Processing $IND..."
    
    # Extract sequences for this individual
    samtools faidx $FASTA $(grep "^${IND}#" hla-f.fa.gz.fai | cut -f1) > $OUTDIR/${IND}.fa 2>/dev/null || continue
    
    # Get sequence length
    SEQ_LEN=$(grep -v "^>" $OUTDIR/${IND}.fa | tr -d '\n' | wc -c)
    
    if [[ $SEQ_LEN -eq 0 ]]; then
        echo "  Skipping $IND - no sequence"
        continue
    fi
    
    # Calculate read pairs needed for target coverage
    READ_PAIRS=$((SEQ_LEN * COVERAGE / 300))  # 150bp paired-end reads
    
    echo "  Sequence length: $SEQ_LEN bp"
    echo "  Simulating $READ_PAIRS read pairs for ${COVERAGE}x coverage"
    
    # Simulate reads
    wgsim -1 150 -2 150 -N $READ_PAIRS -e 0.001 -r 0.001 \
        $OUTDIR/${IND}.fa \
        $OUTDIR/${IND}.1.fq \
        $OUTDIR/${IND}.2.fq 2>/dev/null
    
    # Map to graph paths
    bwa mem -t $THREADS $OUTDIR/paths.fa \
        $OUTDIR/${IND}.1.fq \
        $OUTDIR/${IND}.2.fq 2>/dev/null | \
        samtools view -bS - > $OUTDIR/${IND}.bam
    
    # Inject into graph
    gfainject --gfa $GFA --bam $OUTDIR/${IND}.bam > $OUTDIR/${IND}.gaf 2>/dev/null
    
    # Get coverage
    gafpack --gfa $GFA --gaf $OUTDIR/${IND}.gaf > $OUTDIR/${IND}.coverage.tsv
    
    # Extract coverage values and add to matrix
    tail -n1 $OUTDIR/${IND}.coverage.tsv | cut -f2- | \
        awk -v ind="$IND" '{printf "%s\t%s\n", ind, $0}' >> $OUTDIR/coverage_matrix.tsv
    
    # Clean up intermediate files
    rm -f $OUTDIR/${IND}.{fa,1.fq,2.fq,bam,gaf}
    
    echo "  Done"
done

echo ""
echo "Coverage matrix saved to: $OUTDIR/coverage_matrix.tsv"
echo "This is the REAL reference coverage from ${COVERAGE}x reads"
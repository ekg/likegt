# LikeGT Example Workflow

This document demonstrates a complete workflow using LikeGT for pangenome genotyping validation.

## Complete Workflow Example

### Step 1: Prepare Your Data

```bash
# Assuming you have:
# - hla-f.fa.gz: FASTA file with haplotype sequences
# - hla-f.k51.gfa: Pangenome graph file

# Check your data
zcat hla-f.fa.gz | grep ">" | head -5
# >HG00096#1#haplotype1-0000024:29677236-29701492
# >HG00096#2#haplotype2-0000117:29686034-29710088
# >HG00171#1#haplotype1-0000033:143253257-143277343
# >HG00171#2#haplotype2-0000168:142352281-142376370
# >HG00268#1#haplotype1-0000024:141273481-141297761
```

### Step 2: Validate the Graph

```bash
# Check if the graph is suitable for genotyping
cargo run --release -- check -g hla-f.k51.gfa -v

# Expected output:
# ✓ Graph loaded successfully
# ✓ Graph has 787653 nodes
# ✓ Graph has paths
# ✓ Graph appears suitable for genotyping
```

### Step 3: Run Hold-2-out Validation

#### Single Individual Test
```bash
# Test a single individual with verbose output
cargo run --release -- hold-out \
    -f hla-f.fa.gz \
    -g hla-f.k51.gfa \
    -i HG00096 \
    --verbose

# Output shows pipeline stages:
# 🔍 Stage 1: Extracting 2 sequence(s) for HG00096
# 📊 Stage 2: Loading reference coverage
# 🧪 Stage 4: Simulating reads (wgsim, 150bp, 30x coverage)
# 🎯 Stage 5: Aligning reads with minimap2 (preset: sr)
# 📦 Stage 5: Generating sample coverage from GAF alignments
# 🔬 Stage 6: Running COSIGT genotyping
# ✅ Complete hold-2-out pipeline finished in 1.12s
# 📊 Result: CORRECT (similarity: 0.9991, rank: 1, alignment: 100.0%)
```

#### Batch Processing
```bash
# Test multiple individuals
cargo run --release -- hold-out \
    -f hla-f.fa.gz \
    -g hla-f.k51.gfa \
    -i "HG00096,HG00171,HG00268" \
    --format table

# Output (TSV format):
# sample    true_hap1    true_hap2    called_hap1    called_hap2    similarity    qv    alignment    bias_loss    time
# HG00096   HG00096#1    HG00096#2    HG00733#2      HG01352#1      0.9986        28.6  100.0%       0.0%         1.21s
# HG00171   HG00171#1    HG00171#2    HG00514#2      NA20847#2      0.9992        31.1  100.0%       0.0%         1.10s
# HG00268   HG00268#1    HG00268#2    HG00358#2      HG01114#2      0.9982        27.5  100.0%       0.0%         1.13s
```

#### Process All Individuals
```bash
# Test all individuals in the dataset
cargo run --release -- hold-out \
    -f hla-f.fa.gz \
    -g hla-f.k51.gfa \
    -i all \
    --format tsv \
    -o hold2out_results/all_results.tsv

# Monitor progress
# Processing 67 individuals...
# [1/67] Testing HG00096... PASS
# [2/67] Testing HG00171... PASS
# ...
```

### Step 4: Compute Maximum Attainable QV

```bash
# Analyze the reference coverage matrix
cargo run --release -- max-qv \
    -c ./hold2out_results/reduced_reference_coverage.tsv.gz \
    -v

# Output shows theoretical best performance:
# 🔬 Computing maximum attainable QV for all samples...
# 📊 Loaded 130 haplotypes from ./hold2out_results/reduced_reference_coverage.tsv.gz
# ✅ Computed max QV for 64 individuals
# 
# === MAXIMUM ATTAINABLE QV REPORT ===
# This shows the best possible QV when genotyping against non-self haplotypes
# 
# Total samples analyzed: 64
# Average maximum attainable QV: 60.0
# Average best achievable rank: 2.1
# 
# Sample    Max_QV    Best_Rank    Best_Similarity    Best_Non-Self_Genotype
# HG03520   60.0      2            1.000000           HG03452#1 + NA19129#1
# NA19705   60.0      2            1.000000           HG00732#1 + NA19347#2
# ...
```

### Step 5: Compare Different K-mer Sizes

```bash
# Test multiple k-mer graphs
for k in 25 51 101; do
    echo "Testing k=$k..."
    cargo run --release -- hold-out \
        -f hla-f.fa.gz \
        -g hla-f.k${k}.gfa \
        -i "HG00096,HG00171,HG00268" \
        --format csv \
        > results_k${k}.csv
done

# Compare results
echo "K-value,Mean_Similarity,Mean_QV"
for k in 25 51 101; do
    awk -F'\t' 'NR>1 {sum+=$6; qv+=$7; n++} END {printf "%d,%.4f,%.1f\n", k, sum/n, qv/n}' \
        k=$k results_k${k}.csv
done
```

### Step 6: Advanced Options

#### Use Different Aligners
```bash
# Use BWA-MEM instead of minimap2
cargo run --release -- hold-out \
    -f hla-f.fa.gz \
    -g hla-f.k51.gfa \
    -i HG00096 \
    --aligner bwa-mem \
    --preset short \
    -v
```

#### Adjust Simulation Parameters
```bash
# Higher coverage, longer reads
cargo run --release -- hold-out \
    -f hla-f.fa.gz \
    -g hla-f.k51.gfa \
    -i HG00096 \
    --coverage-depth 50 \
    --read-length 250 \
    --fragment-length 600 \
    --fragment-std 100 \
    -v
```

#### Keep Intermediate Files for Debugging
```bash
cargo run --release -- hold-out \
    -f hla-f.fa.gz \
    -g hla-f.k51.gfa \
    -i HG00096 \
    --keep-files \
    -o debug_output \
    -v

# Check intermediate files
ls -la debug_output/HG00096/
# sequences.fa      # Extracted sequences
# reads_1.fq       # Simulated reads (pair 1)
# reads_2.fq       # Simulated reads (pair 2)
# alignment.sam    # Alignment file
# alignment.gaf    # Graph alignment
# sample.cov.tsv   # Sample coverage
```

### Step 7: Output Formats

#### JSON Output (for programmatic processing)
```bash
cargo run --release -- hold-out \
    -f hla-f.fa.gz \
    -g hla-f.k51.gfa \
    -i HG00096 \
    --format json \
    > result.json

# Parse with jq
cat result.json | jq '.graph_qv, .cosine_similarity'
```

#### Human-Readable Output
```bash
cargo run --release -- hold-out \
    -f hla-f.fa.gz \
    -g hla-f.k51.gfa \
    -i HG00096 \
    --format text

# Output:
# PIPELINE: PASS | Individual: HG00096 | Similarity: 0.9991 | Alignment: 100.0% | QV: 30.7 | Time: 1.12s
```

## Performance Tips

1. **Always use `--release` builds** for production work (10-100x faster)
2. **Batch processing** is more efficient than individual runs
3. **Adjust thread count** based on your system: `-t 8` for 8 threads
4. **Monitor memory usage** with large datasets
5. **Use appropriate output format**: TSV/CSV for analysis, JSON for automation

## Validation Metrics Explained

- **Similarity**: Cosine similarity between called and true genotype coverage vectors (0-1, higher is better)
- **QV (Quality Value)**: Phred-scaled quality score (-10 * log10(error_probability))
  - QV 20 = 99% accuracy
  - QV 30 = 99.9% accuracy
  - QV 40 = 99.99% accuracy
  - QV 60 = Perfect match (maximum)
- **Rank**: Position of true genotype in ranked results (1 = best match was correct)
- **Alignment Rate**: Percentage of reads successfully aligned to graph

## Troubleshooting Common Issues

### Issue: "Command not found"
```bash
# Check if required tools are installed
which minimap2 gfainject gafpack samtools seqtk odgi
# Install missing tools as needed
```

### Issue: High memory usage
```bash
# Process in smaller batches
cargo run --release -- hold-out \
    -f hla-f.fa.gz \
    -g hla-f.k51.gfa \
    -i "HG00096,HG00171" \  # Process 2 at a time
    --format tsv
```

### Issue: Slow performance
```bash
# Use release build and more threads
cargo build --release
./target/release/likegt hold-out \
    -f hla-f.fa.gz \
    -g hla-f.k51.gfa \
    -i all \
    -t 16  # Use 16 threads
```

## Next Steps

1. **Explore different graphs**: Try various k-mer sizes to find optimal parameters
2. **Test different populations**: Validate across diverse genetic backgrounds
3. **Optimize parameters**: Tune coverage depth, read length for your use case
4. **Integrate with pipelines**: Use JSON output for automated workflows
5. **Contribute**: Report issues and submit improvements on GitHub
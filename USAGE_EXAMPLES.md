# LikeGT Usage Examples and Documentation

## Overview

LikeGT is a toolkit for pangenome graph-based genotyping validation. It implements the COSIGT algorithm using cosine similarity between coverage vectors.

## Important Clarification: Coverage-based vs Sequence-based QV

### What the current implementation does:
- **Coverage-based QV**: Computes Quality Values based on cosine similarity between coverage vectors
- The `max-qv` command finds the best matching pair of non-self haplotypes in **coverage space**
- QV is calculated from similarity scores, NOT from actual sequence alignments

### What it does NOT do (yet):
- Does NOT align sequences with allwave to compute edit distances
- Does NOT calculate QV from actual sequence alignment metrics
- The `sequence_qv` module exists but is not integrated with max QV computation

### For true sequence-based max attainable QV, you would need:
1. Extract actual DNA sequences for held-out haplotypes
2. Find best matching non-self haplotype pair
3. Perform 4 alignments (2x2 pairings) with allwave/biWFA
4. Calculate QV from actual edit distances

## Basic Commands

### 1. Show help
```bash
cargo run -- --help
```

### 2. Build a pangenome graph
```bash
# Build graph from FASTA using allwave + seqwish + odgi
cargo run -- build -f input.fa -o output_prefix -k 51,101 -t 8
```

### 3. Check graph suitability
```bash
cargo run -- check -g graph.gfa -v
```

### 4. Run hold-2-out validation
```bash
# Single individual
cargo run -- hold-out -f hla-f.fa.gz -g hla-f.k51.gfa -i HG00096 -v

# All individuals
cargo run -- hold-out -f hla-f.fa.gz -g hla-f.k51.gfa -i all -v

# With sequence QV (uses biWFA, not allwave)
cargo run -- hold-out -f hla-f.fa.gz -g hla-f.k51.gfa -i HG00096 --sequence-qv -v
```

### 5. Compute maximum attainable QV (coverage-based)
```bash
# Compute max QV for all samples in a coverage matrix
cargo run -- max-qv -c reference_coverage.tsv.gz -v

# Save results to file
cargo run -- max-qv -c reference_coverage.tsv.gz -o max_qv_results.txt
```

## Full Pipeline Example

```bash
# 1. Build graph
cargo run --release -- build \
    -f hla-f.fa.gz \
    -o hla-f \
    -k 51 \
    -t 8

# 2. Check the graph
cargo run --release -- check \
    -g hla-f.k51.gfa \
    -v

# 3. Run hold-2-out validation for specific individual
cargo run --release -- hold-out \
    -f hla-f.fa.gz \
    -g hla-f.k51.gfa \
    -i HG00096 \
    --threads 8 \
    --coverage-depth 30 \
    --verbose

# 4. Compute max attainable QV
cargo run --release -- max-qv \
    -c ./hold2out_results/reduced_reference_coverage.tsv.gz \
    -v
```

## Output Formats

### Hold-out validation outputs:
- `text`: Human-readable format
- `json`: JSON format for programmatic processing
- `table`/`tsv`/`csv`: Tab-separated values for analysis

### Example output interpretation:
```
PIPELINE: PASS | Individual: HG00096 | Similarity: 0.9876 | Alignment: 95.2% | QV: 45.3 | Time: 12.34s
```
- **PASS/FAIL**: Whether the correct genotype was called
- **Similarity**: Cosine similarity between called and true genotype (coverage-based)
- **Alignment**: Percentage of reads successfully aligned
- **QV**: Quality Value (higher is better, max 60)
- **Time**: Pipeline execution time

## Testing

```bash
# Run all tests
cargo test

# Run specific test suite
cargo test --lib         # Library tests only
cargo test --test integration_test  # Integration tests

# Run with verbose output
cargo test -- --nocapture
```

## Performance Notes

- Use `--release` flag for production runs (10-100x faster)
- The `max-qv` computation scales as O(n²) for n haplotypes
- Coverage matrices can be large; ensure sufficient memory

## Current Limitations

1. **Coverage-based QV only**: The max QV calculation uses coverage similarity, not actual sequence alignment
2. **Diploid only**: Currently assumes ploidy=2
3. **No allwave integration**: Sequence QV uses biWFA, not allwave
4. **Memory intensive**: Large coverage matrices require significant RAM

## Future Improvements Needed

1. Integrate allwave for true sequence-based max QV computation
2. Support for arbitrary ploidy
3. Batch processing optimizations
4. GPU acceleration for similarity computations
5. Proper integration of sequence alignment QV with max QV calculations
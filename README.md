# LikeGT - Pangenome Graph-based Genotyping Toolkit

A high-performance toolkit for validating and benchmarking pangenome graph-based genotyping using the COSIGT (COSine similarity Genotyping) algorithm.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Commands](#commands)
- [Important Notes](#important-notes)
- [Performance](#performance)
- [Troubleshooting](#troubleshooting)

## Overview

LikeGT implements comprehensive validation pipelines for pangenome genotyping, including:
- Hold-k-out validation (k=0, 1, 2, ..., n)
- Maximum attainable QV computation
- Graph construction and validation
- Reference bias analysis
- Coverage-based genotype calling using cosine similarity

## Features

### Core Functionality
- **Hold-out Validation**: Complete pipeline from read simulation to genotype calling
- **Max QV Analysis**: Compute theoretical upper bounds on genotyping accuracy
- **Graph Building**: Construct pangenome graphs using allwave + seqwish + odgi
- **Graph Validation**: Check graph suitability for genotyping
- **Batch Processing**: Process multiple samples efficiently
- **Multiple Output Formats**: text, JSON, CSV, TSV, table

### Algorithms
- **COSIGT**: Cosine similarity-based genotyping from coverage vectors
- **Coverage-based QV**: Quality values computed from similarity scores
- **Sequence QV**: Optional alignment-based quality values using biWFA

## Installation

### Prerequisites

Required tools (must be in PATH):
```bash
# Core tools
cargo         # Rust compiler
odgi          # Graph manipulation
minimap2      # Read alignment
wgsim         # Read simulation
gfainject     # SAM/BAM to GAF conversion
gafpack       # Coverage calculation
samtools      # SAM/BAM processing
seqtk         # FASTA manipulation

# Optional tools
allwave       # For graph construction
seqwish       # For graph construction
bwa-mem       # Alternative aligner
pbsim3        # Long-read simulation
mason         # Alternative read simulator
```

### Building from Source

```bash
# Clone the repository
git clone https://github.com/yourusername/likegt.git
cd likegt

# Build release version (recommended)
cargo build --release

# Run tests
cargo test

# Install to cargo bin directory
cargo install --path .
```

## Quick Start

### 1. Basic Hold-2-out Validation

```bash
# Single individual
likegt hold-out -f sequences.fa.gz -g graph.gfa -i HG00096 -v

# Multiple individuals
likegt hold-out -f sequences.fa.gz -g graph.gfa -i "HG00096,HG00171,HG00268" -v

# All individuals
likegt hold-out -f sequences.fa.gz -g graph.gfa -i all --format tsv > results.tsv
```

### 2. Compute Maximum Attainable QV

```bash
# Analyze coverage matrix
likegt max-qv -c coverage_matrix.tsv.gz -v

# Save results to file
likegt max-qv -c coverage_matrix.tsv.gz -o max_qv_results.txt
```

### 3. Build Pangenome Graph

```bash
# Build graph with specific k-mer size
likegt build -f sequences.fa -o output_prefix -k 51 -t 8

# Build multiple k-mer graphs
likegt build -f sequences.fa -o output_prefix -k 25,51,101 -t 8
```

## Commands

### `hold-out` - Hold-out Validation Pipeline

Complete pipeline for hold-k-out validation including:
1. Sequence extraction
2. Read simulation
3. Alignment to graph
4. Coverage calculation
5. Genotype calling
6. QV computation

**Options:**
- `-f, --fasta`: Input FASTA file with haplotype sequences
- `-g, --graph`: Pangenome graph (.gfa or .og format)
- `-i, --individual`: Individual(s) to test (name, comma-separated, or "all")
- `--hold`: Number of haplotypes to hold out (default: 2)
- `-t, --threads`: Number of threads (default: 4)
- `--coverage-depth`: Simulated read coverage (default: 30)
- `--read-length`: Read length for simulation (default: 150)
- `--aligner`: Aligner to use (minimap2, bwa-mem, GraphAligner)
- `--simulator`: Read simulator (wgsim, mason, pbsim3)
- `--format`: Output format (text, json, csv, tsv, table)
- `-v, --verbose`: Verbose output

### `max-qv` - Maximum Attainable QV

Computes the best possible QV when genotyping using non-self haplotypes (nearest neighbor approach).

**Options:**
- `-c, --coverage`: Input coverage matrix (.tsv or .tsv.gz)
- `-o, --output`: Output file (optional, defaults to stdout)
- `-p, --ploidy`: Number of haplotypes per individual (default: 2)
- `-v, --verbose`: Verbose output

### `build` - Build Pangenome Graph

Constructs pangenome graphs using allwave + seqwish + odgi pipeline.

**Options:**
- `-f, --fasta`: Input FASTA file
- `-o, --output`: Output prefix for generated files
- `-k, --kmer-sizes`: K-mer sizes (comma-separated)
- `-t, --threads`: Number of threads
- `--keep-intermediates`: Keep intermediate files

### `check` - Validate Graph

Checks if a graph is suitable for genotyping.

**Options:**
- `-g, --graph`: Graph file to check (.gfa or .og)
- `-v, --verbose`: Verbose output

### `validate` - Run Validation Tests

Runs comprehensive validation tests on coverage data.

**Options:**
- `-f, --fasta`: Input FASTA file
- `-o, --output`: Output directory
- `-t, --threads`: Number of threads
- `-k, --kmer-size`: K-mer size for analysis

## Important Notes

### Coverage-based vs Sequence-based QV

**Current Implementation (Coverage-based):**
- QV computed from cosine similarity between coverage vectors
- `max-qv` finds best matching haplotypes in **coverage space**
- Fast and efficient but doesn't consider actual sequence similarity

**Sequence-based QV (Limited Support):**
- Uses biWFA for sequence alignment
- Available via `--sequence-qv` flag in hold-out command
- More accurate but computationally intensive
- NOT integrated with max-qv computation

**For true allwave-based max QV, you would need:**
1. Extract actual DNA sequences
2. Run allwave alignments (4 total for 2x2 pairings)
3. Compute QV from edit distances

### Performance Considerations

- Use `--release` builds for production (10-100x faster)
- Coverage matrices can be large; ensure sufficient RAM
- Batch processing is more efficient than individual runs
- Default threads: 4 (adjust based on your system)

## Performance

### Validated Results on HLA-F Dataset

| K-value | Mean Similarity | Mean QV | Samples | Status |
|---------|----------------|---------|---------|---------|
| k=51    | 99.70%        | 26.1    | 46      | **Optimal** |
| k=25    | 93.97%        | 12.7    | 66      | Acceptable |
| k=101   | 89.46%        | 10.1    | 67      | Poor |

### Benchmarks

- **Hold-2-out validation**: ~1.2s per individual
- **Max QV computation**: <1s for 130 haplotypes
- **Graph building**: Varies with input size
- **Memory usage**: ~2GB for typical HLA dataset

## Troubleshooting

### Common Issues

1. **"Command not found" errors**
   - Ensure all required tools are installed and in PATH
   - Check with: `which minimap2 gfainject gafpack`

2. **High QV values (60.0)**
   - Indicates perfect match - check if hold-out is working
   - Verify reference doesn't contain test individual

3. **Dimension mismatch errors**
   - Coverage matrices must have consistent dimensions
   - Check for metadata columns that need skipping

4. **Memory errors**
   - Large coverage matrices require significant RAM
   - Consider processing in batches or using a machine with more memory

### Debug Mode

Enable verbose output for detailed pipeline information:
```bash
likegt hold-out -f input.fa -g graph.gfa -i HG00096 --verbose
```

### File Formats

- **FASTA**: Gzipped or uncompressed, with haplotype IDs like "HG00096#1"
- **GFA**: Graph format from odgi/vg tools
- **Coverage Matrix**: TSV with haplotypes as rows, nodes as columns
- **Output**: JSON, CSV, TSV, or human-readable text

## Citation

If you use LikeGT in your research, please cite:

```
LikeGT: A toolkit for pangenome graph-based genotyping validation
[Your publication details here]
```

## License

[Specify your license here]

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## Support

For issues, questions, or suggestions:
- Open an issue on GitHub
- Contact: [your contact information]

## Acknowledgments

This toolkit builds upon several excellent tools:
- odgi for graph manipulation
- minimap2 for alignment
- gafpack for coverage calculation
- biWFA for sequence alignment
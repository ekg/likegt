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

LikeGT requires several bioinformatics tools to be installed and available in your PATH.

**Core Dependencies:**
- `cargo` (>= 1.70) - Rust compiler and package manager
- `odgi` (>= 0.8) - Optimized dynamic genome/graph implementation
- `minimap2` (>= 2.24) - Fast sequence alignment program
- `wgsim` (>= 1.0) - Read simulator for next-generation sequencing
- `gfainject` - Tool to project SAM/BAM alignments onto pangenome graphs
- `gafpack` - Coverage calculator for GAF alignments
- `samtools` (>= 1.17) - SAM/BAM file manipulation
- `seqtk` (>= 1.3) - FASTA/FASTQ processing toolkit

**Additional Dependencies:**
- `allwave` - All-vs-all sequence alignment (required for `max-qv` and `build` commands)
- `seqwish` - Graph inducer from alignments (required for `build` command)
- `bwa` (>= 0.7.17) - Alternative aligner (optional, for `--aligner bwa-mem`)
- `zcat`/`gzip` - For handling compressed files

**Python Dependencies (for analysis scripts):**
- Python 3.8+
- Standard library only (no external packages required)

### Quick Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/likegt.git
cd likegt

# Check and install dependencies
./install_dependencies.sh

# Build release version (recommended)
cargo build --release

# Run tests
cargo test

# Install to cargo bin directory
cargo install --path .
```

### Manual Installation

If you prefer to install dependencies manually:

1. **Install Rust** (if not already installed):
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

2. **Install bioinformatics tools**:

Using conda/mamba (recommended):
```bash
conda create -n likegt python=3.10
conda activate likegt
conda install -c bioconda -c conda-forge \
    minimap2 samtools seqtk bwa odgi seqwish wgsim
```

Using apt (Ubuntu/Debian):
```bash
sudo apt update
sudo apt install minimap2 samtools seqtk bwa
```

3. **Install tools from source** (required for some dependencies):
- odgi: https://github.com/pangenome/odgi
- gfainject: https://github.com/ekg/gfainject  
- gafpack: https://github.com/ekg/gafpack
- allwave: https://github.com/ekg/allwave

4. **Build LikeGT**:
```bash
cargo build --release
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
# Single individual - find best non-self match
likegt max-qv -f sequences.fa.gz -i HG00096 -v

# Multiple individuals
likegt max-qv -f sequences.fa.gz -i "HG00096,HG00171,HG00268" -v

# All individuals in the dataset (shows allwave progress with -v)
likegt max-qv -f sequences.fa.gz -i all -v

# Example verbose output with progress:
# 📊 Found 65 individuals to process
# 🧬 Running allwave for all-vs-all alignment...
# [1.0s] 436/1242 (35.1%) 435.9 alignments/sec ETA: 1.8s
# [2.3s] 609/1242 (49.0%) 263.7 alignments/sec ETA: 2.4s
# ✅ Got 1310 alignment records
# [1/65] Processing HG00096...
#   Best match: HG00268#1 + HG04036#2
#   Max attainable QV: 31.4

# Save results to file
likegt max-qv -f sequences.fa.gz -i all -o max_qv_results.tsv
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

The hold-out validation pipeline evaluates genotyping accuracy by simulating a realistic genotyping scenario where the true haplotypes are excluded from the reference panel.

**How the Pipeline Works:**

1. **Sequence Extraction**: Extracts the target individual's haplotype sequences from the input FASTA
2. **Read Simulation**: Simulates paired-end sequencing reads from the extracted sequences using:
   - `wgsim` for short reads (default: 150bp reads, 30x coverage)
   - Configurable error rates and fragment sizes
3. **Read Alignment**: Maps simulated reads to the pangenome graph:
   - First aligns to linear sequences using `minimap2` or `bwa-mem`
   - Produces SAM/BAM alignment file
4. **Graph Projection**: Projects linear alignments onto the graph:
   - Converts SAM to GAF format using `gfainject`
   - Maps linear coordinates to graph node space
5. **Coverage Calculation**: Computes per-node coverage using `gafpack`:
   - Counts read support for each graph node
   - Generates coverage matrix for the sample
6. **Genotype Calling**: Uses COSIGT algorithm to find best matching genotype:
   - Computes cosine similarity between sample and reference coverage vectors
   - Identifies the pair of reference haplotypes with highest similarity
7. **Quality Assessment**: Calculates QV from cosine similarity:
   - QV = -10 * log10(1 - similarity)
   - Reports whether correct genotype was called

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

### `max-qv` - Maximum Attainable QV using Sequence Alignment

Computes the best possible QV achievable when genotyping held-out individuals by finding their most similar non-self haplotypes using actual sequence alignment with allwave.

**Options:**
- `-f, --fasta`: Input FASTA file with haplotype sequences
- `-i, --individual`: Individual to analyze (name, comma-separated list, or "all")
- `-o, --output`: Output file (optional, defaults to stdout)
- `-t, --threads`: Number of threads for allwave (default: 4)
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

### How max-qv Works

The `max-qv` command answers the question: **"What's the best possible genotyping accuracy we could achieve for a held-out individual using only other individuals in the dataset?"**

**Process:**
1. Runs allwave for all-vs-all sequence alignment
2. For each individual, finds the best matching non-self haplotypes
3. Computes QV from alignment identity: `QV = -10 * log10(1 - identity)`
4. Reports the best achievable genotype (pair of haplotypes) and their QV scores

**Understanding the Output:**
- **QV 20** = 99% sequence identity
- **QV 30** = 99.9% sequence identity  
- **QV 40** = 99.99% sequence identity
- **Typical max QV**: 25-35 for human genome regions

**Example Output:**
```
individual  target_hap1  target_hap2  qv_hap1  qv_hap2  avg_qv  identity_hap1  identity_hap2
HG00096     HG00268#1    HG04036#2    29.1     33.8     31.4    0.9988         0.9996
```
This shows that HG00096's best non-self match would be HG00268's haplotype 1 + HG04036's haplotype 2, achieving an average QV of 31.4.

### Coverage-based vs Sequence-based QV

**hold-out command (Coverage-based QV):**
- QV computed from cosine similarity between graph node coverage vectors
- Measures how well genotyping via graph coverage matches the truth
- Fast and reflects actual genotyping algorithm performance

**max-qv command (Sequence-based QV):**
- QV computed from actual DNA sequence alignment using allwave
- Measures theoretical best possible accuracy based on sequence similarity
- Shows upper bound on genotyping performance

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
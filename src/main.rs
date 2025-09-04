use clap::{Parser, Subcommand};
use anyhow::Result;

#[derive(Parser)]
#[command(name = "likegt")]
#[command(about = "A tool for pangenome graph-based genotyping validation")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run validation tests (hold-0-out, hold-2-out, reference bias)
    Validate {
        /// Input FASTA file with sequences
        #[arg(short, long)]
        fasta: String,
        
        /// Output directory for results
        #[arg(short, long, default_value = "validation_results")]
        output: String,
        
        /// Number of threads to use
        #[arg(short, long, default_value = "4")]
        threads: usize,
        
        /// K-mer size for graph construction
        #[arg(short, long, default_value = "51")]
        kmer_size: usize,
    },
    
    /// Check if a graph is suitable for genotyping
    Check {
        /// Input GFA file
        #[arg(short, long)]
        gfa: String,
        
        /// Output report file
        #[arg(short, long, default_value = "genotyping_report.txt")]
        output: String,
    },
    
    /// Build a pangenome graph from FASTA sequences using allwave + seqwish + odgi
    Build {
        /// Input FASTA file (.fa, .fa.gz)
        #[arg(short, long)]
        fasta: String,
        
        /// Output GFA file (or prefix for multiple k-mer sizes)
        #[arg(short, long)]
        output: String,
        
        /// K-mer size for seqwish (or comma-separated list: 0,25,51,101)
        #[arg(short, long, default_value = "51")]
        kmer_sizes: String,
        
        /// Number of threads
        #[arg(short, long, default_value = "8")]
        threads: usize,
        
        /// Allwave pruning mode (none, low, medium, high)
        #[arg(long, default_value = "none")]
        pruning: String,
        
        /// Generate PNG visualizations with odgi viz
        #[arg(long)]
        visualize: bool,
        
        /// Keep intermediate files (PAF, seqwish GFA)
        #[arg(long)]
        keep_intermediates: bool,
    },
    
    /// Run complete hold-out validation pipeline
    HoldOut {
        /// Input FASTA sequences  
        #[arg(short, long)]
        fasta: String,
        
        /// Input graph file (.gfa or .og)
        #[arg(short, long)]
        graph: String,
        
        /// Output directory for results
        #[arg(short, long, default_value = "hold2out_results")]
        output: String,
        
        /// Test individual to hold out (e.g., "HG00096" or "all" for all samples)
        #[arg(short, long, default_value = "all")]
        individual: String,
        
        /// Number of haplotypes to hold out (1 or 2)
        #[arg(long, default_value = "2")]
        hold: usize,
        
        /// Ploidy (number of haplotypes per individual)
        #[arg(short, long, default_value = "2")]
        ploidy: usize,
        
        /// Number of threads
        #[arg(short, long, default_value = "4")]
        threads: usize,
        
        /// K-mer size
        #[arg(short, long, default_value = "51")]
        kmer_size: usize,
        
        /// Read simulator (wgsim, mason, pbsim3)
        #[arg(long, default_value = "wgsim")]
        simulator: String,
        
        /// Read length for simulation
        #[arg(long, default_value = "150")]
        read_length: usize,
        
        /// Coverage depth for read simulation
        #[arg(long, default_value = "30")]
        coverage_depth: usize,
        
        /// Fragment length mean (paired-end reads)
        #[arg(long, default_value = "500")]
        fragment_length: usize,
        
        /// Fragment length std dev
        #[arg(long, default_value = "50")]
        fragment_std: usize,
        
        /// Aligner (minimap2, bwa-mem, strobealign)
        #[arg(long, default_value = "minimap2")]
        aligner: String,
        
        /// Alignment preset (sr for short reads, map-ont for long reads)
        #[arg(long, default_value = "sr")]
        preset: String,
        
        /// Keep intermediate files
        #[arg(long)]
        keep_files: bool,
        
        /// Extract sequences matching prefix from pangenome FASTA for bias filtering
        /// e.g., "CHM13" or "grch38" - filters reads that align to sequences matching prefix  
        #[arg(long)]
        bias_prefix: Option<String>,
        
        /// External reference FASTA file for bias filtering (e.g., grch38_chr6_pansn.fa.gz)
        /// Uses complete reference FASTA file for initial alignment before filtering
        #[arg(long)]
        bias_fasta: Option<String>,
        
        /// Output sequence-level QV validation
        #[arg(long)]
        sequence_qv: bool,
        
        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
        
        /// Output format (text, json, csv, tsv, table)
        #[arg(long, default_value = "text")]
        format: String,
    },
    
    /// Compute maximum attainable QV using allwave sequence alignment
    MaxQv {
        /// Input FASTA file with sequences
        #[arg(short, long)]
        fasta: String,
        
        /// Individual to analyze (name or "all")
        #[arg(short, long)]
        individual: Option<String>,
        
        /// Output file (optional, defaults to stdout)
        #[arg(short, long)]
        output: Option<String>,
        
        /// Number of threads for allwave
        #[arg(short, long, default_value = "4")]
        threads: usize,
        
        /// Sparsification strategy for allwave (default: tree:5:0:0)
        /// Options: "none" for exact all-vs-all, "tree:N:F:R" for neighbor-based
        #[arg(short = 'p', long, default_value = "tree:5:0:0")]
        sparsification: String,
        
        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
    },
}

#[tokio::main]
async fn main() -> Result<()> {
    env_logger::init();
    
    let cli = Cli::parse();
    
    match cli.command {
        Commands::Validate { fasta, output, threads, kmer_size } => {
            likegt::commands::validate::run_validation(
                &fasta,
                &output, 
                threads,
                kmer_size,
            ).await
        }
        
        Commands::Check { gfa, output } => {
            likegt::commands::check::check_graph_genotyping_suitability(&gfa, &output)
        }
        
        Commands::Build { fasta, output, kmer_sizes, threads, pruning, visualize, keep_intermediates } => {
            likegt::pipeline::build::build_graph_allwave_seqwish(
                &fasta,
                &output,
                &kmer_sizes,
                threads,
                &pruning,
                visualize,
                keep_intermediates,
            ).await
        }
        
        Commands::HoldOut { 
            fasta,
            graph,
            output, 
            individual,
            hold,
            ploidy, 
            threads, 
            kmer_size,
            simulator,
            read_length,
            coverage_depth,
            fragment_length,
            fragment_std,
            aligner,
            preset,
            keep_files,
            bias_prefix,
            bias_fasta,
            sequence_qv,
            verbose,
            format 
        } => {
            // Check if batch processing
            if individual == "all" || individual.contains(',') {
                likegt::commands::hold2out::run_batch_hold2out(
                    &fasta,
                    &graph,
                    &output,
                    &individual,
                    hold,
                    ploidy,
                    threads,
                    kmer_size,
                    &simulator,
                    read_length,
                    coverage_depth,
                    fragment_length,
                    fragment_std,
                    &aligner,
                    &preset,
                    keep_files,
                    bias_prefix.as_deref(),
                    bias_fasta.as_deref(),
                    sequence_qv,
                    verbose,
                    &format,
                ).await
            } else {
                likegt::commands::hold2out::run_complete_hold2out_pipeline(
                    &fasta,
                    &graph,
                    &output,
                    &individual,
                    hold,
                    ploidy,
                    threads,
                    kmer_size,
                    &simulator,
                    read_length,
                    coverage_depth,
                    fragment_length,
                    fragment_std,
                    &aligner,
                    &preset,
                    keep_files,
                    bias_prefix.as_deref(),
                    bias_fasta.as_deref(),
                    sequence_qv,
                    verbose,
                    &format,
                ).await
            }
        }
        
        Commands::MaxQv { fasta, individual, output, threads, sparsification, verbose } => {
            likegt::commands::max_qv::run_max_qv_analysis(
                &fasta,
                individual.as_deref(),
                output.as_deref(),
                threads,
                Some(&sparsification),
                verbose,
            ).await
        }
    }
}
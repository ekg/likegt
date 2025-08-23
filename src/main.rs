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
    
    /// Build a pangenome graph from FASTA sequences
    Build {
        /// Input FASTA file
        #[arg(short, long)]
        fasta: String,
        
        /// Output GFA file
        #[arg(short, long)]
        output: String,
        
        /// K-mer size
        #[arg(short, long, default_value = "51")]
        kmer_size: usize,
        
        /// Segment length for pggb
        #[arg(short, long, default_value = "10000")]
        segment_length: usize,
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
        
        Commands::Build { fasta, output, kmer_size, segment_length } => {
            likegt::pipeline::build::build_graph_from_fasta(
                &fasta,
                &output,
                kmer_size,
                segment_length,
            ).await
        }
    }
}
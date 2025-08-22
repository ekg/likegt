use likegt::{genotype, graph, io, validation, hold2out};

use clap::{Parser, Subcommand};
use anyhow::Result;

#[derive(Parser)]
#[command(name = "likegt")]
#[command(about = "Simplified graph-based genotyping", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Genotype a sample using graph coverage
    Genotype {
        /// Path coverage file from odgi paths (gzip TSV)
        #[arg(short, long)]
        paths: String,
        
        /// Sample coverage file from gafpack (gzip TSV)
        #[arg(short, long)]
        gaf: String,
        
        /// Output directory
        #[arg(short, long)]
        output: String,
        
        /// Sample ID
        #[arg(short, long)]
        id: String,
        
        /// Ploidy level
        #[arg(short = 'n', long, default_value = "2")]
        ploidy: usize,
        
        /// Number of threads
        #[arg(short, long, default_value = "8")]
        threads: usize,
    },
    
    /// Build a pangenome graph from sequences
    Build {
        /// Input FASTA file
        #[arg(short, long)]
        input: String,
        
        /// Output directory
        #[arg(short, long)]
        output: String,
        
        /// K-mer sizes for graph construction
        #[arg(short, long, value_delimiter = ',', default_value = "51")]
        k: Vec<usize>,
        
        /// Number of threads
        #[arg(short, long, default_value = "8")]
        threads: usize,
    },
    
    /// Validate genotyping accuracy with hold-0-out
    Validate {
        /// Path to reference coverage TSV (gzipped)
        #[arg(short, long)]
        paths: String,
        
        /// Output file for validation report
        #[arg(short, long, default_value = "validation_report.txt")]
        output: String,
    },
    
    /// Compare hold-0-out vs hold-2-out validation
    CompareHoldout {
        /// Path to reference coverage TSV (gzipped)
        #[arg(short, long)]
        paths: String,
    },
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    
    match cli.command {
        Commands::Genotype { 
            paths, 
            gaf, 
            output, 
            id, 
            ploidy,
            threads 
        } => {
            genotype::run_genotyping(&paths, &gaf, &output, &id, ploidy, threads)?;
        }
        Commands::Build {
            input,
            output,
            k,
            threads,
        } => {
            let builder = graph::GraphBuilder::new(input, output)
                .threads(threads)
                .k_values(k);
            
            let result = builder.build()?;
            log::info!("Built {} graphs for {}", result.graphs.len(), result.base_name);
            
            for graph in &result.graphs {
                println!("Graph k={}:", graph.k);
                println!("  GFA: {}", graph.gfa_path.display());
                println!("  Paths: {}", graph.paths_path.display());
            }
        }
        Commands::Validate { paths, output } => {
            let ref_data = io::read_gzip_tsv(&paths)?;
            let results = validation::validate_all_individuals(&ref_data)?;
            let report = validation::generate_accuracy_report(&results);
            
            std::fs::write(&output, &report)?;
            println!("{}", report);
            println!("\nReport saved to: {}", output);
        }
        Commands::CompareHoldout { paths } => {
            let ref_data = io::read_gzip_tsv(&paths)?;
            hold2out::compare_validation_methods(&ref_data)?;
        }
    }
    
    Ok(())
}

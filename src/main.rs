mod cami;
mod commands;
mod expression;
mod processing;
mod taxonomy;

use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

use commands::{
    fillup::{self as fillup_cmd, FillupConfig},
    filter::{self as filter_cmd, FilterConfig},
    list::{self as list_cmd, ListConfig},
    preview::{self as preview_cmd, PreviewConfig},
    renorm::{self as renorm_cmd, RenormConfig},
};

#[derive(Parser)]
#[command(author, version, about)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Filter CAMI profiling data
    Filter {
        /// Filter expression, e.g. "r==species & a>=0.1"
        expression: String,
        /// Output file (defaults to stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Fill up missing higher ranks using NCBI taxdump (downloads to ~/.cami)
        #[arg(long)]
        fill_up: bool,
        /// Target rank to fill up to (inclusive)
        #[arg(long, default_value = "phylum")]
        to_rank: String,
        /// Renormalize percentages to 100 per rank after filtering/filling
        #[arg(long)]
        renorm: bool,
        /// Input CAMI file
        #[arg(value_name = "INPUT", index = 2)]
        input: Option<PathBuf>,
    },
    /// List samples and per-rank summaries
    List {
        /// Input CAMI file
        #[arg(value_name = "INPUT", index = 1)]
        input: Option<PathBuf>,
    },
    /// Preview first N entries per sample
    Preview {
        /// Number of entries per sample to show
        #[arg(short = 'n', long, default_value_t = 5)]
        n: usize,
        /// Input CAMI file
        #[arg(value_name = "INPUT", index = 1)]
        input: Option<PathBuf>,
    },
    /// Renormalize abundances to 100 per rank for each sample
    Renorm {
        /// Input CAMI file
        #[arg(value_name = "INPUT", index = 1)]
        input: Option<PathBuf>,
        /// Output file (defaults to stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    /// Fill up samples to populate missing ranks using taxonomy lineages
    Fillup {
        /// Input CAMI file
        #[arg(value_name = "INPUT", index = 1)]
        input: Option<PathBuf>,
        /// Output file (defaults to stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Target rank to fill to (inclusive). Defaults to the highest declared rank
        #[arg(long)]
        to_rank: Option<String>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Filter {
            expression,
            output,
            fill_up,
            to_rank,
            renorm,
            input,
        } => {
            let cfg = FilterConfig {
                expression,
                output: output.as_ref(),
                fill_up: *fill_up,
                to_rank,
                renorm: *renorm,
                input: input.as_ref(),
            };
            filter_cmd::run(&cfg)
        }
        Commands::List { input } => {
            let cfg = ListConfig {
                input: input.as_ref(),
            };
            list_cmd::run(&cfg)
        }
        Commands::Preview { n, input } => {
            let cfg = PreviewConfig {
                n: *n,
                input: input.as_ref(),
            };
            preview_cmd::run(&cfg)
        }
        Commands::Renorm { input, output } => {
            let cfg = RenormConfig {
                input: input.as_ref(),
                output: output.as_ref(),
            };
            renorm_cmd::run(&cfg)
        }
        Commands::Fillup {
            input,
            output,
            to_rank,
        } => {
            let cfg = FillupConfig {
                input: input.as_ref(),
                output: output.as_ref(),
                to_rank: to_rank.as_deref(),
            };
            fillup_cmd::run(&cfg)
        }
    }
}

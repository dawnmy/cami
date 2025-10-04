mod cami;
mod commands;
mod expression;
mod processing;
mod taxonomy;

use anyhow::{Result, bail};
use clap::{Parser, Subcommand};
use std::path::PathBuf;

use commands::{
    fillup::{self as fillup_cmd, FillupConfig},
    filter::{self as filter_cmd, FilterConfig},
    list::{self as list_cmd, ListConfig},
    preview::{self as preview_cmd, PreviewConfig},
    renorm::{self as renorm_cmd, RenormConfig},
    sort::{self as sort_cmd, SortConfig, SortMode, TaxPathField},
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
        /// Source rank to aggregate from when filling up (defaults to species if present)
        #[arg(long = "from")]
        from_rank: Option<String>,
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
        /// Source rank to aggregate from when filling up (defaults to species if present)
        #[arg(long = "from")]
        from_rank: Option<String>,
    },
    /// Sort taxa within each rank for every sample
    Sort {
        /// Sort by abundance (descending) and drop zero-abundance taxa
        #[arg(short = 'a', long, conflicts_with = "taxpath")]
        abundance: bool,
        /// Sort by taxonomy path (TAXPATH or TAXPATHSN)
        #[arg(
            short = 't',
            value_enum,
            num_args = 0..=1,
            default_missing_value = "taxpath"
        )]
        taxpath: Option<TaxPathField>,
        /// Input CAMI file
        #[arg(value_name = "INPUT", index = 1)]
        input: Option<PathBuf>,
        /// Output file (defaults to stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Filter {
            expression,
            output,
            fill_up,
            from_rank,
            to_rank,
            renorm,
            input,
        } => {
            let cfg = FilterConfig {
                expression,
                output: output.as_ref(),
                fill_up: *fill_up,
                from_rank: from_rank.as_deref(),
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
            from_rank,
        } => {
            let cfg = FillupConfig {
                input: input.as_ref(),
                output: output.as_ref(),
                to_rank: to_rank.as_deref(),
                from_rank: from_rank.as_deref(),
            };
            fillup_cmd::run(&cfg)
        }
        Commands::Sort {
            abundance,
            taxpath,
            input,
            output,
        } => {
            let mode = if *abundance {
                SortMode::Abundance
            } else if let Some(field) = taxpath {
                SortMode::TaxPath(*field)
            } else {
                bail!("either -a/--abundance or -t must be provided");
            };
            let cfg = SortConfig {
                input: input.as_ref(),
                output: output.as_ref(),
                mode,
            };
            sort_cmd::run(&cfg)
        }
    }
}

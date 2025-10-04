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
#[command(
    author,
    version,
    about = "Explore and post-process CAMI profiling tables",
    long_about = "cami reads CAMI-format abundance tables and provides commands to inspect, filter, and reshape them.",
    disable_help_subcommand = true
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    #[command(
        about = "Filter CAMI profiling data with logical expressions",
        long_about = "Apply boolean filter expressions to CAMI profiling tables. Combine rank, sample, abundance, taxonomy, and cumulative-sum predicates with & (and), | (or), and parentheses. When --fill-up is provided the command completes missing lineages using the NCBI taxdump before writing a filtered CAMI table."
    )]
    Filter {
        #[arg(
            value_name = "EXPR",
            help = "Filter expression combining rank (r), sample (s), abundance (a), taxonomy (t/tax), and cumsum (c) tests.",
            long_help = concat!(
                "Filter expression combining rank (r), sample (s), abundance (a), taxonomy (t/tax), and cumulative-sum tests. Use & (and), | (or), and parentheses.\n\n",
                "Rank selectors: r==rank, r!=rank, r<=rank (current or more specific), r>=rank (current or more general), and range comparisons follow the order declared in @Ranks.\n",
                "Sample selectors: s==id accepts IDs, 1-based indices, comma-separated lists, or inclusive ranges (e.g. s==1:3); use '.' or ':' for all samples; s!= negates; s~'regex' matches sample IDs with a regular expression.\n",
                "Abundance selectors: a>=value, a<=value, a>value, a<value, a==value, and a!=value compare percentages.\n",
                "Taxonomy selectors: t==taxid, t!=taxid, t<=taxid (ancestor or self), t<taxid (strict ancestor); prefix with ! to invert.\n",
                "Cumulative sums: c<=threshold keeps the least-abundant taxa per rank whose cumulative percentage is within the threshold; prefix with ! to discard them instead."
            )
        )]
        expression: String,
        #[arg(short, long, help = "Write output to a file instead of stdout.")]
        output: Option<PathBuf>,
        #[arg(
            long,
            help = "Fill in missing higher ranks using the NCBI taxdump (downloaded to ~/.cami if absent)."
        )]
        fill_up: bool,
        #[arg(
            long = "from",
            help = "Rank to aggregate from when filling (defaults to species when available)."
        )]
        from_rank: Option<String>,
        #[arg(
            long,
            default_value = "phylum",
            help = "Highest rank to build during fill-up (inclusive)."
        )]
        to_rank: String,
        #[arg(
            long,
            help = "Renormalize each rank to 100% after filtering and filling."
        )]
        renorm: bool,
        #[arg(
            value_name = "INPUT",
            index = 2,
            help = "Input CAMI file (defaults to stdin)."
        )]
        input: Option<PathBuf>,
    },
    #[command(
        about = "List samples and per-rank summaries",
        long_about = "Report each sample's declared ranks, total taxa, and the summed abundance at every rank so you can spot coverage gaps before further processing."
    )]
    List {
        #[arg(
            value_name = "INPUT",
            index = 1,
            help = "Input CAMI file (defaults to stdin)."
        )]
        input: Option<PathBuf>,
    },
    #[command(
        about = "Preview the first entries per sample",
        long_about = "Show the first N taxonomic entries for each sample, preserving the CAMI header format so you can quickly verify parsing and taxonomy strings."
    )]
    Preview {
        #[arg(
            short = 'n',
            long,
            default_value_t = 5,
            help = "Number of entries per sample to display.",
            value_name = "N"
        )]
        n: usize,
        #[arg(
            value_name = "INPUT",
            index = 1,
            help = "Input CAMI file (defaults to stdin)."
        )]
        input: Option<PathBuf>,
    },
    #[command(
        about = "Renormalize abundances per rank",
        long_about = "Scale positive abundances within each rank of every sample so they sum to 100 and round the results to five decimal places."
    )]
    Renorm {
        #[arg(
            value_name = "INPUT",
            index = 1,
            help = "Input CAMI file (defaults to stdin)."
        )]
        input: Option<PathBuf>,
        #[arg(short, long, help = "Write output to a file instead of stdout.")]
        output: Option<PathBuf>,
    },
    #[command(
        about = "Fill missing ranks using taxonomy",
        long_about = "Complete partial lineages in each sample by consulting the NCBI taxdump. Newly created abundances are rounded to five decimal places so the output stays tidy."
    )]
    Fillup {
        #[arg(
            value_name = "INPUT",
            index = 1,
            help = "Input CAMI file (defaults to stdin)."
        )]
        input: Option<PathBuf>,
        #[arg(short, long, help = "Write output to a file instead of stdout.")]
        output: Option<PathBuf>,
        #[arg(
            long,
            help = "Highest rank to fill to (inclusive); defaults to the first declared rank."
        )]
        to_rank: Option<String>,
        #[arg(
            long = "from",
            help = "Rank to aggregate from when filling (defaults to species when available)."
        )]
        from_rank: Option<String>,
    },
    #[command(
        about = "Sort taxa within each rank",
        long_about = "Reorder taxa in every sample either by decreasing abundance (dropping zero-abundance entries) or by their taxonomy paths so related lineages stay adjacent."
    )]
    Sort {
        #[arg(
            short = 'a',
            long,
            conflicts_with = "taxpath",
            help = "Sort taxa by descending abundance and drop entries with zero abundance."
        )]
        abundance: bool,
        #[arg(
            short = 't',
            value_enum,
            num_args = 0..=1,
            default_missing_value = "taxpath",
            help = "Sort taxa by their taxonomy path (TAXPATH or TAXPATHSN)."
        )]
        taxpath: Option<TaxPathField>,
        #[arg(
            value_name = "INPUT",
            index = 1,
            help = "Input CAMI file (defaults to stdin)."
        )]
        input: Option<PathBuf>,
        #[arg(short, long, help = "Write output to a file instead of stdout.")]
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

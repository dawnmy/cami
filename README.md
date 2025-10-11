# cami

`cami` is a command-line companion for working with [CAMI](https://cami-challenge.org/file-formats/) taxonomic profiling tables. It helps you inspect samples, clean and reformat abundances, and prepare subsets for downstream analysis without leaving the terminal.

## Feature overview

- Summarize CAMI files to see which samples, ranks, and taxa are present.
- Preview the top entries of each sample before loading the file into another tool.
- Filter taxa with expressive boolean predicates that reference rank (`r`), sample (`s`), abundance (`a`), taxonomy (`t`/`tax`), and cumulative sums (`c`).
- Fill in missing higher ranks by pulling lineage information from the NCBI taxdump and round abundances to five decimal places.
- Renormalize abundances so that every rank in every sample sums to 100.
- Reorder taxa within each rank, either by abundance (dropping zeroes) or by lineage, to make tables easier to scan.
- Benchmark predicted profiles against ground truth tables with precision/recall, abundance error, correlation, diversity, UniFrac, and abundance-rank metrics (ARE and mARE).

The repository includes a small demo table at [`examples/test.cami`](examples/test.cami) that you can use with the examples below.

## Installation

1. Install [Rust](https://www.rust-lang.org/tools/install) if it is not already available.
2. Clone this repository and build the binary:
   ```bash
   git clone https://github.com/dawnmy/cami.git
   cd cami
   cargo install --path .
   ```
3. Run `cami --help` to confirm the command is available.

You can also invoke subcommands directly with `cargo run -- <command>` while developing.

## Command reference

### `cami list`

Print a per-sample summary that counts how many taxa and how much abundance is assigned to each declared rank.

```text
$ cami list examples/test.cami
Sample: s1
  Ranks: superkingdom, phylum, class, order, family, genus, species, strain
  Total taxa: 18
    superkingdom: taxa=1 total=100.000
    phylum: taxa=2 total=100.000
    class: taxa=3 total=100.000
    order: taxa=3 total=100.000
    family: taxa=3 total=100.000
    genus: taxa=3 total=100.000
    species: taxa=3 total=100.000
    strain: taxa=0 total=0.000
...
```

Use this command to ensure the file has the ranks and coverage you expect before diving into more complex operations.

### `cami preview`

Show the first `N` entries per sample (default `5`). This is handy for spot-checking formatting and verifying that the taxonomy paths look correct.

```text
$ cami preview -n 2 examples/test.cami
@SampleID:s1
@Version:0.10.0
@Ranks:superkingdom|phylum|class|order|family|genus|species|strain
@@TAXID  RANK    TAXPATH              TAXPATHSN                                  PERCENTAGE
2        superkingdom 2                Bacteria                                   100
201174   phylum       2|201174         Bacteria|Actinobacteria                    65.67585
...
```

### `cami filter`

Filter taxa with boolean expressions while optionally filling missing ranks and renormalizing abundances. Results are emitted as a valid CAMI table, so you can chain additional commands or redirect to a file. It is recommended to use single quotation marks instead of double quotes. For sample ID matching, you can enclose the sample ID or pattern in double quotes within the single-quoted expression. If you use `!c`, you must use single quotes for the expression.

Common workflow:

```bash
cami filter --fill-up --renorm 's==s1 & r==species & a>=5' examples/test.cami > enriched.cami
```

This keeps species-level entries from sample `s1` that are at least 5% abundant, fills in any missing higher ranks using the NCBI taxonomy, renormalizes each rank to 100%, and writes the output to `enriched.cami`.

#### Expression syntax

Write expressions with `&` (and), `|` (or), and parentheses. Each atom targets one aspect of the data:

| Atom | Purpose | Operators | Notes |
| ---- | ------- | --------- | ----- |
| `r` or `rank` | Match entry ranks | `==`, `!=`, `<=`, `<`, `>=`, `>` | Uses the order declared by `@Ranks`. `r<=class` keeps class and more specific ranks, while `r>class` keeps more general ranks. Comma-separated lists are allowed with `==`/`!=`. |
| `s` or `sample` | Select samples | `==`, `!=`, `~` | `==` accepts sample IDs, 1-based indices, comma-separated lists, and inclusive ranges (`s==1:3`). `.` or `:` match all samples. Use `s~"regex"` to match IDs with a regular expression. |
| `a` or `abundance` | Compare abundances | `==`, `!=`, `>=`, `>`, `<=`, `<` | Values are interpreted as percentages (0–100). |
| `t` or `tax` | Test lineage membership | `==`, `!=`, `<=`, `<` | Compares against TAXID values. With `--fill-up` or when taxonomy data is available, ancestors are resolved through the NCBI taxdump; otherwise the command inspects `TAXPATH`. Prefix with `!` to negate the result. |
| `c` or `cumsum` | Filter by cumulative abundance | `<=` | Keeps the least-abundant taxa within each rank whose cumulative sum is at most the threshold (again using percentage units). Prefix with `!` to discard those instead. |

Examples:

- `r==species & a>=1` keeps species entries that are at least 1% abundant.
- `s==1,3-5 | s~"^gut"` keeps explicit samples plus any whose IDs start with `gut`.
- `t<=562` keeps entries that fall under *Escherichia coli* (taxid 562) or match the taxid exactly.
- `!c<=2` removes the lowest-abundance taxa per rank whose cumulative total is at most 2%.

When `--fill-up` is supplied, the command downloads the NCBI taxdump (stored under `~/.cami`) if necessary. Use `--from <rank>` to specify which rank to aggregate from when filling and `--to <rank>` to control how far up the lineage to build. Combine `--renorm` to ensure each rank sums to 100 after filtering and filling.

### `cami fillup`

Populate missing higher ranks for every sample using the NCBI taxdump. Abundances retain their full precision after the hierarchy is filled.

```bash
cami fillup --to-rank family examples/test.cami > with_family.cami
```

If `--to-rank` is omitted, the command fills to the highest rank declared in each sample. Use `--from <rank>` to choose the source rank used for aggregation (defaults to `species` when available).

### `cami renorm`

Renormalize abundances so that the percentages at each rank sum to 100 for every sample. Entries with zero or negative abundances are ignored during scaling, and positive values keep their full double-precision values.

```bash
cami renorm examples/test.cami > renormalized.cami
```

### `cami sort`

Reorder taxa within each rank for every sample.

- `--abundance` sorts taxa by descending abundance and removes entries whose abundance is exactly zero.
- `--taxpath [taxpath|taxpathsn]` sorts by the lineage strings so related taxa stay together. The default field is `taxpath` when `-t` is passed without a value.

```bash
cami sort --abundance examples/test.cami > sorted.cami
cami sort --taxpathsn examples/test.cami > lineage_sorted.cami
```

### `cami benchmark`

Evaluate one or more predicted profiles against a ground-truth CAMI table. For every sample and rank the command computes detection metrics (TP/FP/FN, precision/purity, recall/completeness, F1, Jaccard), abundance distances (L1 error, Bray–Curtis), diversity summaries (Shannon index and equitability), Pearson/Spearman correlations, weighted/unweighted UniFrac differences, the Abundance Rank Error (ARE), and the mass-weighted Abundance Rank Error (mARE). Results are written to TSV files so you can load them into spreadsheets or plotting notebooks.

```bash
cami benchmark -g truth.cami predictions/profiler1.cami predictions/profiler2.cami \
  -l "profiler1,profiler2" --af 'a>=0.01' -n --by-domain \
  -o benchmark-results -r "phylum,class,order,family,genus,species"
```

- `-g, --ground-truth` selects the reference CAMI table.
- Positional arguments list predicted profiles to score; provide as many as you like.
- `-l, --labels` (optional) supplies comma-separated names used in the output. When omitted the command derives labels from the file names.
- `--af` applies the same expression filter to both the ground truth and predicted profiles before scoring (e.g., `--af 'a>=0.01'`).
- `--gf` filters the ground-truth profile before scoring using the same expression language as `cami filter`.
- `--pf` applies an expression filter to every predicted profile before metrics are computed.
- `-n, --normalize` rescales each sample/rank in every profile so positive abundances sum to 100 prior to computing metrics.
- `--by-domain` produces additional TSV files restricted to Bacteria, Archaea, Eukarya, and Viruses alongside the overall report.
- `-o, --output` points to the directory where reports such as `benchmark.tsv` and `benchmark_bacteria.tsv` are written.
- `-r, --ranks` restricts the evaluation to specific ranks; mix short forms (`p,c,g`) and full names (`phylum,class,genus`).

Each TSV contains one row per profile/sample/rank combination:

```text
profile   sample   rank     tp  fp  fn  precision  recall   f1        jaccard  l1_error  bray_curtis  shannon_pred  shannon_truth  evenness_pred  evenness_truth  pearson  spearman  weighted_unifrac  unweighted_unifrac  abundance_rank_error  mass_weighted_abundance_rank_error
profiler1 s1       species  42  5   3   0.893617   0.933333  0.913043  0.777778 4.210000  0.021053     2.271111      2.318765       0.932842       0.950112        0.981000 0.975000 0.042000           0.018519              0.052632              0.041875
```

### UniFrac normalization in CAMI benchmark

The `cami benchmark` command reports weighted and unweighted UniFrac scores that are always between 0 and 1. Internally the tool builds a taxonomic tree from the lineages present in the ground-truth and predicted profiles, normalizes the mass present at each lineage tip, and then computes branch-wise discrepancies between the two distributions. To avoid ambiguous superkingdom/domain assignments, the UniFrac implementation only considers the canonical ranks `phylum`, `class`, `order`, `family`, `genus`, `species`, and `strain` when building the comparison tree.

For a given evaluation rank the weighted variant sums the absolute differences in relative mass along every branch down to that rank and divides by the maximum possible distance (placing all mass on mismatching leaves whose lowest common ancestor is the root for the depth being evaluated). The unweighted variant measures how much branch length is unique to either profile by counting the number of phylum-to-rank edges that appear exclusively in the ground truth or the prediction and dividing by the maximum number of such edges given the observed support in each profile. Because every branch is treated as having length one, the reported values can be interpreted as the proportion of disagreement in the shared taxonomy. Missing intermediate ranks do not penalize a tool as long as both profiles share the same descendants—the implementation trims and right-aligns the lineages to a common depth before constructing the tree so that absent ancestors do not inflate the distance.

## Working with the filter language

Expressions can be combined freely, allowing complex workflows:

- **Focus on a cohort:** `cami filter 's~"^trial_" & r<=genus' table.cami`
- **Drop rare tails:** `cami filter '!c<=2' table.cami`
- **Isolate a lineage:** `cami filter 't<=1224 & r>=phylum' table.cami`
- **Chain post-processing:** `cami filter --fill-up --renorm 'r==species & a>=1' table.cami | cami sort --abundance`

Remember that each command reads from stdin when no input path is supplied and writes to stdout by default, making it easy to compose multiple steps.

## Taxonomy data cache

Commands that require lineage information (`filter --fill-up`, `fillup`) download the NCBI taxdump once and cache it under `~/.cami`. Subsequent runs reuse the cached files. You can remove the directory to force a refresh.

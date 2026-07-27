# Comprehensive code review

This review covers the shared CAMI parser/writer, taxonomy loader, expression engine,
processing primitives, and every CLI subcommand in `src/commands`. It distinguishes
correctness defects fixed during the review from lower-risk design and performance
follow-ups that should be addressed separately.

## Correctness fixes made during this review

| Area | Finding | Resolution |
| --- | --- | --- |
| All CAMI-reading commands | Malformed percentages were silently converted to `0`, truncated rows were silently skipped, and data before `@SampleID` was ignored. This could silently change profiles and every downstream result. | Parsing now returns a line-numbered error for malformed, non-finite, or negative percentages, truncated rows, and rows before a sample header. |
| `preview` | `@TaxonomyID` was dropped, so preview output did not preserve all sample metadata. | Preserve the taxonomy tag in preview output. |
| `fillup` and `filter --fill-up` | Filling a bounded rank interval replaced the entire entry list and discarded valid entries above and below that interval. | Preserve entries outside the requested interval while rebuilding entries inside it. |
| `benchmark --ranks` | Requested canonical ranks were compared with raw input rank names. For example, selecting `superkingdom` excluded entries whose equivalent CAMI rank was named `domain`. | Canonicalize the entry rank before applying the rank filter. |
| `benchmark` metrics | Duplicate rows for one taxid were overwritten in metric maps, while totals and detail output included every duplicate. Detection, distance, correlation, and abundance-rank results could therefore disagree with one another. | Aggregate duplicate taxids consistently before metric calculation. |
| `convert --norm` | Unknown/deleted/unmappable rows remained in the normalization denominator, allowing normalized output below 100%. | Normalize only successfully mapped rows. |
| `convert` | Negative and non-finite inputs were accepted, duplicate resolution warnings were printed, and completely unmappable input emitted an empty profile. | Validate abundances, retain the taxonomy resolver's precise warning, and reject empty mapped output. |

## Subcommand-by-subcommand assessment

### `list`

The implementation is linear in the number of entries and groups statistics by the
declared rank index. It intentionally counts only positive-abundance taxa in per-rank
statistics, while `Total taxa` reports physical rows. This difference is documented by
the output labels but may surprise users; a future interface could report both row and
positive-taxon counts explicitly.

### `preview`

The command preserves version, ranks, taxonomy tag, modern extended columns, and the
first requested rows. It currently means "first N rows per sample," not "first N rows
per rank". That behavior matches the implementation and should remain explicit in help
text.

### `filter`

Boolean precedence (`&` before `|`), rank aliases, sample selectors, abundance filters,
taxonomy ancestry filters, and per-rank cumulative sums were reviewed. Cumulative sums
are calculated from least abundant to most abundant, matching the documented low-mass
filter behavior. Two follow-ups remain:

1. Numeric selector parsing currently falls back to zero for malformed values, so a
   typo such as `a>abc` behaves like `a>0` instead of producing an error.
2. Invalid regular expressions evaluate to false rather than returning a diagnostic.

Fixing these cleanly requires making expression evaluation fallible and propagating
errors through `apply_filter`; that API change should be isolated in a follow-up.

### `fillup`

Taxonomy lineage lookup is cached per taxid during each sample fill. Mixed-rank input
adds descendant mass to directly assigned ancestor mass, which produces 90% at species
and 100% at genus for the documented 90% species plus 10% genus-only example. Entries
outside a bounded fill interval are now preserved.

There is an inherent ambiguity when an input already contains a fully aggregated
ancestor and all of its descendants: the file format does not indicate whether the
ancestor value is direct-only mass or an existing aggregate. The current additive
behavior is correct for mixed-rank direct assignments but can double-count a profile
whose ancestor rows are already aggregate totals. A future interface should expose an
explicit `direct` versus `already-aggregated` policy rather than guessing.

### `renorm`

Positive values are independently rescaled to 100 within each sample/rank. Zero values
remain zero. The stricter shared parser now prevents negative and non-finite values from
reaching this operation. Floating-point rounding can still leave a displayed sum a few
units in the final decimal place away from exactly 100; this is normal unless an exact
sum-correction policy is added.

### `sort`

Abundance sorting is descending with taxid as a deterministic tie-breaker and drops
zero-abundance rows as documented. Taxonomy-path sorting preserves zeroes. The same
sorting block is repeated for known, remaining, and unknown rank groups; extracting a
single helper would reduce maintenance risk but does not materially change runtime.

### `convert`

The selected taxdump supplies current taxids, ranks, lineage taxids, and scientific
names. Merged identifiers are resolved transitively. Deleted and unknown identifiers
are omitted with warnings. Conversion and mixed-rank behavior now have fixture-backed
tests covering current identifiers/names, normalization, and invalid abundance input.

### `benchmark`

The following metric implementations were checked:

- Detection metrics use taxid presence at positive abundance: TP, FP, FN, precision,
  recall, F1, and Jaccard follow their standard set formulas.
- L1 and Bray-Curtis compare per-profile relative abundance over the union of observed
  taxids. With normalized compositions, Bray-Curtis is one half of L1.
- Shannon diversity uses natural logarithms and ignores zero mass. Pielou-style
  evenness divides Shannon diversity by the log of the number of positive taxa and is
  undefined for fewer than two positive taxa.
- Pearson uses relative abundance over the union of taxa. Spearman converts values to
  average ranks for ties and then applies Pearson correlation to those ranks.
- Weighted and unweighted UniFrac use unit-length canonical taxonomic edges. Each
  profile is normalized to unit mass before weighted edge-flow differences are
  calculated. Unweighted UniFrac divides unique branch length by union branch length.
  The weighted score is divided by twice the root-to-evaluation-rank path length, as
  documented by this project; this is a project-specific normalization rather than a
  claim of equivalence to every external UniFrac implementation.
- ARE is the project's rank-displacement score with explicit penalties for missed and
  extra taxa. mARE weights normalized rank displacement by abundance and assigns full
  penalty to taxa present in only one profile. Both are bounded to `[0, 1]`.

Duplicate taxids are now aggregated consistently, and canonical rank aliases are
filtered correctly. Remaining benchmark follow-ups are primarily scalability and
definition clarity:

1. With `--by-domain`, every prediction file is parsed, taxonomy-updated, and filtered
   once per domain report. Preprocessing each prediction once before the domain loop
   would substantially reduce work on large inputs.
2. The lineage cache is rebuilt for each profile-map construction. A cache shared across
   domain reports could reduce repeated lineage construction.
3. Prediction-only sample IDs and ranks are not reported because iteration is driven by
   ground-truth samples and ranks. This is defensible for a fixed benchmark design but
   should be explicitly documented; users expecting those rows to count as false
   positives may otherwise misinterpret results.
4. ARE and mARE are project-specific metrics. Their formulas should eventually be
   documented with equations and tie-policy examples, not only prose and unit tests.

## Shared taxonomy and I/O assessment

Taxdump parsing is parallelized across nodes, names, merged ids, and deleted ids.
Lineages and ancestry queries are cached. Potential follow-ups include checking HTTP
status before unpacking a downloaded taxdump, validating malformed dump rows instead of
coercing invalid node ids to zero, and avoiding synchronization-backed caches when a
taxonomy is used strictly from one thread.

The CAMI writer preserves optional modern columns and emits grouped rank aliases. The
reader is now deliberately strict for data integrity; callers receive contextual errors
instead of silently altered profiles.

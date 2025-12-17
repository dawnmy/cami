use crate::cami::{Entry, Sample, open_output, write_cami};
use crate::processing::{fill_up_default, round_percentages};
use crate::taxonomy::{Taxonomy, default_taxdump_dir, ensure_taxdump};
use anyhow::{Context, Result, bail};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

const LEGACY_CAMI_VERSION: &str = "0.10.0";
const MODERN_CAMI_VERSION: &str = "0.11.0";

pub struct ConvertConfig<'a> {
    pub input: Option<&'a PathBuf>,
    pub output: Option<&'a PathBuf>,
    pub taxid_column: usize,
    pub abundance_column: usize,
    pub input_is_percent: bool,
    pub normalize: bool,
    pub sample_id: &'a str,
    pub dmp_dir: Option<&'a PathBuf>,
    pub taxonomy_tag: Option<&'a str>,
}

pub fn run(cfg: &ConvertConfig) -> Result<()> {
    ensure_valid_indices(cfg)?;

    let dir = cfg.dmp_dir.cloned().unwrap_or_else(default_taxdump_dir);
    ensure_taxdump(&dir).with_context(|| format!("ensuring taxdump in {}", dir.display()))?;
    let taxonomy = Taxonomy::load(&dir)?;

    let modern = taxonomy.uses_modern_ranks();
    let version = if modern {
        MODERN_CAMI_VERSION
    } else {
        LEGACY_CAMI_VERSION
    };
    let mut sample = Sample {
        id: cfg.sample_id.to_string(),
        version: Some(version.to_string()),
        taxonomy_tag: cfg.taxonomy_tag.map(|tag| tag.to_string()),
        ranks: Vec::new(),
        rank_groups: Vec::new(),
        rank_aliases: HashMap::new(),
        entries: Vec::new(),
    };
    let groups = if modern {
        modern_rank_groups()
    } else {
        legacy_rank_groups()
    };
    sample.set_rank_groups(groups);

    let records = read_records(cfg)?;
    if records.is_empty() {
        bail!("no data rows found in input TSV");
    }

    let scale = if cfg.input_is_percent { 1.0 } else { 100.0 };
    let total: f64 = records
        .iter()
        .map(|(_, _, abundance)| abundance * scale)
        .sum();
    let norm_factor = if cfg.normalize {
        if total == 0.0 {
            bail!("cannot normalize because total abundance is zero");
        }
        100.0 / total
    } else {
        1.0
    };

    let mut aggregated: HashMap<(String, String), f64> = HashMap::new();
    for (taxid_str, taxid_value, abundance) in records {
        let resolved_taxid = match taxonomy.resolve_taxid(taxid_value) {
            Some(value) => value,
            None => {
                eprintln!(
                    "warning: skipping taxid {taxid_str} because it is not present in the taxonomy"
                );
                continue;
            }
        };

        let mut rank = taxonomy
            .rank_of(resolved_taxid)
            .unwrap_or_else(|| "no rank".to_string());
        if sample.rank_index(&rank).is_none() && rank.eq_ignore_ascii_case("no rank") {
            if sample.rank_index("strain").is_some() {
                rank = "strain".to_string();
            }
        }
        let Some((target_taxid, target_rank)) =
            target_rank_for_entry(&sample, &taxonomy, resolved_taxid, &rank)
        else {
            eprintln!(
                "warning: skipping taxid {taxid_str} because none of its ancestor ranks are present in the @Ranks header"
            );
            continue;
        };

        let percentage = abundance * scale * norm_factor;
        aggregated
            .entry((target_taxid, target_rank))
            .and_modify(|value| *value += percentage)
            .or_insert(percentage);
    }

    let mut aggregated_entries: Vec<_> = aggregated.into_iter().collect();
    aggregated_entries.sort_by(|a, b| a.0.0.cmp(&b.0.0).then_with(|| a.0.1.cmp(&b.0.1)));
    for ((taxid, rank), percentage) in aggregated_entries {
        sample.entries.push(Entry {
            taxid,
            rank,
            taxpath: String::new(),
            taxpathsn: String::new(),
            percentage,
            cami_genome_id: None,
            cami_otu: None,
            hosts: None,
        });
    }

    let mut samples = vec![sample];
    fill_up_default(&mut samples, None, &taxonomy);
    round_percentages(&mut samples);

    let mut out = open_output(cfg.output)?;
    write_cami(&samples, &mut *out)?;
    Ok(())
}

fn ensure_valid_indices(cfg: &ConvertConfig) -> Result<()> {
    if cfg.taxid_column == 0 {
        bail!("taxid column indices are 1-based; got 0");
    }
    if cfg.abundance_column == 0 {
        bail!("abundance column indices are 1-based; got 0");
    }
    Ok(())
}

fn read_records(cfg: &ConvertConfig) -> Result<Vec<(String, u32, f64)>> {
    let taxid_idx = cfg.taxid_column - 1;
    let abundance_idx = cfg.abundance_column - 1;
    let reader = open_reader(cfg.input)?;
    let mut records = Vec::new();
    let mut first = true;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("reading line {}", line_no + 1))?;
        if line.trim().is_empty() {
            continue;
        }
        let mut fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= taxid_idx || fields.len() <= abundance_idx {
            fields = line.split_whitespace().collect();
        }
        if fields.len() <= taxid_idx || fields.len() <= abundance_idx {
            continue;
        }
        let taxid_raw = fields[taxid_idx].trim();
        let abundance_raw = fields[abundance_idx].trim();
        if first {
            first = false;
            if taxid_raw.parse::<u32>().is_err() || abundance_raw.parse::<f64>().is_err() {
                continue;
            }
        }
        let taxid_value: u32 = taxid_raw
            .parse()
            .with_context(|| format!("parsing taxid on line {}", line_no + 1))?;
        let abundance: f64 = abundance_raw
            .parse()
            .with_context(|| format!("parsing abundance on line {}", line_no + 1))?;
        records.push((taxid_raw.to_string(), taxid_value, abundance));
    }

    Ok(records)
}

fn target_rank_for_entry(
    sample: &Sample,
    taxonomy: &Taxonomy,
    taxid: u32,
    rank: &str,
) -> Option<(String, String)> {
    if sample.rank_index(rank).is_some() {
        return Some((taxid.to_string(), rank.to_string()));
    }

    let lineage = taxonomy.lineage(taxid);
    for (ancestor_taxid, ancestor_rank, _) in lineage.iter().rev() {
        if *ancestor_taxid == taxid {
            continue;
        }
        if ancestor_rank.eq_ignore_ascii_case("no rank") {
            continue;
        }
        if sample.rank_index(ancestor_rank).is_some() {
            return Some((ancestor_taxid.to_string(), ancestor_rank.clone()));
        }
    }

    None
}

fn open_reader(input: Option<&PathBuf>) -> Result<Box<dyn BufRead>> {
    if let Some(path) = input {
        if path != Path::new("-") {
            let file = File::open(path)
                .with_context(|| format!("opening input TSV {}", path.display()))?;
            return Ok(Box::new(BufReader::new(file)));
        }
    }
    Ok(Box::new(BufReader::new(io::stdin().lock())))
}

fn legacy_rank_groups() -> Vec<Vec<String>> {
    vec![
        vec!["superkingdom".to_string()],
        vec!["phylum".to_string()],
        vec!["class".to_string()],
        vec!["order".to_string()],
        vec!["family".to_string()],
        vec!["genus".to_string()],
        vec!["species".to_string()],
        vec!["strain".to_string()],
    ]
}

fn modern_rank_groups() -> Vec<Vec<String>> {
    vec![
        vec![
            "cellular root".to_string(),
            "acellular root".to_string(),
            "other entries".to_string(),
        ],
        vec!["domain".to_string(), "realm".to_string()],
        vec!["kingdom".to_string()],
        vec!["phylum".to_string()],
        vec!["class".to_string()],
        vec!["order".to_string()],
        vec!["family".to_string()],
        vec!["genus".to_string()],
        vec!["species".to_string()],
        vec!["strain".to_string()],
    ]
}

use crate::cami::{Entry, Sample, open_output, write_cami};
use crate::processing::{fill_up_default, round_percentages};
use crate::taxonomy::{Taxonomy, default_taxdump_dir, ensure_taxdump};
use anyhow::{Context, Result, bail};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

pub struct ConvertConfig<'a> {
    pub input: Option<&'a PathBuf>,
    pub output: Option<&'a PathBuf>,
    pub taxid_column: usize,
    pub abundance_column: usize,
    pub sample_id: &'a str,
    pub dmp_dir: Option<&'a PathBuf>,
}

pub fn run(cfg: &ConvertConfig) -> Result<()> {
    ensure_valid_indices(cfg)?;

    let dir = cfg.dmp_dir.cloned().unwrap_or_else(default_taxdump_dir);
    ensure_taxdump(&dir).with_context(|| format!("ensuring taxdump in {}", dir.display()))?;
    let taxonomy = Taxonomy::load(&dir)?;

    let modern = taxonomy.uses_modern_ranks();
    let mut sample = Sample {
        id: cfg.sample_id.to_string(),
        version: Some(env!("CARGO_PKG_VERSION").to_string()),
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

    for (taxid_str, taxid_value, abundance) in records {
        let mut rank = taxonomy
            .rank_of(taxid_value)
            .unwrap_or_else(|| "no rank".to_string());
        if sample.rank_index(&rank).is_none() && rank.eq_ignore_ascii_case("no rank") {
            if sample.rank_index("strain").is_some() {
                rank = "strain".to_string();
            }
        }
        if sample.rank_index(&rank).is_none() {
            bail!(
                "rank '{}' for taxid {} is not present in the @Ranks header",
                rank,
                taxid_str
            );
        }
        sample.entries.push(Entry {
            taxid: taxid_str,
            rank,
            taxpath: String::new(),
            taxpathsn: String::new(),
            percentage: abundance,
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

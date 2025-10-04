use crate::cami::{Entry, load_samples, open_output, write_cami};
use anyhow::Result;
use clap::ValueEnum;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

#[derive(Clone, Copy)]
pub enum SortMode {
    Abundance,
    TaxPath(TaxPathField),
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum TaxPathField {
    #[value(alias = "taxpath")]
    Taxpath,
    #[value(alias = "taxpathsn")]
    Taxpathsn,
}

impl TaxPathField {
    fn key<'a>(&self, entry: &'a Entry) -> &'a str {
        match self {
            TaxPathField::Taxpath => &entry.taxpath,
            TaxPathField::Taxpathsn => &entry.taxpathsn,
        }
    }
}

pub struct SortConfig<'a> {
    pub input: Option<&'a PathBuf>,
    pub output: Option<&'a PathBuf>,
    pub mode: SortMode,
}

pub fn run(cfg: &SortConfig) -> Result<()> {
    let mut samples = load_samples(cfg.input)?;

    for sample in samples.iter_mut() {
        let entries = std::mem::take(&mut sample.entries);
        let mut by_rank: HashMap<String, Vec<Entry>> = HashMap::new();
        for entry in entries {
            by_rank.entry(entry.rank.clone()).or_default().push(entry);
        }

        let mut ordered_ranks = Vec::new();
        let mut seen: HashSet<String> = HashSet::new();
        for rank in &sample.ranks {
            if seen.insert(rank.clone()) {
                ordered_ranks.push(rank.clone());
            }
        }
        for rank in by_rank.keys() {
            if seen.insert(rank.clone()) {
                ordered_ranks.push(rank.clone());
            }
        }

        let mut new_entries = Vec::new();
        for rank in ordered_ranks {
            if let Some(mut entries) = by_rank.remove(&rank) {
                match cfg.mode {
                    SortMode::Abundance => {
                        entries.retain(|e| e.percentage > 0.0);
                        entries.sort_by(|a, b| {
                            match b
                                .percentage
                                .partial_cmp(&a.percentage)
                                .unwrap_or(Ordering::Equal)
                            {
                                Ordering::Equal => a.taxid.cmp(&b.taxid),
                                other => other,
                            }
                        });
                    }
                    SortMode::TaxPath(field) => {
                        entries.sort_by(|a, b| {
                            let key_a = field.key(a);
                            let key_b = field.key(b);
                            key_a.cmp(key_b).then_with(|| a.taxid.cmp(&b.taxid))
                        });
                    }
                }
                new_entries.extend(entries);
            }
        }
        sample.entries = new_entries;
    }

    let mut out = open_output(cfg.output)?;
    write_cami(&samples, &mut *out)?;
    Ok(())
}

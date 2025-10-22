use crate::cami::{Entry, load_samples, open_output, write_cami};
use anyhow::Result;
use clap::ValueEnum;
use std::cmp::Ordering;
use std::collections::HashMap;
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
        let mut by_rank: HashMap<usize, Vec<Entry>> = HashMap::new();
        let mut unknown_order: Vec<String> = Vec::new();
        let mut unknown_map: HashMap<String, Vec<Entry>> = HashMap::new();
        for entry in entries {
            if let Some(idx) = sample.rank_index(&entry.rank) {
                by_rank.entry(idx).or_default().push(entry);
            } else {
                let key = entry.rank.clone();
                if !unknown_map.contains_key(&key) {
                    unknown_order.push(key.clone());
                }
                unknown_map.entry(key).or_default().push(entry);
            }
        }

        let mut new_entries = Vec::new();
        for idx in 0..sample.ranks.len() {
            if let Some(mut entries) = by_rank.remove(&idx) {
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
        let mut remaining_known: Vec<(usize, Vec<Entry>)> = by_rank.into_iter().collect();
        remaining_known.sort_by_key(|(idx, _)| *idx);
        for (_idx, mut entries) in remaining_known {
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
        for rank in unknown_order {
            if let Some(mut entries) = unknown_map.remove(&rank) {
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

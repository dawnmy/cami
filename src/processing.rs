use crate::cami::{Entry, Sample};
use crate::taxonomy::{Taxonomy, parse_taxid};
use std::collections::HashMap;

pub fn renormalize(samples: &mut [Sample]) {
    for sample in samples.iter_mut() {
        let mut by_rank: HashMap<String, Vec<usize>> = HashMap::new();
        for (idx, entry) in sample.entries.iter().enumerate() {
            by_rank.entry(entry.rank.clone()).or_default().push(idx);
        }
        for idxs in by_rank.values() {
            let sum: f64 = idxs
                .iter()
                .map(|&i| {
                    let v = sample.entries[i].percentage;
                    if v > 0.0 { v } else { 0.0 }
                })
                .sum();
            if sum <= 0.0 {
                continue;
            }
            for &i in idxs {
                if sample.entries[i].percentage > 0.0 {
                    sample.entries[i].percentage = sample.entries[i].percentage / sum * 100.0;
                }
            }
        }
    }
}

pub fn fill_up_to(samples: &mut [Sample], to_rank: &str, taxonomy: &Taxonomy) {
    for sample in samples.iter_mut() {
        if sample.entries.is_empty() {
            continue;
        }
        let Some(target_idx) = sample.rank_index(to_rank) else {
            continue;
        };
        let mut accum: HashMap<(String, String), f64> = HashMap::new();
        let mut cache: HashMap<String, HashMap<String, (String, String)>> = HashMap::new();

        for entry in &sample.entries {
            let Some(entry_rank_idx) = sample.rank_index(&entry.rank) else {
                continue;
            };
            let Some(rank_map) = rank_map_for(sample, taxonomy, &entry.taxid, &mut cache) else {
                continue;
            };
            let start = entry_rank_idx.min(target_idx);
            let end = entry_rank_idx.max(target_idx);
            let mut targets = Vec::new();
            for rank in &sample.ranks[start..=end] {
                if let Some((tid, _name)) = rank_map.get(rank) {
                    targets.push((rank.clone(), tid.clone()));
                }
            }
            for (rank_name, taxid_for_rank) in targets {
                let key = (taxid_for_rank.clone(), rank_name.clone());
                *accum.entry(key).or_insert(0.0) += entry.percentage;
                let _ = rank_map_for(sample, taxonomy, taxid_for_rank.as_str(), &mut cache);
            }
        }

        let mut new_entries: Vec<Entry> = Vec::new();
        for rank in &sample.ranks {
            let mut rank_entries: Vec<(String, f64)> = accum
                .iter()
                .filter_map(|((taxid, r), pct)| {
                    if r == rank {
                        Some((taxid.clone(), *pct))
                    } else {
                        None
                    }
                })
                .collect();
            rank_entries.sort_by(|a, b| a.0.cmp(&b.0));
            for (taxid, pct) in rank_entries {
                if pct == 0.0 {
                    continue;
                }
                let Some(rank_map) = rank_map_for(sample, taxonomy, &taxid, &mut cache) else {
                    continue;
                };
                let Some((_tid, _name)) = rank_map.get(rank) else {
                    continue;
                };
                let (taxpath, taxpathsn) = build_paths(sample, &rank_map, rank);
                new_entries.push(Entry {
                    taxid,
                    rank: rank.clone(),
                    taxpath,
                    taxpathsn,
                    percentage: pct,
                });
            }
        }
        sample.entries = new_entries;
    }
}

pub fn fill_up_default(samples: &mut [Sample], taxonomy: &Taxonomy) {
    for sample in samples.iter_mut() {
        if sample.entries.is_empty() {
            continue;
        }
        if let Some(target_rank) = sample.ranks.first().cloned() {
            fill_up_to(std::slice::from_mut(sample), &target_rank, taxonomy);
        }
    }
}

fn rank_map_for<'a>(
    sample: &Sample,
    taxonomy: &Taxonomy,
    taxid: &str,
    cache: &'a mut HashMap<String, HashMap<String, (String, String)>>,
) -> Option<&'a HashMap<String, (String, String)>> {
    if !cache.contains_key(taxid) {
        if let Some(map) = build_rank_map(sample, taxonomy, taxid) {
            cache.insert(taxid.to_string(), map);
        } else {
            cache.insert(taxid.to_string(), HashMap::new());
        }
    }
    let map = cache.get(taxid)?;
    if map.is_empty() { None } else { Some(map) }
}

fn build_rank_map(
    sample: &Sample,
    taxonomy: &Taxonomy,
    taxid: &str,
) -> Option<HashMap<String, (String, String)>> {
    let tid = parse_taxid(taxid)?;
    let lineage = taxonomy.lineage(tid);
    let mut map = HashMap::new();
    for (tid_u32, rank, name) in lineage {
        if sample.ranks.iter().any(|r| r == &rank) {
            map.insert(rank, (tid_u32.to_string(), name));
        }
    }
    if map.is_empty() { None } else { Some(map) }
}

fn build_paths(
    sample: &Sample,
    rank_map: &HashMap<String, (String, String)>,
    upto_rank: &str,
) -> (String, String) {
    let mut taxids = Vec::new();
    let mut names = Vec::new();
    for rank in &sample.ranks {
        if let Some((tid, name)) = rank_map.get(rank) {
            taxids.push(tid.clone());
            names.push(name.clone());
        }
        if rank == upto_rank {
            break;
        }
    }
    (taxids.join("|"), names.join("|"))
}

use crate::cami::{Entry, Sample};
use crate::taxonomy::{Taxonomy, parse_taxid};
use std::collections::{BTreeSet, HashMap};

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

pub fn fill_up_to(
    samples: &mut [Sample],
    from_rank: Option<&str>,
    to_rank: &str,
    taxonomy: &Taxonomy,
) {
    for sample in samples.iter_mut() {
        if sample.entries.is_empty() {
            continue;
        }
        let Some(target_idx) = sample.rank_index(to_rank) else {
            continue;
        };
        let Some(base_idx) = select_base_rank(sample, from_rank) else {
            continue;
        };
        if base_idx < target_idx {
            continue;
        }

        let start_idx = target_idx.min(base_idx);
        let end_idx = target_idx.max(base_idx);

        let mut existing_by_key: HashMap<(usize, String), Entry> = HashMap::new();
        let mut existing_by_rank: HashMap<usize, Vec<String>> = HashMap::new();
        for entry in &sample.entries {
            if let Some(idx) = sample.rank_index(&entry.rank) {
                existing_by_key
                    .entry((idx, entry.taxid.clone()))
                    .or_insert_with(|| entry.clone());
                existing_by_rank
                    .entry(idx)
                    .or_default()
                    .push(entry.taxid.clone());
            }
        }

        let mut sums: HashMap<(usize, String), f64> = HashMap::new();
        let mut cache: HashMap<String, HashMap<String, (String, String)>> = HashMap::new();

        for entry in &sample.entries {
            let Some(entry_rank_idx) = sample.rank_index(&entry.rank) else {
                continue;
            };
            if entry_rank_idx != base_idx {
                continue;
            }
            let fallback_entry = existing_by_key.get(&(entry_rank_idx, entry.taxid.clone()));
            let Some(rank_map) = rank_map_for(sample, taxonomy, &entry.taxid, &mut cache, || {
                fallback_entry.and_then(|e| fallback_rank_map(e, sample))
            }) else {
                continue;
            };

            for rank_idx in start_idx..=end_idx {
                if rank_idx == entry_rank_idx {
                    continue;
                }
                let rank = &sample.ranks[rank_idx];
                if let Some((tid, _name)) = rank_map.get(rank) {
                    if existing_by_key.contains_key(&(rank_idx, tid.clone())) {
                        continue;
                    }
                    *sums.entry((rank_idx, tid.clone())).or_insert(0.0) += entry.percentage;
                }
            }
        }

        let mut new_entries: Vec<Entry> = Vec::new();
        for (idx, rank) in sample.ranks.iter().enumerate() {
            let mut taxids: BTreeSet<String> = BTreeSet::new();
            if let Some(existing) = existing_by_rank.get(&idx) {
                for taxid in existing {
                    taxids.insert(taxid.clone());
                }
            }
            for ((s_idx, taxid), pct) in &sums {
                if *s_idx == idx && *pct > 0.0 {
                    taxids.insert(taxid.clone());
                }
            }

            for taxid in taxids {
                let existing_entry = existing_by_key.get(&(idx, taxid.clone()));
                let percentage = existing_entry
                    .map(|e| e.percentage)
                    .unwrap_or_else(|| *sums.get(&(idx, taxid.clone())).unwrap_or(&0.0));
                if percentage <= 0.0 {
                    continue;
                }

                let Some(rank_map) = rank_map_for(sample, taxonomy, &taxid, &mut cache, || {
                    existing_entry.and_then(|e| fallback_rank_map(e, sample))
                }) else {
                    if let Some(entry) = existing_entry {
                        new_entries.push(entry.clone());
                    }
                    continue;
                };
                if !rank_map.contains_key(rank) {
                    if let Some(entry) = existing_entry {
                        new_entries.push(entry.clone());
                    }
                    continue;
                }
                let (taxpath, taxpathsn) = build_paths(sample, &rank_map, rank);
                new_entries.push(Entry {
                    taxid: taxid.clone(),
                    rank: rank.clone(),
                    taxpath,
                    taxpathsn,
                    percentage,
                });
            }
        }

        sample.entries = new_entries;
    }
}

pub fn fill_up_default(samples: &mut [Sample], from_rank: Option<&str>, taxonomy: &Taxonomy) {
    for sample in samples.iter_mut() {
        if sample.entries.is_empty() {
            continue;
        }
        if let Some(target_rank) = sample.ranks.first().cloned() {
            fill_up_to(
                std::slice::from_mut(sample),
                from_rank,
                &target_rank,
                taxonomy,
            );
        }
    }
}

fn select_base_rank(sample: &Sample, from_rank: Option<&str>) -> Option<usize> {
    if let Some(rank) = from_rank {
        if let Some(idx) = sample.rank_index(rank) {
            if has_entries_at(sample, &sample.ranks[idx]) {
                return Some(idx);
            }
        }
    }
    if let Some(idx) = sample.rank_index("species") {
        if has_entries_at(sample, &sample.ranks[idx]) {
            return Some(idx);
        }
    }
    sample
        .ranks
        .iter()
        .enumerate()
        .rev()
        .find(|(_, rank)| has_entries_at(sample, rank))
        .map(|(idx, _)| idx)
}

fn has_entries_at(sample: &Sample, rank: &str) -> bool {
    sample.entries.iter().any(|e| e.rank == rank)
}

fn rank_map_for<'a, F>(
    sample: &Sample,
    taxonomy: &Taxonomy,
    taxid: &str,
    cache: &'a mut HashMap<String, HashMap<String, (String, String)>>,
    mut fallback: F,
) -> Option<&'a HashMap<String, (String, String)>>
where
    F: FnMut() -> Option<HashMap<String, (String, String)>>,
{
    if !cache.contains_key(taxid) {
        if let Some(map) = build_rank_map(sample, taxonomy, taxid) {
            cache.insert(taxid.to_string(), map);
        } else if let Some(map) = fallback() {
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

fn fallback_rank_map(entry: &Entry, sample: &Sample) -> Option<HashMap<String, (String, String)>> {
    let ids: Vec<&str> = entry.taxpath.split('|').collect();
    let names: Vec<&str> = entry.taxpathsn.split('|').collect();
    if ids.is_empty() || names.is_empty() {
        return None;
    }
    let mut map = HashMap::new();
    let upto = ids.len().min(names.len()).min(sample.ranks.len());
    for idx in 0..upto {
        let rank = &sample.ranks[idx];
        let taxid = ids[idx].trim();
        let name = names[idx].trim();
        if !taxid.is_empty() {
            map.insert(rank.clone(), (taxid.to_string(), name.to_string()));
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

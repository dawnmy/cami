use std::cmp::Ordering;
use std::collections::{BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

use anyhow::{Context, Result, bail, ensure};

use crate::cami;
use crate::cami::{Entry, Sample};
use crate::expression::{apply_filter, expr_needs_taxdump, parse_expression};
use crate::taxonomy::{Taxonomy, ensure_taxdump, parse_taxid};

#[derive(Clone)]
pub struct BenchmarkConfig {
    pub ground_truth: PathBuf,
    pub predictions: Vec<PathBuf>,
    pub labels: Vec<String>,
    pub all_filter: Option<String>,
    pub ground_filter: Option<String>,
    pub pred_filter: Option<String>,
    pub normalize: bool,
    pub by_domain: bool,
    pub output: PathBuf,
    pub ranks: Option<Vec<String>>,
}

#[derive(Clone)]
struct ProfileEntry {
    taxid: String,
    rank: String,
    percentage: f64,
    taxpath: String,
    taxpathsn: String,
}

pub fn run(cfg: &BenchmarkConfig) -> Result<()> {
    if cfg.predictions.is_empty() {
        bail!("at least one predicted profile must be provided");
    }

    if !cfg.labels.is_empty() && cfg.labels.len() != cfg.predictions.len() {
        bail!("labels count must match number of predicted profiles");
    }

    let all_expr = cfg
        .all_filter
        .as_ref()
        .map(|s| parse_expression(s))
        .transpose()
        .context("parsing global filter expression")?;
    let ground_expr = cfg
        .ground_filter
        .as_ref()
        .map(|s| parse_expression(s))
        .transpose()
        .context("parsing ground-truth filter expression")?;
    let pred_expr = cfg
        .pred_filter
        .as_ref()
        .map(|s| parse_expression(s))
        .transpose()
        .context("parsing predicted-profile filter expression")?;

    let needs_taxdump = cfg.by_domain
        || all_expr
            .as_ref()
            .is_some_and(|expr| expr_needs_taxdump(expr))
        || ground_expr
            .as_ref()
            .is_some_and(|expr| expr_needs_taxdump(expr))
        || pred_expr
            .as_ref()
            .is_some_and(|expr| expr_needs_taxdump(expr));

    let taxonomy = if needs_taxdump {
        let dir = taxonomy_dir();
        ensure_taxdump(&dir).with_context(|| format!("ensuring taxdump in {}", dir.display()))?;
        Some(Taxonomy::load(&dir)?)
    } else {
        None
    };

    let mut gt_samples = cami::parse_cami(&cfg.ground_truth)?;
    if gt_samples.is_empty() {
        bail!("ground truth profile has no samples");
    }

    if let Some(expr) = all_expr.as_ref() {
        gt_samples = apply_filter(&gt_samples, expr, taxonomy.as_ref());
    }

    if let Some(expr) = ground_expr.as_ref() {
        gt_samples = apply_filter(&gt_samples, expr, taxonomy.as_ref());
        if gt_samples.is_empty() {
            bail!("ground truth profile has no samples after applying filters");
        }
    } else if gt_samples.is_empty() {
        bail!("ground truth profile has no samples after applying filters");
    }

    fs::create_dir_all(&cfg.output)
        .with_context(|| format!("creating output directory {}", cfg.output.display()))?;

    let ranks = cfg.ranks.as_ref().map(|v| canonical_ranks(v)).transpose()?;

    let labels: Vec<String> = if cfg.labels.is_empty() {
        cfg.predictions
            .iter()
            .map(|p| {
                p.file_stem()
                    .and_then(|s| s.to_str())
                    .map(|s| s.to_string())
                    .unwrap_or_else(|| p.to_string_lossy().into_owned())
            })
            .collect()
    } else {
        cfg.labels.clone()
    };

    let domains = if cfg.by_domain {
        vec![
            None,
            Some("Bacteria".to_string()),
            Some("Archaea".to_string()),
            Some("Eukarya".to_string()),
            Some("Viruses".to_string()),
        ]
    } else {
        vec![None]
    };

    for domain in domains {
        let suffix = domain
            .as_ref()
            .map(|d| format!("_{}", d.to_lowercase()))
            .unwrap_or_else(|| "".to_string());
        let path = cfg.output.join(format!("benchmark{}.tsv", suffix));
        let mut writer = File::create(&path)
            .with_context(|| format!("creating benchmark report {}", path.display()))?;
        writeln!(
            writer,
            "profile\tsample\trank\ttp\tfp\tfn\tprecision\trecall\tf1\tjaccard\tl1_error\tbray_curtis\tshannon_pred\tshannon_truth\tevenness_pred\tevenness_truth\tpearson\tspearman\tweighted_unifrac\tunweighted_unifrac\tabundance_rank_error\tmass_weighted_abundance_rank_error"
        )?;

        let gt_map = build_profile_map(
            &gt_samples,
            ranks.as_ref(),
            cfg.normalize,
            domain.as_deref(),
            taxonomy.as_ref(),
        );

        let mut sample_ids: Vec<_> = gt_map.keys().cloned().collect();
        sample_ids.sort();

        for (pred_path, label) in cfg.predictions.iter().zip(labels.iter()) {
            let mut pred_samples = cami::parse_cami(pred_path)?;
            if let Some(expr) = all_expr.as_ref() {
                pred_samples = apply_filter(&pred_samples, expr, taxonomy.as_ref());
            }
            if let Some(expr) = pred_expr.as_ref() {
                pred_samples = apply_filter(&pred_samples, expr, taxonomy.as_ref());
            }
            let pred_map = build_profile_map(
                &pred_samples,
                ranks.as_ref(),
                cfg.normalize,
                domain.as_deref(),
                taxonomy.as_ref(),
            );

            for sample_id in &sample_ids {
                let Some(gt_ranks) = gt_map.get(sample_id) else {
                    continue;
                };
                let mut rank_names: Vec<_> = gt_ranks.keys().cloned().collect();
                rank_names.sort();
                let pred_ranks = pred_map.get(sample_id);
                for rank in rank_names {
                    if rank_is_above_phylum(&rank) {
                        continue;
                    }
                    let gt_entries = gt_ranks.get(&rank).unwrap();
                    let pred_entries = pred_ranks
                        .and_then(|map| map.get(&rank))
                        .map(|v| v.as_slice())
                        .unwrap_or(&[]);
                    let metrics = compute_metrics(&rank, gt_entries, pred_entries)?;
                    write_metrics(&mut writer, label, sample_id, &rank, &metrics)?;
                }
            }
        }
    }

    Ok(())
}

fn canonical_ranks(ranks: &[String]) -> Result<Vec<String>> {
    ranks.iter().map(|r| canonical_rank(r)).collect()
}

fn rank_is_above_phylum(rank: &str) -> bool {
    matches!(
        rank.trim().to_lowercase().as_str(),
        "superkingdom" | "domain" | "kingdom" | "realm" | "cellular root" | "acellular root"
    )
}

fn canonical_rank(rank: &str) -> Result<String> {
    let trimmed = rank.trim();
    if trimmed.is_empty() {
        bail!("rank names must not be empty");
    }
    let lower = trimmed.to_lowercase();
    let canonical = match lower.as_str() {
        "d" | "domain" | "superkingdom" | "realm" | "cellular root" | "acellular root" => {
            "superkingdom"
        }
        "k" | "kingdom" => "kingdom",
        "p" | "phylum" => "phylum",
        "c" | "class" => "class",
        "o" | "order" => "order",
        "f" | "family" => "family",
        "g" | "genus" => "genus",
        "s" | "species" => "species",
        "t" | "strain" => "strain",
        other => other,
    };
    Ok(canonical.to_string())
}

fn build_profile_map(
    samples: &[Sample],
    ranks: Option<&Vec<String>>,
    normalize: bool,
    domain: Option<&str>,
    taxonomy: Option<&Taxonomy>,
) -> HashMap<String, HashMap<String, Vec<ProfileEntry>>> {
    let mut map: HashMap<String, HashMap<String, Vec<ProfileEntry>>> = HashMap::new();
    let rank_filter: Option<BTreeSet<String>> = ranks.map(|r| r.iter().cloned().collect());
    let domain_lower = domain.map(|d| d.to_lowercase());

    for sample in samples {
        let mut sample_map: HashMap<String, Vec<ProfileEntry>> = HashMap::new();
        for entry in &sample.entries {
            if entry.percentage <= 0.0 {
                continue;
            }
            if let Some(filter) = &rank_filter {
                if !filter.contains(&entry.rank) {
                    continue;
                }
            }
            if let Some(domain) = &domain_lower {
                if !entry_belongs_to_domain(entry, domain, taxonomy) {
                    continue;
                }
            }
            let rank_name = canonical_rank(&entry.rank).unwrap_or_else(|_| entry.rank.clone());
            sample_map
                .entry(rank_name.clone())
                .or_default()
                .push(ProfileEntry {
                    taxid: entry.taxid.clone(),
                    rank: rank_name,
                    percentage: entry.percentage,
                    taxpath: entry.taxpath.clone(),
                    taxpathsn: entry.taxpathsn.clone(),
                });
        }

        if normalize {
            for entries in sample_map.values_mut() {
                let sum: f64 = entries.iter().map(|e| e.percentage).sum();
                if sum > 0.0 {
                    for entry in entries.iter_mut() {
                        entry.percentage = entry.percentage / sum * 100.0;
                    }
                }
            }
        }

        if !sample_map.is_empty() {
            map.insert(sample.id.clone(), sample_map);
        }
    }

    map
}

fn entry_belongs_to_domain(entry: &Entry, domain: &str, taxonomy: Option<&Taxonomy>) -> bool {
    if entry.taxpathsn.split('|').any(|name| {
        let trimmed = name.trim();
        !trimmed.is_empty() && trimmed.eq_ignore_ascii_case(domain)
    }) {
        return true;
    }

    if let Some(tax) = taxonomy {
        if let Some(tid) = highest_taxid_in_taxpath(entry) {
            if let Some(name) = tax.domain_of(tid) {
                if name.eq_ignore_ascii_case(domain) {
                    return true;
                }
            }
        }
        if let Some(tid) = parse_taxid(&entry.taxid) {
            if let Some(name) = tax.domain_of(tid) {
                if name.eq_ignore_ascii_case(domain) {
                    return true;
                }
            }
        }
    }

    false
}

fn highest_taxid_in_taxpath(entry: &Entry) -> Option<u32> {
    entry
        .taxpath
        .split('|')
        .map(|tid| tid.trim())
        .find_map(|tid| {
            if tid.is_empty() {
                None
            } else {
                parse_taxid(tid).filter(|id| *id > 1)
            }
        })
        .or_else(|| parse_taxid(&entry.taxid))
}

#[derive(Default)]
struct Metrics {
    tp: usize,
    fp: usize,
    fn_: usize,
    precision: Option<f64>,
    recall: Option<f64>,
    f1: Option<f64>,
    jaccard: Option<f64>,
    l1_error: f64,
    bray_curtis: Option<f64>,
    shannon_pred: f64,
    shannon_truth: f64,
    evenness_pred: Option<f64>,
    evenness_truth: Option<f64>,
    pearson: Option<f64>,
    spearman: Option<f64>,
    weighted_unifrac: Option<f64>,
    unweighted_unifrac: Option<f64>,
    abundance_rank_error: Option<f64>,
    mass_weighted_abundance_rank_error: Option<f64>,
}

fn compute_metrics(
    rank: &str,
    gt_entries: &[ProfileEntry],
    pred_entries: &[ProfileEntry],
) -> Result<Metrics> {
    let mut metrics = Metrics::default();
    let mut gt_map: HashMap<String, f64> = HashMap::new();
    for entry in gt_entries {
        ensure!(
            !entry.percentage.is_nan(),
            "ground truth abundance for taxid {} is NaN",
            entry.taxid
        );
        ensure!(
            entry.percentage >= 0.0,
            "ground truth abundance for taxid {} is negative",
            entry.taxid
        );
        if entry.percentage > 0.0 {
            gt_map.insert(entry.taxid.clone(), entry.percentage);
        }
    }
    let mut pred_map: HashMap<String, f64> = HashMap::new();
    for entry in pred_entries {
        ensure!(
            !entry.percentage.is_nan(),
            "predicted abundance for taxid {} is NaN",
            entry.taxid
        );
        ensure!(
            entry.percentage >= 0.0,
            "predicted abundance for taxid {} is negative",
            entry.taxid
        );
        if entry.percentage > 0.0 {
            pred_map.insert(entry.taxid.clone(), entry.percentage);
        }
    }

    let gt_total: f64 = gt_entries.iter().map(|e| e.percentage).sum();
    let pred_total: f64 = pred_entries.iter().map(|e| e.percentage).sum();

    let mut union: BTreeSet<String> = BTreeSet::new();
    for taxid in gt_map.keys() {
        union.insert(taxid.clone());
    }
    for taxid in pred_map.keys() {
        union.insert(taxid.clone());
    }

    let mut gt_values: Vec<f64> = Vec::new();
    let mut pred_values: Vec<f64> = Vec::new();

    let mut sum_abs = 0.0;
    let mut sum_tot = 0.0;

    for taxid in &union {
        let g = gt_map.get(taxid).copied().unwrap_or(0.0);
        let p = pred_map.get(taxid).copied().unwrap_or(0.0);
        if g > 0.0 && p > 0.0 {
            metrics.tp += 1;
        } else if p > 0.0 {
            metrics.fp += 1;
        } else if g > 0.0 {
            metrics.fn_ += 1;
        }

        let g_norm = if gt_total > 0.0 { g / gt_total } else { 0.0 };
        let p_norm = if pred_total > 0.0 {
            p / pred_total
        } else {
            0.0
        };

        sum_abs += (g_norm - p_norm).abs();
        sum_tot += g_norm + p_norm;
        gt_values.push(g_norm);
        pred_values.push(p_norm);
    }

    metrics.l1_error = sum_abs;
    metrics.precision = ratio(metrics.tp as f64, (metrics.tp + metrics.fp) as f64);
    metrics.recall = ratio(metrics.tp as f64, (metrics.tp + metrics.fn_) as f64);
    metrics.f1 = match (metrics.precision, metrics.recall) {
        (Some(p), Some(r)) if p + r > 0.0 => Some(2.0 * p * r / (p + r)),
        _ => None,
    };
    metrics.jaccard = ratio(
        metrics.tp as f64,
        (metrics.tp + metrics.fp + metrics.fn_) as f64,
    );
    metrics.bray_curtis = if sum_tot > 0.0 {
        Some(sum_abs / sum_tot)
    } else {
        None
    };

    metrics.shannon_truth = shannon(&gt_values);
    metrics.shannon_pred = shannon(&pred_values);
    metrics.evenness_truth = evenness(&gt_values, metrics.shannon_truth);
    metrics.evenness_pred = evenness(&pred_values, metrics.shannon_pred);
    metrics.pearson = pearson(&gt_values, &pred_values);
    metrics.spearman = spearman(&gt_values, &pred_values);

    let unifrac = unifrac(rank, gt_entries, pred_entries);
    metrics.weighted_unifrac = unifrac.weighted;
    metrics.unweighted_unifrac = unifrac.unweighted;
    metrics.abundance_rank_error = Some(abundance_rank_error(&gt_map, &pred_map)?);
    metrics.mass_weighted_abundance_rank_error =
        Some(mass_weighted_abundance_rank_error(&gt_map, &pred_map, 1.0)?);

    Ok(metrics)
}

fn ratio(num: f64, denom: f64) -> Option<f64> {
    if denom > 0.0 { Some(num / denom) } else { None }
}

fn shannon(values: &[f64]) -> f64 {
    let sum: f64 = values.iter().copied().sum();
    if sum <= 0.0 {
        return 0.0;
    }
    let mut h = 0.0;
    for v in values {
        if *v <= 0.0 {
            continue;
        }
        let p = v / sum;
        h -= p * p.ln();
    }
    h
}

fn evenness(values: &[f64], shannon: f64) -> Option<f64> {
    let count = values.iter().filter(|v| **v > 0.0).count();
    if count <= 1 || shannon <= 0.0 {
        return None;
    }
    let denom = (count as f64).ln();
    if denom <= 0.0 {
        None
    } else {
        Some(shannon / denom)
    }
}

fn pearson(x: &[f64], y: &[f64]) -> Option<f64> {
    if x.len() != y.len() || x.len() < 2 {
        return None;
    }
    let mean_x = x.iter().copied().sum::<f64>() / x.len() as f64;
    let mean_y = y.iter().copied().sum::<f64>() / y.len() as f64;
    let mut num = 0.0;
    let mut denom_x = 0.0;
    let mut denom_y = 0.0;
    for (&xi, &yi) in x.iter().zip(y.iter()) {
        let dx = xi - mean_x;
        let dy = yi - mean_y;
        num += dx * dy;
        denom_x += dx * dx;
        denom_y += dy * dy;
    }
    if denom_x <= 0.0 || denom_y <= 0.0 {
        None
    } else {
        Some(num / (denom_x.sqrt() * denom_y.sqrt()))
    }
}

fn spearman(x: &[f64], y: &[f64]) -> Option<f64> {
    if x.len() != y.len() || x.len() < 2 {
        return None;
    }
    let rx = rank_values(x);
    let ry = rank_values(y);
    pearson(&rx, &ry)
}

fn rank_values(values: &[f64]) -> Vec<f64> {
    let mut pairs: Vec<(usize, f64)> = values.iter().copied().enumerate().collect();
    pairs.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let mut ranks = vec![0.0; values.len()];
    let mut i = 0;
    while i < pairs.len() {
        let value = pairs[i].1;
        let mut j = i + 1;
        while j < pairs.len() && (pairs[j].1 - value).abs() < f64::EPSILON {
            j += 1;
        }
        let rank = (i + j - 1) as f64 / 2.0 + 1.0;
        for k in i..j {
            ranks[pairs[k].0] = rank;
        }
        i = j;
    }
    ranks
}

#[derive(Clone, Copy, Debug, Default)]
pub struct UniFracResult {
    pub weighted: Option<f64>,
    pub unweighted: Option<f64>,
}

struct UnifracComponents {
    weighted_raw: f64,
    unweighted_raw: f64,
    gt_unweighted_max: f64,
    max_weighted: f64,
}

const CANONICAL_RANKS: [&str; 8] = [
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "strain",
];

const BRANCH_LENGTHS: [f64; CANONICAL_RANKS.len()] = [1.0; CANONICAL_RANKS.len()];
const MASS_EPSILON: f64 = 1e-12;

fn branch_length_for_level(level: usize) -> f64 {
    BRANCH_LENGTHS
        .get(level)
        .copied()
        .unwrap_or_else(|| BRANCH_LENGTHS.last().copied().unwrap_or(1.0))
}

fn unifrac(
    rank: &str,
    gt_entries: &[ProfileEntry],
    pred_entries: &[ProfileEntry],
) -> UniFracResult {
    let Some(rank_index) = canonical_rank_index(rank) else {
        return UniFracResult::default();
    };

    match unifrac_components(rank_index, gt_entries, pred_entries) {
        Some(components) => {
            let weighted = if components.max_weighted > 0.0 {
                Some((components.weighted_raw / components.max_weighted).clamp(0.0, 1.0))
            } else {
                None
            };
            let unweighted = if components.gt_unweighted_max > MASS_EPSILON {
                Some((components.unweighted_raw / components.gt_unweighted_max).clamp(0.0, 1.0))
            } else {
                None
            };
            UniFracResult {
                weighted,
                unweighted,
            }
        }
        None => UniFracResult::default(),
    }
}

fn unifrac_components(
    rank_index: usize,
    gt_entries: &[ProfileEntry],
    pred_entries: &[ProfileEntry],
) -> Option<UnifracComponents> {
    let gt_total: f64 = gt_entries
        .iter()
        .filter(|e| e.percentage > 0.0)
        .map(|e| e.percentage)
        .sum();
    let pred_total: f64 = pred_entries
        .iter()
        .filter(|e| e.percentage > 0.0)
        .map(|e| e.percentage)
        .sum();

    if gt_total <= 0.0 && pred_total <= 0.0 {
        return None;
    }

    let mut tree = TaxTree::new();

    if gt_total > 0.0 {
        for entry in gt_entries.iter().filter(|e| e.percentage > 0.0) {
            let mass = entry.percentage / gt_total;
            if !mass.is_finite() || mass <= 0.0 {
                continue;
            }
            let path = parse_taxpath(entry);
            tree.add_mass(&path, rank_index, mass, true);
        }
    }

    if pred_total > 0.0 {
        for entry in pred_entries.iter().filter(|e| e.percentage > 0.0) {
            let mass = entry.percentage / pred_total;
            if !mass.is_finite() || mass <= 0.0 {
                continue;
            }
            let path = parse_taxpath(entry);
            tree.add_mass(&path, rank_index, mass, false);
        }
    }

    let (weighted_raw, unweighted_raw, gt_unweighted_max) = tree.compute_edge_flows();
    let path_length = BRANCH_LENGTHS
        .iter()
        .take(rank_index + 1)
        .copied()
        .sum::<f64>();
    let max_weighted = 2.0 * path_length;
    let clamped_weighted = if max_weighted > 0.0 {
        weighted_raw.max(0.0).min(max_weighted)
    } else {
        weighted_raw.max(0.0)
    };

    Some(UnifracComponents {
        weighted_raw: clamped_weighted,
        unweighted_raw: unweighted_raw.max(0.0),
        gt_unweighted_max,
        max_weighted,
    })
}

fn canonical_rank_index(rank: &str) -> Option<usize> {
    let canonical = canonical_rank(rank).ok()?;
    match canonical.as_str() {
        "superkingdom" => Some(0),
        "phylum" => Some(1),
        "class" => Some(2),
        "order" => Some(3),
        "family" => Some(4),
        "genus" => Some(5),
        "species" => Some(6),
        "strain" => Some(7),
        _ => None,
    }
}

fn parse_taxpath(entry: &ProfileEntry) -> Vec<Option<String>> {
    let entry_rank_index = canonical_rank_index(&entry.rank)
        .unwrap_or_else(|| CANONICAL_RANKS.len().saturating_sub(1));
    let mut result = vec![None; CANONICAL_RANKS.len()];

    let mut parts = if !entry.taxpathsn.trim().is_empty() {
        split_taxpath(&entry.taxpathsn)
    } else {
        split_taxpath(&entry.taxpath)
    };

    if parts.is_empty() {
        parts = split_taxpath(&entry.taxpath);
    }

    let mut normalized: Vec<String> = parts
        .iter()
        .filter_map(|part| normalize_taxpath_component(part))
        .collect();

    if normalized.is_empty() {
        if let Some(value) = normalize_taxpath_component(&entry.taxid) {
            normalized.push(value);
        }
    }

    if normalized.len() > entry_rank_index + 1 {
        normalized = normalized.split_off(normalized.len() - (entry_rank_index + 1));
    }

    let offset = entry_rank_index + 1 - normalized.len();
    for (idx, name) in normalized.into_iter().enumerate() {
        result[offset + idx] = Some(name);
    }

    if result[0].is_none() {
        if let Some(superkingdom) = infer_superkingdom(entry) {
            result[0] = Some(superkingdom);
        }
    }

    result
}

fn split_taxpath(path: &str) -> Vec<String> {
    path.split('|')
        .map(|part| part.trim())
        .filter(|part| !part.is_empty())
        .map(|part| part.to_string())
        .collect()
}

fn normalize_taxpath_component(raw: &str) -> Option<String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    let lower = trimmed.to_lowercase();
    if matches!(
        lower.as_str(),
        "root" | "cellular organisms" | "cellular root" | "acellular organisms" | "acellular root"
    ) {
        return None;
    }
    if lower.contains("viruses") && !lower.eq("viruses") {
        return Some("Viruses".to_string());
    }
    Some(trimmed.to_string())
}

fn infer_superkingdom(entry: &ProfileEntry) -> Option<String> {
    for component in split_taxpath(&entry.taxpathsn) {
        if let Some(label) =
            normalize_taxpath_component(&component).and_then(|value| detect_superkingdom(&value))
        {
            return Some(label);
        }
        if let Some(label) = detect_superkingdom(&component) {
            return Some(label);
        }
    }

    for component in split_taxpath(&entry.taxpath) {
        if let Some(label) =
            normalize_taxpath_component(&component).and_then(|value| detect_superkingdom(&value))
        {
            return Some(label);
        }
        if let Some(label) = detect_superkingdom(&component) {
            return Some(label);
        }
        if let Ok(taxid) = component.parse::<u32>() {
            if let Some(name) = superkingdom_from_taxid(taxid) {
                return Some(name.to_string());
            }
        }
    }

    detect_superkingdom(&entry.taxid)
}

fn detect_superkingdom(raw: &str) -> Option<String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    let canonical = trimmed
        .split("__")
        .last()
        .unwrap_or(trimmed)
        .split(':')
        .last()
        .unwrap_or(trimmed)
        .trim();
    let lower = canonical.to_lowercase();
    if lower.contains("viruses") {
        Some("Viruses".to_string())
    } else if lower.contains("bacteria") {
        Some("Bacteria".to_string())
    } else if lower.contains("archaea") || lower.contains("archaeota") || lower.contains("archaeon")
    {
        Some("Archaea".to_string())
    } else if lower.contains("eukary") {
        Some("Eukaryota".to_string())
    } else {
        None
    }
}

fn superkingdom_from_taxid(taxid: u32) -> Option<&'static str> {
    match taxid {
        2 => Some("Bacteria"),
        2157 => Some("Archaea"),
        2759 => Some("Eukaryota"),
        10239 => Some("Viruses"),
        _ => None,
    }
}

struct TaxTree {
    nodes: Vec<TaxNode>,
    index: HashMap<(usize, String), usize>,
}

impl TaxTree {
    fn new() -> Self {
        let root = TaxNode {
            name: "root".to_string(),
            parent: None,
            children: Vec::new(),
            gt_mass: 0.0,
            pred_mass: 0.0,
            branch_length: 0.0,
        };
        Self {
            nodes: vec![root],
            index: HashMap::new(),
        }
    }

    fn add_mass(&mut self, path: &[Option<String>], rank_index: usize, mass: f64, is_gt: bool) {
        if mass <= 0.0 {
            return;
        }

        let mut parent = 0usize;
        let mut parent_name = self.nodes[parent].name.clone();
        let mut target_node = parent;

        for level in 0..CANONICAL_RANKS.len() {
            let node_name = match path.get(level).and_then(|p| p.clone()) {
                Some(name) => name,
                None => format!("__missing_{}_under_{}", CANONICAL_RANKS[level], parent_name),
            };
            let child = self.get_or_create_child(parent, level, node_name.clone());
            parent = child;
            if level == rank_index {
                target_node = child;
            }
            parent_name = node_name;
        }

        if is_gt {
            self.nodes[target_node].gt_mass += mass;
        } else {
            self.nodes[target_node].pred_mass += mass;
        }
    }

    fn get_or_create_child(&mut self, parent: usize, rank_index: usize, name: String) -> usize {
        let key = (parent, name.clone());
        if let Some(id) = self.index.get(&key) {
            return *id;
        }

        let id = self.nodes.len();
        self.index.insert(key, id);
        self.nodes[parent].children.push(id);
        self.nodes.push(TaxNode {
            name,
            parent: Some(parent),
            children: Vec::new(),
            gt_mass: 0.0,
            pred_mass: 0.0,
            branch_length: branch_length_for_level(rank_index),
        });
        id
    }

    fn compute_edge_flows(&self) -> (f64, f64, f64) {
        let mut gt_subtree: Vec<f64> = self.nodes.iter().map(|node| node.gt_mass).collect();
        let mut pred_subtree: Vec<f64> = self.nodes.iter().map(|node| node.pred_mass).collect();

        for node_id in (1..self.nodes.len()).rev() {
            if let Some(parent) = self.nodes[node_id].parent {
                gt_subtree[parent] += gt_subtree[node_id];
                pred_subtree[parent] += pred_subtree[node_id];
            }
        }

        let mut weighted = 0.0;
        let mut unweighted = 0.0;
        let mut gt_unweighted = 0.0;

        for node_id in 1..self.nodes.len() {
            let branch_length = self.nodes[node_id].branch_length;
            if branch_length <= 0.0 {
                continue;
            }

            let gt_mass = gt_subtree[node_id];
            let pred_mass = pred_subtree[node_id];
            let diff = (gt_mass - pred_mass).abs();
            if diff > MASS_EPSILON {
                weighted += branch_length * diff;
            }

            let gt_present = gt_mass > MASS_EPSILON;
            let pred_present = pred_mass > MASS_EPSILON;
            let presence_diff = match (gt_present, pred_present) {
                (true, false) | (false, true) => 1.0,
                _ => 0.0,
            };
            if presence_diff > 0.0 {
                unweighted += branch_length * presence_diff;
            }
            if gt_present {
                gt_unweighted += branch_length;
            }
        }

        (weighted, unweighted, gt_unweighted)
    }
}

struct TaxNode {
    name: String,
    parent: Option<usize>,
    children: Vec<usize>,
    gt_mass: f64,
    pred_mass: f64,
    branch_length: f64,
}

fn abundance_rank_error(
    gt_map: &HashMap<String, f64>,
    pred_map: &HashMap<String, f64>,
) -> Result<f64> {
    let num_gt = gt_map.len();
    let num_pred = pred_map.len();

    if num_gt == 0 && num_pred == 0 {
        return Ok(0.0);
    }

    let mut gt_sorted: Vec<(String, f64)> = gt_map
        .iter()
        .map(|(taxon, value)| (taxon.clone(), *value))
        .collect();
    gt_sorted.sort_by(
        |a, b| match b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal) {
            Ordering::Equal => a.0.cmp(&b.0),
            other => other,
        },
    );

    let mut pred_sorted: Vec<(String, f64)> = pred_map
        .iter()
        .map(|(taxon, value)| (taxon.clone(), *value))
        .collect();
    pred_sorted.sort_by(
        |a, b| match b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal) {
            Ordering::Equal => a.0.cmp(&b.0),
            other => other,
        },
    );

    let mut gt_ranks: HashMap<String, usize> = HashMap::new();
    for (idx, (taxon, _)) in gt_sorted.iter().enumerate() {
        gt_ranks.insert(taxon.clone(), idx + 1);
    }

    let mut pred_ranks: HashMap<String, usize> = HashMap::new();
    for (idx, (taxon, _)) in pred_sorted.iter().enumerate() {
        pred_ranks.insert(taxon.clone(), idx + 1);
    }

    let mut total_error = 0.0;

    if num_gt > 0 {
        let num_gt_f = num_gt as f64;
        for (taxon, &gt_rank) in &gt_ranks {
            if let Some(&pred_rank) = pred_ranks.get(taxon) {
                let weight = (num_gt - gt_rank + 1) as f64 / num_gt_f;
                let rank_diff = if gt_rank > pred_rank {
                    (gt_rank - pred_rank) as f64
                } else {
                    (pred_rank - gt_rank) as f64
                };
                total_error += rank_diff * weight;
            }
        }

        for (taxon, &gt_rank) in &gt_ranks {
            if !pred_ranks.contains_key(taxon) {
                let offset = (num_gt - gt_rank + 1) as f64;
                let weight = offset / num_gt_f;
                total_error += offset * weight;
            }
        }
    }

    if num_pred > 0 {
        let num_pred_f = num_pred as f64;
        for (taxon, &pred_rank) in &pred_ranks {
            if !gt_ranks.contains_key(taxon) {
                let offset = (num_pred - pred_rank + 1) as f64;
                let weight = offset / num_pred_f;
                total_error += offset * weight;
            }
        }
    }

    let mut normalization = 0.0;
    if num_gt > 0 {
        let num_gt_f = num_gt as f64;
        for j in 1..=num_gt {
            let term = (num_gt + 1 - j) as f64;
            normalization += (term * term) / num_gt_f;
        }
    }
    if num_pred > 0 {
        let num_pred_f = num_pred as f64;
        for k in 1..=num_pred {
            let term = (num_pred + 1 - k) as f64;
            normalization += (term * term) / num_pred_f;
        }
    }

    if normalization <= 0.0 {
        return Ok(0.0);
    }

    let ratio = total_error / normalization;
    if ratio < 0.0 {
        Ok(0.0)
    } else if ratio > 1.0 {
        Ok(1.0)
    } else {
        Ok(ratio)
    }
}

fn mass_weighted_abundance_rank_error(
    gt_map: &HashMap<String, f64>,
    pred_map: &HashMap<String, f64>,
    p: f64,
) -> Result<f64> {
    ensure!(p > 0.0, "mARE exponent must be positive");

    let (gt_weights, gt_ranks, gt_count, gt_mass) = prepare_profile(gt_map);
    let (pred_weights, pred_ranks, pred_count, pred_mass) = prepare_profile(pred_map);

    let denominator = gt_mass + pred_mass;
    if denominator <= 0.0 {
        return Ok(0.0);
    }

    let mut union: BTreeSet<String> = BTreeSet::new();
    for taxon in gt_weights.keys() {
        union.insert(taxon.clone());
    }
    for taxon in pred_weights.keys() {
        union.insert(taxon.clone());
    }

    if union.is_empty() {
        return Ok(0.0);
    }

    let max_rank_diff = (gt_count.max(pred_count).saturating_sub(1)) as f64;

    let mut numerator = 0.0;
    for taxon in union {
        match (gt_weights.get(&taxon), pred_weights.get(&taxon)) {
            (Some(&p_gt), Some(&p_pred)) => {
                let rank_penalty = if max_rank_diff > 0.0 {
                    let gt_rank = gt_ranks.get(&taxon).copied().unwrap_or(1);
                    let pred_rank = pred_ranks.get(&taxon).copied().unwrap_or(1);
                    let diff = (gt_rank as i64 - pred_rank as i64).abs() as f64;
                    (diff / max_rank_diff).powf(p)
                } else {
                    0.0
                };
                numerator += rank_penalty * (p_gt + p_pred);
            }
            (Some(&p_gt), None) => {
                numerator += p_gt;
            }
            (None, Some(&p_pred)) => {
                numerator += p_pred;
            }
            (None, None) => {}
        }
    }

    let ratio = numerator / denominator;
    Ok(ratio.max(0.0).min(1.0))
}

fn prepare_profile(
    abundances: &HashMap<String, f64>,
) -> (HashMap<String, f64>, HashMap<String, usize>, usize, f64) {
    let mut entries: Vec<(String, f64)> = abundances
        .iter()
        .map(|(taxon, value)| (taxon.clone(), *value))
        .collect();

    let total: f64 = entries.iter().map(|(_, v)| *v).sum();
    if !total.is_finite() || total <= 0.0 {
        return (HashMap::new(), HashMap::new(), 0, 0.0);
    }

    for (_, value) in entries.iter_mut() {
        *value /= total;
    }

    entries.sort_by(
        |a, b| match b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal) {
            Ordering::Equal => a.0.cmp(&b.0),
            other => other,
        },
    );

    let mut ranks = HashMap::new();
    let mut weights = HashMap::new();
    for (idx, (taxon, value)) in entries.into_iter().enumerate() {
        ranks.insert(taxon.clone(), idx + 1);
        weights.insert(taxon, value);
    }

    let sum_norm = weights.values().copied().sum();
    let count = ranks.len();
    (weights, ranks, count, sum_norm)
}

fn write_metrics<W: Write>(
    writer: &mut W,
    label: &str,
    sample: &str,
    rank: &str,
    metrics: &Metrics,
) -> Result<()> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        label,
        sample,
        rank,
        metrics.tp,
        metrics.fp,
        metrics.fn_,
        format_opt(metrics.precision),
        format_opt(metrics.recall),
        format_opt(metrics.f1),
        format_opt(metrics.jaccard),
        format_float(metrics.l1_error),
        format_opt(metrics.bray_curtis),
        format_float(metrics.shannon_pred),
        format_float(metrics.shannon_truth),
        format_opt(metrics.evenness_pred),
        format_opt(metrics.evenness_truth),
        format_opt(metrics.pearson),
        format_opt(metrics.spearman),
        format_opt(metrics.weighted_unifrac),
        format_opt(metrics.unweighted_unifrac),
        format_opt(metrics.abundance_rank_error),
        format_opt(metrics.mass_weighted_abundance_rank_error)
    )?;
    Ok(())
}

fn format_opt(value: Option<f64>) -> String {
    value
        .map(|v| format_float(v))
        .unwrap_or_else(|| "NA".to_string())
}

fn format_float(value: f64) -> String {
    format!("{:.5}", value)
}

fn taxonomy_dir() -> PathBuf {
    dirs::home_dir()
        .map(|p| p.join(".cami"))
        .unwrap_or_else(|| PathBuf::from(".cami"))
}

#[cfg(test)]
mod tests {
    use super::{
        CANONICAL_RANKS, ProfileEntry, abundance_rank_error, canonical_rank_index,
        highest_taxid_in_taxpath, mass_weighted_abundance_rank_error, unifrac, unifrac_components,
    };
    use crate::cami::Entry;
    use std::collections::HashMap;

    fn map(entries: &[(&str, f64)]) -> HashMap<String, f64> {
        entries
            .iter()
            .map(|(name, value)| ((*name).to_string(), *value))
            .collect()
    }

    fn profile_entry_for_test(taxpath: &str, percentage: f64) -> ProfileEntry {
        let parts: Vec<&str> = taxpath.split('|').filter(|s| !s.is_empty()).collect();
        let taxid = parts.last().copied().unwrap_or(taxpath);
        let rank_index = parts
            .len()
            .checked_sub(1)
            .unwrap_or(0)
            .min(CANONICAL_RANKS.len() - 1);
        let rank = CANONICAL_RANKS[rank_index].to_string();
        ProfileEntry {
            taxid: taxid.to_string(),
            rank,
            percentage,
            taxpath: taxpath.to_string(),
            taxpathsn: taxpath.to_string(),
        }
    }

    #[test]
    fn mass_weighted_abundance_rank_error_perfect_match() {
        let gt = map(&[("A", 0.5), ("B", 0.3), ("C", 0.2)]);
        let pred = map(&[("A", 0.5), ("B", 0.3), ("C", 0.2)]);
        let score = mass_weighted_abundance_rank_error(&gt, &pred, 1.0).unwrap();
        assert!((score - 0.0).abs() < 1e-9);
    }

    #[test]
    fn mass_weighted_abundance_rank_error_only_missed_taxa() {
        let gt = map(&[("A", 0.6), ("B", 0.4)]);
        let pred = map(&[]);
        let score = mass_weighted_abundance_rank_error(&gt, &pred, 1.0).unwrap();
        assert!((score - 1.0).abs() < 1e-9);
    }

    #[test]
    fn unifrac_identical_profiles_are_zero() {
        let rank = "species";
        let gt = vec![
            profile_entry_for_test("sk|p|c|o|f|g|s", 60.0),
            profile_entry_for_test("sk|p|c|o|f|h|i", 40.0),
        ];
        let pred = vec![
            profile_entry_for_test("sk|p|c|o|f|g|s", 60.0),
            profile_entry_for_test("sk|p|c|o|f|h|i", 40.0),
        ];
        let result = unifrac(rank, &gt, &pred);
        assert!(result.weighted.unwrap_or(0.0) < 1e-12);
        assert!(result.unweighted.unwrap_or(0.0) < 1e-12);
    }

    #[test]
    fn unifrac_strain_singletons_on_opposite_branches_reach_max() {
        let rank = "strain";
        let gt = vec![profile_entry_for_test("a|b|c|d|e|f|g|t1", 100.0)];
        let pred = vec![profile_entry_for_test("z|y|x|w|v|u|s|t2", 100.0)];
        let result = unifrac(rank, &gt, &pred);
        assert!((result.weighted.unwrap() - 1.0).abs() < 1e-9);
        assert!((result.unweighted.unwrap() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn unweighted_unifrac_max_follows_opal_formula() {
        let rank = "genus";
        let gt = vec![
            profile_entry_for_test("sk|p|c|o|f|g1", 60.0),
            profile_entry_for_test("sk|p|c|o|f|g2", 40.0),
        ];
        let pred = vec![
            profile_entry_for_test("sk|p|c|o|f|g1", 60.0),
            profile_entry_for_test("sk|p|c|o|f|g3", 40.0),
        ];
        let result = unifrac(rank, &gt, &pred);
        let rank_index = canonical_rank_index(rank).unwrap();
        let components = unifrac_components(rank_index, &gt, &pred).unwrap();
        let baseline = unifrac_components(rank_index, &gt, &[]).unwrap();
        let expected_normalized = components.unweighted_raw / baseline.gt_unweighted_max;
        assert!((result.unweighted.unwrap() - expected_normalized).abs() < 1e-9);
    }

    #[test]
    fn weighted_unifrac_is_monotonic_with_deeper_mismatch() {
        let rank = "strain";
        let gt = vec![profile_entry_for_test("sk|p|c|o|f|g1|s1|t1", 100.0)];
        let pred_genus = vec![profile_entry_for_test("sk|p|c|o|f|g2", 100.0)];
        let pred_strain = vec![profile_entry_for_test("sk|p|c|o|f|g2|s2|t2", 100.0)];

        let rank_index = canonical_rank_index(rank).unwrap();
        let genus_components = unifrac_components(rank_index, &gt, &pred_genus).unwrap();
        let strain_components = unifrac_components(rank_index, &gt, &pred_strain).unwrap();

        assert!(genus_components.weighted_raw > 0.0);
        assert!(strain_components.weighted_raw >= genus_components.weighted_raw - 1e-9);
        assert!(strain_components.weighted_raw <= genus_components.max_weighted + 1e-9);
    }

    #[test]
    fn mass_weighted_abundance_rank_error_only_false_positives() {
        let gt = map(&[]);
        let pred = map(&[("X", 0.7), ("Y", 0.3)]);
        let score = mass_weighted_abundance_rank_error(&gt, &pred, 1.0).unwrap();
        assert!((score - 1.0).abs() < 1e-9);
    }

    #[test]
    fn mass_weighted_abundance_rank_error_mixed_case() {
        let gt = map(&[("A", 0.6), ("B", 0.3), ("C", 0.1)]);
        let pred = map(&[("A", 0.55), ("C", 0.35), ("D", 0.10)]);
        let score = mass_weighted_abundance_rank_error(&gt, &pred, 1.0).unwrap();
        assert!((score - 0.3125).abs() < 1e-9);
    }

    #[test]
    fn mass_weighted_abundance_rank_error_rank_reversal() {
        let gt = map(&[("A", 0.5), ("B", 0.3), ("C", 0.2)]);
        let pred = map(&[("A", 0.2), ("B", 0.3), ("C", 0.5)]);
        let score = mass_weighted_abundance_rank_error(&gt, &pred, 1.0).unwrap();
        assert!((score - 0.7).abs() < 1e-9);
    }

    #[test]
    fn abundance_rank_error_perfect_match() {
        let gt = map(&[("A", 0.5), ("B", 0.3), ("C", 0.2)]);
        let pred = map(&[("A", 0.5), ("B", 0.3), ("C", 0.2)]);
        let score = abundance_rank_error(&gt, &pred).unwrap();
        assert!((score - 0.0).abs() < 1e-9);
    }

    #[test]
    fn abundance_rank_error_rank_swap() {
        let gt = map(&[("A", 0.5), ("B", 0.3), ("C", 0.2)]);
        let pred = map(&[("B", 0.5), ("A", 0.3), ("C", 0.2)]);
        let score = abundance_rank_error(&gt, &pred).unwrap();
        assert!(score > 0.0 && score < 1.0);
    }

    #[test]
    fn abundance_rank_error_only_missed_taxa() {
        let gt = map(&[("A", 0.6), ("B", 0.4)]);
        let pred = map(&[]);
        let score = abundance_rank_error(&gt, &pred).unwrap();
        assert!((score - 1.0).abs() < 1e-9);
    }

    #[test]
    fn abundance_rank_error_only_false_positives() {
        let gt = map(&[]);
        let pred = map(&[("X", 0.7), ("Y", 0.3)]);
        let score = abundance_rank_error(&gt, &pred).unwrap();
        assert!((score - 1.0).abs() < 1e-9);
    }

    #[test]
    fn abundance_rank_error_mixed_case() {
        let gt = map(&[("A", 0.6), ("B", 0.3), ("C", 0.1)]);
        let pred = map(&[("A", 0.55), ("C", 0.35), ("D", 0.10)]);
        let score = abundance_rank_error(&gt, &pred).unwrap();
        assert!(score > 0.0 && score < 1.0);
    }

    #[test]
    fn highest_taxid_prefers_first_taxpath_entry() {
        let entry = Entry {
            taxid: "2956277".to_string(),
            rank: "species".to_string(),
            taxpath: "2731618|2731619|3424659|1198980|2956277".to_string(),
            taxpathsn: "".to_string(),
            percentage: 1.0,
        };

        assert_eq!(highest_taxid_in_taxpath(&entry), Some(2_731_618));
    }

    #[test]
    fn highest_taxid_falls_back_to_entry_taxid() {
        let entry = Entry {
            taxid: "12345".to_string(),
            rank: "species".to_string(),
            taxpath: "||abc||".to_string(),
            taxpathsn: "".to_string(),
            percentage: 1.0,
        };

        assert_eq!(highest_taxid_in_taxpath(&entry), Some(12_345));
    }
}

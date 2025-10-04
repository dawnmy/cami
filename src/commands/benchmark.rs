use std::collections::{BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

use anyhow::{Context, Result, bail};

use crate::cami;
use crate::cami::{Entry, Sample};

#[derive(Clone)]
pub struct BenchmarkConfig {
    pub ground_truth: PathBuf,
    pub predictions: Vec<PathBuf>,
    pub labels: Vec<String>,
    pub gmin: Option<f64>,
    pub normalize: bool,
    pub by_domain: bool,
    pub output: PathBuf,
    pub ranks: Option<Vec<String>>,
}

#[derive(Clone)]
struct ProfileEntry {
    taxid: String,
    percentage: f64,
    taxpath: String,
}

pub fn run(cfg: &BenchmarkConfig) -> Result<()> {
    if cfg.predictions.is_empty() {
        bail!("at least one predicted profile must be provided");
    }

    if !cfg.labels.is_empty() && cfg.labels.len() != cfg.predictions.len() {
        bail!("labels count must match number of predicted profiles");
    }

    let gt_samples = cami::parse_cami(&cfg.ground_truth)?;
    if gt_samples.is_empty() {
        bail!("ground truth profile has no samples");
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
            "profile\tsample\trank\ttp\tfp\tfn\tprecision\trecall\tf1\tjaccard\tl1_error\tbray_curtis\tshannon_pred\tshannon_truth\tevenness_pred\tevenness_truth\tpearson\tspearman\tweighted_unifrac\tunweighted_unifrac"
        )?;

        let gt_map = build_profile_map(
            &gt_samples,
            ranks.as_ref(),
            cfg.gmin,
            cfg.normalize,
            domain.as_deref(),
        );

        let mut sample_ids: Vec<_> = gt_map.keys().cloned().collect();
        sample_ids.sort();

        for (pred_path, label) in cfg.predictions.iter().zip(labels.iter()) {
            let pred_samples = cami::parse_cami(pred_path)?;
            let pred_map = build_profile_map(
                &pred_samples,
                ranks.as_ref(),
                None,
                cfg.normalize,
                domain.as_deref(),
            );

            for sample_id in &sample_ids {
                let Some(gt_ranks) = gt_map.get(sample_id) else {
                    continue;
                };
                let mut rank_names: Vec<_> = gt_ranks.keys().cloned().collect();
                rank_names.sort();
                let pred_ranks = pred_map.get(sample_id);
                for rank in rank_names {
                    let gt_entries = gt_ranks.get(&rank).unwrap();
                    let pred_entries = pred_ranks
                        .and_then(|map| map.get(&rank))
                        .map(|v| v.as_slice())
                        .unwrap_or(&[]);
                    let metrics = compute_metrics(gt_entries, pred_entries);
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

fn canonical_rank(rank: &str) -> Result<String> {
    let trimmed = rank.trim();
    if trimmed.is_empty() {
        bail!("rank names must not be empty");
    }
    let lower = trimmed.to_lowercase();
    let canonical = match lower.as_str() {
        "d" | "domain" | "superkingdom" => "superkingdom",
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
    min_threshold: Option<f64>,
    normalize: bool,
    domain: Option<&str>,
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
            if let Some(threshold) = min_threshold {
                if entry.percentage < threshold {
                    continue;
                }
            }
            if let Some(filter) = &rank_filter {
                if !filter.contains(&entry.rank) {
                    continue;
                }
            }
            if let Some(domain) = &domain_lower {
                if !entry_belongs_to_domain(entry, domain) {
                    continue;
                }
            }
            sample_map
                .entry(entry.rank.clone())
                .or_default()
                .push(ProfileEntry {
                    taxid: entry.taxid.clone(),
                    percentage: entry.percentage,
                    taxpath: entry.taxpath.clone(),
                });
        }

        if normalize {
            for entries in sample_map.values_mut() {
                let sum: f64 = entries.iter().map(|e| e.percentage).sum();
                if sum > 0.0 {
                    for entry in entries.iter_mut() {
                        entry.percentage = entry.percentage / sum * 100.0;
                        entry.percentage = (entry.percentage * 100000.0).round() / 100000.0;
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

fn entry_belongs_to_domain(entry: &Entry, domain: &str) -> bool {
    entry.taxpathsn.split('|').any(|name| {
        let trimmed = name.trim();
        !trimmed.is_empty() && trimmed.eq_ignore_ascii_case(domain)
    })
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
}

fn compute_metrics(gt_entries: &[ProfileEntry], pred_entries: &[ProfileEntry]) -> Metrics {
    let mut metrics = Metrics::default();
    let mut gt_map: HashMap<String, &ProfileEntry> = HashMap::new();
    for entry in gt_entries {
        if entry.percentage > 0.0 {
            gt_map.insert(entry.taxid.clone(), entry);
        }
    }
    let mut pred_map: HashMap<String, &ProfileEntry> = HashMap::new();
    for entry in pred_entries {
        if entry.percentage > 0.0 {
            pred_map.insert(entry.taxid.clone(), entry);
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
        let g = gt_map.get(taxid).map(|e| e.percentage).unwrap_or(0.0);
        let p = pred_map.get(taxid).map(|e| e.percentage).unwrap_or(0.0);
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

    let (weighted, unweighted) = unifrac(&gt_map, &pred_map, gt_total, pred_total);
    metrics.weighted_unifrac = weighted;
    metrics.unweighted_unifrac = unweighted;

    metrics
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

fn unifrac(
    gt_map: &HashMap<String, &ProfileEntry>,
    pred_map: &HashMap<String, &ProfileEntry>,
    gt_total: f64,
    pred_total: f64,
) -> (Option<f64>, Option<f64>) {
    if gt_map.is_empty() && pred_map.is_empty() {
        return (None, None);
    }

    let mut root = TreeNode::default();

    for entry in gt_map.values() {
        let mass = if gt_total > 0.0 {
            entry.percentage / gt_total
        } else {
            0.0
        };
        if mass <= 0.0 {
            continue;
        }
        let parts = split_taxpath(&entry.taxpath);
        add_mass(&mut root, &parts, 0, mass, true);
    }

    for entry in pred_map.values() {
        let mass = if pred_total > 0.0 {
            entry.percentage / pred_total
        } else {
            0.0
        };
        if mass <= 0.0 {
            continue;
        }
        let parts = split_taxpath(&entry.taxpath);
        add_mass(&mut root, &parts, 0, mass, false);
    }

    let (weighted_num, weighted_denom) = weighted_unifrac_contribution(&root);
    let weighted = if weighted_denom > 0.0 {
        Some(weighted_num / weighted_denom)
    } else {
        None
    };

    let (unweighted_num, unweighted_denom) = unweighted_unifrac_contribution(&root);
    let unweighted = if unweighted_denom > 0.0 {
        Some(unweighted_num / unweighted_denom)
    } else {
        None
    };

    (weighted, unweighted)
}

#[derive(Default)]
struct TreeNode {
    children: HashMap<String, TreeNode>,
    gt_mass: f64,
    pred_mass: f64,
}

fn split_taxpath(path: &str) -> Vec<String> {
    path.split('|')
        .map(|part| part.trim())
        .filter(|part| !part.is_empty())
        .map(|part| part.to_string())
        .collect()
}

fn add_mass(node: &mut TreeNode, parts: &[String], idx: usize, mass: f64, is_gt: bool) {
    if is_gt {
        node.gt_mass += mass;
    } else {
        node.pred_mass += mass;
    }

    if idx >= parts.len() {
        return;
    }

    let child = node
        .children
        .entry(parts[idx].clone())
        .or_insert_with(TreeNode::default);
    add_mass(child, parts, idx + 1, mass, is_gt);
}

fn weighted_unifrac_contribution(node: &TreeNode) -> (f64, f64) {
    let mut num = 0.0;
    let mut denom = 0.0;
    for child in node.children.values() {
        num += (child.gt_mass - child.pred_mass).abs();
        denom += child.gt_mass + child.pred_mass;
        let (child_num, child_denom) = weighted_unifrac_contribution(child);
        num += child_num;
        denom += child_denom;
    }
    (num, denom)
}

fn unweighted_unifrac_contribution(node: &TreeNode) -> (f64, f64) {
    let mut num = 0.0;
    let mut denom = 0.0;
    for child in node.children.values() {
        let gt_present = child.gt_mass > 0.0;
        let pred_present = child.pred_mass > 0.0;
        if gt_present || pred_present {
            denom += 1.0;
            if gt_present ^ pred_present {
                num += 1.0;
            }
        }
        let (child_num, child_denom) = unweighted_unifrac_contribution(child);
        num += child_num;
        denom += child_denom;
    }
    (num, denom)
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
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
        format_opt(metrics.unweighted_unifrac)
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

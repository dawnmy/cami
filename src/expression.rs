use crate::cami::{Entry, Sample};
use crate::taxonomy::{Taxonomy, parse_taxid};
use anyhow::{Result, anyhow};

#[derive(Debug, Clone)]
pub enum Expr {
    And(Box<Expr>, Box<Expr>),
    Or(Box<Expr>, Box<Expr>),
    Atom(String),
}

pub fn parse_expression(s: &str) -> Result<Expr> {
    let normalized = s.replace("&&", "&").replace("||", "|");
    let chars: Vec<char> = normalized.chars().collect();
    let mut pos = 0;

    fn skip_ws(chars: &[char], pos: &mut usize) {
        while *pos < chars.len() && chars[*pos].is_whitespace() {
            *pos += 1;
        }
    }

    fn parse_primary(chars: &[char], pos: &mut usize) -> Result<Expr> {
        skip_ws(chars, pos);
        if *pos >= chars.len() {
            return Err(anyhow!("unexpected end of expression"));
        }
        if chars[*pos] == '(' {
            *pos += 1;
            let expr = parse_expr(chars, pos)?;
            skip_ws(chars, pos);
            if *pos >= chars.len() || chars[*pos] != ')' {
                return Err(anyhow!("missing closing parenthesis"));
            }
            *pos += 1;
            Ok(expr)
        } else {
            let start = *pos;
            while *pos < chars.len()
                && chars[*pos] != '&'
                && chars[*pos] != '|'
                && chars[*pos] != ')'
            {
                *pos += 1;
            }
            let atom: String = chars[start..*pos].iter().collect();
            Ok(Expr::Atom(atom.trim().to_string()))
        }
    }

    fn parse_term(chars: &[char], pos: &mut usize) -> Result<Expr> {
        let mut left = parse_primary(chars, pos)?;
        loop {
            skip_ws(chars, pos);
            if *pos + 1 < chars.len() && chars[*pos] == '&' && chars[*pos + 1] == '&' {
                *pos += 1;
            }
            if *pos >= chars.len() || chars[*pos] != '&' {
                break;
            }
            *pos += 1;
            let right = parse_primary(chars, pos)?;
            left = Expr::And(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn parse_expr(chars: &[char], pos: &mut usize) -> Result<Expr> {
        let mut left = parse_term(chars, pos)?;
        loop {
            skip_ws(chars, pos);
            if *pos + 1 < chars.len() && chars[*pos] == '|' && chars[*pos + 1] == '|' {
                *pos += 1;
            }
            if *pos >= chars.len() || chars[*pos] != '|' {
                break;
            }
            *pos += 1;
            let right = parse_term(chars, pos)?;
            left = Expr::Or(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    let expr = parse_expr(&chars, &mut pos)?;
    skip_ws(&chars, &mut pos);
    if pos != chars.len() {
        return Err(anyhow!("unexpected characters after expression"));
    }
    Ok(expr)
}

pub fn expr_needs_taxdump(expr: &Expr) -> bool {
    match expr {
        Expr::And(a, b) | Expr::Or(a, b) => expr_needs_taxdump(a) || expr_needs_taxdump(b),
        Expr::Atom(s) => {
            let trimmed = s.trim().trim_start_matches('!').trim_start();
            trimmed.starts_with("tax") || trimmed.starts_with('t')
        }
    }
}

pub fn apply_filter(samples: &[Sample], expr: &Expr, taxonomy: Option<&Taxonomy>) -> Vec<Sample> {
    let mut filtered = Vec::new();
    for (sample_index, sample) in samples.iter().enumerate() {
        let mut ns = sample.clone();
        ns.entries = ns
            .entries
            .iter()
            .cloned()
            .filter(|entry| eval_expr(expr, samples, sample, entry, sample_index, taxonomy))
            .collect();
        if !ns.entries.is_empty() {
            filtered.push(ns);
        }
    }
    filtered
}

fn eval_expr(
    e: &Expr,
    samples: &[Sample],
    sample: &Sample,
    entry: &Entry,
    sample_index: usize,
    taxonomy: Option<&Taxonomy>,
) -> bool {
    match e {
        Expr::And(a, b) => {
            eval_expr(a, samples, sample, entry, sample_index, taxonomy)
                && eval_expr(b, samples, sample, entry, sample_index, taxonomy)
        }
        Expr::Or(a, b) => {
            eval_expr(a, samples, sample, entry, sample_index, taxonomy)
                || eval_expr(b, samples, sample, entry, sample_index, taxonomy)
        }
        Expr::Atom(s) => eval_atom(s, samples, sample, entry, sample_index, taxonomy),
    }
}

fn eval_atom(
    atom: &str,
    samples: &[Sample],
    sample: &Sample,
    entry: &Entry,
    sample_index: usize,
    taxonomy: Option<&Taxonomy>,
) -> bool {
    let atom = atom.trim();

    if let Some(res) = eval_rank(atom, sample, entry) {
        return res;
    }
    if let Some(res) = eval_sample(atom, samples, sample, sample_index) {
        return res;
    }
    if let Some(res) = eval_abundance(atom, entry) {
        return res;
    }
    if let Some(res) = eval_tax(atom, entry, taxonomy) {
        return res;
    }

    false
}

fn eval_rank(atom: &str, sample: &Sample, entry: &Entry) -> Option<bool> {
    let rest = if let Some(r) = atom.strip_prefix("rank") {
        r
    } else if let Some(r) = atom.strip_prefix('r') {
        r
    } else {
        return None;
    };
    let rest = rest.trim_start();
    if let Some(v) = rest.strip_prefix("==") {
        return Some(entry.rank == v.trim());
    }
    if let Some(v) = rest.strip_prefix("<=") {
        return Some(rank_compare(sample, &entry.rank, v.trim(), |a, b| a >= b));
    }
    if let Some(v) = rest.strip_prefix("<") {
        return Some(rank_compare(sample, &entry.rank, v.trim(), |a, b| a > b));
    }
    if let Some(v) = rest.strip_prefix(">=") {
        return Some(rank_compare(sample, &entry.rank, v.trim(), |a, b| a <= b));
    }
    if let Some(v) = rest.strip_prefix('>') {
        return Some(rank_compare(sample, &entry.rank, v.trim(), |a, b| a < b));
    }
    None
}

fn rank_compare<F>(sample: &Sample, entry_rank: &str, other: &str, cmp: F) -> bool
where
    F: Fn(usize, usize) -> bool,
{
    let Some(entry_idx) = sample.rank_index(entry_rank) else {
        return false;
    };
    let Some(other_idx) = sample.rank_index(other) else {
        return false;
    };
    cmp(entry_idx, other_idx)
}

fn eval_sample(
    atom: &str,
    samples: &[Sample],
    sample: &Sample,
    sample_index: usize,
) -> Option<bool> {
    let rest = if let Some(r) = atom.strip_prefix("sample") {
        r
    } else if let Some(r) = atom.strip_prefix('s') {
        r
    } else {
        return None;
    };
    let rest = rest.trim_start();
    if let Some(v) = rest.strip_prefix("==") {
        return Some(match_sample(v.trim(), samples, sample, sample_index));
    }
    None
}

fn match_sample(selector: &str, samples: &[Sample], sample: &Sample, sample_index: usize) -> bool {
    if selector.is_empty() || selector == ":" {
        return true;
    }
    let mut matched = false;
    for part in selector
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
    {
        if part.contains(':') {
            if match_sample_range(part, samples, sample_index) {
                matched = true;
                break;
            }
        } else if match_sample_value(part, samples, sample, sample_index) {
            matched = true;
            break;
        }
    }
    matched
}

fn match_sample_value(
    value: &str,
    samples: &[Sample],
    sample: &Sample,
    sample_index: usize,
) -> bool {
    if let Ok(idx) = value.parse::<usize>() {
        return idx == sample_index + 1;
    }
    if let Some(idx) = samples.iter().position(|s| s.id == value) {
        return idx == sample_index;
    }
    value == sample.id
}

fn match_sample_range(range: &str, samples: &[Sample], sample_index: usize) -> bool {
    let mut parts = range.split(':');
    let start_raw = parts.next().unwrap_or("").trim();
    let end_raw = parts.next().unwrap_or("").trim();

    if start_raw.is_empty() && end_raw.is_empty() {
        return true;
    }

    let start_idx = start_raw
        .parse::<usize>()
        .ok()
        .map(|n| n.saturating_sub(1))
        .or_else(|| samples.iter().position(|s| s.id == start_raw));
    if !start_raw.is_empty() && start_idx.is_none() {
        return false;
    }
    let end_idx = end_raw
        .parse::<usize>()
        .ok()
        .map(|n| n.saturating_sub(1))
        .or_else(|| samples.iter().position(|s| s.id == end_raw));
    if !end_raw.is_empty() && end_idx.is_none() {
        return false;
    }

    let start = start_idx.unwrap_or(0);
    let end = end_idx.unwrap_or_else(|| samples.len().saturating_sub(1));

    let (low, high) = if start <= end {
        (start, end)
    } else {
        (end, start)
    };
    (low..=high).contains(&sample_index)
}

fn eval_abundance(atom: &str, entry: &Entry) -> Option<bool> {
    let rest = if let Some(r) = atom.strip_prefix("abundance") {
        r
    } else if let Some(r) = atom.strip_prefix('a') {
        r
    } else {
        return None;
    };
    let rest = rest.trim_start();
    if let Some(v) = rest.strip_prefix("<=") {
        return Some(entry.percentage <= parse_f64(v.trim()));
    }
    if let Some(v) = rest.strip_prefix(">=") {
        return Some(entry.percentage >= parse_f64(v.trim()));
    }
    if let Some(v) = rest.strip_prefix("==") {
        let target = parse_f64(v.trim());
        return Some((entry.percentage - target).abs() < 1e-9);
    }
    if let Some(v) = rest.strip_prefix('>') {
        return Some(entry.percentage > parse_f64(v.trim()));
    }
    if let Some(v) = rest.strip_prefix('<') {
        return Some(entry.percentage < parse_f64(v.trim()));
    }
    None
}

fn parse_f64(s: &str) -> f64 {
    s.parse().unwrap_or(0.0)
}

fn eval_tax(atom: &str, entry: &Entry, taxonomy: Option<&Taxonomy>) -> Option<bool> {
    let mut working = atom.trim();
    let mut negate = false;
    if let Some(rest) = working.strip_prefix('!') {
        negate = true;
        working = rest.trim_start();
    }

    let rest = if let Some(r) = working.strip_prefix("tax") {
        r
    } else if let Some(r) = working.strip_prefix('t') {
        r
    } else {
        return None;
    };

    let rest = rest.trim_start();
    let (op, value) = if let Some(v) = rest.strip_prefix("==") {
        ("==", v.trim())
    } else if let Some(v) = rest.strip_prefix("<=") {
        ("<=", v.trim())
    } else if let Some(v) = rest.strip_prefix('<') {
        ("<", v.trim())
    } else {
        return Some(false);
    };

    let result = if let Some(taxonomy) = taxonomy {
        eval_taxonomy(entry, taxonomy, op, value)
    } else {
        eval_taxpath(entry, op, value)
    };
    Some(if negate { !result } else { result })
}

fn eval_taxonomy(entry: &Entry, taxonomy: &Taxonomy, op: &str, value: &str) -> bool {
    let entry_taxid = parse_taxid(&entry.taxid);
    let target_taxid = parse_taxid(value);

    match op {
        "==" => match (entry_taxid, target_taxid) {
            (Some(e), Some(t)) => e == t,
            _ => entry.taxid == value,
        },
        "<=" => {
            if let (Some(e), Some(t)) = (entry_taxid, target_taxid) {
                if e == t {
                    return true;
                }
                taxonomy
                    .ancestors_of(e)
                    .iter()
                    .any(|ancestor| *ancestor == t)
            } else {
                entry.taxpath.split('|').any(|tid| tid == value)
            }
        }
        "<" => {
            if let (Some(e), Some(t)) = (entry_taxid, target_taxid) {
                if e == t {
                    return false;
                }
                taxonomy
                    .ancestors_of(e)
                    .iter()
                    .any(|ancestor| *ancestor == t)
            } else {
                entry
                    .taxpath
                    .split('|')
                    .any(|tid| tid == value && tid != entry.taxid)
            }
        }
        _ => false,
    }
}

fn eval_taxpath(entry: &Entry, op: &str, value: &str) -> bool {
    match op {
        "==" => entry
            .taxpath
            .split('|')
            .last()
            .map(|tid| tid == value)
            .unwrap_or(false),
        "<=" => entry.taxpath.split('|').any(|tid| tid == value),
        "<" => {
            let mut parts: Vec<&str> = entry.taxpath.split('|').collect();
            if let Some(last) = parts.pop() {
                parts.contains(&value) && last != value
            } else {
                false
            }
        }
        _ => false,
    }
}

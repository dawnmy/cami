use anyhow::{Context, Result, anyhow};
use clap::{Parser, Subcommand};
use std::cmp::{max, min};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;
// std::io::Write not directly needed in this scope
use flate2::read::GzDecoder;
use reqwest::blocking::get;
use std::fs;
use tar::Archive;

#[derive(Parser)]
#[command(author, version, about)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Filter CAMI profiling data
    Filter {
        /// Filter expression, e.g. "rank==species & abundance>=0.1"
        expression: String,
        /// Output file (defaults to stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Fill up missing higher ranks using NCBI taxdump (downloads to ~/.cami)
        #[arg(long)]
        fill_up: bool,
        /// Target rank to fill up to (inclusive), e.g. phylum
        #[arg(long, default_value = "phylum")]
        to_rank: String,
        /// Renormalize percentages to 100 per rank after filtering/filling
        #[arg(long)]
        renorm: bool,
        /// Input CAMI file (positional, optional). If omitted, defaults to examples/text.cami
        #[arg(value_name = "INPUT", index = 2)]
        input: Option<PathBuf>,
    },
    /// List samples and summary statistics
    List {
        /// Input CAMI file (positional, optional). If omitted, defaults to examples/text.cami
        #[arg(value_name = "INPUT", index = 1)]
        input: Option<PathBuf>,
    },
    /// Preview first N entries per sample
    Preview {
        /// Number of entries per sample to show
        #[arg(short = 'n', long, default_value_t = 5)]
        n: usize,
        /// Input CAMI file (positional, optional). If omitted, defaults to examples/text.cami
        #[arg(value_name = "INPUT", index = 1)]
        input: Option<PathBuf>,
    },
}

#[derive(Debug, Clone)]
struct Sample {
    id: String,
    version: Option<String>,
    ranks: Vec<String>,
    entries: Vec<Entry>,
}

#[derive(Debug, Clone)]
struct Entry {
    taxid: String,
    rank: String,
    taxpath: String,
    taxpathsn: String,
    percentage: f64,
}

#[derive(Debug, Clone)]
struct TaxNode {
    parent: String,
    rank: String,
}

#[derive(Debug, Clone)]
struct Taxonomy {
    nodes: HashMap<String, TaxNode>,
    names: HashMap<String, String>,
    ancestors: HashMap<String, Vec<String>>,
}

fn parse_cami(path: &PathBuf) -> Result<Vec<Sample>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut samples = Vec::new();
    let mut current: Option<Sample> = None;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        if line.starts_with("@SampleID:") {
            if let Some(s) = current.take() {
                samples.push(s);
            }
            let id = line[10..].trim().to_string();
            current = Some(Sample {
                id,
                version: None,
                ranks: Vec::new(),
                entries: Vec::new(),
            });
            continue;
        }
        if line.starts_with("@Version:") {
            if let Some(s) = current.as_mut() {
                s.version = Some(line[9..].trim().to_string());
            }
            continue;
        }
        if line.starts_with("@Ranks:") {
            if let Some(s) = current.as_mut() {
                let ranks = line[7..].trim().split('|').map(|s| s.to_string()).collect();
                s.ranks = ranks;
            }
            continue;
        }
        if line.starts_with("@@TAXID") {
            // header for table, skip
            continue;
        }

        // Data row: TAXID<TAB>RANK<TAB>TAXPATH<TAB>TAXPATHSN<TAB>PERCENTAGE
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            // some files may use multiple spaces; try splitting on whitespace
            let fields_ws: Vec<&str> = line.split_whitespace().collect();
            if fields_ws.len() < 5 {
                continue;
            }
            if let Some(s) = current.as_mut() {
                let entry = Entry {
                    taxid: fields_ws[0].to_string(),
                    rank: fields_ws[1].to_string(),
                    taxpath: fields_ws[2].to_string(),
                    taxpathsn: fields_ws[3].to_string(),
                    percentage: fields_ws[4].parse().unwrap_or(0.0),
                };
                s.entries.push(entry);
            }
            continue;
        }

        if let Some(s) = current.as_mut() {
            let entry = Entry {
                taxid: fields[0].to_string(),
                rank: fields[1].to_string(),
                taxpath: fields[2].to_string(),
                taxpathsn: fields[3].to_string(),
                percentage: fields[4].parse().unwrap_or(0.0),
            };
            s.entries.push(entry);
        }
    }

    if let Some(s) = current.take() {
        samples.push(s);
    }

    Ok(samples)
}

#[derive(Debug)]
enum Expr {
    And(Box<Expr>, Box<Expr>),
    Or(Box<Expr>, Box<Expr>),
    Atom(String),
}

// Tokenize expression into atoms, parentheses and operators & and |
fn parse_expression(s: &str) -> Result<Expr> {
    let s = s.replace("&&", "&").replace("||", "|");
    let chars: Vec<char> = s.chars().collect();

    let mut pos = 0;

    fn skip_ws(chars: &Vec<char>, pos: &mut usize) {
        while *pos < chars.len() && chars[*pos].is_whitespace() {
            *pos += 1;
        }
    }

    fn parse_primary(chars: &Vec<char>, pos: &mut usize) -> Result<Expr> {
        skip_ws(chars, pos);
        if *pos >= chars.len() {
            return Err(anyhow::anyhow!("unexpected end"));
        }
        if chars[*pos] == '(' {
            *pos += 1;
            let e = parse_or(chars, pos)?;
            skip_ws(chars, pos);
            if *pos < chars.len() && chars[*pos] == ')' {
                *pos += 1;
                Ok(e)
            } else {
                Err(anyhow::anyhow!("missing )"))
            }
        } else {
            // read until whitespace or & or |
            let start = *pos;
            while *pos < chars.len()
                && !chars[*pos].is_whitespace()
                && chars[*pos] != '&'
                && chars[*pos] != '|'
                && chars[*pos] != ')'
            {
                *pos += 1;
            }
            let atom: String = chars[start..*pos]
                .iter()
                .collect::<String>()
                .trim()
                .to_string();
            Ok(Expr::Atom(atom))
        }
    }

    fn parse_and(chars: &Vec<char>, pos: &mut usize) -> Result<Expr> {
        let mut left = parse_primary(chars, pos)?;
        skip_ws(chars, pos);
        while *pos < chars.len() && chars[*pos] == '&' {
            *pos += 1;
            let right = parse_primary(chars, pos)?;
            left = Expr::And(Box::new(left), Box::new(right));
            skip_ws(chars, pos);
        }
        Ok(left)
    }

    fn parse_or(chars: &Vec<char>, pos: &mut usize) -> Result<Expr> {
        let mut left = parse_and(chars, pos)?;
        skip_ws(chars, pos);
        while *pos < chars.len() && chars[*pos] == '|' {
            *pos += 1;
            let right = parse_and(chars, pos)?;
            left = Expr::Or(Box::new(left), Box::new(right));
            skip_ws(chars, pos);
        }
        Ok(left)
    }

    let e = parse_or(&chars, &mut pos)?;
    Ok(e)
}

fn eval_expr(
    e: &Expr,
    samples: &Vec<Sample>,
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
    s: &str,
    samples: &Vec<Sample>,
    sample: &Sample,
    entry: &Entry,
    sample_index: usize,
    taxonomy: Option<&Taxonomy>,
) -> bool {
    let s = s.trim();
    if s.starts_with("rank==") {
        let v = s[6..].trim();
        return entry.rank == v;
    }
    if s.starts_with("rank<=") {
        let v = s[6..].trim();
        if let Some(pos_v) = sample.ranks.iter().position(|r| r == v) {
            if let Some(pos_e) = sample.ranks.iter().position(|r| r == &entry.rank) {
                return pos_e >= pos_v;
            }
        }
        return false;
    }
    if s.starts_with("rank<") {
        let v = s[5..].trim();
        if let Some(pos_v) = sample.ranks.iter().position(|r| r == v) {
            if let Some(pos_e) = sample.ranks.iter().position(|r| r == &entry.rank) {
                return pos_e > pos_v;
            }
        }
        return false;
    }
    if s.starts_with("rank>=") {
        let v = s[7..].trim();
        if let Some(pos_v) = sample.ranks.iter().position(|r| r == v) {
            if let Some(pos_e) = sample.ranks.iter().position(|r| r == &entry.rank) {
                return pos_e <= pos_v;
            }
        }
        return false;
    }
    if s.starts_with("rank>") {
        let v = s[5..].trim();
        if let Some(pos_v) = sample.ranks.iter().position(|r| r == v) {
            if let Some(pos_e) = sample.ranks.iter().position(|r| r == &entry.rank) {
                return pos_e < pos_v;
            }
        }
        return false;
    }

    if s.starts_with("sample==") {
        let v = s[8..].trim();
        if v.contains(':') {
            let parts: Vec<&str> = v.split(':').map(|p| p.trim()).collect();
            if parts.len() == 2 {
                let start = parts[0];
                let end = parts[1];
                if let (Ok(sa), Ok(ea)) = (start.parse::<usize>(), end.parse::<usize>()) {
                    let low = min(sa, ea);
                    let high = max(sa, ea);
                    return (low..=high).contains(&(sample_index + 1));
                }
                let mut start_idx: Option<usize> = None;
                let mut end_idx: Option<usize> = None;
                for (idx, smp) in samples.iter().enumerate() {
                    if start_idx.is_none() && smp.id == start {
                        start_idx = Some(idx);
                    }
                    if smp.id == end {
                        end_idx = Some(idx);
                    }
                }
                if let (Some(a), Some(b)) = (start_idx, end_idx) {
                    let low = min(a, b);
                    let high = max(a, b);
                    return (low..=high).contains(&sample_index);
                }
            }
            return false;
        }
        if v.contains(',') {
            let vals: Vec<&str> = v.split(',').map(|p| p.trim()).collect();
            return vals.iter().any(|val| {
                if let Ok(idx) = val.parse::<usize>() {
                    idx == sample_index + 1
                } else {
                    *val == sample.id
                }
            });
        }
        return if let Ok(idx) = v.parse::<usize>() {
            idx == sample_index + 1
        } else {
            v == sample.id
        };
    }

    if s.starts_with("abundance") {
        if let Some(idx) = s.find("<=") {
            let val: f64 = s[idx + 2..].trim().parse().unwrap_or(0.0);
            return entry.percentage <= val;
        }
        if let Some(idx) = s.find(">=") {
            let val: f64 = s[idx + 2..].trim().parse().unwrap_or(0.0);
            return entry.percentage >= val;
        }
        if let Some(idx) = s.find("==") {
            let val: f64 = s[idx + 2..].trim().parse().unwrap_or(0.0);
            return (entry.percentage - val).abs() < 1e-9;
        }
        if let Some(idx) = s.find('>') {
            let val: f64 = s[idx + 1..].trim().parse().unwrap_or(0.0);
            return entry.percentage > val;
        }
        if let Some(idx) = s.find('<') {
            let val: f64 = s[idx + 1..].trim().parse().unwrap_or(0.0);
            return entry.percentage < val;
        }
    }

    // tax operations: tax==, tax<=, tax<, with optional leading ! for negation
    let mut negate = false;
    let mut s2 = s;
    if s2.starts_with('!') {
        negate = true;
        s2 = &s2[1..];
    }
    if s2.starts_with("tax") {
        if taxonomy.is_none() {
            // without taxdump, fall back to simple taxpath string checks (less robust)
            if s2.starts_with("tax==") {
                let v = s2[5..].trim();
                let contains = entry
                    .taxpath
                    .split('|')
                    .last()
                    .map(|t| t == v)
                    .unwrap_or(false);
                return if negate { !contains } else { contains };
            }
            if s2.starts_with("tax<=") {
                let v = s2[5..].trim();
                let contains = entry.taxpath.split('|').any(|t| t == v);
                return if negate { !contains } else { contains };
            }
            if s2.starts_with("tax<") {
                let v = s2[4..].trim();
                let contains = entry.taxpath.split('|').any(|t| t == v);
                let eq = entry
                    .taxpath
                    .split('|')
                    .last()
                    .map(|t| t == v)
                    .unwrap_or(false);
                let res = contains && !eq;
                return if negate { !res } else { res };
            }
            return false;
        }

        let taxonomy = taxonomy.unwrap();
        let v = if s2.starts_with("tax==") {
            Some(s2[5..].trim())
        } else if s2.starts_with("tax<=") {
            Some(s2[5..].trim())
        } else if s2.starts_with("tax<") {
            Some(s2[4..].trim())
        } else {
            None
        };
        if v.is_none() {
            return false;
        }
        let target = v.unwrap();
        if s2.starts_with("tax==") {
            let res = entry.taxid == target;
            return if negate { !res } else { res };
        }
        let ancestors = taxonomy.ancestors.get(&entry.taxid);
        if s2.starts_with("tax<=") {
            let res = ancestors
                .map(|anc| anc.iter().any(|t| t == target))
                .unwrap_or(false)
                || entry.taxid == target;
            return if negate { !res } else { res };
        }
        if s2.starts_with("tax<") {
            let is_ancestor = ancestors
                .map(|anc| anc.iter().any(|t| t == target))
                .unwrap_or(false);
            let res = (is_ancestor || entry.taxid == target) && entry.taxid != target;
            return if negate { !res } else { res };
        }
    }

    false
}

fn expr_needs_taxdump(e: &Expr) -> bool {
    match e {
        Expr::And(a, b) | Expr::Or(a, b) => expr_needs_taxdump(a) || expr_needs_taxdump(b),
        Expr::Atom(s) => s.contains("tax"),
    }
}

fn write_cami(samples: &Vec<Sample>, out: &mut dyn io::Write) -> Result<()> {
    for s in samples {
        writeln!(out, "@SampleID:{}", s.id)?;
        if let Some(v) = &s.version {
            writeln!(out, "@Version:{}", v)?;
        }
        if !s.ranks.is_empty() {
            writeln!(out, "@Ranks:{}", s.ranks.join("|"))?;
        }
        writeln!(out, "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE")?;
        for e in &s.entries {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}",
                e.taxid, e.rank, e.taxpath, e.taxpathsn, e.percentage
            )?;
        }
    }
    Ok(())
}

fn ensure_taxdump(dir: &PathBuf) -> Result<()> {
    fs::create_dir_all(dir)?;
    let nodes = dir.join("nodes.dmp");
    let names = dir.join("names.dmp");
    if nodes.exists() && names.exists() {
        return Ok(());
    }
    // download taxdump.tar.gz
    let url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
    let resp = get(url)?;
    let bytes = resp.bytes()?;
    let gz = GzDecoder::new(&bytes[..]);
    let mut ar = Archive::new(gz);
    ar.unpack(dir)?;
    Ok(())
}

fn build_ancestors(nodes: &HashMap<String, TaxNode>) -> HashMap<String, Vec<String>> {
    fn helper(
        taxid: &str,
        nodes: &HashMap<String, TaxNode>,
        cache: &mut HashMap<String, Vec<String>>,
    ) -> Vec<String> {
        if let Some(existing) = cache.get(taxid) {
            return existing.clone();
        }
        let mut res = Vec::new();
        if let Some(node) = nodes.get(taxid) {
            if node.parent != taxid {
                res.push(node.parent.clone());
                let mut parent_line = helper(&node.parent, nodes, cache);
                res.append(&mut parent_line);
            }
        }
        cache.insert(taxid.to_string(), res.clone());
        res
    }

    let mut cache: HashMap<String, Vec<String>> = HashMap::new();
    let keys: Vec<String> = nodes.keys().cloned().collect();
    for taxid in keys {
        if !cache.contains_key(&taxid) {
            helper(&taxid, nodes, &mut cache);
        }
    }
    cache
}

fn parse_taxdump(dir: &PathBuf) -> Result<Taxonomy> {
    let mut nodes: HashMap<String, TaxNode> = HashMap::new();
    let mut names = HashMap::new();
    let nodes_path = dir.join("nodes.dmp");
    let names_path = dir.join("names.dmp");
    let f = File::open(&names_path).with_context(|| "opening names.dmp")?;
    for line in BufReader::new(f).lines() {
        let l = line?;
        let parts: Vec<&str> = l.split("\t|\t").collect();
        if parts.len() >= 2 {
            let taxid = parts[0].trim().to_string();
            let name = parts[1].trim().to_string();
            names.insert(taxid, name);
        }
    }
    let f2 = File::open(&nodes_path).with_context(|| "opening nodes.dmp")?;
    for line in BufReader::new(f2).lines() {
        let l = line?;
        let parts: Vec<&str> = l.split("\t|\t").collect();
        if parts.len() >= 3 {
            let taxid = parts[0].trim().to_string();
            let parent = parts[1].trim().to_string();
            let rank = parts[2].trim().to_string();
            nodes.insert(taxid, TaxNode { parent, rank });
        }
    }
    let ancestors = build_ancestors(&nodes);
    Ok(Taxonomy {
        nodes,
        names,
        ancestors,
    })
}

impl Taxonomy {
    fn lineage(&self, taxid: &str) -> Vec<(String, String, String)> {
        let mut res = Vec::new();
        let mut cur = taxid.to_string();
        let mut seen = std::collections::HashSet::new();
        while !seen.contains(&cur) {
            seen.insert(cur.clone());
            if let Some(node) = self.nodes.get(&cur) {
                let name = self.names.get(&cur).cloned().unwrap_or_else(|| cur.clone());
                res.push((cur.clone(), node.rank.clone(), name));
                if node.parent == cur {
                    break;
                }
                cur = node.parent.clone();
            } else {
                let name = self.names.get(&cur).cloned().unwrap_or_else(|| cur.clone());
                res.push((cur.clone(), "no_rank".to_string(), name));
                break;
            }
        }
        res.reverse();
        res
    }
}

fn fill_up_samples(samples: &mut Vec<Sample>, to_rank: &str, taxonomy: &Taxonomy) {
    for s in samples.iter_mut() {
        if s.entries.is_empty() {
            continue;
        }
        let Some(target_pos) = s.ranks.iter().position(|r| r == to_rank) else {
            continue;
        };
        let mut acc: HashMap<(String, String), f64> = HashMap::new();
        for e in &s.entries {
            let lineage = taxonomy.lineage(&e.taxid);
            if lineage.is_empty() {
                continue;
            }
            let mut rank_to_taxid: HashMap<String, String> = HashMap::new();
            for (tid, rnk, _name) in &lineage {
                rank_to_taxid.insert(rnk.clone(), tid.clone());
            }
            if let Some(entry_pos) = s.ranks.iter().position(|r| r == &e.rank) {
                let start = min(entry_pos, target_pos);
                let end = max(entry_pos, target_pos);
                for pos in start..=end {
                    let target_rank = &s.ranks[pos];
                    if let Some(taxid_for_rank) = rank_to_taxid.get(target_rank) {
                        let key = (taxid_for_rank.clone(), target_rank.clone());
                        *acc.entry(key).or_insert(0.0) += e.percentage;
                    }
                }
            }
        }
        if acc.is_empty() {
            continue;
        }
        let mut new_entries: Vec<Entry> = acc
            .into_iter()
            .map(|((tid, rank), pct)| {
                let lineage = taxonomy.lineage(&tid);
                let taxpath = lineage
                    .iter()
                    .map(|(t, _, _)| t.clone())
                    .collect::<Vec<_>>()
                    .join("|");
                let taxpathsn = lineage
                    .iter()
                    .map(|(_, _, n)| n.clone())
                    .collect::<Vec<_>>()
                    .join("|");
                Entry {
                    taxid: tid,
                    rank,
                    taxpath,
                    taxpathsn,
                    percentage: pct,
                }
            })
            .collect();
        new_entries.sort_by(|a, b| {
            let pa = s
                .ranks
                .iter()
                .position(|r| r == &a.rank)
                .unwrap_or(usize::MAX);
            let pb = s
                .ranks
                .iter()
                .position(|r| r == &b.rank)
                .unwrap_or(usize::MAX);
            if pa == pb {
                a.taxpath.cmp(&b.taxpath)
            } else {
                pa.cmp(&pb)
            }
        });
        s.entries = new_entries;
    }
}

fn renormalize(samples: &mut Vec<Sample>) {
    for s in samples.iter_mut() {
        let mut by_rank: HashMap<String, Vec<usize>> = HashMap::new();
        for (i, e) in s.entries.iter().enumerate() {
            by_rank.entry(e.rank.clone()).or_default().push(i);
        }
        for (_r, idxs) in by_rank {
            let sum: f64 = idxs.iter().map(|&i| s.entries[i].percentage).sum();
            if sum == 0.0 {
                continue;
            }
            for &i in &idxs {
                s.entries[i].percentage = s.entries[i].percentage / sum * 100.0;
            }
        }
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Filter {
            expression,
            output,
            fill_up,
            to_rank,
            renorm,
            input,
        } => {
            let input = input
                .as_ref()
                .cloned()
                .unwrap_or_else(|| PathBuf::from("examples/text.cami"));
            let samples = parse_cami(&input)?;
            let expr = parse_expression(expression).with_context(|| "parsing expression")?;
            let mut out: Box<dyn io::Write> = if let Some(p) = output {
                Box::new(File::create(p)?)
            } else {
                Box::new(io::stdout())
            };
            let mut filtered: Vec<Sample> = Vec::new();
            let needs_taxdump = expr_needs_taxdump(&expr) || *fill_up;
            let taxonomy = if needs_taxdump {
                let cami_dir = dirs::home_dir()
                    .map(|p| p.join(".cami"))
                    .unwrap_or_else(|| PathBuf::from(".cami"));
                let need_download =
                    !(cami_dir.join("nodes.dmp").exists() && cami_dir.join("names.dmp").exists());
                if need_download {
                    eprintln!("downloading taxdump to {}", cami_dir.display());
                }
                ensure_taxdump(&cami_dir)
                    .with_context(|| format!("ensuring taxdump in {}", cami_dir.display()))?;
                Some(parse_taxdump(&cami_dir)?)
            } else {
                None
            };

            for (i, s) in samples.iter().enumerate() {
                let mut ns = s.clone();
                ns.entries = ns
                    .entries
                    .into_iter()
                    .filter(|e| eval_expr(&expr, &samples, s, e, i, taxonomy.as_ref()))
                    .collect();
                if !ns.entries.is_empty() {
                    filtered.push(ns);
                }
            }
            if *fill_up {
                let tax = taxonomy
                    .as_ref()
                    .ok_or_else(|| anyhow!("fill-up requires taxonomy data"))?;
                fill_up_samples(&mut filtered, to_rank, tax);
                if *renorm {
                    renormalize(&mut filtered);
                }
                write_cami(&filtered, &mut *out)?;
            } else if *renorm {
                let mut renormed = filtered.clone();
                renormalize(&mut renormed);
                write_cami(&renormed, &mut *out)?;
            } else {
                write_cami(&filtered, &mut *out)?;
            }
        }
        Commands::Preview { n, input } => {
            let input = input
                .as_ref()
                .cloned()
                .unwrap_or_else(|| PathBuf::from("examples/text.cami"));
            let samples = parse_cami(&input)?;
            let mut out = io::stdout();
            for s in &samples {
                writeln!(out, "@SampleID:{}", s.id)?;
                if let Some(v) = &s.version {
                    writeln!(out, "@Version:{}", v)?;
                }
                if !s.ranks.is_empty() {
                    writeln!(out, "@Ranks:{}", s.ranks.join("|"))?;
                }
                writeln!(out, "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE")?;
                for e in s.entries.iter().take(*n) {
                    writeln!(
                        out,
                        "{}\t{}\t{}\t{}\t{}",
                        e.taxid, e.rank, e.taxpath, e.taxpathsn, e.percentage
                    )?;
                }
            }
        }
        Commands::List { input } => {
            let input = input
                .as_ref()
                .cloned()
                .unwrap_or_else(|| PathBuf::from("examples/text.cami"));
            let samples = parse_cami(&input)?;
            let mut out = io::stdout();
            for s in &samples {
                writeln!(out, "Sample: {}", s.id)?;
                writeln!(out, "  Ranks: {}", s.ranks.join(", "))?;
                writeln!(out, "  Taxa: {}", s.entries.len())?;
                // totals per rank
                let mut totals: HashMap<String, f64> = HashMap::new();
                for e in &s.entries {
                    *totals.entry(e.rank.clone()).or_insert(0.0) += e.percentage;
                }
                for r in &s.ranks {
                    let t = totals.get(r).cloned().unwrap_or(0.0);
                    writeln!(out, "    {}: {:.6}", r, t)?;
                }
            }
        }
    }

    Ok(())
}

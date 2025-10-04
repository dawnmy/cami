use anyhow::{Context, Result};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct Sample {
    pub id: String,
    pub version: Option<String>,
    pub ranks: Vec<String>,
    pub entries: Vec<Entry>,
}

#[derive(Debug, Clone)]
pub struct Entry {
    pub taxid: String,
    pub rank: String,
    pub taxpath: String,
    pub taxpathsn: String,
    pub percentage: f64,
}

impl Sample {
    pub fn rank_index(&self, rank: &str) -> Option<usize> {
        self.ranks.iter().position(|r| r == rank)
    }
}

pub fn parse_cami_reader<R: BufRead>(reader: R) -> Result<Vec<Sample>> {
    let mut samples = Vec::new();
    let mut current: Option<Sample> = None;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim_end();
        if line.is_empty() || line.starts_with('#') {
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
                let ranks = line[7..]
                    .trim()
                    .split('|')
                    .map(|r| r.trim().to_string())
                    .collect();
                s.ranks = ranks;
            }
            continue;
        }
        if line.starts_with("@@TAXID") {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        let (taxid, rank, taxpath, taxpathsn, percentage) = if fields.len() >= 5 {
            (fields[0], fields[1], fields[2], fields[3], fields[4])
        } else {
            let fields_ws: Vec<&str> = line.split_whitespace().collect();
            if fields_ws.len() < 5 {
                continue;
            }
            (
                fields_ws[0],
                fields_ws[1],
                fields_ws[2],
                fields_ws[3],
                fields_ws[4],
            )
        };

        if let Some(s) = current.as_mut() {
            let entry = Entry {
                taxid: taxid.to_string(),
                rank: rank.to_string(),
                taxpath: taxpath.to_string(),
                taxpathsn: taxpathsn.to_string(),
                percentage: percentage.parse().unwrap_or(0.0),
            };
            s.entries.push(entry);
        }
    }

    if let Some(s) = current.take() {
        samples.push(s);
    }
    Ok(samples)
}

pub fn parse_cami(path: &Path) -> Result<Vec<Sample>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader = BufReader::new(file);
    parse_cami_reader(reader)
}

pub fn load_samples(input: Option<&PathBuf>) -> Result<Vec<Sample>> {
    match input {
        Some(path) if path != Path::new("-") => parse_cami(path),
        _ => {
            let stdin = io::stdin();
            let handle = stdin.lock();
            parse_cami_reader(handle)
        }
    }
}

pub fn write_cami(samples: &[Sample], out: &mut dyn Write) -> Result<()> {
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

pub fn open_output(path: Option<&PathBuf>) -> Result<Box<dyn Write>> {
    let writer: Box<dyn Write> = match path {
        Some(p) if p != Path::new("-") => Box::new(File::create(p)?),
        _ => Box::new(io::stdout()),
    };
    Ok(writer)
}

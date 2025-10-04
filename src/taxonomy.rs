use anyhow::{Context, Result, anyhow};
use flate2::read::GzDecoder;
use reqwest::blocking::get;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::RwLock;
use std::thread;
use tar::Archive;

#[derive(Debug, Clone)]
pub struct TaxNode {
    pub parent: u32,
    pub rank: String,
}

#[derive(Debug)]
pub struct Taxonomy {
    nodes: HashMap<u32, TaxNode>,
    names: HashMap<u32, String>,
    ancestors: RwLock<HashMap<u32, Vec<u32>>>,
}

impl Taxonomy {
    pub fn load(dir: &Path) -> Result<Self> {
        let nodes_path = dir.join("nodes.dmp");
        let names_path = dir.join("names.dmp");

        let nodes_handle = {
            let path = nodes_path.clone();
            thread::spawn(move || parse_nodes(&path))
        };
        let names_handle = {
            let path = names_path.clone();
            thread::spawn(move || parse_names(&path))
        };

        let nodes = nodes_handle
            .join()
            .map_err(|_| anyhow!("failed to parse nodes.dmp"))??;
        let names = names_handle
            .join()
            .map_err(|_| anyhow!("failed to parse names.dmp"))??;

        Ok(Self {
            nodes,
            names,
            ancestors: RwLock::new(HashMap::new()),
        })
    }

    pub fn ancestors_of(&self, taxid: u32) -> Vec<u32> {
        if let Some(cached) = self.ancestors.read().unwrap().get(&taxid) {
            return cached.clone();
        }

        let mut lineage = Vec::new();
        let mut current = taxid;
        let mut seen = HashSet::new();
        while let Some(node) = self.nodes.get(&current) {
            if node.parent == current || seen.contains(&current) {
                break;
            }
            seen.insert(current);
            lineage.push(node.parent);
            current = node.parent;
        }

        self.ancestors
            .write()
            .unwrap()
            .insert(taxid, lineage.clone());
        lineage
    }

    pub fn lineage(&self, taxid: u32) -> Vec<(u32, String, String)> {
        let mut stack = Vec::new();
        let mut current = Some(taxid);
        let mut visited = HashSet::new();
        while let Some(tid) = current {
            if visited.contains(&tid) {
                break;
            }
            visited.insert(tid);
            let rank = self
                .nodes
                .get(&tid)
                .map(|n| n.rank.clone())
                .unwrap_or_else(|| "no_rank".to_string());
            let name = self
                .names
                .get(&tid)
                .cloned()
                .unwrap_or_else(|| tid.to_string());
            stack.push((tid, rank, name));
            current = self.nodes.get(&tid).and_then(|n| {
                if n.parent == tid {
                    None
                } else {
                    Some(n.parent)
                }
            });
        }
        stack.reverse();
        stack
    }
}

pub fn ensure_taxdump(dir: &Path) -> Result<()> {
    fs::create_dir_all(dir)?;
    let nodes = dir.join("nodes.dmp");
    let names = dir.join("names.dmp");
    if nodes.exists() && names.exists() {
        return Ok(());
    }
    let url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
    let resp = get(url)?;
    let bytes = resp.bytes()?;
    let gz = GzDecoder::new(&bytes[..]);
    let mut ar = Archive::new(gz);
    ar.unpack(dir)?;
    Ok(())
}

fn parse_nodes(path: &Path) -> Result<HashMap<u32, TaxNode>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut nodes = HashMap::new();
    for line in reader.lines() {
        let l = line?;
        let parts: Vec<&str> = l.split('|').collect();
        if parts.len() < 3 {
            continue;
        }
        let taxid: u32 = parts[0]
            .trim_matches(|c: char| c.is_whitespace())
            .parse()
            .unwrap_or(0);
        let parent: u32 = parts[1]
            .trim_matches(|c: char| c.is_whitespace())
            .parse()
            .unwrap_or(taxid);
        let rank = parts[2]
            .trim_matches(|c: char| c.is_whitespace())
            .to_string();
        nodes.insert(taxid, TaxNode { parent, rank });
    }
    Ok(nodes)
}

fn parse_names(path: &Path) -> Result<HashMap<u32, String>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut names = HashMap::new();
    for line in reader.lines() {
        let l = line?;
        let parts: Vec<&str> = l.split('|').collect();
        if parts.len() < 4 {
            continue;
        }
        let taxid: u32 = parts[0]
            .trim_matches(|c: char| c.is_whitespace())
            .parse()
            .unwrap_or(0);
        let name = parts[1]
            .trim_matches(|c: char| c.is_whitespace())
            .to_string();
        let class = parts[3].trim_matches(|c: char| c.is_whitespace());
        if class != "scientific name" {
            continue;
        }
        names.entry(taxid).or_insert(name);
    }
    Ok(names)
}

pub fn parse_taxid(taxid: &str) -> Option<u32> {
    taxid.parse().ok()
}

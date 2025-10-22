use anyhow::{Context, Result, anyhow};
use flate2::read::GzDecoder;
use reqwest::blocking::get;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::{Arc, RwLock};
use std::thread;
use tar::Archive;

#[derive(Debug, Clone)]
pub struct TaxNode {
    pub parent: u32,
    pub rank: String,
}

type LineageEntry = (u32, String, String);

#[derive(Debug)]
pub struct Taxonomy {
    nodes: HashMap<u32, TaxNode>,
    names: HashMap<u32, String>,
    ancestors: RwLock<HashMap<u32, Arc<[u32]>>>,
    lineages: RwLock<HashMap<u32, Arc<[LineageEntry]>>>,
    domains: RwLock<HashMap<u32, Option<String>>>,
    is_modern: bool,
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

        let is_modern = nodes.values().any(|node| {
            matches!(
                node.rank.to_ascii_lowercase().as_str(),
                "acellular root" | "cellular root" | "realm"
            )
        });

        Ok(Self {
            nodes,
            names,
            ancestors: RwLock::new(HashMap::new()),
            lineages: RwLock::new(HashMap::new()),
            domains: RwLock::new(HashMap::new()),
            is_modern,
        })
    }

    pub fn ancestors_of(&self, taxid: u32) -> Arc<[u32]> {
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

        let cached: Arc<[u32]> = lineage.into_boxed_slice().into();
        self.ancestors
            .write()
            .unwrap()
            .insert(taxid, cached.clone());
        cached
    }

    pub fn lineage(&self, taxid: u32) -> Arc<[LineageEntry]> {
        if let Some(cached) = self.lineages.read().unwrap().get(&taxid) {
            return cached.clone();
        }
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
        let cached: Arc<[LineageEntry]> = stack.into_boxed_slice().into();
        self.lineages.write().unwrap().insert(taxid, cached.clone());
        cached
    }

    pub fn name_of(&self, taxid: u32) -> Option<String> {
        self.names.get(&taxid).cloned()
    }

    pub fn rank_of(&self, taxid: u32) -> Option<String> {
        self.nodes.get(&taxid).map(|n| n.rank.clone())
    }

    pub fn domain_of(&self, taxid: u32) -> Option<String> {
        if let Some(cached) = self.domains.read().unwrap().get(&taxid) {
            return cached.clone();
        }

        let lineage = self.lineage(taxid);
        let result = lineage.iter().find_map(|(_, rank, name)| {
            if rank.eq_ignore_ascii_case("superkingdom")
                || rank.eq_ignore_ascii_case("domain")
                || rank.eq_ignore_ascii_case("acellular root")
            {
                Some(name.clone())
            } else {
                None
            }
        });

        let final_result = result.or_else(|| self.name_of(taxid));
        self.domains
            .write()
            .unwrap()
            .insert(taxid, final_result.clone());
        final_result
    }

    pub fn uses_modern_ranks(&self) -> bool {
        self.is_modern
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

pub fn default_taxdump_dir() -> PathBuf {
    dirs::home_dir()
        .map(|p| p.join(".cami"))
        .unwrap_or_else(|| PathBuf::from(".cami"))
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
        let rank = parts[2].trim_matches(|c: char| c.is_whitespace());
        nodes.insert(
            taxid,
            TaxNode {
                parent,
                rank: canonicalize_rank(rank),
            },
        );
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

fn canonicalize_rank(rank: &str) -> String {
    rank.trim().to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn domain_of_includes_acellular_root() {
        let mut nodes = HashMap::new();
        nodes.insert(
            1,
            TaxNode {
                parent: 1,
                rank: "no rank".to_string(),
            },
        );
        nodes.insert(
            10239,
            TaxNode {
                parent: 1,
                rank: "acellular root".to_string(),
            },
        );

        let mut names = HashMap::new();
        names.insert(1, "root".to_string());
        names.insert(10239, "Viruses".to_string());

        let taxonomy = Taxonomy {
            nodes,
            names,
            ancestors: RwLock::new(HashMap::new()),
            lineages: RwLock::new(HashMap::new()),
            domains: RwLock::new(HashMap::new()),
            is_modern: true,
        };

        assert_eq!(taxonomy.domain_of(10239), Some("Viruses".to_string()));
    }
}

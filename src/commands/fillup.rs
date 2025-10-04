use crate::cami::{load_samples, open_output, write_cami};
use crate::processing::{fill_up_default, fill_up_to};
use crate::taxonomy::{Taxonomy, ensure_taxdump};
use anyhow::{Context, Result};
use std::path::PathBuf;

pub struct FillupConfig<'a> {
    pub input: Option<&'a PathBuf>,
    pub output: Option<&'a PathBuf>,
    pub to_rank: Option<&'a str>,
    pub from_rank: Option<&'a str>,
}

pub fn run(cfg: &FillupConfig) -> Result<()> {
    let mut samples = load_samples(cfg.input)?;
    let dir = taxonomy_dir();
    ensure_taxdump(&dir).with_context(|| format!("ensuring taxdump in {}", dir.display()))?;
    let taxonomy = Taxonomy::load(&dir)?;

    if let Some(rank) = cfg.to_rank {
        fill_up_to(&mut samples, cfg.from_rank, rank, &taxonomy);
    } else {
        fill_up_default(&mut samples, cfg.from_rank, &taxonomy);
    }

    let mut out = open_output(cfg.output)?;
    write_cami(&samples, &mut *out)?;
    Ok(())
}

fn taxonomy_dir() -> PathBuf {
    dirs::home_dir()
        .map(|p| p.join(".cami"))
        .unwrap_or_else(|| PathBuf::from(".cami"))
}

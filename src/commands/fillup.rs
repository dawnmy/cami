use crate::cami::{load_samples, open_output, write_cami};
use crate::processing::{fill_up_default, fill_up_to, round_percentages};
use crate::taxonomy::{Taxonomy, default_taxdump_dir, ensure_taxdump};
use anyhow::{Context, Result};
use std::path::PathBuf;

pub struct FillupConfig<'a> {
    pub input: Option<&'a PathBuf>,
    pub output: Option<&'a PathBuf>,
    pub to_rank: Option<&'a str>,
    pub from_rank: Option<&'a str>,
    pub dmp_dir: Option<&'a PathBuf>,
}

pub fn run(cfg: &FillupConfig) -> Result<()> {
    let mut samples = load_samples(cfg.input)?;
    let dir = cfg.dmp_dir.cloned().unwrap_or_else(default_taxdump_dir);
    ensure_taxdump(&dir).with_context(|| format!("ensuring taxdump in {}", dir.display()))?;
    let taxonomy = Taxonomy::load(&dir)?;

    if let Some(rank) = cfg.to_rank {
        fill_up_to(&mut samples, cfg.from_rank, rank, &taxonomy);
    } else {
        fill_up_default(&mut samples, cfg.from_rank, &taxonomy);
    }

    round_percentages(&mut samples);

    let mut out = open_output(cfg.output)?;
    write_cami(&samples, &mut *out)?;
    Ok(())
}

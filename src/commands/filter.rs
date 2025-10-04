use crate::cami::{load_samples, open_output, write_cami};
use crate::expression::{apply_filter, expr_needs_taxdump, parse_expression};
use crate::processing::{fill_up_to, renormalize};
use crate::taxonomy::{Taxonomy, ensure_taxdump};
use anyhow::{Context, Result, anyhow};
use std::path::PathBuf;

pub struct FilterConfig<'a> {
    pub expression: &'a str,
    pub output: Option<&'a PathBuf>,
    pub fill_up: bool,
    pub from_rank: Option<&'a str>,
    pub to_rank: &'a str,
    pub renorm: bool,
    pub input: Option<&'a PathBuf>,
}

pub fn run(cfg: &FilterConfig) -> Result<()> {
    let samples = load_samples(cfg.input)?;
    let expr = parse_expression(cfg.expression).context("parsing expression")?;
    let needs_taxdump = expr_needs_taxdump(&expr) || cfg.fill_up;

    let taxonomy = if needs_taxdump {
        let dir = taxonomy_dir();
        ensure_taxdump(&dir).with_context(|| format!("ensuring taxdump in {}", dir.display()))?;
        Some(Taxonomy::load(&dir)?)
    } else {
        None
    };

    let mut filtered = apply_filter(&samples, &expr, taxonomy.as_ref());

    if cfg.fill_up {
        let tax = taxonomy
            .as_ref()
            .ok_or_else(|| anyhow!("fill-up requires taxonomy data"))?;
        fill_up_to(&mut filtered, cfg.from_rank, cfg.to_rank, tax);
    }
    if cfg.renorm {
        renormalize(&mut filtered);
    }

    let mut out = open_output(cfg.output)?;
    write_cami(&filtered, &mut *out)?;
    Ok(())
}

fn taxonomy_dir() -> PathBuf {
    dirs::home_dir()
        .map(|p| p.join(".cami"))
        .unwrap_or_else(|| PathBuf::from(".cami"))
}

use crate::cami::{load_samples, open_output, write_cami};
use crate::processing::renormalize;
use anyhow::Result;
use std::path::PathBuf;

pub struct RenormConfig<'a> {
    pub input: Option<&'a PathBuf>,
    pub output: Option<&'a PathBuf>,
}

pub fn run(cfg: &RenormConfig) -> Result<()> {
    let mut samples = load_samples(cfg.input)?;
    renormalize(&mut samples);
    let mut out = open_output(cfg.output)?;
    write_cami(&samples, &mut *out)?;
    Ok(())
}

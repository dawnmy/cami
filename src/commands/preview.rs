use crate::cami::load_samples;
use anyhow::Result;
use std::io::{self, Write};
use std::path::PathBuf;

pub struct PreviewConfig<'a> {
    pub n: usize,
    pub input: Option<&'a PathBuf>,
}

pub fn run(cfg: &PreviewConfig) -> Result<()> {
    let samples = load_samples(cfg.input)?;
    let mut out = io::stdout();
    for sample in &samples {
        writeln!(out, "@SampleID:{}", sample.id)?;
        if let Some(version) = &sample.version {
            writeln!(out, "@Version:{}", version)?;
        }
        if !sample.ranks.is_empty() {
            writeln!(out, "@Ranks:{}", sample.ranks.join("|"))?;
        }
        writeln!(out, "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE")?;
        for entry in sample.entries.iter().take(cfg.n) {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}",
                entry.taxid, entry.rank, entry.taxpath, entry.taxpathsn, entry.percentage
            )?;
        }
    }
    Ok(())
}

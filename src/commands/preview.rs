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
        let rank_tokens = sample.header_rank_tokens();
        if !rank_tokens.is_empty() {
            writeln!(out, "@Ranks:{}", rank_tokens.join("|"))?;
        }
        let extended = sample.is_modern_format()
            || sample
                .entries
                .iter()
                .any(|e| e.cami_genome_id.is_some() || e.cami_otu.is_some() || e.hosts.is_some());
        if extended {
            writeln!(
                out,
                "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\t_CAMI_genomeID\t_CAMI_OTU\tHOSTS"
            )?;
        } else {
            writeln!(out, "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE")?;
        }
        for entry in sample.entries.iter().take(cfg.n) {
            if extended {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    entry.taxid,
                    entry.rank,
                    entry.taxpath,
                    entry.taxpathsn,
                    entry.percentage,
                    entry.cami_genome_id.as_deref().unwrap_or(""),
                    entry.cami_otu.as_deref().unwrap_or(""),
                    entry.hosts.as_deref().unwrap_or("")
                )?;
            } else {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}",
                    entry.taxid, entry.rank, entry.taxpath, entry.taxpathsn, entry.percentage
                )?;
            }
        }
    }
    Ok(())
}

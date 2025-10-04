use crate::cami::load_samples;
use anyhow::Result;
use std::collections::HashMap;
use std::io::{self, Write};
use std::path::PathBuf;

pub struct ListConfig<'a> {
    pub input: Option<&'a PathBuf>,
}

pub fn run(cfg: &ListConfig) -> Result<()> {
    let samples = load_samples(cfg.input)?;
    let mut out = io::stdout();

    for sample in &samples {
        writeln!(out, "Sample: {}", sample.id)?;
        writeln!(out, "  Ranks: {}", sample.ranks.join(", "))?;
        writeln!(out, "  Total taxa: {}", sample.entries.len())?;

        let mut stats: HashMap<&str, (usize, f64)> = HashMap::new();
        for entry in &sample.entries {
            let stat = stats.entry(&entry.rank).or_insert((0, 0.0));
            if entry.percentage > 0.0 {
                stat.0 += 1;
                stat.1 += entry.percentage;
            }
        }

        for rank in &sample.ranks {
            let (count, total) = stats.get(rank.as_str()).cloned().unwrap_or((0, 0.0));
            writeln!(out, "    {}: taxa={} total={:.3}", rank, count, total)?;

        }
    }

    Ok(())
}

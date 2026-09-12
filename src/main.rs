mod hash;
mod run;
mod scan;
mod seq;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "kstrata", about = "exact-match k-mer strata across genomes")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Uniform random genome (control for chance k-mer sharing) -> <prefix>.seq + <prefix>.idx
    Random { prefix: PathBuf, #[arg(long)] length: u64, #[arg(long, default_value_t = 1)] seed: u64, #[arg(long, default_value_t = 100_000_000)] contig: u64 },
    /// FASTA(.gz) -> <prefix>.seq + <prefix>.idx
    Prep { fasta: PathBuf, prefix: PathBuf, /// keep only contigs whose name contains this string
        #[arg(long)] include: Option<String> },
    /// Presence of every query k-mer in each subject, and its multiplicity in the query
    Run {
        #[arg(long)] query: PathBuf,
        #[arg(long, default_value = "query")] query_name: String,
        /// name=prefix, repeatable (up to 8)
        #[arg(long = "subject")] subjects: Vec<String>,
        #[arg(short, long)] k: usize,
        #[arg(long, default_value_t = 1)] stride: u64,
        #[arg(long)] query_parts: Option<u32>,
        #[arg(long, default_value_t = 3.0)] budget_gb: f64,
        #[arg(long, default_value_t = 8.0)] chunk_mb: f64,
        #[arg(long)] no_self: bool,
        #[arg(long)] threads: Option<usize>,
        #[arg(long)] out: PathBuf,
    },
}

fn main() -> Result<()> {
    match Cli::parse().cmd {
        Cmd::Random { prefix, length, seed, contig } => seq::random(&prefix, length, seed, contig),
        Cmd::Prep { fasta, prefix, include } => seq::prep(&fasta, &prefix, include.as_deref()),
        Cmd::Run { query, query_name, subjects, k, stride, query_parts, budget_gb, chunk_mb, no_self, threads, out } => {
            if let Some(t) = threads { rayon::ThreadPoolBuilder::new().num_threads(t).build_global()?; }
            if stride == 0 { bail!("stride must be >= 1"); }
            let q = seq::load(&query, &query_name)?;
            let mut subs = Vec::new();
            for s in &subjects {
                let (name, prefix) = s.split_once('=').with_context(|| format!("subject must be name=prefix: {}", s))?;
                subs.push(seq::load(&PathBuf::from(prefix), name)?);
            }
            let opts = run::RunOpts {
                k, stride, query_parts, budget_bytes: (budget_gb * 1e9) as u64,
                chunk_len: (chunk_mb * 1e6) as u64, out: out.clone(), no_self,
            };
            let meta = run::run(&q, &subs, &opts)?;
            std::fs::write(out.with_extension("json"), serde_json::to_string_pretty(&meta)?)?;
            eprintln!("done in {:.1}s, max rss {:.2} GB", meta.seconds, meta.max_rss_bytes as f64 / 1e9);
            Ok(())
        }
    }
}

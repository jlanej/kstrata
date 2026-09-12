//! One run = one k. For each partition of the key space, the query's (key, position) pairs are
//! collected and sorted; the query's own keys give the self multiplicity of every entry; each
//! subject genome is scanned (in as many sub-partitions as the memory budget requires), sorted,
//! deduplicated and merged to set one presence bit per subject. Results are written to two
//! byte arrays indexed by position/stride: presence bits and self multiplicity (0 = no k-mer).
use crate::hash::{Exact, Poly};
use crate::scan::scan;
use crate::seq::Genome;
use anyhow::{bail, Result};
use memmap2::MmapMut;
use rayon::prelude::*;
use serde::Serialize;
use std::fs::OpenOptions;
use std::path::PathBuf;
use std::time::Instant;

pub const HMAX: usize = 1 << 20;

pub struct RunOpts {
    pub k: usize,
    pub stride: u64,
    pub query_parts: Option<u32>,
    pub budget_bytes: u64,
    pub chunk_len: u64,
    pub out: PathBuf,
    pub no_self: bool,
}

#[derive(Serialize)]
pub struct SubjectStats { pub name: String, pub length: u64, pub sub_parts: u32, pub distinct: u64, pub shared_distinct: u64, pub shared_entries: u64 }

#[derive(Serialize)]
pub struct Meta {
    pub k: usize, pub stride: u64, pub mode: String, pub query: String, pub query_length: u64,
    pub n_entries: u64, pub valid_entries: u64, pub query_parts: u32, pub budget_bytes: u64,
    pub query_distinct: u64, pub subjects: Vec<SubjectStats>,
    /// (multiplicity, number of distinct k-mers with that multiplicity); the last bin is >= HMAX
    pub self_histogram: Vec<(u64, u64)>,
    pub seconds: f64, pub max_rss_bytes: u64,
}

fn next_pow2(x: u64) -> u64 { let mut p = 1; while p < x { p <<= 1; } p }
fn log2(x: u64) -> u32 { x.trailing_zeros() }

fn range(total_bits: u32, index: u64) -> (u64, u64) {
    if total_bits == 0 { return (0, u64::MAX); }
    let lo = index << (64 - total_bits);
    let hi = lo | ((1u64 << (64 - total_bits)) - 1);
    (lo, hi)
}

fn max_rss() -> u64 {
    let mut ru: libc::rusage = unsafe { std::mem::zeroed() };
    unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut ru) };
    ru.ru_maxrss as u64 // bytes on macOS, kilobytes on Linux
}

fn log(t0: &Instant, msg: &str) {
    eprintln!("[{:8.1}s rss {:5.1}G] {}", t0.elapsed().as_secs_f64(), max_rss() as f64 / 1e9, msg);
}

fn scan_keys(exact: bool, g: &Genome, k: usize, lo: u64, hi: u64, chunk: u64) -> Vec<u64> {
    if exact { scan::<Exact, u64, _>(g, k, lo, hi, 1, chunk, |key, _| key) }
    else { scan::<Poly, u64, _>(g, k, lo, hi, 1, chunk, |key, _| key) }
}

fn scan_pairs(exact: bool, g: &Genome, k: usize, lo: u64, hi: u64, stride: u64, chunk: u64) -> Vec<(u64, u64)> {
    if exact { scan::<Exact, (u64, u64), _>(g, k, lo, hi, stride, chunk, |key, pos| (key, pos)) }
    else { scan::<Poly, (u64, u64), _>(g, k, lo, hi, stride, chunk, |key, pos| (key, pos)) }
}

pub fn run(query: &Genome, subjects: &[Genome], o: &RunOpts) -> Result<Meta> {
    let t0 = Instant::now();
    if subjects.len() > 8 { bail!("at most 8 subjects (one presence bit each)"); }
    if o.k < 1 || o.k > (1 << 24) { bail!("k out of range"); }
    let exact = o.k <= 32;
    let n_total = query.len();
    let n_entries = (n_total + o.stride - 1) / o.stride;
    let query_parts = match o.query_parts {
        Some(p) => next_pow2(p as u64) as u32,
        None => next_pow2((n_entries * 16 + o.budget_bytes - 1) / o.budget_bytes) as u32,
    };
    let pbits = log2(query_parts as u64);
    let sub_parts = |len: u64| -> u32 {
        let per_part = len / query_parts as u64 + 1;
        next_pow2((per_part * 8 + o.budget_bytes - 1) / o.budget_bytes).max(1) as u32
    };
    log(&t0, &format!("k={} stride={} mode={} query={} ({} bp, {} entries) parts={} subjects={}",
        o.k, o.stride, if exact { "exact" } else { "hash61" }, query.name, n_total, n_entries, query_parts, subjects.len()));

    let open = |suffix: &str| -> Result<MmapMut> {
        let path = o.out.with_extension(suffix);
        let f = OpenOptions::new().read(true).write(true).create(true).truncate(true).open(&path)?;
        f.set_len(n_entries)?;
        Ok(unsafe { MmapMut::map_mut(&f)? })
    };
    let mut pres = open("pres.u8")?;
    let mut mult_out = open("mult.u8")?;

    let mut hist = vec![0u64; HMAX + 1];
    let mut query_distinct = 0u64;
    let mut valid_entries = 0u64;
    let mut stats: Vec<SubjectStats> = subjects.iter().map(|s| SubjectStats {
        name: s.name.clone(), length: s.len(), sub_parts: sub_parts(s.len()), distinct: 0, shared_distinct: 0, shared_entries: 0 }).collect();
    let self_sub = sub_parts(n_total);

    for part in 0..query_parts as u64 {
        let (lo, hi) = range(pbits, part);
        let mut q = scan_pairs(exact, query, o.k, lo, hi, o.stride, o.chunk_len);
        q.par_sort_unstable_by_key(|e| e.0);
        valid_entries += q.len() as u64;
        let mut flags = vec![0u8; q.len()];
        let mut mult = vec![0u8; q.len()];
        log(&t0, &format!("part {}/{}: {} query entries sorted", part + 1, query_parts, q.len()));

        if !o.no_self {
            if o.stride == 1 {
                let mut i = 0;
                while i < q.len() {
                    let key = q[i].0;
                    let mut j = i;
                    while j < q.len() && q[j].0 == key { j += 1; }
                    let c = j - i;
                    hist[c.min(HMAX)] += 1;
                    query_distinct += 1;
                    let m = c.min(255) as u8;
                    for e in &mut mult[i..j] { *e = m; }
                    i = j;
                }
            } else {
                let mbits = log2(self_sub as u64);
                for sub in 0..self_sub as u64 {
                    let (slo, shi) = range(pbits + mbits, (part << mbits) | sub);
                    let mut keys = scan_keys(exact, query, o.k, slo, shi, o.chunk_len);
                    keys.par_sort_unstable();
                    let (mut i, mut a) = (0usize, 0usize);
                    while a < keys.len() {
                        let key = keys[a];
                        let mut b = a;
                        while b < keys.len() && keys[b] == key { b += 1; }
                        let c = b - a;
                        hist[c.min(HMAX)] += 1;
                        query_distinct += 1;
                        while i < q.len() && q[i].0 < key { i += 1; }
                        while i < q.len() && q[i].0 == key { mult[i] = c.min(255) as u8; i += 1; }
                        a = b;
                    }
                }
            }
            log(&t0, "self multiplicity done");
        }

        for (si, s) in subjects.iter().enumerate() {
            let bit = 1u8 << si;
            let mbits = log2(stats[si].sub_parts as u64);
            for sub in 0..stats[si].sub_parts as u64 {
                let (slo, shi) = range(pbits + mbits, (part << mbits) | sub);
                let mut keys = scan_keys(exact, s, o.k, slo, shi, o.chunk_len);
                keys.par_sort_unstable();
                keys.dedup();
                stats[si].distinct += keys.len() as u64;
                let (mut i, mut j) = (0usize, 0usize);
                while i < q.len() {
                    let key = q[i].0;
                    while j < keys.len() && keys[j] < key { j += 1; }
                    let hit = j < keys.len() && keys[j] == key;
                    let mut i2 = i;
                    while i2 < q.len() && q[i2].0 == key {
                        if hit { flags[i2] |= bit; }
                        i2 += 1;
                    }
                    if hit { stats[si].shared_distinct += 1; stats[si].shared_entries += (i2 - i) as u64; }
                    i = i2;
                }
            }
            log(&t0, &format!("subject {} merged ({} sub-parts)", s.name, stats[si].sub_parts));
        }

        // write: reorder by position, then sequential writes into the output arrays
        q.par_iter_mut().zip(flags.par_iter()).zip(mult.par_iter()).for_each(|((e, f), m)| {
            *e = (e.1, (*f as u64) | ((*m as u64) << 8));
        });
        drop(flags); drop(mult);
        q.par_sort_unstable_by_key(|e| e.0);
        for (pos, packed) in &q {
            let idx = (pos / o.stride) as usize;
            pres[idx] = (*packed & 0xff) as u8;
            mult_out[idx] = ((*packed >> 8) & 0xff) as u8;
        }
        log(&t0, "written");
    }
    pres.flush()?; mult_out.flush()?;

    let self_histogram: Vec<(u64, u64)> = hist.iter().enumerate().filter(|(_, &n)| n > 0).map(|(c, &n)| (c as u64, n)).collect();
    Ok(Meta {
        k: o.k, stride: o.stride, mode: if exact { "exact".into() } else { "hash61".into() },
        query: query.name.clone(), query_length: n_total, n_entries, valid_entries, query_parts,
        budget_bytes: o.budget_bytes, query_distinct, subjects: stats, self_histogram,
        seconds: t0.elapsed().as_secs_f64(), max_rss_bytes: max_rss(),
    })
}

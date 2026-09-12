//! Parallel scan of a prepared genome: every valid k-mer (no N, inside one contig) whose
//! canonical key falls in [lo, hi] and whose global start is a multiple of `stride` is emitted.
use crate::hash::Roller;
use crate::seq::Genome;
use rayon::prelude::*;

pub struct Chunk { pub contig_end: u64, pub a: u64, pub b: u64 }

pub fn chunks(g: &Genome, chunk_len: u64) -> Vec<Chunk> {
    let mut v = Vec::new();
    for c in &g.contigs {
        let mut a = c.start;
        let end = c.start + c.len;
        while a < end {
            let b = (a + chunk_len).min(end);
            v.push(Chunk { contig_end: end, a, b });
            a = b;
        }
    }
    v
}

pub fn scan<R: Roller, T: Send, F: Fn(u64, u64) -> Option<T> + Sync>(
    g: &Genome, k: usize, lo: u64, hi: u64, stride: u64, chunk_len: u64, make: F,
) -> Vec<T> {
    let _ = g.seq.advise(memmap2::Advice::Sequential);
    let seq: &[u8] = &g.seq;
    let parts: Vec<Vec<T>> = chunks(g, chunk_len)
        .par_iter()
        .map(|ch| {
            let mut out: Vec<T> = Vec::new();
            let k64 = k as u64;
            if ch.contig_end < ch.a + k64 { return out; }
            let feed_end = (ch.b + k64 - 1).min(ch.contig_end);
            let mut r = R::new(k);
            let mut valid_from = ch.a;
            let mut q = ch.a;
            while q < feed_end {
                let raw = seq[q as usize];
                let (c, is_n) = if raw > 3 { (0u8, true) } else { (raw, false) };
                if is_n { valid_from = q + 1; }
                let full = q >= ch.a + k64;
                let c_old = if full { let o = seq[(q - k64) as usize]; if o > 3 { 0 } else { o } } else { 0 };
                r.push(c, c_old, full);
                if q + 1 >= ch.a + k64 {
                    let s = q + 1 - k64;
                    if s >= valid_from && (stride == 1 || s % stride == 0) {
                        let key = r.key();
                        if key >= lo && key <= hi { if let Some(v) = make(key, s) { out.push(v); } }
                    }
                }
                q += 1;
            }
            out
        })
        .collect();
    // drop the genome's file pages from the resident set: the next scan re-reads them from disk/cache
    let _ = unsafe { g.seq.unchecked_advise(memmap2::UncheckedAdvice::DontNeed) };
    let total: usize = parts.iter().map(|p| p.len()).sum();
    let mut all = Vec::with_capacity(total);
    for p in parts { all.extend(p); }
    all
}

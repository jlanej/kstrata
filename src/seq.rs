//! Genome preparation and loading. A prepared genome is a flat byte array of base codes
//! (A=0, C=1, G=2, T=3, anything else=4) over all contigs concatenated, plus a TSV index.
use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use memmap2::Mmap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

pub const N_CODE: u8 = 4;

static CODE: [u8; 256] = {
    let mut t = [N_CODE; 256];
    t[b'A' as usize] = 0; t[b'a' as usize] = 0;
    t[b'C' as usize] = 1; t[b'c' as usize] = 1;
    t[b'G' as usize] = 2; t[b'g' as usize] = 2;
    t[b'T' as usize] = 3; t[b't' as usize] = 3;
    t
};

#[derive(Clone, Debug)]
pub struct Contig { pub name: String, pub start: u64, pub len: u64 }

pub struct Genome { pub name: String, pub seq: Mmap, pub contigs: Vec<Contig> }

impl Genome {
    pub fn len(&self) -> u64 { self.seq.len() as u64 }
}

/// FASTA (plain or gzip) -> <prefix>.seq (codes) + <prefix>.idx (name, start, len).
pub fn prep(fasta: &Path, prefix: &Path, include: Option<&str>) -> Result<()> {
    let f = File::open(fasta).with_context(|| format!("open {}", fasta.display()))?;
    let reader: Box<dyn Read> = if fasta.extension().map(|e| e == "gz").unwrap_or(false) {
        Box::new(MultiGzDecoder::new(f))
    } else { Box::new(f) };
    let mut rd = BufReader::with_capacity(1 << 22, reader);
    let mut seq = BufWriter::with_capacity(1 << 22, File::create(prefix.with_extension("seq"))?);
    let mut idx = BufWriter::new(File::create(prefix.with_extension("idx"))?);
    let mut line = Vec::with_capacity(1 << 16);
    let mut offset: u64 = 0;
    let mut cur: Option<(String, u64)> = None;
    let mut buf = Vec::with_capacity(1 << 16);
    loop {
        line.clear();
        let n = rd.read_until(b'\n', &mut line)?;
        if n == 0 { break; }
        if line[0] == b'>' {
            if let Some((name, start)) = cur.take() {
                writeln!(idx, "{}\t{}\t{}", name, start, offset - start)?;
            }
            let hdr = std::str::from_utf8(&line[1..])?.trim();
            let name = hdr.split_whitespace().next().unwrap_or("").to_string();
            let keep = include.map(|s| name.contains(s)).unwrap_or(true);
            cur = if keep { Some((name, offset)) } else { None };
        } else if cur.is_none() {
            continue;
        } else {
            buf.clear();
            buf.extend(line.iter().filter(|&&b| b > 32).map(|&b| CODE[b as usize]));
            seq.write_all(&buf)?;
            offset += buf.len() as u64;
        }
    }
    if let Some((name, start)) = cur.take() {
        writeln!(idx, "{}\t{}\t{}", name, start, offset - start)?;
    }
    seq.flush()?; idx.flush()?;
    Ok(())
}

pub fn load(prefix: &Path, name: &str) -> Result<Genome> {
    let f = File::open(prefix.with_extension("seq")).with_context(|| format!("open {}.seq", prefix.display()))?;
    let seq = unsafe { Mmap::map(&f)? };
    let idx = std::fs::read_to_string(prefix.with_extension("idx"))?;
    let mut contigs = Vec::new();
    for l in idx.lines() {
        let mut it = l.split('\t');
        let name = it.next().context("idx name")?.to_string();
        let start: u64 = it.next().context("idx start")?.parse()?;
        let len: u64 = it.next().context("idx len")?.parse()?;
        contigs.push(Contig { name, start, len });
    }
    Ok(Genome { name: name.to_string(), seq, contigs })
}

/// Uniform random genome of `length` bases in contigs of `contig` bases (splitmix64 stream).
pub fn random(prefix: &Path, length: u64, seed: u64, contig: u64) -> Result<()> {
    let mut seq = BufWriter::with_capacity(1 << 22, File::create(prefix.with_extension("seq"))?);
    let mut idx = BufWriter::new(File::create(prefix.with_extension("idx"))?);
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(0x1234_5678);
    let mut buf = vec![0u8; 1 << 20];
    let mut written = 0u64;
    let mut ci = 0;
    while written < length {
        let len = contig.min(length - written);
        writeln!(idx, "rand{}\t{}\t{}", ci, written, len)?;
        let mut left = len;
        while left > 0 {
            let n = (buf.len() as u64).min(left) as usize;
            for i in (0..n).step_by(32) {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
                z ^= z >> 31;
                for j in 0..32.min(n - i) { buf[i + j] = ((z >> (2 * j)) & 3) as u8; }
            }
            seq.write_all(&buf[..n])?;
            left -= n as u64;
        }
        written += len; ci += 1;
    }
    seq.flush()?; idx.flush()?;
    Ok(())
}

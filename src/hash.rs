//! Canonical k-mer keys. k <= 32: the exact 2-bit packing of min(kmer, revcomp) passed through a
//! bijective 64-bit mixer (no collisions). k > 32: polynomial rolling hash modulo 2^61-1 of the
//! k-mer and of its reverse complement, canonical = min, mixed the same way. The mixer makes the
//! top bits uniform so the key space can be partitioned by prefix.
pub const P61: u64 = (1u64 << 61) - 1;
/// Fixed random base for the polynomial hash (reproducible). Any value in (4, P61).
pub const BASE: u64 = 0x1A2B_3C4D_5E6F_7081 % P61;

#[inline(always)]
pub fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

#[inline(always)]
pub fn mulmod(a: u64, b: u64) -> u64 {
    let x = (a as u128) * (b as u128);
    let mut r = ((x as u64) & P61) + ((x >> 61) as u64);
    r = (r & P61) + (r >> 61);
    if r >= P61 { r - P61 } else { r }
}

#[inline(always)]
fn addmod(a: u64, b: u64) -> u64 { let s = a + b; if s >= P61 { s - P61 } else { s } }

pub fn powmod(mut b: u64, mut e: u64) -> u64 {
    let mut r = 1u64;
    while e > 0 { if e & 1 == 1 { r = mulmod(r, b); } b = mulmod(b, b); e >>= 1; }
    r
}

pub trait Roller {
    fn new(k: usize) -> Self;
    /// Feed the next base code (0..3). `c_old` is the code leaving the window when `full`.
    fn push(&mut self, c_new: u8, c_old: u8, full: bool);
    fn key(&self) -> u64;
}

pub struct Exact { fwd: u64, rc: u64, mask: u64, shift: u32 }

impl Roller for Exact {
    fn new(k: usize) -> Self {
        assert!(k >= 1 && k <= 32);
        let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
        Exact { fwd: 0, rc: 0, mask, shift: (2 * (k - 1)) as u32 }
    }
    #[inline(always)]
    fn push(&mut self, c_new: u8, _c_old: u8, _full: bool) {
        self.fwd = ((self.fwd << 2) | c_new as u64) & self.mask;
        self.rc = (self.rc >> 2) | (((3 - c_new) as u64) << self.shift);
    }
    #[inline(always)]
    fn key(&self) -> u64 { mix64(self.fwd.min(self.rc)) }
}

pub struct Poly { fwd: u64, rc: u64, binv: u64, powj: u64, d_bk1: [u64; 4], dc_bk1: [u64; 4] }

impl Roller for Poly {
    fn new(k: usize) -> Self {
        let bk1 = powmod(BASE, (k - 1) as u64);
        let binv = powmod(BASE, P61 - 2);
        let mut d_bk1 = [0u64; 4];
        let mut dc_bk1 = [0u64; 4];
        for c in 0..4u64 {
            d_bk1[c as usize] = mulmod(c + 1, bk1);
            dc_bk1[c as usize] = mulmod((3 - c) + 1, bk1);
        }
        Poly { fwd: 0, rc: 0, binv, powj: 1, d_bk1, dc_bk1 }
    }
    #[inline(always)]
    fn push(&mut self, c_new: u8, c_old: u8, full: bool) {
        let d = c_new as u64 + 1;
        let dc = (3 - c_new) as u64 + 1;
        if !full {
            self.fwd = addmod(mulmod(self.fwd, BASE), d);
            self.rc = addmod(self.rc, mulmod(dc, self.powj));
            self.powj = mulmod(self.powj, BASE);
        } else {
            let f = addmod(self.fwd, P61 - self.d_bk1[c_old as usize]);
            self.fwd = addmod(mulmod(f, BASE), d);
            let dco = (3 - c_old) as u64 + 1;
            let r = addmod(self.rc, P61 - dco);
            self.rc = addmod(mulmod(r, self.binv), self.dc_bk1[c_new as usize]);
        }
    }
    #[inline(always)]
    fn key(&self) -> u64 { mix64(self.fwd.min(self.rc)) }
}

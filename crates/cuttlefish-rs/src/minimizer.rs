//! Reference minimizer implementation, exercised only by tests.
//!
//! Production partitioning uses the fused rolling scan in
//! [`crate::partition`] (`for_each_weak_superkmer_impl`), which never calls
//! into this module. Kept as the readable specification the fused scan is
//! checked against.

use crate::dna::Base;
use crate::hash::wyhash_u64;
use crate::kmer::Kmer;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Minimizer {
    pub lmer: u64,
    pub hash: u64,
    pub offset: usize,
}

impl Minimizer {
    #[inline]
    fn candidate(lmer: u64, offset: usize, seed: u64) -> Self {
        Self {
            lmer,
            hash: wyhash_u64(lmer, seed),
            offset,
        }
    }
}

impl Ord for Minimizer {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.hash, self.lmer, self.offset).cmp(&(other.hash, other.lmer, other.offset))
    }
}

impl PartialOrd for Minimizer {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

pub fn minimizer<const K: usize>(kmer: Kmer<K>, l: usize, seed: u64) -> Minimizer {
    assert!(l > 0 && l <= K && l <= 32);
    let mut lmer = 0u64;
    for idx in 0..l {
        lmer = (lmer << 2) | kmer.get(idx).bits() as u64;
    }

    let mut min = Minimizer::candidate(lmer, 0, seed);
    let mask = if l == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * l)) - 1
    };

    for offset in 1..=K - l {
        lmer = ((lmer << 2) | kmer.get(offset + l - 1).bits() as u64) & mask;
        let cand = Minimizer::candidate(lmer, offset, seed);
        if cand < min {
            min = cand;
        }
    }

    min
}

pub fn canonical_minimizer<const K: usize>(kmer: Kmer<K>, l: usize, seed: u64) -> Minimizer {
    let fwd = minimizer(kmer, l, seed);
    let rev = minimizer(kmer.reverse_complement(), l, seed);
    fwd.min(rev)
}

pub fn encode_lmer_ascii(seq: &[u8]) -> Option<u64> {
    if seq.len() > 32 {
        return None;
    }

    let mut value = 0u64;
    for &ch in seq {
        let base = Base::from_ascii(ch);
        if !base.is_dna() {
            return None;
        }
        value = (value << 2) | base.bits() as u64;
    }
    Some(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minimizer_chooses_lowest_hash_then_lmer_then_offset() {
        let k = Kmer::<8>::from_ascii(b"ACGTACGT").unwrap();
        let m = minimizer(k, 3, 0);
        for off in 0..=5 {
            let lmer = encode_lmer_ascii(&b"ACGTACGT"[off..off + 3]).unwrap();
            let cand = Minimizer {
                lmer,
                hash: wyhash_u64(lmer, 0),
                offset: off,
            };
            assert!(m <= cand);
        }
    }

    #[test]
    fn canonical_is_symmetric_under_reverse_complement() {
        let k = Kmer::<15>::from_ascii(b"ACGTTTTACGTACGA").unwrap();
        assert_eq!(
            canonical_minimizer(k, 5, 0),
            canonical_minimizer(k.reverse_complement(), 5, 0)
        );
    }
}

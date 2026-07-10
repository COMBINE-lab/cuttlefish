use crate::dna::Base;
use crate::hash::{hash_bytes, hash_u64};
use std::hash::{Hash, Hasher};

#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Kmer<const K: usize> {
    words: [u64; 2],
}

impl<const K: usize> Kmer<K> {
    #[inline]
    pub fn zero() -> Self {
        Self { words: [0, 0] }
    }

    pub fn from_ascii(seq: &[u8]) -> Result<Self, KmerError> {
        Self::check_k()?;
        if seq.len() < K {
            return Err(KmerError::TooShort {
                expected: K,
                got: seq.len(),
            });
        }

        let mut value = 0u128;
        for &ch in &seq[..K] {
            let base = Base::from_ascii(ch);
            if !base.is_dna() {
                return Err(KmerError::InvalidBase(ch));
            }
            value = (value << 2) | base.bits() as u128;
        }

        Ok(Self::from_u128(value))
    }

    #[inline]
    pub fn as_u128(self) -> u128 {
        self.words[0] as u128 | ((self.words[1] as u128) << 64)
    }

    #[inline]
    pub fn words(self) -> [u64; 2] {
        self.words
    }

    #[inline]
    pub fn get(self, idx: usize) -> Base {
        assert!(idx < K);
        let shift = 2 * (K - 1 - idx);
        match ((self.as_u128() >> shift) & 0b11) as u8 {
            0 => Base::A,
            1 => Base::C,
            2 => Base::G,
            3 => Base::T,
            _ => unreachable!(),
        }
    }

    #[inline]
    pub fn front(self) -> Base {
        self.get(0)
    }

    #[inline]
    pub fn back(self) -> Base {
        self.get(K - 1)
    }

    pub fn reverse_complement(self) -> Self {
        let mut v = 0u128;
        let mut x = self.as_u128();
        for _ in 0..K {
            v = (v << 2) | ((!x) & 0b11);
            x >>= 2;
        }
        Self::from_u128(v)
    }

    #[inline]
    pub fn canonical(self) -> Self {
        self.min(self.reverse_complement())
    }

    #[inline]
    pub fn is_canonical(self) -> bool {
        self <= self.reverse_complement()
    }

    pub fn roll_forward(self, next: Base) -> Self {
        assert!(next.is_dna());
        let mask = if K == 64 {
            u128::MAX
        } else {
            (1u128 << (2 * K)) - 1
        };
        Self::from_u128(((self.as_u128() << 2) | next.bits() as u128) & mask)
    }

    pub fn roll_backward(self, prev: Base) -> Self {
        assert!(prev.is_dna());
        let high = (prev.bits() as u128) << (2 * (K - 1));
        Self::from_u128(high | (self.as_u128() >> 2))
    }

    pub fn to_ascii_string(self) -> String {
        let mut s = Vec::with_capacity(K);
        for idx in 0..K {
            s.push(self.get(idx).to_ascii());
        }
        String::from_utf8(s).unwrap()
    }

    pub fn hash64(self, seed: u64) -> u64 {
        let byte_len = (2 * K).div_ceil(8);
        if byte_len <= 8 {
            return hash_u64(self.words[0], seed);
        }
        let bytes = self.as_u128().to_le_bytes();
        hash_bytes(&bytes[..byte_len], seed)
    }

    #[inline]
    fn from_u128(value: u128) -> Self {
        Self {
            words: [value as u64, (value >> 64) as u64],
        }
    }

    #[inline]
    pub(crate) fn from_bits(value: u128) -> Self {
        Self::from_u128(value)
    }

    #[inline]
    fn check_k() -> Result<(), KmerError> {
        if K == 0 || K > 63 {
            Err(KmerError::UnsupportedK(K))
        } else {
            Ok(())
        }
    }
}

impl<const K: usize> Ord for Kmer<K> {
    #[inline]
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.as_u128().cmp(&other.as_u128())
    }
}

impl<const K: usize> PartialOrd for Kmer<K> {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<const K: usize> Hash for Kmer<K> {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u64(self.words[0]);
        state.write_u64(self.words[1]);
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KmerError {
    UnsupportedK(usize),
    TooShort { expected: usize, got: usize },
    InvalidBase(u8),
}

impl std::fmt::Display for KmerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedK(k) => write!(f, "unsupported k-mer length {k}; expected 1..=63"),
            Self::TooShort { expected, got } => {
                write!(
                    f,
                    "sequence too short for k-mer: expected {expected}, got {got}"
                )
            }
            Self::InvalidBase(b) => write!(f, "invalid DNA base '{}'", *b as char),
        }
    }
}

impl std::error::Error for KmerError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_ascii() {
        let k = Kmer::<31>::from_ascii(b"ACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
        assert_eq!(k.to_ascii_string(), "ACGTACGTACGTACGTACGTACGTACGTACG");
        assert_eq!(k.front(), Base::A);
        assert_eq!(k.back(), Base::G);
    }

    #[test]
    fn reverse_complement_round_trips() {
        let k = Kmer::<7>::from_ascii(b"ACGTAGC").unwrap();
        assert_eq!(k.reverse_complement().to_ascii_string(), "GCTACGT");
        assert_eq!(k.reverse_complement().reverse_complement(), k);
    }

    #[test]
    fn rolling_matches_reparse() {
        let k = Kmer::<5>::from_ascii(b"ACGTA").unwrap();
        assert_eq!(
            k.roll_forward(Base::C),
            Kmer::<5>::from_ascii(b"CGTAC").unwrap()
        );
        assert_eq!(
            k.roll_backward(Base::T),
            Kmer::<5>::from_ascii(b"TACGT").unwrap()
        );
    }

    #[test]
    fn one_word_hash_matches_byte_hash() {
        let k = Kmer::<31>::from_ascii(b"ACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
        for seed in [0, 1, 17, u64::MAX] {
            assert_eq!(
                k.hash64(seed),
                hash_bytes(&k.words()[0].to_le_bytes(), seed)
            );
        }
    }
}

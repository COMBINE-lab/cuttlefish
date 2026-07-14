use crate::dna::Base;
use crate::hash::{hash_bytes, wyhash_u64};
use std::hash::{Hash, Hasher};

const fn ascii_quad_table() -> [u32; 256] {
    let mut table = [0u32; 256];
    let mut byte = 0usize;
    while byte < table.len() {
        let mut encoded = 0u32;
        let mut base = 0usize;
        while base < 4 {
            let bits = (byte >> (2 * (3 - base))) & 0b11;
            let ascii = match bits {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            encoded |= (ascii as u32) << (8 * base);
            base += 1;
        }
        table[byte] = encoded;
        byte += 1;
    }
    table
}

const ASCII_QUADS: [u32; 256] = ascii_quad_table();

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

    #[inline]
    pub fn reverse_complement(self) -> Self {
        if K <= 32 {
            let mut x = self.words[0];
            x = ((x & 0x3333_3333_3333_3333) << 2) | ((x >> 2) & 0x3333_3333_3333_3333);
            x = ((x & 0x0f0f_0f0f_0f0f_0f0f) << 4) | ((x >> 4) & 0x0f0f_0f0f_0f0f_0f0f);
            x = x.swap_bytes();
            let shift = 64 - 2 * K;
            let mask = if K == 32 {
                u64::MAX
            } else {
                (1u64 << (2 * K)) - 1
            };
            return Self::from_u128(((!x >> shift) & mask) as u128);
        }

        let mut x = self.as_u128();
        x = ((x & 0x3333_3333_3333_3333_3333_3333_3333_3333) << 2)
            | ((x >> 2) & 0x3333_3333_3333_3333_3333_3333_3333_3333);
        x = ((x & 0x0f0f_0f0f_0f0f_0f0f_0f0f_0f0f_0f0f_0f0f) << 4)
            | ((x >> 4) & 0x0f0f_0f0f_0f0f_0f0f_0f0f_0f0f_0f0f_0f0f);
        x = x.swap_bytes();
        let shift = 128 - 2 * K;
        let mask = (1u128 << (2 * K)) - 1;
        Self::from_u128((!x >> shift) & mask)
    }

    #[inline]
    pub fn canonical(self) -> Self {
        self.min(self.reverse_complement())
    }

    #[inline]
    pub fn is_canonical(self) -> bool {
        self <= self.reverse_complement()
    }

    #[inline]
    pub fn roll_forward(self, next: Base) -> Self {
        assert!(next.is_dna());
        if K <= 32 {
            let mask = if K == 32 {
                u64::MAX
            } else {
                (1u64 << (2 * K)) - 1
            };
            return Self::from_u128(
                (((self.words[0] << 2) | u64::from(next.bits())) & mask) as u128,
            );
        }
        let mask = if K == 64 {
            u128::MAX
        } else {
            (1u128 << (2 * K)) - 1
        };
        Self::from_u128(((self.as_u128() << 2) | next.bits() as u128) & mask)
    }

    #[inline]
    pub fn roll_backward(self, prev: Base) -> Self {
        assert!(prev.is_dna());
        if K <= 32 {
            let high = u64::from(prev.bits()) << (2 * (K - 1));
            return Self::from_u128((high | (self.words[0] >> 2)) as u128);
        }
        let high = (prev.bits() as u128) << (2 * (K - 1));
        Self::from_u128(high | (self.as_u128() >> 2))
    }

    pub fn to_ascii_string(self) -> String {
        let mut s = Vec::with_capacity(K);
        self.append_ascii(&mut s);
        String::from_utf8(s).unwrap()
    }

    #[inline]
    pub(crate) fn append_ascii(self, output: &mut Vec<u8>) {
        if K <= 32 {
            output.reserve(K);
            let leading = K & 3;
            for index in 0..leading {
                output.push(self.get(index).to_ascii());
            }
            let quads = K / 4;
            for index in 0..quads {
                let shift = 8 * (quads - 1 - index);
                let byte = ((self.words[0] >> shift) & 0xff) as usize;
                output.extend_from_slice(&ASCII_QUADS[byte].to_le_bytes());
            }
            return;
        }
        output.reserve(K);
        for index in 0..K {
            output.push(self.get(index).to_ascii());
        }
    }

    #[inline(always)]
    pub fn hash64(self, seed: u64) -> u64 {
        let byte_len = (2 * K).div_ceil(8);
        if byte_len <= 8 {
            return wyhash_u64(self.words[0], seed);
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
        if K <= 32 {
            self.words[0].cmp(&other.words[0])
        } else {
            self.words[1]
                .cmp(&other.words[1])
                .then_with(|| self.words[0].cmp(&other.words[0]))
        }
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

    fn slow_reverse_complement<const K: usize>(kmer: Kmer<K>) -> Kmer<K> {
        let mut value = 0u128;
        let mut input = kmer.as_u128();
        for _ in 0..K {
            value = (value << 2) | ((!input) & 0b11);
            input >>= 2;
        }
        Kmer::from_bits(value)
    }

    #[test]
    fn packed_reverse_complement_matches_basewise_reference() {
        fn check<const K: usize>() {
            let mask = (1u128 << (2 * K)) - 1;
            for value in [
                0,
                1,
                mask,
                0x0123_4567_89ab_cdef_7654_3210_fedc_ba98 & mask,
                0xaaaa_5555_3333_cccc_0f0f_f0f0_9696_6969 & mask,
            ] {
                let kmer = Kmer::<K>::from_bits(value);
                assert_eq!(kmer.reverse_complement(), slow_reverse_complement(kmer));
            }
        }

        check::<1>();
        check::<4>();
        check::<31>();
        check::<32>();
        check::<33>();
        check::<63>();
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
    fn one_word_hash_matches_cpp_wyhash() {
        let k = Kmer::<31>::from_ascii(b"ACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
        for seed in [0, 1, 17, u64::MAX] {
            assert_eq!(k.hash64(seed), wyhash_u64(k.words()[0], seed));
        }
    }
}

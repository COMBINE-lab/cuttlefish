#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Base {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    N = 4,
    E = 5,
}

pub const INVALID_BASE_BITS: u8 = 4;

pub const ASCII_BASE_BITS: [u8; 256] = {
    let mut bits = [INVALID_BASE_BITS; 256];
    bits[b'A' as usize] = 0;
    bits[b'a' as usize] = 0;
    bits[b'C' as usize] = 1;
    bits[b'c' as usize] = 1;
    bits[b'G' as usize] = 2;
    bits[b'g' as usize] = 2;
    bits[b'T' as usize] = 3;
    bits[b't' as usize] = 3;
    bits[b'U' as usize] = 3;
    bits[b'u' as usize] = 3;
    bits
};

impl Base {
    #[inline]
    pub const fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
            Self::N => Self::N,
            Self::E => Self::E,
        }
    }

    #[inline]
    pub const fn to_ascii(self) -> u8 {
        match self {
            Self::A => b'A',
            Self::C => b'C',
            Self::G => b'G',
            Self::T => b'T',
            Self::N => b'N',
            Self::E => b'$',
        }
    }

    #[inline]
    pub const fn from_ascii(byte: u8) -> Self {
        match byte {
            b'A' | b'a' => Self::A,
            b'C' | b'c' => Self::C,
            b'G' | b'g' => Self::G,
            b'T' | b't' | b'U' | b'u' => Self::T,
            _ => Self::N,
        }
    }

    #[inline]
    pub const fn is_dna(self) -> bool {
        matches!(self, Self::A | Self::C | Self::G | Self::T)
    }

    #[inline]
    pub const fn bits(self) -> u8 {
        self as u8
    }
}

#[inline]
pub const fn is_dna_ascii(byte: u8) -> bool {
    ASCII_BASE_BITS[byte as usize] != INVALID_BASE_BITS
}

/// ASCII complement for every byte, built from [`Base::from_ascii`] so the table
/// is exactly equivalent to computing it, including for non-DNA bytes.
///
/// Reverse-strand label assembly complements one base at a time in its inner
/// loop; going through the enum compiles to range checks and a jump table, while
/// this is a single indexed load.
pub const ASCII_COMPLEMENT: [u8; 256] = {
    let mut table = [0u8; 256];
    let mut byte = 0usize;
    while byte < 256 {
        table[byte] = Base::from_ascii(byte as u8).complement().to_ascii();
        byte += 1;
    }
    table
};

#[inline]
pub const fn complement_ascii(byte: u8) -> u8 {
    ASCII_COMPLEMENT[byte as usize]
}

#[inline]
pub const fn ascii_base_bits(byte: u8) -> Option<u8> {
    let bits = ASCII_BASE_BITS[byte as usize];
    if bits == INVALID_BASE_BITS {
        None
    } else {
        Some(bits)
    }
}

#[inline]
pub const fn valid_ascii_base_bits(byte: u8) -> u8 {
    ((byte >> 2) ^ (byte >> 1)) & 0b11
}

#[inline]
pub const fn ascii_complement_bits(byte: u8) -> Option<u8> {
    match ascii_base_bits(byte) {
        Some(bits) => Some(bits ^ 0b11),
        None => None,
    }
}

/// Reverse complement of an ASCII label, tolerating non-ACGT bytes.
#[inline]
pub(crate) fn reverse_complement_label(label: &[u8]) -> Vec<u8> {
    label
        .iter()
        .rev()
        .map(|&base| complement_ascii(base))
        .collect()
}

/// Lexicographically least rotation of `label`, as a fresh buffer.
pub(crate) fn minimal_rotation(label: &[u8]) -> Vec<u8> {
    let start = least_rotation_start(label);
    label[start..]
        .iter()
        .chain(label[..start].iter())
        .copied()
        .collect()
}

/// Booth's algorithm: start index of the lexicographically least rotation.
pub(crate) fn least_rotation_start(s: &[u8]) -> usize {
    let n = s.len();
    if n <= 1 {
        return 0;
    }

    let mut i = 0;
    let mut j = 1;
    let mut k = 0;
    while i < n && j < n && k < n {
        let a = s[(i + k) % n];
        let b = s[(j + k) % n];
        if a == b {
            k += 1;
        } else if a > b {
            i += k + 1;
            if i <= j {
                i = j + 1;
            }
            k = 0;
        } else {
            j += k + 1;
            if j <= i {
                j = i + 1;
            }
            k = 0;
        }
    }

    i.min(j)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ascii_mapping_matches_cpp_encoding() {
        assert_eq!(Base::from_ascii(b'A').bits(), 0);
        assert_eq!(Base::from_ascii(b'C').bits(), 1);
        assert_eq!(Base::from_ascii(b'G').bits(), 2);
        assert_eq!(Base::from_ascii(b'T').bits(), 3);
        assert_eq!(Base::from_ascii(b'N'), Base::N);
    }

    #[test]
    fn complement_table_matches_enum_path() {
        for byte in 0..=u8::MAX {
            assert_eq!(
                complement_ascii(byte),
                Base::from_ascii(byte).complement().to_ascii(),
                "complement table diverges at byte {byte}"
            );
        }
    }

    #[test]
    fn complements_are_involutions() {
        for b in [Base::A, Base::C, Base::G, Base::T, Base::N, Base::E] {
            assert_eq!(b.complement().complement(), b);
        }
    }
}

/// Bytes a label of `bases` takes packed 2 bits per base.
#[inline]
pub(crate) const fn packed_2bit_len(bases: usize) -> usize {
    bases.div_ceil(4)
}

/// Appends an ACGT label packed 2 bits per base, first base in the high bits
/// of the first byte, the last byte zero-padded. Each call starts on a byte
/// boundary, so labels packed one after another stay independently
/// addressable by byte offset.
pub(crate) fn pack_2bit_extend(dst: &mut Vec<u8>, ascii: &[u8]) {
    debug_assert!(
        ascii
            .iter()
            .all(|&base| matches!(base, b'A' | b'C' | b'G' | b'T'))
    );
    dst.reserve(packed_2bit_len(ascii.len()));
    let (quads, tail) = ascii.as_chunks::<4>();
    for quad in quads {
        dst.push(
            valid_ascii_base_bits(quad[0]) << 6
                | valid_ascii_base_bits(quad[1]) << 4
                | valid_ascii_base_bits(quad[2]) << 2
                | valid_ascii_base_bits(quad[3]),
        );
    }
    if !tail.is_empty() {
        let mut byte = 0u8;
        for (index, &base) in tail.iter().enumerate() {
            byte |= valid_ascii_base_bits(base) << (6 - 2 * index);
        }
        dst.push(byte);
    }
}

/// Appends the `len` ASCII bases of a label packed by [`pack_2bit_extend`].
pub(crate) fn unpack_2bit_extend(dst: &mut Vec<u8>, packed: &[u8], len: usize) {
    debug_assert!(packed.len() >= packed_2bit_len(len));
    let start = dst.len();
    dst.reserve(packed.len() * 4);
    for &byte in &packed[..packed_2bit_len(len)] {
        dst.extend_from_slice(&crate::kmer::ASCII_QUADS[usize::from(byte)].to_le_bytes());
    }
    dst.truncate(start + len);
}

#[cfg(test)]
mod packed_label_tests {
    use super::*;

    #[test]
    fn packed_labels_round_trip_at_every_length() {
        let mut state = 0x243f_6a88_85a3_08d3u64;
        for len in 0..=130usize {
            let label: Vec<u8> = (0..len)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    b"ACGT"[(state % 4) as usize]
                })
                .collect();
            let mut packed = vec![0xAAu8];
            pack_2bit_extend(&mut packed, &label);
            assert_eq!(packed.len(), 1 + packed_2bit_len(len));
            let mut out = b"prefix".to_vec();
            unpack_2bit_extend(&mut out, &packed[1..], len);
            assert_eq!(&out[..6], b"prefix");
            assert_eq!(&out[6..], &label[..], "length {len}");
        }
    }
}

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

#[inline]
pub const fn complement_ascii(byte: u8) -> u8 {
    Base::from_ascii(byte).complement().to_ascii()
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
    ASCII_BASE_BITS[byte as usize] & 0b11
}

#[inline]
pub const fn ascii_complement_bits(byte: u8) -> Option<u8> {
    match ascii_base_bits(byte) {
        Some(bits) => Some(bits ^ 0b11),
        None => None,
    }
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
    fn complements_are_involutions() {
        for b in [Base::A, Base::C, Base::G, Base::T, Base::N, Base::E] {
            assert_eq!(b.complement().complement(), b);
        }
    }
}

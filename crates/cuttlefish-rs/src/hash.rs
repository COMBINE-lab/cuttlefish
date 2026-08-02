use std::hash::{BuildHasherDefault, Hasher};

const WY_SALT: [u64; 4] = [
    4167021922371662411,
    7320285940802167691,
    14307255741305819987,
    10859488101230029397,
];

#[inline]
fn mum(a: u64, b: u64) -> u64 {
    let r = (a as u128).wrapping_mul(b as u128);
    (r as u64) ^ (r >> 64) as u64
}

#[inline]
fn mum_parts(a: u64, b: u64) -> (u64, u64) {
    let product = (a as u128).wrapping_mul(b as u128);
    (product as u64, (product >> 64) as u64)
}

#[inline]
fn mix(a: u64, b: u64) -> u64 {
    mum(a ^ WY_SALT[0], b ^ WY_SALT[1])
}

#[inline]
fn read_u64_le(bytes: &[u8]) -> u64 {
    let mut buf = [0u8; 8];
    buf[..bytes.len()].copy_from_slice(bytes);
    u64::from_le_bytes(buf)
}

pub fn hash_bytes(bytes: &[u8], seed: u64) -> u64 {
    let mut h = seed ^ WY_SALT[2] ^ bytes.len() as u64;
    let mut chunks = bytes.chunks_exact(16);
    for chunk in &mut chunks {
        let a = u64::from_le_bytes(chunk[..8].try_into().unwrap());
        let b = u64::from_le_bytes(chunk[8..16].try_into().unwrap());
        h = mix(h ^ a, b);
    }

    let rem = chunks.remainder();
    if rem.len() > 8 {
        h = mix(h ^ read_u64_le(&rem[..8]), read_u64_le(&rem[8..]));
    } else if !rem.is_empty() {
        h = mix(h, read_u64_le(rem));
    }

    mum(h ^ WY_SALT[3], h.rotate_left(32) ^ WY_SALT[0])
}

pub type FastBuildHasher = BuildHasherDefault<FastHasher>;

#[derive(Clone, Default)]
pub struct FastHasher {
    state: u64,
}

impl Hasher for FastHasher {
    #[inline]
    fn finish(&self) -> u64 {
        mum(
            self.state ^ WY_SALT[3],
            self.state.rotate_left(32) ^ WY_SALT[0],
        )
    }

    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut chunks = bytes.chunks_exact(8);
        for chunk in &mut chunks {
            self.write_u64(u64::from_le_bytes(chunk.try_into().unwrap()));
        }
        let rem = chunks.remainder();
        if !rem.is_empty() {
            self.write_u64(read_u64_le(rem));
        }
    }

    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.write_u64(i as u64);
    }

    #[inline]
    fn write_u16(&mut self, i: u16) {
        self.write_u64(i as u64);
    }

    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.write_u64(i as u64);
    }

    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.state = mix(self.state ^ i, i.rotate_left(32));
    }

    #[inline]
    fn write_u128(&mut self, i: u128) {
        self.write_u64(i as u64);
        self.write_u64((i >> 64) as u64);
    }

    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.write_u64(i as u64);
    }
}

#[inline]
pub fn hash_u64(value: u64, seed: u64) -> u64 {
    let h = seed ^ WY_SALT[2] ^ 8;
    let h = mix(h, value);
    mum(h ^ WY_SALT[3], h.rotate_left(32) ^ WY_SALT[0])
}

// Exact len=8 specialization of the wyhash version used by C++ minimizers.
#[inline]
pub fn wyhash_u64(value: u64, mut seed: u64) -> u64 {
    seed ^= mum(seed ^ WY_SALT[0], WY_SALT[1]);
    let a = value.rotate_left(32) ^ WY_SALT[1];
    let b = value ^ seed;
    let (a, b) = mum_parts(a, b);
    mum(a ^ WY_SALT[0] ^ 8, b ^ WY_SALT[1])
}

#[inline]
pub fn hash_two_u64(first: u64, second: u64) -> u64 {
    let h = mix(first, first.rotate_left(32));
    let h = mix(h ^ second, second.rotate_left(32));
    mum(h ^ WY_SALT[3], h.rotate_left(32) ^ WY_SALT[0])
}

#[inline]
pub fn fast_u64_hash(value: u64, seed: u64) -> u64 {
    let value = value ^ seed;
    let mixed = (value ^ (value >> 23)).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    mixed ^ (mixed >> 29)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hash_is_seeded_and_stable() {
        assert_eq!(hash_u64(42, 0), hash_u64(42, 0));
        assert_ne!(hash_u64(42, 0), hash_u64(42, 1));
    }

    #[test]
    fn hash_u64_matches_byte_hash() {
        for seed in [0, 1, 17, u64::MAX] {
            for value in [0, 1, 42, u32::MAX as u64, u64::MAX] {
                assert_eq!(
                    hash_u64(value, seed),
                    hash_bytes(&value.to_le_bytes(), seed)
                );
            }
        }
    }

    #[test]
    fn hash_two_u64_matches_fast_hasher_writes() {
        for first in [0, 1, 42, u64::MAX] {
            for second in [0, 7, u32::MAX as u64, u64::MAX] {
                let mut hasher = FastHasher::default();
                hasher.write_u64(first);
                hasher.write_u64(second);
                assert_eq!(hash_two_u64(first, second), hasher.finish());
            }
        }
    }

    #[test]
    fn wyhash_u64_matches_cpp_minimizer_hash() {
        assert_eq!(wyhash_u64(0, 0), 7_824_564_342_066_666_581);
        assert_eq!(wyhash_u64(1, 0), 17_550_070_827_293_468_892);
        assert_eq!(wyhash_u64(42, 0), 7_038_757_716_059_073_700);
        assert_eq!(
            wyhash_u64(0x1234_5678_9abc_def0, 0),
            6_800_848_065_093_387_582
        );
    }
}

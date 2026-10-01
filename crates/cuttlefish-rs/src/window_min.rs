//! The partition's minimizer scan: the minimum canonical l-mer hash of every
//! `window`-base window of a fragment.
//!
//! A window's minimum hash picks its subgraph, so this runs over every base of
//! the input. The vectorized form follows simd-minimizers -- R. Groot Koerkamp
//! and I. Martayan, "SimdMinimizers: Computing random minimizers, fast", SEA
//! 2025, doi:10.4230/LIPIcs.SEA.2025.20, <https://crates.io/crates/simd-minimizers>
//! -- in three respects:
//!
//! * The windows are split into eight contiguous chunks, one per 32-bit lane
//!   (one AVX2 register, or two NEON registers), and the lanes scan their
//!   chunks in lock step. Each lane reads its own bases -- with a gather on
//!   AVX2, with plain loads on NEON, which has none -- so the sequence is
//!   never transposed.
//! * The sliding minimum is the two-stacks scheme: a running prefix minimum
//!   of the current block of `w` hashes, plus suffix minima of the previous
//!   block, recomputed once every `w` steps. That is about three minimum
//!   operations per l-mer, with no data-dependent branches.
//! * l-mers are hashed with a cheap function that vectorizes (a 32-bit
//!   multiply-xorshift) rather than 64-bit wyhash, whose 128-bit products have
//!   no AVX2 counterpart.
//!
//! It departs from simd-minimizers in one respect. simd-minimizers compares
//! only the upper 16 bits of each hash, with the position in the lower 16,
//! so ties are broken by position -- differently on the two strands. Here the
//! subgraph is taken from the minimum hash *value*, which must be the same
//! for a window and its reverse complement, so whole hashes are compared.
//! The hash is a bijection on canonical l-mers (`l <= 16`), so equal hashes
//! mean equal l-mers and the minimum value is strand-invariant by
//! construction.
//!
//! The scalar scan computes exactly the same values. The AVX2 scan is chosen
//! at run time and the NEON scan is always available on aarch64, so output
//! never depends on the CPU.

/// Longest l-mer this scan hashes: a canonical l-mer must fit 32 bits.
pub(crate) const MAX_LMER_LEN: usize = 16;

/// Fewest windows per lane worth vectorizing. Each lane first spends
/// `window - 1` steps filling its window, so very short chunks mostly warm up.
#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
const MIN_WINDOWS_PER_LANE: usize = 32;

/// Mixed into each l-mer before hashing. The finalizer maps 0 to 0, which
/// would make the poly-A l-mer -- canonical value 0, and common in
/// assemblies -- the minimum of every window holding it, piling all of them
/// into one subgraph.
const LMER_HASH_SEED: u32 = 0x9e37_79b9;

/// Hashes a canonical l-mer: the MurmurHash3 32-bit finalizer of the seeded
/// l-mer, a bijection on `u32` whose low bits -- which pick the subgraph --
/// mix every input bit.
#[inline(always)]
pub(crate) fn lmer_hash(canonical: u32) -> u32 {
    let mut h = canonical ^ LMER_HASH_SEED;
    h ^= h >> 16;
    h = h.wrapping_mul(0x85eb_ca6b);
    h ^= h >> 13;
    h = h.wrapping_mul(0xc2b2_ae35);
    h ^ (h >> 16)
}

#[inline(always)]
fn base_bits(byte: u8) -> u32 {
    // A, C, G, T -> 0, 1, 2, 3; the input is known to be ACGT.
    u32::from(((byte >> 2) ^ (byte >> 1)) & 0b11)
}

/// Fills `out[j]` with the minimum canonical l-mer hash of the `window`-base
/// window starting at base `first + j` of `seq`.
///
/// `seq` must be ACGT only, `1 <= l <= MAX_LMER_LEN`, `l <= window`, and the
/// windows must lie inside `seq`.
pub(crate) fn fill_window_mins(seq: &[u8], l: usize, window: usize, first: usize, out: &mut [u32]) {
    assert!((1..=MAX_LMER_LEN).contains(&l) && l <= window && window - l < 64);
    assert!(first + out.len() + window - 1 <= seq.len());
    let done = fill_window_mins_simd(seq, l, window, first, out);
    window_mins_scalar(seq, l, window, first + done, &mut out[done..]);
}

/// Fills a prefix of `out` with the vectorized scan where the CPU allows it,
/// returning how many windows it covered; the scalar scan does the rest.
#[cfg(target_arch = "x86_64")]
fn fill_window_mins_simd(
    seq: &[u8],
    l: usize,
    window: usize,
    first: usize,
    out: &mut [u32],
) -> usize {
    if !std::arch::is_x86_feature_detected!("avx2") || seq.len() >= i32::MAX as usize {
        return 0;
    }
    // The last lane's final gather reads four bytes from up to
    // `8 * lane_len + window + 1` past `first`; keep that inside `seq`.
    let available = seq.len() - first - window + 1;
    let lane_len = out.len().min(available.saturating_sub(3)) / 8;
    if lane_len < MIN_WINDOWS_PER_LANE {
        return 0;
    }
    // SAFETY: AVX2 was just detected, and the bound above keeps every gather
    // inside `seq`.
    unsafe { window_mins_avx2(seq, l, window, first, lane_len, &mut out[..8 * lane_len]) };
    8 * lane_len
}

/// The NEON form of the eight-lane scan; NEON is part of the aarch64 baseline.
#[cfg(target_arch = "aarch64")]
fn fill_window_mins_simd(
    seq: &[u8],
    l: usize,
    window: usize,
    first: usize,
    out: &mut [u32],
) -> usize {
    // The same bound as AVX2: the last lane's final four-byte load stays
    // inside `seq`.
    let available = seq.len() - first - window + 1;
    let lane_len = out.len().min(available.saturating_sub(3)) / 8;
    if lane_len < MIN_WINDOWS_PER_LANE {
        return 0;
    }
    // SAFETY: NEON is always present on aarch64.
    unsafe { window_mins_neon(seq, l, window, first, lane_len, &mut out[..8 * lane_len]) };
    8 * lane_len
}

/// Other architectures have no vectorized scan: the scalar scan covers every
/// window.
#[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
fn fill_window_mins_simd(
    _seq: &[u8],
    _l: usize,
    _window: usize,
    _first: usize,
    _out: &mut [u32],
) -> usize {
    0
}

fn window_mins_scalar(seq: &[u8], l: usize, window: usize, first: usize, out: &mut [u32]) {
    if out.is_empty() {
        return;
    }
    let w = window - l + 1;
    let mask = if l == 16 {
        u32::MAX
    } else {
        (1 << (2 * l)) - 1
    };
    let rev_shift = 2 * (l - 1);
    let mut ring = [u32::MAX; 64];
    let mut prefix = u32::MAX;
    let mut slot = 0;
    let (mut fwd, mut rev) = (0u32, 0u32);
    for (t, &byte) in seq[first..first + out.len() + window - 1]
        .iter()
        .enumerate()
    {
        let bits = base_bits(byte);
        fwd = ((fwd << 2) | bits) & mask;
        rev = (rev >> 2) | ((bits ^ 0b11) << rev_shift);
        if t + 1 < l {
            continue;
        }
        let hash = lmer_hash(fwd.min(rev));
        ring[slot] = hash;
        prefix = prefix.min(hash);
        slot += 1;
        if slot == w {
            slot = 0;
            for j in (0..w - 1).rev() {
                ring[j] = ring[j].min(ring[j + 1]);
            }
            prefix = u32::MAX;
        }
        if t + 1 >= window {
            out[t + 1 - window] = prefix.min(ring[slot]);
        }
    }
}

/// The scalar scan run across eight lanes: lane `i` covers the `lane_len`
/// windows starting at `first + i * lane_len`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn window_mins_avx2(
    seq: &[u8],
    l: usize,
    window: usize,
    first: usize,
    lane_len: usize,
    out: &mut [u32],
) {
    use std::arch::x86_64::*;

    let w = window - l + 1;
    let mask = _mm256_set1_epi32(if l == 16 {
        -1
    } else {
        ((1u32 << (2 * l)) - 1) as i32
    });
    let rev_shift = _mm_cvtsi32_si128(2 * (l as i32 - 1));
    let three = _mm256_set1_epi32(0b11);
    let byte_mask = _mm256_set1_epi32(0xff);
    let max = _mm256_set1_epi32(-1);
    let mut ring = [max; 64];
    let mut prefix = max;
    let mut slot = 0;
    let (mut fwd, mut rev) = (_mm256_setzero_si256(), _mm256_setzero_si256());
    let lane = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    let mut index = _mm256_add_epi32(
        _mm256_set1_epi32(first as i32),
        _mm256_mullo_epi32(lane, _mm256_set1_epi32(lane_len as i32)),
    );
    let four = _mm256_set1_epi32(4);
    let seed = _mm256_set1_epi32(LMER_HASH_SEED as i32);
    let mut lanes = [0u32; 8];
    let steps = lane_len + window - 1;
    let mut t = 0;
    while t < steps {
        // Four bases per lane per gather.
        // SAFETY: `fill_window_mins` bounds every index so all four bytes
        // lie inside `seq`.
        let mut bytes = unsafe { _mm256_i32gather_epi32::<1>(seq.as_ptr().cast(), index) };
        index = _mm256_add_epi32(index, four);
        for _ in 0..4.min(steps - t) {
            let byte = _mm256_and_si256(bytes, byte_mask);
            bytes = _mm256_srli_epi32::<8>(bytes);
            let bits = _mm256_and_si256(
                _mm256_xor_si256(_mm256_srli_epi32::<2>(byte), _mm256_srli_epi32::<1>(byte)),
                three,
            );
            fwd = _mm256_and_si256(_mm256_or_si256(_mm256_slli_epi32::<2>(fwd), bits), mask);
            rev = _mm256_or_si256(
                _mm256_srli_epi32::<2>(rev),
                _mm256_sll_epi32(_mm256_xor_si256(bits, three), rev_shift),
            );
            t += 1;
            if t < l {
                continue;
            }
            let mut h = _mm256_xor_si256(_mm256_min_epu32(fwd, rev), seed);
            h = _mm256_xor_si256(h, _mm256_srli_epi32::<16>(h));
            h = _mm256_mullo_epi32(h, _mm256_set1_epi32(0x85eb_ca6bu32 as i32));
            h = _mm256_xor_si256(h, _mm256_srli_epi32::<13>(h));
            h = _mm256_mullo_epi32(h, _mm256_set1_epi32(0xc2b2_ae35u32 as i32));
            h = _mm256_xor_si256(h, _mm256_srli_epi32::<16>(h));
            ring[slot] = h;
            prefix = _mm256_min_epu32(prefix, h);
            slot += 1;
            if slot == w {
                slot = 0;
                for j in (0..w - 1).rev() {
                    ring[j] = _mm256_min_epu32(ring[j], ring[j + 1]);
                }
                prefix = max;
            }
            if t >= window {
                let min = _mm256_min_epu32(prefix, ring[slot]);
                // SAFETY: `lanes` is eight u32s.
                unsafe { _mm256_storeu_si256(lanes.as_mut_ptr().cast(), min) };
                let at = t - window;
                for (i, &value) in lanes.iter().enumerate() {
                    out[i * lane_len + at] = value;
                }
            }
        }
    }
}

/// The AVX2 scan's eight lanes as two four-lane NEON registers, `[lanes
/// 0-3, lanes 4-7]`, so lane `i` covers the same windows as it does there.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn window_mins_neon(
    seq: &[u8],
    l: usize,
    window: usize,
    first: usize,
    lane_len: usize,
    out: &mut [u32],
) {
    use std::arch::aarch64::*;

    let w = window - l + 1;
    let mask = vdupq_n_u32(if l == 16 {
        u32::MAX
    } else {
        (1u32 << (2 * l)) - 1
    });
    let rev_shift = vdupq_n_s32(2 * (l as i32 - 1));
    let three = vdupq_n_u32(0b11);
    let byte_mask = vdupq_n_u32(0xff);
    let max = vdupq_n_u32(u32::MAX);
    let seed = vdupq_n_u32(LMER_HASH_SEED);
    let (mul1, mul2) = (vdupq_n_u32(0x85eb_ca6b), vdupq_n_u32(0xc2b2_ae35));
    let mut ring = [[max; 2]; 64];
    let mut prefix = [max; 2];
    let mut slot = 0;
    let (mut fwd, mut rev) = ([vdupq_n_u32(0); 2], [vdupq_n_u32(0); 2]);
    let starts: [usize; 8] = std::array::from_fn(|lane| first + lane * lane_len);
    let mut words = [0u32; 8];
    let mut lanes = [0u32; 8];
    let steps = lane_len + window - 1;
    let mut t = 0;
    while t < steps {
        // Four bases per lane per load; the indexing is bounds-checked, and
        // `fill_window_mins_simd` keeps every load inside `seq`.
        for (word, &start) in words.iter_mut().zip(&starts) {
            let at = start + t;
            *word = u32::from_le_bytes([seq[at], seq[at + 1], seq[at + 2], seq[at + 3]]);
        }
        // SAFETY: `words` holds eight u32s.
        let mut bytes = unsafe { [vld1q_u32(words.as_ptr()), vld1q_u32(words.as_ptr().add(4))] };
        for _ in 0..4.min(steps - t) {
            for half in 0..2 {
                let byte = vandq_u32(bytes[half], byte_mask);
                bytes[half] = vshrq_n_u32::<8>(bytes[half]);
                let bits = vandq_u32(
                    veorq_u32(vshrq_n_u32::<2>(byte), vshrq_n_u32::<1>(byte)),
                    three,
                );
                fwd[half] = vandq_u32(vorrq_u32(vshlq_n_u32::<2>(fwd[half]), bits), mask);
                rev[half] = vorrq_u32(
                    vshrq_n_u32::<2>(rev[half]),
                    vshlq_u32(veorq_u32(bits, three), rev_shift),
                );
            }
            t += 1;
            if t < l {
                continue;
            }
            for half in 0..2 {
                let mut h = veorq_u32(vminq_u32(fwd[half], rev[half]), seed);
                h = veorq_u32(h, vshrq_n_u32::<16>(h));
                h = vmulq_u32(h, mul1);
                h = veorq_u32(h, vshrq_n_u32::<13>(h));
                h = vmulq_u32(h, mul2);
                h = veorq_u32(h, vshrq_n_u32::<16>(h));
                ring[slot][half] = h;
                prefix[half] = vminq_u32(prefix[half], h);
            }
            slot += 1;
            if slot == w {
                slot = 0;
                for j in (0..w - 1).rev() {
                    let next = ring[j + 1];
                    for (current, next) in ring[j].iter_mut().zip(next) {
                        *current = vminq_u32(*current, next);
                    }
                }
                prefix = [max; 2];
            }
            if t >= window {
                // SAFETY: `lanes` holds eight u32s.
                unsafe {
                    vst1q_u32(lanes.as_mut_ptr(), vminq_u32(prefix[0], ring[slot][0]));
                    vst1q_u32(
                        lanes.as_mut_ptr().add(4),
                        vminq_u32(prefix[1], ring[slot][1]),
                    );
                }
                let at = t - window;
                for (lane, &value) in lanes.iter().enumerate() {
                    out[lane * lane_len + at] = value;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_seq(len: usize, seed: u64) -> Vec<u8> {
        let mut state = seed | 1;
        (0..len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state >> 32) as usize & 3]
            })
            .collect()
    }

    fn reverse_complement(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect()
    }

    /// The specification: every window's minimum, computed from scratch.
    fn naive(seq: &[u8], l: usize, window: usize) -> Vec<u32> {
        let encode = |s: &[u8]| s.iter().fold(0u32, |v, &b| (v << 2) | base_bits(b));
        (0..=seq.len() - window)
            .map(|start| {
                (start..=start + window - l)
                    .map(|at| {
                        let lmer = &seq[at..at + l];
                        lmer_hash(encode(lmer).min(encode(&reverse_complement(lmer))))
                    })
                    .min()
                    .unwrap()
            })
            .collect()
    }

    #[test]
    fn scans_match_the_specification() {
        for (l, window) in [
            (12, 30),
            (1, 4),
            (5, 5),
            (16, 16),
            (16, 62),
            (11, 32),
            (7, 20),
        ] {
            for len in [window, window + 1, window + 40, 300, 1_000, 4_099] {
                let seq = random_seq(len, (l * 1_000 + len) as u64);
                let expected = naive(&seq, l, window);
                let mut dispatched = vec![0; expected.len()];
                fill_window_mins(&seq, l, window, 0, &mut dispatched);
                assert_eq!(dispatched, expected, "l={l} window={window} len={len}");
                let mut scalar = vec![0; expected.len()];
                window_mins_scalar(&seq, l, window, 0, &mut scalar);
                assert_eq!(scalar, expected, "scalar l={l} window={window} len={len}");
                // A window range starting mid-sequence.
                if expected.len() > 10 {
                    let mut tail = vec![0; expected.len() - 7];
                    fill_window_mins(&seq, l, window, 7, &mut tail);
                    assert_eq!(
                        tail,
                        expected[7..],
                        "offset l={l} window={window} len={len}"
                    );
                }
            }
        }
    }

    #[test]
    fn window_minima_are_strand_invariant() {
        let seq = random_seq(10_000, 42);
        let rc = reverse_complement(&seq);
        let windows = seq.len() - 30 + 1;
        let (mut forward, mut reverse) = (vec![0; windows], vec![0; windows]);
        fill_window_mins(&seq, 12, 30, 0, &mut forward);
        fill_window_mins(&rc, 12, 30, 0, &mut reverse);
        reverse.reverse();
        assert_eq!(forward, reverse);
    }

    /// Poly-A must not hash to the smallest possible value, or it would be
    /// the minimum of every window that holds it.
    #[test]
    fn poly_a_is_not_a_fixed_minimum() {
        assert_ne!(lmer_hash(0), 0);
    }

    #[test]
    fn lmer_hash_is_a_bijection_on_small_inputs() {
        let mut seen = std::collections::HashSet::new();
        for value in 0..1 << 16 {
            assert!(seen.insert(lmer_hash(value)));
        }
    }
}

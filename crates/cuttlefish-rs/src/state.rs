//! Compact vertex state and colored-coordinate representations.
//!
//! These types are used in the hottest local-contraction tables. Their sizes
//! and bit allocations are deliberate; avoid adding fields without measuring
//! memory bandwidth and peak RSS at scale.

use crate::Side;
use crate::dna::Base;
use xxhash_rust::xxh3::xxh3_64;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct EdgeFrequency {
    packed: u32,
}

impl EdgeFrequency {
    const MAX: u32 = 0xF;

    #[inline]
    pub fn add_edge(&mut self, side: Side, edge: Base) {
        assert!(matches!(edge, Base::A | Base::C | Base::G | Base::T));
        let off = Self::offset(side, edge);
        let mask = Self::MAX << off;
        let cur = (self.packed & mask) >> off;
        if cur < Self::MAX {
            self.packed = (self.packed & !mask) | ((cur + 1) << off);
        }
    }

    #[inline]
    pub fn edge_count(&self, side: Side, cutoff: u32) -> u32 {
        let side_off = Self::side_offset(side);
        let packed = (self.packed >> side_off) & 0xffff;
        u32::from((packed & 0x000f) >= cutoff)
            + u32::from(((packed >> 4) & 0x000f) >= cutoff)
            + u32::from(((packed >> 8) & 0x000f) >= cutoff)
            + u32::from(((packed >> 12) & 0x000f) >= cutoff)
    }

    #[inline]
    pub fn edge_at(&self, side: Side, cutoff: u32) -> Base {
        let side_off = Self::side_offset(side);
        let packed = (self.packed >> side_off) & 0xffff;
        let edge_a = u32::from((packed & 0x000f) >= cutoff);
        let edge_c = u32::from(((packed >> 4) & 0x000f) >= cutoff);
        let edge_g = u32::from(((packed >> 8) & 0x000f) >= cutoff);
        let edge_t = u32::from(((packed >> 12) & 0x000f) >= cutoff);
        match edge_a + edge_c + edge_g + edge_t {
            0 => Base::E,
            1 => match edge_c + 2 * edge_g + 3 * edge_t {
                0 => Base::A,
                1 => Base::C,
                2 => Base::G,
                3 => Base::T,
                _ => unreachable!(),
            },
            _ => Base::N,
        }
    }

    #[inline]
    pub fn frequency(&self, side: Side, base_bits: u32) -> u32 {
        let off = Self::side_offset(side) + 4 * base_bits;
        (self.packed >> off) & Self::MAX
    }

    #[inline]
    fn side_offset(side: Side) -> u32 {
        match side {
            Side::Front => 0,
            Side::Back => 16,
        }
    }

    #[inline]
    fn offset(side: Side, edge: Base) -> u32 {
        Self::side_offset(side) + 4 * edge.bits() as u32
    }
}

/// What a vertex stores about colors, which is nothing unless the build is
/// colored.
///
/// This is the Rust spelling of C++'s `State_Config<bool Colored_>`
/// specialization, extended to edges. Local contraction is
/// memory-bandwidth-bound, so every byte of a vertex costs: dropping the
/// colour field from uncolored states was worth 14% of the phase (the R1
/// profile in the performance record), and keeping edges as presence bits in
/// the flags at cutoff 1 took the colored slot from 24 to 20 bytes and the
/// uncolored from 16 to 12.
///
/// | slot | colours | edges | state |
/// | --- | --- | --- | ---: |
/// | [`Presence`] | none | presence (cutoff 1) | 4 B |
/// | [`ColoredPresence`] | set hash | presence (cutoff 1) | 12 B |
/// | [`Counted`] | none | 4-bit counts | 8 B |
/// | [`ColoredCounted`] | set hash | 4-bit counts | 16 B |
pub trait ColorSlot: Copy + Default + std::fmt::Debug + PartialEq + Eq {
    /// The empty slot, usable in constants.
    const ZERO: Self;

    /// Whether the slot keeps saturating edge counts, which any cutoff can
    /// read. Otherwise edges are presence bits in the state's flags, which
    /// answer only for cutoff 1.
    const COUNTS_EDGES: bool;

    /// Padding for a K > 31 flat-map slot, whose 16-byte key and state
    /// otherwise pack without any. A slot of an odd size straddles cache
    /// lines, so a colored presence slot is padded from 28 to 32 bytes; the
    /// others (20 and 24) are small enough that their density wins.
    type WidePad: Copy + Default + std::fmt::Debug;

    /// The colour-set hash, or zero when the build carries no colours.
    fn hash(self) -> u64;

    /// Folds another source's hash in. A no-op without a slot to fold into,
    /// which is why the shared contraction code can stay generic.
    fn combine(&mut self, source_hash: u64);

    /// The edge counts of a counting slot. Presence slots keep none, so
    /// asking one is a caller bug.
    #[inline(always)]
    fn edge_frequency(self) -> EdgeFrequency {
        debug_assert!(Self::COUNTS_EDGES, "presence slots keep no edge counts");
        EdgeFrequency::default()
    }

    /// The slot with its edge counts replaced; unchanged for presence slots.
    #[inline(always)]
    fn with_edge_frequency(self, _edges: EdgeFrequency) -> Self {
        self
    }
}

/// The uncolored cutoff-1 slot: zero-sized, so it costs a vertex nothing.
/// Edges live as presence bits in the state's flags.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Presence;

impl ColorSlot for Presence {
    const ZERO: Self = Presence;
    const COUNTS_EDGES: bool = false;
    type WidePad = ();

    #[inline(always)]
    fn hash(self) -> u64 {
        0
    }

    #[inline(always)]
    fn combine(&mut self, _source_hash: u64) {}
}

/// The colored cutoff-1 slot: the running hash of the vertex's colour set.
/// Edges live as presence bits in the state's flags.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(C, packed(4))]
pub struct ColoredPresence {
    hash: u64,
}

impl ColorSlot for ColoredPresence {
    const ZERO: Self = Self { hash: 0 };
    const COUNTS_EDGES: bool = false;
    type WidePad = u32;

    #[inline(always)]
    fn hash(self) -> u64 {
        self.hash
    }

    #[inline(always)]
    fn combine(&mut self, source_hash: u64) {
        self.hash = hash_combine(self.hash, source_hash);
    }
}

/// The uncolored slot for any cutoff: saturating edge counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Counted {
    edges: EdgeFrequency,
}

impl ColorSlot for Counted {
    const ZERO: Self = Self {
        edges: EdgeFrequency { packed: 0 },
    };
    const COUNTS_EDGES: bool = true;
    type WidePad = ();

    #[inline(always)]
    fn hash(self) -> u64 {
        0
    }

    #[inline(always)]
    fn combine(&mut self, _source_hash: u64) {}

    #[inline(always)]
    fn edge_frequency(self) -> EdgeFrequency {
        self.edges
    }

    #[inline(always)]
    fn with_edge_frequency(self, edges: EdgeFrequency) -> Self {
        Self { edges }
    }
}

/// The colored slot for any cutoff: edge counts and the colour-set hash,
/// packed to 12 bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(C, packed(4))]
pub struct ColoredCounted {
    edges: EdgeFrequency,
    hash: u64,
}

impl ColorSlot for ColoredCounted {
    const ZERO: Self = Self {
        edges: EdgeFrequency { packed: 0 },
        hash: 0,
    };
    const COUNTS_EDGES: bool = true;
    type WidePad = ();

    #[inline(always)]
    fn hash(self) -> u64 {
        self.hash
    }

    #[inline(always)]
    fn combine(&mut self, source_hash: u64) {
        self.hash = hash_combine(self.hash, source_hash);
    }

    #[inline(always)]
    fn edge_frequency(self) -> EdgeFrequency {
        self.edges
    }

    #[inline(always)]
    fn with_edge_frequency(self, edges: EdgeFrequency) -> Self {
        Self {
            edges,
            hash: self.hash,
        }
    }
}

/// Largest source ID representable in the packed last-source field.
///
/// Colored builds track the previously seen source per vertex in 21 bits of
/// `flags`; partitioning rejects larger source sets before contraction so this
/// bound is never reached at run time. It lives outside [`VertexState`] so a
/// caller can check it without naming a colour slot.
pub const MAX_SOURCE_ID: u32 = 0x1F_FFFF;

/// A vertex's edges, flags and colours.
///
/// `flags` holds, from bit 0: visited, the two discontinuity marks, eight
/// edge-presence bits (front A..T, back A..T; unused when the slot counts
/// edges), and from bit 11 the last source seen. Packed to 4-byte alignment
/// so a flat-map slot is the key and state with no padding.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(C, packed(4))]
pub struct VertexState<C: ColorSlot = Counted> {
    flags: u32,
    slot: C,
}

impl<C: ColorSlot> VertexState<C> {
    pub const EMPTY: Self = Self {
        flags: 0,
        slot: C::ZERO,
    };

    const VISITED: u32 = 1 << 0;
    const DISC_FRONT: u32 = 1 << 1;
    const DISC_BACK: u32 = 1 << 2;
    const SOURCE_SHIFT: u32 = 11;
    const SOURCE_MASK: u32 = 0x1F_FFFF << Self::SOURCE_SHIFT;

    const EDGE_SHIFT: u32 = 3;

    #[inline(always)]
    fn edge_bits(&self, side: Side) -> u32 {
        let off = Self::EDGE_SHIFT + if side == Side::Front { 0 } else { 4 };
        (self.flags >> off) & 0xf
    }

    #[inline(always)]
    pub fn update_edges(&mut self, front: Base, back: Base) {
        if C::COUNTS_EDGES {
            let slot = self.slot;
            let mut edges = slot.edge_frequency();
            if front != Base::E {
                edges.add_edge(Side::Front, front);
            }
            if back != Base::E {
                edges.add_edge(Side::Back, back);
            }
            self.slot = slot.with_edge_frequency(edges);
            return;
        }
        // An N here would set a bit outside the eight presence bits, in the
        // last-source field. Callers skip N-adjacent edges, as for counts.
        debug_assert!(
            front.bits() < 4 || front == Base::E,
            "front edge {front:?} is not ACGT"
        );
        debug_assert!(
            back.bits() < 4 || back == Base::E,
            "back edge {back:?} is not ACGT"
        );
        if front != Base::E {
            self.flags |= 1 << (Self::EDGE_SHIFT + front.bits() as u32);
        }
        if back != Base::E {
            self.flags |= 1 << (Self::EDGE_SHIFT + 4 + back.bits() as u32);
        }
    }

    #[inline]
    pub fn edge_at(&self, side: Side, cutoff: u32) -> Base {
        if C::COUNTS_EDGES {
            let slot = self.slot;
            return slot.edge_frequency().edge_at(side, cutoff);
        }
        debug_assert!(cutoff <= 1, "presence bits answer only for cutoff 1");
        let bits = self.edge_bits(side);
        match bits.count_ones() {
            0 => Base::E,
            1 => match bits.trailing_zeros() {
                0 => Base::A,
                1 => Base::C,
                2 => Base::G,
                _ => Base::T,
            },
            _ => Base::N,
        }
    }

    #[inline]
    pub fn is_branching_side(&self, side: Side, cutoff: u32) -> bool {
        if C::COUNTS_EDGES {
            let slot = self.slot;
            return slot.edge_frequency().edge_count(side, cutoff) > 1;
        }
        debug_assert!(cutoff <= 1, "presence bits answer only for cutoff 1");
        self.edge_bits(side).count_ones() > 1
    }

    #[inline]
    pub fn is_empty_side(&self, side: Side, cutoff: u32) -> bool {
        if C::COUNTS_EDGES {
            let slot = self.slot;
            return slot.edge_frequency().edge_count(side, cutoff) == 0;
        }
        debug_assert!(cutoff <= 1, "presence bits answer only for cutoff 1");
        self.edge_bits(side) == 0
    }

    #[inline]
    pub fn is_isolated(&self, cutoff: u32) -> bool {
        self.is_empty_side(Side::Front, cutoff) && self.is_empty_side(Side::Back, cutoff)
    }

    #[inline]
    pub fn mark_visited(&mut self) {
        self.flags |= Self::VISITED;
    }

    #[inline]
    pub fn is_visited(&self) -> bool {
        self.flags & Self::VISITED != 0
    }

    #[inline]
    pub fn mark_discontinuous(&mut self, side: Side) {
        self.flags |= match side {
            Side::Front => Self::DISC_FRONT,
            Side::Back => Self::DISC_BACK,
        };
    }

    #[inline]
    pub fn is_discontinuous(&self, side: Side) -> bool {
        self.flags
            & match side {
                Side::Front => Self::DISC_FRONT,
                Side::Back => Self::DISC_BACK,
            }
            != 0
    }

    #[inline(always)]
    pub fn add_source(&mut self, source: u32) {
        self.add_source_hashed(source, source_hash(source));
    }

    #[inline(always)]
    pub fn add_source_hashed(&mut self, source: u32, source_hash: u64) {
        debug_assert!(
            source <= MAX_SOURCE_ID,
            "source IDs are bounded during partitioning"
        );
        let last = (self.flags & Self::SOURCE_MASK) >> Self::SOURCE_SHIFT;
        if source != last {
            let mut color = self.slot;
            color.combine(source_hash);
            self.slot = color;
            self.flags = (self.flags & !Self::SOURCE_MASK) | (source << Self::SOURCE_SHIFT);
        }
    }

    #[inline]
    pub fn color_hash(&self) -> u64 {
        let color = self.slot;
        color.hash()
    }
}

/// The point of the colour-slot parameter: at cutoff 1 an uncolored vertex
/// costs 4 bytes and a colored one 12 (8 and 16 when the slot counts edges),
/// so a flat-map slot with its 8-byte key is 12 or 20 bytes.
const _: () = assert!(std::mem::size_of::<VertexState<Presence>>() == 4);
const _: () = assert!(std::mem::size_of::<VertexState<ColoredPresence>>() == 12);
const _: () = assert!(std::mem::size_of::<VertexState<Counted>>() == 8);
const _: () = assert!(std::mem::size_of::<VertexState<ColoredCounted>>() == 16);

#[inline]
pub fn source_hash(source: u32) -> u64 {
    debug_assert!(source > 0 && source < (1 << 21));
    xxh3_64(&source.to_le_bytes()[..3])
}

#[inline]
pub fn hash_combine(lhs: u64, rhs: u64) -> u64 {
    lhs ^ rhs
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
/// Compact location of a deduplicated source set in the color repository.
///
/// Published coordinates use 8 worker bits and 32 worker-local index bits.
/// The high bit is reserved for concurrent insertion state.
pub struct ColorCoordinate(u64);

impl ColorCoordinate {
    const IN_PROCESS: u64 = 1u64 << 63;
    const INDEX_SHIFT: u32 = 8;

    pub fn discovered(worker: u64, index: u64) -> Self {
        assert!(worker < (1u64 << Self::INDEX_SHIFT));
        assert!(index < (1u64 << 32));
        Self(worker | (index << Self::INDEX_SHIFT))
    }

    pub fn from_u40(value: u64) -> Self {
        assert!(value < (1u64 << 40));
        Self(value)
    }

    #[inline]
    pub fn is_in_process(self) -> bool {
        self.0 & Self::IN_PROCESS != 0
    }

    #[inline]
    pub fn as_u40(self) -> u64 {
        assert!(self.0 < (1u64 << 40));
        self.0
    }

    #[inline]
    pub fn worker(self) -> usize {
        assert!(!self.is_in_process());
        (self.0 & 0xff) as usize
    }

    #[inline]
    pub fn index(self) -> u32 {
        assert!(!self.is_in_process());
        (self.0 >> Self::INDEX_SHIFT) as u32
    }
}

#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// A packed positional color run.
///
/// The low 24 bits hold a unitig vertex offset and the upper 40 bits hold a
/// [`ColorCoordinate`]. This transparent 64-bit layout is written directly to
/// private intermediate streams.
pub struct UnitigColor(u64);

impl UnitigColor {
    /// Packs a run starting at `offset` and referring to `coord`.
    pub fn new(offset: u32, coord: ColorCoordinate) -> Self {
        assert!(offset <= 0xFF_FFFF);
        Self((coord.as_u40() << 24) | offset as u64)
    }

    /// Returns the zero-based vertex offset where this color run begins.
    #[inline]
    pub fn offset(self) -> u32 {
        (self.0 & 0xFF_FFFF) as u32
    }

    /// Returns the raw 40-bit color-repository coordinate.
    #[inline]
    pub fn coordinate(self) -> u64 {
        self.0 >> 24
    }

    /// Returns the complete packed representation.
    #[inline]
    pub fn raw(self) -> u64 {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// At cutoff 1 the presence slots must answer every edge question as the
    /// counting slots do, and the colored ones must hash sources alike.
    #[test]
    fn presence_slots_match_counting_slots_at_cutoff_one() {
        fn check<P: ColorSlot, Q: ColorSlot>(seed: u64) {
            let mut state = seed | 1;
            let mut next = |n: u64| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                state % n
            };
            let bases = [Base::A, Base::C, Base::G, Base::T, Base::E];
            let (mut presence, mut counted) =
                (VertexState::<P>::default(), VertexState::<Q>::default());
            for _ in 0..next(40) {
                let (front, back) = (bases[next(5) as usize], bases[next(5) as usize]);
                presence.update_edges(front, back);
                counted.update_edges(front, back);
                let source = 1 + next(4) as u32;
                presence.add_source(source);
                counted.add_source(source);
                if next(7) == 0 {
                    let side = if next(2) == 0 {
                        Side::Front
                    } else {
                        Side::Back
                    };
                    presence.mark_discontinuous(side);
                    counted.mark_discontinuous(side);
                }
                if next(11) == 0 {
                    presence.mark_visited();
                    counted.mark_visited();
                }
                assert_eq!(presence.is_visited(), counted.is_visited());
                for side in [Side::Front, Side::Back] {
                    assert_eq!(presence.edge_at(side, 1), counted.edge_at(side, 1));
                    assert_eq!(
                        presence.is_branching_side(side, 1),
                        counted.is_branching_side(side, 1)
                    );
                    assert_eq!(
                        presence.is_empty_side(side, 1),
                        counted.is_empty_side(side, 1)
                    );
                    assert_eq!(
                        presence.is_discontinuous(side),
                        counted.is_discontinuous(side)
                    );
                }
                assert_eq!(presence.is_isolated(1), counted.is_isolated(1));
                assert_eq!(presence.color_hash(), counted.color_hash());
            }
        }
        for seed in 0..2000 {
            check::<Presence, Counted>(seed);
            check::<ColoredPresence, ColoredCounted>(seed);
        }
    }

    #[test]
    fn edge_frequency_saturates_and_respects_cutoff() {
        let mut f = EdgeFrequency::default();
        f.add_edge(Side::Back, Base::A);
        assert_eq!(f.edge_at(Side::Back, 1), Base::A);
        assert_eq!(f.edge_at(Side::Back, 2), Base::E);
        f.add_edge(Side::Back, Base::C);
        assert_eq!(f.edge_at(Side::Back, 1), Base::N);
        for _ in 0..20 {
            f.add_edge(Side::Back, Base::A);
        }
        assert_eq!(f.frequency(Side::Back, Base::A.bits() as u32), 15);
    }

    #[test]
    fn color_coordinate_packing_matches_limits() {
        let c = ColorCoordinate::discovered(7, 42);
        assert_eq!(c.as_u40(), 7 | (42 << 8));
        let u = UnitigColor::new(123, c);
        assert_eq!(u.offset(), 123);
        assert_eq!(u.coordinate(), c.as_u40());
    }
}

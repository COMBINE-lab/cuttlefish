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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct VertexState {
    edges: EdgeFrequency,
    flags: u32,
    color_hash: u64,
}

impl VertexState {
    const VISITED: u32 = 1 << 0;
    const DISC_FRONT: u32 = 1 << 1;
    const DISC_BACK: u32 = 1 << 2;
    const SOURCE_SHIFT: u32 = 11;
    const SOURCE_MASK: u32 = 0x1F_FFFF << Self::SOURCE_SHIFT;

    #[inline(always)]
    pub fn update_edges(&mut self, front: Base, back: Base) {
        if front != Base::E {
            self.edges.add_edge(Side::Front, front);
        }
        if back != Base::E {
            self.edges.add_edge(Side::Back, back);
        }
    }

    #[inline]
    pub fn edge_at(&self, side: Side, cutoff: u32) -> Base {
        self.edges.edge_at(side, cutoff)
    }

    #[inline]
    pub fn is_branching_side(&self, side: Side, cutoff: u32) -> bool {
        self.edges.edge_count(side, cutoff) > 1
    }

    #[inline]
    pub fn is_empty_side(&self, side: Side, cutoff: u32) -> bool {
        self.edges.edge_count(side, cutoff) == 0
    }

    #[inline]
    pub fn is_isolated(&self, cutoff: u32) -> bool {
        self.is_empty_side(Side::Front, cutoff) && self.is_empty_side(Side::Back, cutoff)
    }

    #[inline]
    pub fn is_discontinuity(&self) -> bool {
        self.is_discontinuous(Side::Front) || self.is_discontinuous(Side::Back)
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
        assert!(source <= 0x1F_FFFF);
        let last = (self.flags & Self::SOURCE_MASK) >> Self::SOURCE_SHIFT;
        if source != last {
            self.color_hash = hash_combine(self.color_hash, source_hash);
            self.flags = (self.flags & !Self::SOURCE_MASK) | (source << Self::SOURCE_SHIFT);
        }
    }

    #[inline]
    pub fn color_hash(&self) -> u64 {
        self.color_hash
    }
}

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
pub struct ColorCoordinate(u64);

impl ColorCoordinate {
    const IN_PROCESS: u64 = 1u64 << 63;
    const INDEX_SHIFT: u32 = 8;

    pub fn in_process(worker: u64) -> Self {
        assert!(worker < (1u64 << Self::INDEX_SHIFT));
        Self(Self::IN_PROCESS | worker)
    }

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
    pub fn processing_worker(self) -> u64 {
        assert!(self.is_in_process());
        self.0 & !Self::IN_PROCESS
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
pub struct UnitigColor(u64);

impl UnitigColor {
    pub fn new(offset: u32, coord: ColorCoordinate) -> Self {
        assert!(offset <= 0xFF_FFFF);
        Self((coord.as_u40() << 24) | offset as u64)
    }

    #[inline]
    pub fn offset(self) -> u32 {
        (self.0 & 0xFF_FFFF) as u32
    }

    #[inline]
    pub fn coordinate(self) -> u64 {
        self.0 >> 24
    }

    #[inline]
    pub fn raw(self) -> u64 {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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

//! Local subgraph construction and contraction.
//!
//! Each weak-super-k-mer bucket is converted to a dense vertex-state table and
//! contracted independently. The map is reused between buckets to amortize
//! allocation. Emitted unitigs retain only discontinuity endpoints, labels,
//! and optional positional colors needed by the global external pipeline.

use crate::Side;
#[cfg(test)]
use crate::buckets::BucketRecord;
use crate::buckets::{BorrowedBucketPackedRecord, BucketError, BucketManifestEntry, BucketStore};
use crate::color::{ColorError, ConcurrentColorRepository};
use crate::dna::{Base, minimal_rotation};
use crate::hash::{FastBuildHasher, fast_u64_hash, hash_two_u64};
use crate::kmer::{Kmer, KmerError};
use crate::state::{ColorCoordinate, UnitigColor, VertexState, source_hash};
use hashbrown::HashMap;
use hashbrown::HashTable;
use std::collections::HashSet;

/// Open-addressed vertex map for k-mers that fit 62 bits (K <= 31), with a
/// probe loop this crate owns so record processing can software-prefetch.
///
/// The R1 profile put ~39% of uncolored local contraction in the vertex-map
/// probe chain, and the inline-entry rearrangement measured *worse* (+16%
/// local_s): the cost is one DRAM miss per probe, not the number of lines
/// touched, so the only lever left is hiding that miss. `prefetch` computes
/// the home slot for a key about to be probed and pulls its line while the
/// caller still has register work to do; `add_packed_parts` pre-rolls each
/// record and prefetches every vertex before the state-update loop begins.
///
/// `u64::MAX` marks an empty slot, unreachable while keys carry at most 62
/// bits. K = 32 would collide with the marker, but k must be odd, so the
/// next reachable size after 31 is 33 -- served by the two-word map.
#[derive(Debug, Clone)]
pub struct FlatVertexMap {
    slots: Vec<(u64, VertexState)>,
    /// Insertion order, for the dispatch iteration.
    keys: Vec<u64>,
    mask: usize,
}

const FLAT_EMPTY_KEY: u64 = u64::MAX;

impl FlatVertexMap {
    /// Slot count giving <= 4/5 load for `capacity` live keys.
    fn slot_count(capacity: usize) -> usize {
        (capacity.max(8) * 5 / 4 + 1).next_power_of_two()
    }

    fn with_capacity(capacity: usize) -> Self {
        let slots = Self::slot_count(capacity);
        Self {
            slots: vec![(FLAT_EMPTY_KEY, VertexState::default()); slots],
            keys: Vec::with_capacity(capacity),
            mask: slots - 1,
        }
    }

    fn clear_and_reserve(&mut self, capacity: usize) {
        let needed = Self::slot_count(capacity.max(self.keys.len()));
        if self.slots.len() < needed {
            self.slots = vec![(FLAT_EMPTY_KEY, VertexState::default()); needed];
            self.mask = needed - 1;
        } else {
            self.slots.fill((FLAT_EMPTY_KEY, VertexState::default()));
        }
        self.keys.clear();
        if self.keys.capacity() < capacity {
            self.keys.reserve(capacity - self.keys.capacity());
        }
    }

    #[inline(always)]
    fn prefetch(&self, key: u64) {
        #[cfg(target_arch = "x86_64")]
        unsafe {
            use std::arch::x86_64::{_MM_HINT_T0, _mm_prefetch};
            let index = (local_u64_hash(key) as usize) & self.mask;
            // SAFETY: `index` is masked into `slots`, and prefetching any
            // mapped address is side-effect free.
            _mm_prefetch(self.slots.as_ptr().add(index).cast::<i8>(), _MM_HINT_T0);
        }
        #[cfg(not(target_arch = "x86_64"))]
        let _ = key;
    }

    #[inline(always)]
    fn probe(&self, key: u64) -> (usize, bool) {
        let mut index = (local_u64_hash(key) as usize) & self.mask;
        loop {
            let stored = self.slots[index].0;
            if stored == key {
                return (index, true);
            }
            if stored == FLAT_EMPTY_KEY {
                return (index, false);
            }
            index = (index + 1) & self.mask;
        }
    }

    #[inline(always)]
    fn get(&self, key: u64) -> Option<&VertexState> {
        let (index, hit) = self.probe(key);
        hit.then(|| &self.slots[index].1)
    }

    #[inline(always)]
    fn get_mut(&mut self, key: u64) -> Option<&mut VertexState> {
        let (index, hit) = self.probe(key);
        hit.then(|| &mut self.slots[index].1)
    }

    #[inline(always)]
    fn state_or_default(&mut self, key: u64) -> &mut VertexState {
        let (index, hit) = self.probe(key);
        if hit {
            return &mut self.slots[index].1;
        }
        if (self.keys.len() + 1) * 5 > self.slots.len() * 4 {
            self.grow();
            let (index, _) = self.probe(key);
            self.keys.push(key);
            self.slots[index] = (key, VertexState::default());
            return &mut self.slots[index].1;
        }
        self.keys.push(key);
        self.slots[index] = (key, VertexState::default());
        &mut self.slots[index].1
    }

    #[cold]
    fn grow(&mut self) {
        let new_len = self.slots.len() * 2;
        let old = std::mem::replace(
            &mut self.slots,
            vec![(FLAT_EMPTY_KEY, VertexState::default()); new_len],
        );
        self.mask = new_len - 1;
        for (key, state) in old {
            if key != FLAT_EMPTY_KEY {
                let (index, hit) = self.probe(key);
                debug_assert!(!hit);
                self.slots[index] = (key, state);
            }
        }
    }
}

impl PartialEq for FlatVertexMap {
    fn eq(&self, other: &Self) -> bool {
        self.keys.len() == other.keys.len()
            && self.keys.iter().all(|&key| other.get(key) == self.get(key))
    }
}

impl Eq for FlatVertexMap {}

/// The K > 31 twin of [`FlatVertexMap`].
///
/// Same open-addressed layout and the same reason for existing -- one DRAM
/// miss per probe, hidden by prefetching the home slot -- but keyed on the
/// whole two-word k-mer. `u128::MAX` marks an empty slot: k is capped at 63,
/// so a key carries at most 126 bits and can never reach the marker.
///
/// Before this, K > 31 used a `HashMap`, which cost the probe loop its
/// prefetch and cost the contraction loops their dense iteration order (the
/// hash map cannot index its keys, so every subgraph materialised a key
/// vector instead).
#[derive(Debug, Clone)]
pub struct WideFlatVertexMap {
    slots: Vec<(u128, VertexState)>,
    /// Insertion order, for the dispatch iteration.
    keys: Vec<u128>,
    mask: usize,
}

const WIDE_FLAT_EMPTY_KEY: u128 = u128::MAX;

impl WideFlatVertexMap {
    fn slot_count(capacity: usize) -> usize {
        (capacity.max(8) * 5 / 4 + 1).next_power_of_two()
    }

    fn with_capacity(capacity: usize) -> Self {
        let slots = Self::slot_count(capacity);
        Self {
            slots: vec![(WIDE_FLAT_EMPTY_KEY, VertexState::default()); slots],
            keys: Vec::with_capacity(capacity),
            mask: slots - 1,
        }
    }

    fn clear_and_reserve(&mut self, capacity: usize) {
        let needed = Self::slot_count(capacity.max(self.keys.len()));
        if self.slots.len() < needed {
            self.slots = vec![(WIDE_FLAT_EMPTY_KEY, VertexState::default()); needed];
            self.mask = needed - 1;
        } else {
            self.slots
                .fill((WIDE_FLAT_EMPTY_KEY, VertexState::default()));
        }
        self.keys.clear();
        if self.keys.capacity() < capacity {
            self.keys.reserve(capacity - self.keys.capacity());
        }
    }

    #[inline(always)]
    fn hash(key: u128) -> u64 {
        hash_two_u64(key as u64, (key >> 64) as u64)
    }

    #[inline(always)]
    fn prefetch(&self, key: u128) {
        #[cfg(target_arch = "x86_64")]
        unsafe {
            use std::arch::x86_64::{_MM_HINT_T0, _mm_prefetch};
            let index = (Self::hash(key) as usize) & self.mask;
            // SAFETY: `index` is masked into `slots`, and prefetching any
            // mapped address is side-effect free.
            _mm_prefetch(self.slots.as_ptr().add(index).cast::<i8>(), _MM_HINT_T0);
        }
        #[cfg(not(target_arch = "x86_64"))]
        let _ = key;
    }

    #[inline(always)]
    fn probe(&self, key: u128) -> (usize, bool) {
        let mut index = (Self::hash(key) as usize) & self.mask;
        loop {
            let stored = self.slots[index].0;
            if stored == key {
                return (index, true);
            }
            if stored == WIDE_FLAT_EMPTY_KEY {
                return (index, false);
            }
            index = (index + 1) & self.mask;
        }
    }

    #[inline(always)]
    fn get(&self, key: u128) -> Option<&VertexState> {
        let (index, hit) = self.probe(key);
        hit.then(|| &self.slots[index].1)
    }

    #[inline(always)]
    fn get_mut(&mut self, key: u128) -> Option<&mut VertexState> {
        let (index, hit) = self.probe(key);
        hit.then(|| &mut self.slots[index].1)
    }

    #[inline(always)]
    fn state_or_default(&mut self, key: u128) -> &mut VertexState {
        let (index, hit) = self.probe(key);
        if hit {
            return &mut self.slots[index].1;
        }
        if (self.keys.len() + 1) * 5 > self.slots.len() * 4 {
            self.grow();
            let (index, _) = self.probe(key);
            self.keys.push(key);
            self.slots[index] = (key, VertexState::default());
            return &mut self.slots[index].1;
        }
        self.keys.push(key);
        self.slots[index] = (key, VertexState::default());
        &mut self.slots[index].1
    }

    #[cold]
    fn grow(&mut self) {
        let new_len = self.slots.len() * 2;
        let old = std::mem::replace(
            &mut self.slots,
            vec![(WIDE_FLAT_EMPTY_KEY, VertexState::default()); new_len],
        );
        self.mask = new_len - 1;
        for (key, state) in old {
            if key != WIDE_FLAT_EMPTY_KEY {
                let (index, hit) = self.probe(key);
                debug_assert!(!hit);
                self.slots[index] = (key, state);
            }
        }
    }
}

impl PartialEq for WideFlatVertexMap {
    fn eq(&self, other: &Self) -> bool {
        self.keys.len() == other.keys.len()
            && self.keys.iter().all(|&key| other.get(key) == self.get(key))
    }
}

impl Eq for WideFlatVertexMap {}

#[derive(Default)]
struct DenseWantedColorMap {
    indices: HashTable<u32>,
    entries: Vec<(u64, usize)>,
}

impl DenseWantedColorMap {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            indices: HashTable::with_capacity(capacity),
            entries: Vec::with_capacity(capacity),
        }
    }

    #[inline(always)]
    fn get(&self, key: u64) -> Option<usize> {
        let hash = local_u64_hash(key);
        let index = *self
            .indices
            .find(hash, |&index| self.entries[index as usize].0 == key)?;
        Some(self.entries[index as usize].1)
    }

    fn insert(&mut self, key: u64, value: usize) {
        let hash = local_u64_hash(key);
        if let Some(&index) = self
            .indices
            .find(hash, |&index| self.entries[index as usize].0 == key)
        {
            self.entries[index as usize].1 = value;
            return;
        }
        let index = self.entries.len() as u32;
        self.entries.push((key, value));
        let entries = &self.entries;
        self.indices.insert_unique(hash, index, |&stored| {
            local_u64_hash(entries[stored as usize].0)
        });
    }
}

enum WantedColorMap<const K: usize> {
    OneWord(DenseWantedColorMap),
    TwoWord(HashMap<Kmer<K>, usize, FastBuildHasher>),
}

impl<const K: usize> WantedColorMap<K> {
    fn with_capacity(capacity: usize) -> Self {
        if K <= 32 {
            Self::OneWord(DenseWantedColorMap::with_capacity(capacity))
        } else {
            Self::TwoWord(HashMap::with_capacity_and_hasher(
                capacity,
                FastBuildHasher::default(),
            ))
        }
    }

    fn insert(&mut self, key: Kmer<K>, value: usize) {
        match self {
            Self::OneWord(map) => map.insert(key.as_u128() as u64, value),
            Self::TwoWord(map) => {
                map.insert(key, value);
            }
        }
    }

    #[inline(always)]
    fn get(&self, key: Kmer<K>) -> Option<usize> {
        match self {
            Self::OneWord(map) => map.get(key.as_u128() as u64),
            Self::TwoWord(map) => map.get(&key).copied(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LocalVertexMap<const K: usize> {
    Flat(FlatVertexMap),
    WideFlat(WideFlatVertexMap),
    #[doc(hidden)]
    Unused(std::marker::PhantomData<[(); K]>),
}

impl<const K: usize> LocalVertexMap<K> {
    #[inline(always)]
    fn with_capacity(capacity: usize) -> Self {
        if K <= 31 {
            // 62-bit keys leave u64::MAX free as the flat map's empty marker.
            // K = 32 would collide with it, but k must be odd, so the next
            // reachable size after 31 is 33 -- straight to the two-word map.
            Self::Flat(FlatVertexMap::with_capacity(capacity))
        } else {
            Self::WideFlat(WideFlatVertexMap::with_capacity(capacity))
        }
    }

    fn clear_and_reserve(&mut self, capacity: usize) {
        match self {
            Self::Flat(map) => map.clear_and_reserve(capacity),
            Self::WideFlat(map) => map.clear_and_reserve(capacity),
            Self::Unused(_) => unreachable!("placeholder variant is never built"),
        }
    }

    /// Pulls the home cache line of every vertex in a packed record before
    /// the state-update loop probes them. Flat map only: it is the one table
    /// whose slot addresses this crate can compute.
    #[inline]
    fn prefetch_packed_record(&self, words: &[u64], len: usize) {
        let last_vertex_offset = len - K;
        let mut observed = Kmer::<K>::from_bits(packed_prefix_bits::<K>(words));
        let mut reverse = observed.reverse_complement();
        for offset in 0..=last_vertex_offset {
            let canonical = if observed <= reverse {
                observed
            } else {
                reverse
            };
            match self {
                Self::Flat(map) => map.prefetch(canonical.as_u128() as u64),
                Self::WideFlat(map) => map.prefetch(canonical.as_u128()),
                Self::Unused(_) => {}
            }
            if offset < last_vertex_offset {
                let succ_base = packed_base(words, offset + K);
                observed = observed.roll_forward(succ_base);
                reverse = reverse.roll_backward(succ_base.complement());
            }
        }
    }

    pub fn len(&self) -> usize {
        match self {
            Self::Flat(map) => map.keys.len(),
            Self::WideFlat(map) => map.keys.len(),
            Self::Unused(_) => 0,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Slots the map can hold without growing.
    ///
    /// Clearing walks this rather than [`len`](Self::len), so it is what decides
    /// whether carrying the map to a smaller subgraph is worth the cost.
    pub fn capacity(&self) -> usize {
        match self {
            Self::Flat(map) => map.slots.len() * 4 / 5,
            Self::WideFlat(map) => map.slots.len() * 4 / 5,
            Self::Unused(_) => 0,
        }
    }

    #[inline(always)]
    fn get(&self, kmer: &Kmer<K>) -> Option<&VertexState> {
        match self {
            Self::Flat(map) => map.get(kmer.as_u128() as u64),
            Self::WideFlat(map) => map.get(kmer.as_u128()),
            Self::Unused(_) => None,
        }
    }

    #[inline(always)]
    fn get_mut(&mut self, kmer: &Kmer<K>) -> Option<&mut VertexState> {
        match self {
            Self::Flat(map) => map.get_mut(kmer.as_u128() as u64),
            Self::WideFlat(map) => map.get_mut(kmer.as_u128()),
            Self::Unused(_) => None,
        }
    }

    fn keys_vec(&self) -> Vec<Kmer<K>> {
        match self {
            Self::Flat(map) => map
                .keys
                .iter()
                .map(|&key| Kmer::<K>::from_bits(key as u128))
                .collect(),
            Self::WideFlat(map) => map
                .keys
                .iter()
                .map(|&key| Kmer::<K>::from_bits(key))
                .collect(),
            Self::Unused(_) => Vec::new(),
        }
    }

    fn dense_key_state(&self, index: usize) -> Option<(Kmer<K>, VertexState)> {
        match self {
            Self::Flat(map) => {
                let &key = map.keys.get(index)?;
                map.get(key)
                    .map(|&state| (Kmer::<K>::from_bits(key as u128), state))
            }
            Self::WideFlat(map) => {
                let &key = map.keys.get(index)?;
                map.get(key)
                    .map(|&state| (Kmer::<K>::from_bits(key), state))
            }
            Self::Unused(_) => None,
        }
    }

    fn dense_len(&self) -> Option<usize> {
        match self {
            Self::Flat(map) => Some(map.keys.len()),
            Self::WideFlat(map) => Some(map.keys.len()),
            Self::Unused(_) => None,
        }
    }

    #[inline(always)]
    fn state_or_default(&mut self, kmer: Kmer<K>) -> &mut VertexState {
        match self {
            Self::Flat(map) => map.state_or_default(kmer.as_u128() as u64),
            Self::WideFlat(map) => map.state_or_default(kmer.as_u128()),
            Self::Unused(_) => unreachable!("placeholder variant is never built"),
        }
    }
}
type LocalEdgeSet<const K: usize> = HashSet<LocalEdge<K>, FastBuildHasher>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalSubgraph<const K: usize> {
    pub graph_id: usize,
    pub colored: bool,
    pub cutoff: u32,
    pub vertices: LocalVertexMap<K>,
    pub edges: LocalEdgeSet<K>,
    pub stats: LocalSubgraphStats,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LocalEdge<const K: usize> {
    pub from: Kmer<K>,
    pub to: Kmer<K>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct LocalSubgraphStats {
    pub weak_superkmers: u64,
    pub weak_superkmer_bases: u64,
    pub observed_vertices: u64,
    pub unique_vertices: u64,
    pub observed_edges: u64,
    pub unique_edges: u64,
    pub discontinuity_fronts: u64,
    pub discontinuity_backs: u64,
    pub isolated_vertices: u64,
    pub unitigs: u64,
    pub trivial_unitigs: u64,
    pub cyclic_unitigs: u64,
    pub discontinuity_exits: u64,
    pub unitig_bases: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalUnitig<const K: usize> {
    pub label: Vec<u8>,
    pub vertices: Vec<Kmer<K>>,
    color_hashes: Vec<u64>,
    pub left_exit: Option<(Kmer<K>, Side)>,
    pub right_exit: Option<(Kmer<K>, Side)>,
    pub is_cycle: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalUnitigColorRun {
    pub offset: u32,
    pub color_hash: u64,
    pub coordinate: Option<ColorCoordinate>,
    pub sources: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ColoredLocalUnitig<const K: usize> {
    pub unitig: LocalUnitig<K>,
    pub colors: Vec<LocalUnitigColorRun>,
}

/// A maximal local unitig whose label lives in a caller-owned scratch buffer.
///
/// The production emit paths hand one of these per unitig; the label bytes are
/// valid only for the duration of the callback. This is what removes the
/// per-unitig allocation the owned [`LocalUnitig`] requires -- roughly 1.1
/// billion of them on the full Salmonella corpus.
#[derive(Debug, Clone, Copy)]
pub(crate) struct LocalUnitigRef<'a, const K: usize> {
    pub label: &'a [u8],
    pub left_exit: Option<(Kmer<K>, Side)>,
    pub right_exit: Option<(Kmer<K>, Side)>,
    pub is_cycle: bool,
}

/// Walk endpoints of a compact extraction; the label went to the scratch.
struct CompactUnitigEnds<const K: usize> {
    left_exit: Option<(Kmer<K>, Side)>,
    right_exit: Option<(Kmer<K>, Side)>,
    is_cycle: bool,
}

type PendingColorRun = (u32, u64, Option<ColorCoordinate>);

struct PendingColoredContraction<const K: usize> {
    unitigs: Vec<LocalUnitig<K>>,
    runs: Vec<Vec<PendingColorRun>>,
    representative_indices: HashMap<u64, usize, FastBuildHasher>,
    source_sets: Vec<Vec<u32>>,
}

struct PendingColoredData {
    runs: Vec<Vec<PendingColorRun>>,
    representative_indices: HashMap<u64, usize, FastBuildHasher>,
    source_sets: Vec<Vec<u32>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct UnitigWalk<const K: usize> {
    label: Vec<u8>,
    vertices: WalkVertices<K>,
    color_hashes: Vec<u64>,
    is_cycle: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct WalkVertices<const K: usize> {
    low: Vec<u64>,
    high: Vec<u64>,
}

impl<const K: usize> WalkVertices<K> {
    #[inline]
    fn clear(&mut self) {
        self.low.clear();
        self.high.clear();
    }

    #[inline]
    fn push(&mut self, vertex: Kmer<K>) {
        let words = vertex.words();
        self.low.push(words[0]);
        if K > 32 {
            self.high.push(words[1]);
        }
    }

    #[inline]
    fn len(&self) -> usize {
        self.low.len()
    }

    #[inline]
    fn get(&self, index: usize) -> Kmer<K> {
        let high = if K > 32 { self.high[index] } else { 0 };
        Kmer::from_bits(self.low[index] as u128 | ((high as u128) << 64))
    }
}

impl<const K: usize> Default for WalkVertices<K> {
    fn default() -> Self {
        Self {
            low: Vec::new(),
            high: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum WalkTermination<const K: usize> {
    Branched,
    Crossed,
    DeadEnded,
    Exited(DirectedKmer<K>),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct DirectedKmer<const K: usize> {
    observed: Kmer<K>,
    reverse: Kmer<K>,
}

impl<const K: usize> DirectedKmer<K> {
    #[inline]
    fn canonical(self) -> Kmer<K> {
        self.observed.min(self.reverse)
    }

    #[inline]
    fn in_canonical_form(self) -> bool {
        self.observed <= self.reverse
    }

    #[inline]
    fn entrance_side(self) -> Side {
        if self.in_canonical_form() {
            Side::Front
        } else {
            Side::Back
        }
    }

    #[inline]
    fn roll_forward(self, base: Base) -> Self {
        Self {
            observed: self.observed.roll_forward(base),
            reverse: self.reverse.roll_backward(base.complement()),
        }
    }
}

impl<const K: usize> UnitigWalk<K> {
    fn reset(&mut self, v: DirectedKmer<K>, collect_vertices: bool, color_hash: u64) {
        self.label.clear();
        v.observed.append_ascii(&mut self.label);
        self.vertices.clear();
        self.color_hashes.clear();
        if collect_vertices {
            self.vertices.push(v.canonical());
            self.color_hashes.push(color_hash);
        }
        self.is_cycle = false;
    }

    fn extend(
        &mut self,
        v: DirectedKmer<K>,
        base: Base,
        anchor: Kmer<K>,
        collect_vertices: bool,
        color_hash: u64,
    ) -> bool {
        if v.canonical() == anchor {
            self.is_cycle = true;
            return false;
        }

        self.label.push(base.to_ascii());
        if collect_vertices {
            self.vertices.push(v.canonical());
            self.color_hashes.push(color_hash);
        }
        true
    }
}

impl<const K: usize> Default for UnitigWalk<K> {
    fn default() -> Self {
        Self {
            label: Vec::new(),
            vertices: WalkVertices::default(),
            color_hashes: Vec::new(),
            is_cycle: false,
        }
    }
}

fn reverse_complement_label(label: &[u8]) -> Vec<u8> {
    label
        .iter()
        .rev()
        .map(|&b| complement_valid_ascii(b))
        .collect()
}

#[inline(always)]
fn complement_valid_ascii(base: u8) -> u8 {
    const COMPLEMENT: [u8; 8] = *b"NTNGANNC";
    COMPLEMENT[(base & 7) as usize]
}

fn canonical_cycle_label<const K: usize>(label: Vec<u8>) -> Vec<u8> {
    if label.len() <= K {
        return label;
    }

    let cycle_len = label.len() - K + 1;
    let forward = minimal_rotation(&label[..cycle_len]);
    let reverse = minimal_rotation(&reverse_complement_label(&label[..cycle_len]));
    let canonical_body = if reverse < forward {
        reverse.as_slice()
    } else {
        forward.as_slice()
    };
    linearize_cycle_body::<K>(canonical_body)
}

fn linearize_cycle_body<const K: usize>(body: &[u8]) -> Vec<u8> {
    let mut label = Vec::with_capacity(body.len() + K - 1);
    label.extend_from_slice(body);
    for i in 0..K - 1 {
        label.push(body[i % body.len()]);
    }
    label
}

/// Iteration order over a subgraph's vertex table, held apart from `self` so
/// the contraction loops can mutate vertex states while walking it.
///
/// One plan serves the three contraction loops (owned, compact, colored),
/// which previously each duplicated the dense-vs-sorted dispatch.
enum VertexOrder<const K: usize> {
    Dense(usize),
    Keys(Vec<Kmer<K>>),
}

impl<const K: usize> VertexOrder<K> {
    fn plan(subgraph: &LocalSubgraph<K>) -> Self {
        if let Some(count) = subgraph.vertices.dense_len()
            && !sort_local_vertices_diagnostic()
        {
            Self::Dense(count)
        } else {
            let mut keys = subgraph.vertices.keys_vec();
            if sort_local_vertices_diagnostic() {
                keys.sort_unstable();
            }
            Self::Keys(keys)
        }
    }

    fn len(&self) -> usize {
        match self {
            Self::Dense(count) => *count,
            Self::Keys(keys) => keys.len(),
        }
    }

    fn get(&self, subgraph: &LocalSubgraph<K>, index: usize) -> Option<(Kmer<K>, VertexState)> {
        match self {
            Self::Dense(_) => subgraph.vertices.dense_key_state(index),
            Self::Keys(keys) => {
                let v_hat = keys[index];
                subgraph.vertices.get(&v_hat).copied().map(|s| (v_hat, s))
            }
        }
    }
}

/// Whether to contract vertices in sorted order instead of storage order
/// (diagnostic; makes unitig emission order deterministic for debugging).
fn sort_local_vertices_diagnostic() -> bool {
    static SORT: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *SORT.get_or_init(|| std::env::var_os("CF3_RS_SORT_LOCAL_VERTICES").is_some())
}

/// Whether to decode bucket records without building the graph (diagnostic).
fn decode_only_diagnostic() -> bool {
    static DECODE_ONLY: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *DECODE_ONLY.get_or_init(|| std::env::var_os("CF3_RS_DECODE_ONLY").is_some())
}

impl<const K: usize> LocalSubgraph<K> {
    pub fn from_manifest_entries(
        store: &BucketStore,
        entries: &[BucketManifestEntry],
        cutoff: u32,
    ) -> Result<Self, LocalSubgraphError> {
        Self::from_entries_with_capacity(store, entries, cutoff, 0, None)
    }

    pub(crate) fn from_manifest_entries_reusing(
        store: &BucketStore,
        entries: &[BucketManifestEntry],
        cutoff: u32,
        vertices: Option<LocalVertexMap<K>>,
    ) -> Result<Self, LocalSubgraphError> {
        Self::from_entries_with_capacity(store, entries, cutoff, 0, vertices)
    }

    fn from_entries_with_capacity(
        store: &BucketStore,
        entries: &[BucketManifestEntry],
        cutoff: u32,
        vertex_capacity: usize,
        reusable_vertices: Option<LocalVertexMap<K>>,
    ) -> Result<Self, LocalSubgraphError> {
        if cutoff == 0 {
            return Err(LocalSubgraphError::InvalidCutoff);
        }
        let Some(first_entry) = entries.first() else {
            return Err(LocalSubgraphError::EmptyBucketGroup);
        };
        let mut reader = store.reader(first_entry)?;
        if reader.header().k as usize != K {
            return Err(LocalSubgraphError::KMismatch {
                expected: K,
                got: reader.header().k as usize,
            });
        }

        let graph_id = reader.header().graph_id;
        let colored = reader.header().colored;
        let mut vertices =
            reusable_vertices.unwrap_or_else(|| LocalVertexMap::with_capacity(vertex_capacity));
        vertices.clear_and_reserve(vertex_capacity);
        let mut subgraph = Self {
            graph_id,
            colored,
            cutoff,
            vertices,
            edges: HashSet::with_hasher(FastBuildHasher::default()),
            stats: LocalSubgraphStats::default(),
        };

        // Diagnostic: `CF3_RS_DECODE_ONLY` reads and decodes bucket records
        // without inserting vertices, separating decode cost from table cost.
        let decode_only = decode_only_diagnostic();
        reader.try_for_each_borrowed_packed_record(|record| {
            if decode_only {
                std::hint::black_box(record.words.first());
                return Ok(());
            }
            subgraph.add_borrowed_packed_record(record)
        })?;
        for entry in &entries[1..] {
            let mut reader = store.reader(entry)?;
            if reader.header().k as usize != K {
                return Err(LocalSubgraphError::KMismatch {
                    expected: K,
                    got: reader.header().k as usize,
                });
            }
            if reader.header().graph_id != subgraph.graph_id {
                return Err(LocalSubgraphError::GraphMismatch {
                    expected: subgraph.graph_id,
                    got: reader.header().graph_id,
                });
            }
            if reader.header().colored != subgraph.colored {
                return Err(LocalSubgraphError::MalformedRecord);
            }
            reader.try_for_each_borrowed_packed_record(|record| {
                if decode_only {
                    std::hint::black_box(record.words.first());
                    return Ok(());
                }
                subgraph.add_borrowed_packed_record(record)
            })?;
        }

        subgraph.stats.unique_vertices = subgraph.vertices.len() as u64;
        subgraph.stats.unique_edges = subgraph.edges.len() as u64;

        Ok(subgraph)
    }

    pub(crate) fn into_vertex_map(self) -> LocalVertexMap<K> {
        self.vertices
    }

    pub fn vertex_state(&self, kmer: Kmer<K>) -> Option<&VertexState> {
        self.vertices.get(&kmer)
    }

    pub fn contract(&mut self) -> Result<Vec<LocalUnitig<K>>, LocalSubgraphError> {
        self.contract_internal(true)
    }

    pub(crate) fn contract_compact_with<F>(&mut self, mut emit: F) -> Result<(), LocalSubgraphError>
    where
        F: FnMut(LocalUnitigRef<'_, K>),
    {
        let mut back_walk = UnitigWalk::default();
        let mut front_walk = UnitigWalk::default();
        let mut label = Vec::new();
        let order = VertexOrder::plan(self);
        for index in 0..order.len() {
            let Some((v_hat, state)) = order.get(self, index) else {
                continue;
            };
            self.contract_vertex_compact(
                v_hat,
                state,
                &mut emit,
                &mut back_walk,
                &mut front_walk,
                &mut label,
            )?;
        }
        Ok(())
    }

    pub fn contract_colored(
        &mut self,
        store: &BucketStore,
        entries: &[BucketManifestEntry],
    ) -> Result<Vec<ColoredLocalUnitig<K>>, LocalSubgraphError> {
        self.contract_colored_with_known(store, entries, |_| None)
    }

    pub fn contract_colored_with_known<F>(
        &mut self,
        store: &BucketStore,
        entries: &[BucketManifestEntry],
        color_is_known: F,
    ) -> Result<Vec<ColoredLocalUnitig<K>>, LocalSubgraphError>
    where
        F: Fn(u64) -> Option<ColorCoordinate>,
    {
        let pending = self.contract_colored_impl(store, entries, color_is_known)?;
        pending
            .unitigs
            .into_iter()
            .zip(pending.runs)
            .map(|(unitig, runs)| {
                let colors = runs
                    .into_iter()
                    .map(|(offset, color_hash, coordinate)| {
                        let sources = match coordinate {
                            Some(_) => Vec::new(),
                            None => {
                                let &index = pending
                                    .representative_indices
                                    .get(&color_hash)
                                    .ok_or(LocalSubgraphError::MissingVertex)?;
                                pending.source_sets[index].clone()
                            }
                        };
                        Ok(LocalUnitigColorRun {
                            offset,
                            color_hash,
                            coordinate,
                            sources,
                        })
                    })
                    .collect::<Result<Vec<_>, LocalSubgraphError>>()?;
                Ok(ColoredLocalUnitig { unitig, colors })
            })
            .collect()
    }

    pub(crate) fn contract_colored_resolved_with<F>(
        &mut self,
        store: &BucketStore,
        entries: &[BucketManifestEntry],
        repository: &ConcurrentColorRepository,
        worker: usize,
        emit: F,
    ) -> Result<Vec<Vec<UnitigColor>>, LocalSubgraphError>
    where
        F: for<'u> FnMut(LocalUnitigRef<'u, K>),
    {
        let pending = self.contract_colored_impl_with(
            store,
            entries,
            |color_hash| repository.get(color_hash),
            emit,
        )?;
        let mut color_runs = Vec::with_capacity(pending.runs.len());
        for runs in pending.runs {
            let mut packed = Vec::with_capacity(runs.len());
            for (offset, color_hash, coordinate) in runs {
                let coordinate = match coordinate {
                    Some(coordinate) => coordinate,
                    None => {
                        let &index = pending
                            .representative_indices
                            .get(&color_hash)
                            .ok_or(LocalSubgraphError::MissingVertex)?;
                        repository
                            .resolve_or_insert(color_hash, &pending.source_sets[index], worker)
                            .map_err(LocalSubgraphError::Color)?
                    }
                };
                packed.push(UnitigColor::new(offset, coordinate));
            }
            color_runs.push(packed);
        }
        Ok(color_runs)
    }

    fn contract_colored_impl<F>(
        &mut self,
        store: &BucketStore,
        entries: &[BucketManifestEntry],
        color_is_known: F,
    ) -> Result<PendingColoredContraction<K>, LocalSubgraphError>
    where
        F: Fn(u64) -> Option<ColorCoordinate>,
    {
        let mut unitigs = Vec::new();
        let pending =
            self.contract_colored_impl_with(store, entries, color_is_known, |unitig| {
                unitigs.push(LocalUnitig {
                    label: unitig.label.to_vec(),
                    vertices: Vec::new(),
                    color_hashes: Vec::new(),
                    left_exit: unitig.left_exit,
                    right_exit: unitig.right_exit,
                    is_cycle: unitig.is_cycle,
                })
            })?;
        Ok(PendingColoredContraction {
            unitigs,
            runs: pending.runs,
            representative_indices: pending.representative_indices,
            source_sets: pending.source_sets,
        })
    }

    fn contract_colored_impl_with<F, G>(
        &mut self,
        store: &BucketStore,
        entries: &[BucketManifestEntry],
        color_is_known: F,
        mut emit: G,
    ) -> Result<PendingColoredData, LocalSubgraphError>
    where
        F: Fn(u64) -> Option<ColorCoordinate>,
        G: for<'u> FnMut(LocalUnitigRef<'u, K>),
    {
        if !self.colored {
            return Err(LocalSubgraphError::MalformedRecord);
        }
        let mut representative_indices = HashMap::<u64, usize, FastBuildHasher>::default();
        let mut coordinate_cache =
            HashMap::<u64, Option<ColorCoordinate>, FastBuildHasher>::default();
        let mut representatives = Vec::<Kmer<K>>::new();
        let mut unitig_hash_runs = Vec::new();
        let mut back_walk = UnitigWalk::default();
        let mut front_walk = UnitigWalk::default();
        let mut label = Vec::new();
        let order = VertexOrder::plan(self);
        for index in 0..order.len() {
            let Some((v_hat, state)) = order.get(self, index) else {
                continue;
            };
            self.contract_colored_vertex(
                v_hat,
                state,
                &color_is_known,
                &mut emit,
                &mut unitig_hash_runs,
                &mut representative_indices,
                &mut coordinate_cache,
                &mut representatives,
                &mut back_walk,
                &mut front_walk,
                &mut label,
            )?;
        }

        let mut wanted = WantedColorMap::<K>::with_capacity(representatives.len());
        for (index, &vertex) in representatives.iter().enumerate() {
            wanted.insert(vertex, index);
        }
        let mut source_sets = vec![Vec::<u32>::new(); representatives.len()];
        for entry in entries {
            let mut reader = store.reader(entry)?;
            reader.try_for_each_borrowed_packed_record(|record| {
                collect_wanted_color_relations::<K>(record, &wanted, &mut source_sets)
            })?;
        }
        normalize_source_sets(&mut source_sets);

        Ok(PendingColoredData {
            runs: unitig_hash_runs,
            representative_indices,
            source_sets,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn contract_colored_vertex<F, G>(
        &mut self,
        v_hat: Kmer<K>,
        state: VertexState,
        color_is_known: &F,
        emit: &mut G,
        unitig_hash_runs: &mut Vec<Vec<PendingColorRun>>,
        representative_indices: &mut HashMap<u64, usize, FastBuildHasher>,
        coordinate_cache: &mut HashMap<u64, Option<ColorCoordinate>, FastBuildHasher>,
        representatives: &mut Vec<Kmer<K>>,
        back_walk: &mut UnitigWalk<K>,
        front_walk: &mut UnitigWalk<K>,
        label: &mut Vec<u8>,
    ) -> Result<(), LocalSubgraphError>
    where
        F: Fn(u64) -> Option<ColorCoordinate>,
        G: for<'u> FnMut(LocalUnitigRef<'u, K>),
    {
        if state.is_visited() {
            return Ok(());
        }
        if state.is_isolated(self.cutoff) {
            if let Some(state) = self.vertices.get_mut(&v_hat) {
                state.mark_visited();
            }
            self.stats.isolated_vertices += 1;
            return Ok(());
        }

        let ends =
            self.extract_maximal_unitig_compact(v_hat, back_walk, front_walk, true, label)?;
        self.stats.unitigs += 1;
        self.stats.unitig_bases += label.len() as u64;
        if ends.is_cycle {
            self.stats.cyclic_unitigs += 1;
        }
        if ends.left_exit.is_none() && ends.right_exit.is_none() {
            self.stats.trivial_unitigs += 1;
        } else {
            self.stats.discontinuity_exits +=
                u64::from(ends.left_exit.is_some()) + u64::from(ends.right_exit.is_some());
        }

        let mut runs = Vec::<(u32, u64, Option<ColorCoordinate>)>::new();
        let mut record_hash = |offset: usize, vertex: Kmer<K>, color_hash: u64| {
            if runs
                .last()
                .is_none_or(|&(_, previous, _)| previous != color_hash)
            {
                let coordinate = *coordinate_cache
                    .entry(color_hash)
                    .or_insert_with(|| color_is_known(color_hash));
                runs.push((offset as u32, color_hash, coordinate));
                if coordinate.is_none() {
                    representative_indices.entry(color_hash).or_insert_with(|| {
                        let index = representatives.len();
                        representatives.push(vertex);
                        index
                    });
                }
            }
        };
        if ends.is_cycle {
            let mut vertex = Kmer::<K>::from_ascii(&label[..K])?;
            for offset in 0..=label.len() - K {
                let canonical = vertex.canonical();
                let color_hash = self
                    .vertices
                    .get(&canonical)
                    .ok_or(LocalSubgraphError::MissingVertex)?
                    .color_hash();
                record_hash(offset, canonical, color_hash);
                if offset < label.len() - K {
                    vertex = vertex.roll_forward(Base::from_ascii(label[offset + K]));
                }
            }
        } else {
            debug_assert_eq!(front_walk.vertices.len(), front_walk.color_hashes.len());
            debug_assert_eq!(back_walk.vertices.len(), back_walk.color_hashes.len());
            let mut offset = 0;
            for index in (0..front_walk.vertices.len()).rev() {
                record_hash(
                    offset,
                    front_walk.vertices.get(index),
                    front_walk.color_hashes[index],
                );
                offset += 1;
            }
            for index in 1..back_walk.vertices.len() {
                record_hash(
                    offset,
                    back_walk.vertices.get(index),
                    back_walk.color_hashes[index],
                );
                offset += 1;
            }
        }
        emit(LocalUnitigRef {
            label,
            left_exit: ends.left_exit,
            right_exit: ends.right_exit,
            is_cycle: ends.is_cycle,
        });
        unitig_hash_runs.push(runs);
        Ok(())
    }

    fn contract_internal(
        &mut self,
        collect_vertices: bool,
    ) -> Result<Vec<LocalUnitig<K>>, LocalSubgraphError> {
        let mut unitigs = Vec::new();
        let mut back_walk = UnitigWalk::default();
        let mut front_walk = UnitigWalk::default();
        let order = VertexOrder::plan(self);
        for index in 0..order.len() {
            let Some((v_hat, state)) = order.get(self, index) else {
                continue;
            };
            self.contract_vertex(
                v_hat,
                state,
                collect_vertices,
                &mut |unitig| unitigs.push(unitig),
                &mut back_walk,
                &mut front_walk,
            )?;
        }
        Ok(unitigs)
    }

    fn contract_vertex_compact<F>(
        &mut self,
        v_hat: Kmer<K>,
        state: VertexState,
        emit: &mut F,
        back_walk: &mut UnitigWalk<K>,
        front_walk: &mut UnitigWalk<K>,
        label: &mut Vec<u8>,
    ) -> Result<(), LocalSubgraphError>
    where
        F: for<'u> FnMut(LocalUnitigRef<'u, K>),
    {
        if state.is_visited() {
            return Ok(());
        }
        if state.is_isolated(self.cutoff) {
            if let Some(st) = self.vertices.get_mut(&v_hat) {
                st.mark_visited();
            }
            self.stats.isolated_vertices += 1;
            return Ok(());
        }

        let ends =
            self.extract_maximal_unitig_compact(v_hat, back_walk, front_walk, false, label)?;
        self.stats.unitigs += 1;
        self.stats.unitig_bases += label.len() as u64;
        if ends.is_cycle {
            self.stats.cyclic_unitigs += 1;
        }
        if ends.left_exit.is_none() && ends.right_exit.is_none() {
            self.stats.trivial_unitigs += 1;
        } else {
            self.stats.discontinuity_exits +=
                u64::from(ends.left_exit.is_some()) + u64::from(ends.right_exit.is_some());
        }
        emit(LocalUnitigRef {
            label,
            left_exit: ends.left_exit,
            right_exit: ends.right_exit,
            is_cycle: ends.is_cycle,
        });
        Ok(())
    }

    fn contract_vertex<F>(
        &mut self,
        v_hat: Kmer<K>,
        state: VertexState,
        collect_vertices: bool,
        emit: &mut F,
        back_walk: &mut UnitigWalk<K>,
        front_walk: &mut UnitigWalk<K>,
    ) -> Result<(), LocalSubgraphError>
    where
        F: FnMut(LocalUnitig<K>),
    {
        if state.is_visited() {
            return Ok(());
        }
        if state.is_isolated(self.cutoff) {
            if let Some(st) = self.vertices.get_mut(&v_hat) {
                st.mark_visited();
            }
            self.stats.isolated_vertices += 1;
            return Ok(());
        }

        let unitig = self.extract_maximal_unitig(v_hat, collect_vertices, back_walk, front_walk)?;
        self.stats.unitigs += 1;
        self.stats.unitig_bases += unitig.label.len() as u64;
        if unitig.is_cycle {
            self.stats.cyclic_unitigs += 1;
        }
        if unitig.left_exit.is_none() && unitig.right_exit.is_none() {
            self.stats.trivial_unitigs += 1;
        } else {
            self.stats.discontinuity_exits +=
                u64::from(unitig.left_exit.is_some()) + u64::from(unitig.right_exit.is_some());
        }
        emit(unitig);
        Ok(())
    }

    fn extract_maximal_unitig(
        &mut self,
        v_hat: Kmer<K>,
        collect_vertices: bool,
        back_walk: &mut UnitigWalk<K>,
        front_walk: &mut UnitigWalk<K>,
    ) -> Result<LocalUnitig<K>, LocalSubgraphError> {
        let back_term = self.walk_unitig(v_hat, Side::Back, collect_vertices, back_walk)?;

        if back_walk.is_cycle {
            return Ok(LocalUnitig {
                label: canonical_cycle_label::<K>(back_walk.label.clone()),
                vertices: (0..back_walk.vertices.len())
                    .map(|index| back_walk.vertices.get(index))
                    .collect(),
                color_hashes: back_walk.color_hashes.clone(),
                left_exit: None,
                right_exit: None,
                is_cycle: true,
            });
        }

        let front_term = self.walk_unitig(v_hat, Side::Front, collect_vertices, front_walk)?;
        let vertices = if collect_vertices {
            (0..front_walk.vertices.len())
                .rev()
                .map(|index| front_walk.vertices.get(index))
                .chain((1..back_walk.vertices.len()).map(|index| back_walk.vertices.get(index)))
                .collect()
        } else {
            Vec::new()
        };
        let color_hashes = if collect_vertices {
            front_walk
                .color_hashes
                .iter()
                .rev()
                .chain(back_walk.color_hashes.iter().skip(1))
                .copied()
                .collect()
        } else {
            Vec::new()
        };

        let mut label = Vec::with_capacity(front_walk.label.len() + back_walk.label.len() - K);
        label.extend(
            front_walk
                .label
                .iter()
                .rev()
                .map(|&base| complement_valid_ascii(base)),
        );
        label.extend_from_slice(&back_walk.label[K..]);
        let left_exit = match front_term {
            WalkTermination::Exited(v) => Some((v.canonical(), v.entrance_side())),
            _ => None,
        };
        let right_exit = match back_term {
            WalkTermination::Exited(v) => Some((v.canonical(), v.entrance_side())),
            _ => None,
        };

        Ok(LocalUnitig {
            label,
            vertices,
            color_hashes,
            left_exit,
            right_exit,
            is_cycle: false,
        })
    }

    /// As [`Self::extract_maximal_unitig`], but the label is written into
    /// `label_out` instead of a fresh allocation, and only the walk endpoints
    /// are returned. `collect_walk_vertices` keeps per-vertex color hashes in
    /// the walk buffers for the colored run recording; the uncolored path
    /// skips that work.
    fn extract_maximal_unitig_compact(
        &mut self,
        v_hat: Kmer<K>,
        back_walk: &mut UnitigWalk<K>,
        front_walk: &mut UnitigWalk<K>,
        collect_walk_vertices: bool,
        label_out: &mut Vec<u8>,
    ) -> Result<CompactUnitigEnds<K>, LocalSubgraphError> {
        label_out.clear();
        let back_term = self.walk_unitig(v_hat, Side::Back, collect_walk_vertices, back_walk)?;

        if back_walk.is_cycle {
            // Cycles are rare enough that the rotation's own allocations stay.
            label_out.extend_from_slice(&canonical_cycle_label::<K>(back_walk.label.clone()));
            return Ok(CompactUnitigEnds {
                left_exit: None,
                right_exit: None,
                is_cycle: true,
            });
        }

        let front_term = self.walk_unitig(v_hat, Side::Front, collect_walk_vertices, front_walk)?;
        label_out.reserve(front_walk.label.len() + back_walk.label.len() - K);
        label_out.extend(
            front_walk
                .label
                .iter()
                .rev()
                .map(|&base| complement_valid_ascii(base)),
        );
        label_out.extend_from_slice(&back_walk.label[K..]);
        let left_exit = match front_term {
            WalkTermination::Exited(v) => Some((v.canonical(), v.entrance_side())),
            _ => None,
        };
        let right_exit = match back_term {
            WalkTermination::Exited(v) => Some((v.canonical(), v.entrance_side())),
            _ => None,
        };

        Ok(CompactUnitigEnds {
            left_exit,
            right_exit,
            is_cycle: false,
        })
    }

    fn walk_unitig(
        &mut self,
        v_hat: Kmer<K>,
        start_side: Side,
        collect_vertices: bool,
        walk: &mut UnitigWalk<K>,
    ) -> Result<WalkTermination<K>, LocalSubgraphError> {
        let icc_return_side = start_side.inverse();
        let v_hat_reverse = v_hat.reverse_complement();
        let mut v = if start_side == Side::Back {
            DirectedKmer {
                observed: v_hat,
                reverse: v_hat_reverse,
            }
        } else {
            DirectedKmer {
                observed: v_hat_reverse,
                reverse: v_hat,
            }
        };
        let mut side = start_side;
        let mut state = {
            let state = self
                .vertices
                .get_mut(&v.canonical())
                .ok_or(LocalSubgraphError::MissingVertex)?;
            let copied = *state;
            state.mark_visited();
            copied
        };
        walk.reset(v, collect_vertices, state.color_hash());

        loop {
            let mut edge = state.edge_at(side, self.cutoff);
            if edge == Base::N {
                return Ok(WalkTermination::Branched);
            }
            if edge == Base::E {
                if !state.is_discontinuous(side) {
                    return Ok(WalkTermination::DeadEnded);
                }
                return Ok(WalkTermination::Exited(v));
            }

            if side == Side::Front {
                edge = edge.complement();
            }
            v = v.roll_forward(edge);

            let next_state = self
                .vertices
                .get_mut(&v.canonical())
                .ok_or(LocalSubgraphError::MissingVertex)?;
            let next_state_copy = *next_state;
            side = v.entrance_side();
            if next_state_copy.is_branching_side(side, self.cutoff) {
                return Ok(WalkTermination::Crossed);
            }
            if next_state_copy.is_visited() {
                if v.canonical() == v_hat && side == icc_return_side {
                    walk.is_cycle = true;
                }
                return Ok(WalkTermination::Crossed);
            }

            next_state.mark_visited();
            if !walk.extend(
                v,
                edge,
                v_hat,
                collect_vertices,
                next_state_copy.color_hash(),
            ) {
                return Ok(WalkTermination::Crossed);
            }
            state = next_state_copy;
            side = side.inverse();
        }
    }

    #[cfg(test)]
    fn add_record(&mut self, record: &BucketRecord) -> Result<(), LocalSubgraphError> {
        if record.graph_id != self.graph_id {
            return Err(LocalSubgraphError::GraphMismatch {
                expected: self.graph_id,
                got: record.graph_id,
            });
        }
        if record.len != record.label.len() || record.label.len() < K {
            return Err(LocalSubgraphError::MalformedRecord);
        }

        self.stats.weak_superkmers += 1;
        self.stats.weak_superkmer_bases += record.label.len() as u64;

        let last_vertex_offset = record.label.len() - K;
        let mut observed = Kmer::<K>::from_ascii(&record.label[..K])?;
        let mut reverse = observed.reverse_complement();
        let mut prev = None;
        for offset in 0..=last_vertex_offset {
            let in_canonical_form = observed <= reverse;
            let canonical = if in_canonical_form { observed } else { reverse };
            let pred_base = if offset == 0 {
                Base::E
            } else {
                Base::from_ascii(record.label[offset - 1])
            };
            let succ_base = if offset == last_vertex_offset {
                Base::E
            } else {
                Base::from_ascii(record.label[offset + K])
            };
            let mut front = if in_canonical_form {
                pred_base
            } else {
                succ_base.complement()
            };
            let mut back = if in_canonical_form {
                succ_base
            } else {
                pred_base.complement()
            };

            if offset > 0 && Some(canonical) == prev {
                if in_canonical_form {
                    front = Base::E;
                } else {
                    back = Base::E;
                }
            }

            let mut discontinuity_fronts = 0;
            let mut discontinuity_backs = 0;
            {
                let state = self.vertex_state_or_default(canonical);
                state.update_edges(front, back);
                if let Some(source_id) = record.source_id {
                    state.add_source(source_id);
                }
                if offset == 0 && record.left_discontinuous {
                    let side = if in_canonical_form {
                        Side::Front
                    } else {
                        Side::Back
                    };
                    state.mark_discontinuous(side);
                    match side {
                        Side::Front => discontinuity_fronts += 1,
                        Side::Back => discontinuity_backs += 1,
                    }
                }
                if offset == last_vertex_offset && record.right_discontinuous {
                    let side = if in_canonical_form {
                        Side::Back
                    } else {
                        Side::Front
                    };
                    state.mark_discontinuous(side);
                    match side {
                        Side::Front => discontinuity_fronts += 1,
                        Side::Back => discontinuity_backs += 1,
                    }
                }
            }
            self.stats.discontinuity_fronts += discontinuity_fronts;
            self.stats.discontinuity_backs += discontinuity_backs;
            self.stats.observed_vertices += 1;

            if let Some(from) = prev {
                #[cfg(debug_assertions)]
                {
                    self.edges.insert(LocalEdge {
                        from,
                        to: canonical,
                    });
                }
                #[cfg(not(debug_assertions))]
                let _ = from;
                self.stats.observed_edges += 1;
            }
            prev = Some(canonical);
            if offset < last_vertex_offset {
                observed = observed.roll_forward(succ_base);
                reverse = reverse.roll_backward(succ_base.complement());
            }
        }

        Ok(())
    }

    fn add_borrowed_packed_record(
        &mut self,
        record: BorrowedBucketPackedRecord<'_>,
    ) -> Result<(), LocalSubgraphError> {
        self.add_packed_parts(
            record.graph_id,
            record.len,
            record.source_id,
            record.left_discontinuous,
            record.right_discontinuous,
            record.words,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn add_packed_parts(
        &mut self,
        graph_id: usize,
        len: usize,
        source_id: Option<u32>,
        left_discontinuous: bool,
        right_discontinuous: bool,
        words: &[u64],
    ) -> Result<(), LocalSubgraphError> {
        if graph_id != self.graph_id {
            return Err(LocalSubgraphError::GraphMismatch {
                expected: self.graph_id,
                got: graph_id,
            });
        }
        if len < K || len > words.len() * 32 {
            return Err(LocalSubgraphError::MalformedRecord);
        }

        self.stats.weak_superkmers += 1;
        self.stats.weak_superkmer_bases += len as u64;

        self.vertices.prefetch_packed_record(words, len);
        let last_vertex_offset = len - K;
        let observed_bits = packed_prefix_bits::<K>(words);
        let mut observed = Kmer::<K>::from_bits(observed_bits);
        let mut reverse = observed.reverse_complement();
        let mut prev = None;
        let source = source_id.map(|source| (source, source_hash(source)));
        for offset in 0..=last_vertex_offset {
            let in_canonical_form = observed <= reverse;
            let canonical = if in_canonical_form { observed } else { reverse };
            let pred_base = if offset == 0 {
                Base::E
            } else {
                packed_base(words, offset - 1)
            };
            let succ_base = if offset == last_vertex_offset {
                Base::E
            } else {
                packed_base(words, offset + K)
            };
            let mut front = if in_canonical_form {
                pred_base
            } else {
                succ_base.complement()
            };
            let mut back = if in_canonical_form {
                succ_base
            } else {
                pred_base.complement()
            };

            if offset > 0 && Some(canonical) == prev {
                if in_canonical_form {
                    front = Base::E;
                } else {
                    back = Base::E;
                }
            }

            let mut discontinuity_fronts = 0;
            let mut discontinuity_backs = 0;
            {
                let state = self.vertex_state_or_default(canonical);
                state.update_edges(front, back);
                if let Some((source_id, hash)) = source {
                    state.add_source_hashed(source_id, hash);
                }
                if offset == 0 && left_discontinuous {
                    let side = if in_canonical_form {
                        Side::Front
                    } else {
                        Side::Back
                    };
                    state.mark_discontinuous(side);
                    match side {
                        Side::Front => discontinuity_fronts += 1,
                        Side::Back => discontinuity_backs += 1,
                    }
                }
                if offset == last_vertex_offset && right_discontinuous {
                    let side = if in_canonical_form {
                        Side::Back
                    } else {
                        Side::Front
                    };
                    state.mark_discontinuous(side);
                    match side {
                        Side::Front => discontinuity_fronts += 1,
                        Side::Back => discontinuity_backs += 1,
                    }
                }
            }
            self.stats.discontinuity_fronts += discontinuity_fronts;
            self.stats.discontinuity_backs += discontinuity_backs;
            self.stats.observed_vertices += 1;

            if let Some(from) = prev {
                #[cfg(debug_assertions)]
                {
                    self.edges.insert(LocalEdge {
                        from,
                        to: canonical,
                    });
                }
                #[cfg(not(debug_assertions))]
                let _ = from;
                self.stats.observed_edges += 1;
            }
            prev = Some(canonical);
            if offset < last_vertex_offset {
                observed = observed.roll_forward(succ_base);
                reverse = reverse.roll_backward(succ_base.complement());
            }
        }

        Ok(())
    }

    #[inline]
    fn vertex_state_or_default(&mut self, kmer: Kmer<K>) -> &mut VertexState {
        self.vertices.state_or_default(kmer)
    }
}

fn normalize_source_sets(source_sets: &mut [Vec<u32>]) {
    let max_source = source_sets
        .iter()
        .flat_map(|sources| sources.iter().copied())
        .max()
        .unwrap_or(0) as usize;
    let mut words = vec![0u64; max_source / 64 + 1];
    let mut touched = Vec::<usize>::new();
    for sources in source_sets {
        if sources.len() < 2 {
            continue;
        }
        touched.clear();
        for source in sources.drain(..) {
            let word_index = source as usize / 64;
            if words[word_index] == 0 {
                touched.push(word_index);
            }
            words[word_index] |= 1u64 << (source & 63);
        }
        touched.sort_unstable();
        for &word_index in &touched {
            let mut word = words[word_index];
            while word != 0 {
                let bit = word.trailing_zeros();
                sources.push((word_index as u32) * 64 + bit);
                word &= word - 1;
            }
            words[word_index] = 0;
        }
    }
}

fn collect_wanted_color_relations<const K: usize>(
    record: BorrowedBucketPackedRecord<'_>,
    wanted: &WantedColorMap<K>,
    source_sets: &mut [Vec<u32>],
) -> Result<(), LocalSubgraphError> {
    let source = record
        .source_id
        .ok_or(LocalSubgraphError::MalformedRecord)?;
    if record.len < K || record.len > record.words.len() * 32 {
        return Err(LocalSubgraphError::MalformedRecord);
    }
    let mut vertex = Kmer::<K>::from_bits(packed_prefix_bits::<K>(record.words));
    let mut reverse = vertex.reverse_complement();
    for offset in 0..=record.len - K {
        let canonical = vertex.min(reverse);
        if let Some(color_index) = wanted.get(canonical) {
            let sources = &mut source_sets[color_index];
            if sources.last().copied() != Some(source) {
                sources.push(source);
            }
        }
        if offset < record.len - K {
            let next = packed_base(record.words, offset + K);
            vertex = vertex.roll_forward(next);
            reverse = reverse.roll_backward(next.complement());
        }
    }
    Ok(())
}

#[inline]
fn packed_base(words: &[u64], idx: usize) -> Base {
    let word_idx = idx / 32;
    let shift = 2 * (31 - (idx % 32));
    match ((words[word_idx] >> shift) & 0b11) as u8 {
        0 => Base::A,
        1 => Base::C,
        2 => Base::G,
        3 => Base::T,
        _ => unreachable!(),
    }
}

#[inline]
fn packed_prefix_bits<const K: usize>(words: &[u64]) -> u128 {
    debug_assert!((1..=63).contains(&K));
    debug_assert!(words.len() * 32 >= K);
    if K <= 32 {
        (words[0] >> (2 * (32 - K))) as u128
    } else {
        let tail_bases = K - 32;
        ((words[0] as u128) << (2 * tail_bases)) | (words[1] >> (2 * (32 - tail_bases))) as u128
    }
}

#[inline]
fn local_u64_hash(key: u64) -> u64 {
    fast_u64_hash(key, 0xAAAA_AAAA_5555_5555)
}

#[derive(Debug)]
pub enum LocalSubgraphError {
    Bucket(BucketError),
    Color(ColorError),
    Kmer(KmerError),
    InvalidCutoff,
    EmptyBucketGroup,
    KMismatch { expected: usize, got: usize },
    GraphMismatch { expected: usize, got: usize },
    MalformedRecord,
    MissingVertex,
}

impl From<BucketError> for LocalSubgraphError {
    fn from(value: BucketError) -> Self {
        Self::Bucket(value)
    }
}

impl From<KmerError> for LocalSubgraphError {
    fn from(value: KmerError) -> Self {
        Self::Kmer(value)
    }
}

impl std::fmt::Display for LocalSubgraphError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Bucket(err) => write!(f, "{err}"),
            Self::Color(err) => write!(f, "{err}"),
            Self::Kmer(err) => write!(f, "{err}"),
            Self::InvalidCutoff => write!(f, "local subgraph cutoff must be at least 1"),
            Self::EmptyBucketGroup => write!(f, "local subgraph bucket group is empty"),
            Self::KMismatch { expected, got } => {
                write!(f, "bucket k mismatch: expected {expected}, got {got}")
            }
            Self::GraphMismatch { expected, got } => {
                write!(
                    f,
                    "bucket record graph mismatch: expected {expected}, got {got}"
                )
            }
            Self::MalformedRecord => write!(f, "malformed bucket record for local subgraph"),
            Self::MissingVertex => write!(f, "local subgraph edge references a missing vertex"),
        }
    }
}

impl std::error::Error for LocalSubgraphError {}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::GraphInput;
    use crate::buckets::BucketEmitter;
    use crate::params::BuildParams;
    use crate::partition::WeakSuperKmer;
    use std::fs;

    fn assert_packed_prefix<const K: usize>(seq: &[u8]) {
        let mut words = vec![0u64; seq.len().div_ceil(32)];
        for (idx, &base) in seq.iter().enumerate() {
            words[idx / 32] |= (Base::from_ascii(base).bits() as u64) << (2 * (31 - (idx % 32)));
        }
        assert_eq!(
            packed_prefix_bits::<K>(&words),
            Kmer::<K>::from_ascii(&seq[..K]).unwrap().as_u128()
        );
    }

    #[test]
    fn packed_prefix_decode_spans_word_boundary() {
        let seq = b"ACGTTGCATGTCGCATACGATCGTAGCTAGCTTGCATGACCTAGGCTAACGTTCGATGCATAC";
        assert_packed_prefix::<31>(seq);
        assert_packed_prefix::<33>(seq);
        assert_packed_prefix::<63>(seq);
    }

    #[test]
    fn builds_state_from_record_label() {
        let record = BucketRecord {
            graph_id: 3,
            len: 5,
            source_id: Some(1),
            left_discontinuous: true,
            right_discontinuous: true,
            label: b"ACGTT".to_vec(),
        };
        let mut subgraph = LocalSubgraph::<3> {
            graph_id: 3,
            colored: true,
            cutoff: 1,
            vertices: LocalVertexMap::with_capacity(0),
            edges: HashSet::with_hasher(FastBuildHasher::default()),
            stats: LocalSubgraphStats::default(),
        };

        subgraph.add_record(&record).unwrap();
        subgraph.stats.unique_vertices = subgraph.vertices.len() as u64;
        #[cfg(debug_assertions)]
        {
            subgraph.stats.unique_edges = subgraph.edges.len() as u64;
        }

        assert_eq!(subgraph.stats.observed_vertices, 3);
        assert_eq!(subgraph.stats.observed_edges, 2);
        assert_eq!(
            subgraph.stats.discontinuity_fronts + subgraph.stats.discontinuity_backs,
            2
        );
        #[cfg(debug_assertions)]
        assert_eq!(subgraph.stats.unique_edges, 2);

        let acg = Kmer::<3>::from_ascii(b"ACG").unwrap();
        let state = subgraph.vertex_state(acg).unwrap();
        assert!(state.is_discontinuous(Side::Front));
        assert_eq!(state.edge_at(Side::Back, 1), Base::T);
        assert_ne!(state.color_hash(), 0);
    }

    #[test]
    fn canonical_cycle_label_normalizes_rotation_and_strand() {
        let expected = b"AAACCGTTAAAC".to_vec();

        assert_eq!(
            canonical_cycle_label::<5>(b"CGTTAAACCGTT".to_vec()),
            expected
        );
        assert_eq!(
            canonical_cycle_label::<5>(b"GTTTAACGGTTT".to_vec()),
            expected
        );
    }

    #[test]
    fn source_sorted_buckets_make_color_hashes_exact() {
        // A vertex whose records arrive interleaved (source 1, then 2, then 1
        // again) would hash h(1)^h(2)^h(1) = h(2) if buckets were written in
        // arrival order: exactly the hash of a genuine {2}-only class, which
        // the reuse path would silently conflate. Source-sorted buckets keep
        // the classes distinct and make the hash order-independent.
        let build = |orders: &[(u32, &[u8; 5])], tag: &str| {
            let dir = std::env::temp_dir().join(format!(
                "cf3-colored-exact-{}-{:?}-{tag}",
                std::process::id(),
                std::thread::current().id(),
            ));
            let mut params = BuildParams::new(GraphInput::References, "colored-exact".into());
            params.k = 3;
            params.minimizer_len = 2;
            params.color = true;
            let mut emitter = BucketEmitter::create_in_dir(&params, 1, dir.clone()).unwrap();
            for &(source_id, seq) in orders {
                emitter
                    .add(
                        &WeakSuperKmer {
                            graph_id: 0,
                            offset: 0,
                            len: 5,
                            source_id: Some(source_id),
                            left_discontinuous: false,
                            right_discontinuous: false,
                        },
                        seq,
                    )
                    .unwrap();
            }
            emitter.finish().unwrap();
            let (store, entries) = BucketStore::open_dir(&dir).unwrap();
            let mut subgraph =
                LocalSubgraph::<3>::from_manifest_entries(&store, &entries, 1).unwrap();
            let mut unitigs = subgraph.contract_colored(&store, &entries).unwrap();
            unitigs.sort_by(|a, b| a.unitig.label.cmp(&b.unitig.label));
            fs::remove_dir_all(dir).unwrap();
            unitigs
        };

        // Interleaved arrival order: {1,2}-vertex records split around the
        // {2}-only sequence.
        let interleaved = build(
            &[(1, b"AACGT"), (2, b"TGCAA"), (2, b"AACGT"), (1, b"AACGT")],
            "interleaved",
        );
        // The same records already grouped by source.
        let grouped = build(
            &[(1, b"AACGT"), (1, b"AACGT"), (2, b"AACGT"), (2, b"TGCAA")],
            "grouped",
        );

        assert_eq!(interleaved.len(), 2);
        let hashes = |unitigs: &[ColoredLocalUnitig<3>]| {
            unitigs
                .iter()
                .map(|u| (u.unitig.label.clone(), u.colors[0].color_hash))
                .collect::<Vec<_>>()
        };
        // Order-independence: both arrival orders hash identically.
        assert_eq!(hashes(&interleaved), hashes(&grouped));
        // Distinctness: the {1,2} class must not alias the {2} class.
        assert_ne!(
            interleaved[0].colors[0].color_hash,
            interleaved[1].colors[0].color_hash
        );
        let sets = interleaved
            .iter()
            .map(|u| u.colors[0].sources.clone())
            .collect::<Vec<_>>();
        assert!(sets.contains(&vec![1, 2]) && sets.contains(&vec![2]));
    }

    #[test]
    fn colored_contraction_extracts_source_sets_at_color_transitions() {
        let dir = std::env::temp_dir().join(format!(
            "cf3-colored-local-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        let mut params = BuildParams::new(GraphInput::References, "colored-local".into());
        params.k = 3;
        params.minimizer_len = 2;
        params.color = true;
        let mut emitter = BucketEmitter::create_in_dir(&params, 1, dir.clone()).unwrap();
        for source_id in [2, 1, 2] {
            emitter
                .add(
                    &WeakSuperKmer {
                        graph_id: 0,
                        offset: 0,
                        len: 5,
                        source_id: Some(source_id),
                        left_discontinuous: false,
                        right_discontinuous: false,
                    },
                    b"AACGT",
                )
                .unwrap();
        }
        emitter.finish().unwrap();
        let (store, entries) = BucketStore::open_dir(&dir).unwrap();
        let mut subgraph = LocalSubgraph::<3>::from_manifest_entries(&store, &entries, 1).unwrap();
        let unitigs = subgraph.contract_colored(&store, &entries).unwrap();
        assert_eq!(unitigs.len(), 1);
        assert_eq!(unitigs[0].colors.len(), 1);
        assert_eq!(unitigs[0].colors[0].offset, 0);
        assert_eq!(unitigs[0].colors[0].sources, [1, 2]);
        fs::remove_dir_all(dir).unwrap();
    }
}

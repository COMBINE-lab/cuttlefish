use crate::Side;
#[cfg(test)]
use crate::buckets::BucketRecord;
use crate::buckets::{BorrowedBucketPackedRecord, BucketError, BucketManifestEntry, BucketReader};
use crate::color::{ColorError, ConcurrentColorRepository};
use crate::dna::Base;
use crate::hash::{FastBuildHasher, fast_u64_hash, hash_two_u64};
use crate::kmer::{Kmer, KmerError};
use crate::state::{ColorCoordinate, VertexState, source_hash};
use hashbrown::HashMap;
use hashbrown::HashTable;
use hashbrown::hash_map::RawEntryMut;
use std::collections::HashSet;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct DenseLocalVertexMap {
    indices: HashTable<u32>,
    entries: Vec<(u64, VertexState)>,
}

impl DenseLocalVertexMap {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            indices: HashTable::with_capacity(capacity),
            entries: Vec::with_capacity(capacity),
        }
    }

    fn clear_and_reserve(&mut self, capacity: usize) {
        self.indices.clear();
        self.entries.clear();
        if self.entries.capacity() < capacity {
            self.entries.reserve(capacity - self.entries.capacity());
        }
        if self.indices.capacity() < capacity {
            let entries = &self.entries;
            self.indices
                .reserve(capacity - self.indices.capacity(), |&index| {
                    local_u64_hash(entries[index as usize].0)
                });
        }
    }

    #[inline(always)]
    fn get(&self, key: u64) -> Option<&VertexState> {
        let hash = local_u64_hash(key);
        let index = *self
            .indices
            .find(hash, |&index| self.entries[index as usize].0 == key)?;
        Some(&self.entries[index as usize].1)
    }

    #[inline(always)]
    fn get_mut(&mut self, key: u64) -> Option<&mut VertexState> {
        let hash = local_u64_hash(key);
        let index = *self
            .indices
            .find(hash, |&index| self.entries[index as usize].0 == key)?;
        Some(&mut self.entries[index as usize].1)
    }

    #[inline(always)]
    fn state_or_default(&mut self, key: u64) -> &mut VertexState {
        let hash = local_u64_hash(key);
        if let Some(&index) = self
            .indices
            .find(hash, |&index| self.entries[index as usize].0 == key)
        {
            return &mut self.entries[index as usize].1;
        }
        let index = self.entries.len() as u32;
        self.entries.push((key, VertexState::default()));
        let entries = &self.entries;
        self.indices.insert_unique(hash, index, |&stored| {
            local_u64_hash(entries[stored as usize].0)
        });
        &mut self.entries[index as usize].1
    }
}

impl PartialEq for DenseLocalVertexMap {
    fn eq(&self, other: &Self) -> bool {
        self.entries.len() == other.entries.len()
            && self
                .entries
                .iter()
                .all(|&(key, state)| other.get(key) == Some(&state))
    }
}

impl Eq for DenseLocalVertexMap {}

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
    OneWord(DenseLocalVertexMap),
    TwoWord(HashMap<Kmer<K>, VertexState, FastBuildHasher>),
}

impl<const K: usize> LocalVertexMap<K> {
    #[inline(always)]
    fn with_capacity(capacity: usize) -> Self {
        if K <= 32 {
            Self::OneWord(DenseLocalVertexMap::with_capacity(capacity))
        } else {
            Self::TwoWord(HashMap::with_capacity_and_hasher(
                capacity,
                FastBuildHasher::default(),
            ))
        }
    }

    fn clear_and_reserve(&mut self, capacity: usize) {
        match self {
            Self::OneWord(map) => {
                map.clear_and_reserve(capacity);
            }
            Self::TwoWord(map) => {
                map.clear();
                if map.capacity() < capacity {
                    map.reserve(capacity);
                }
            }
        }
    }

    pub fn len(&self) -> usize {
        match self {
            Self::OneWord(map) => map.entries.len(),
            Self::TwoWord(map) => map.len(),
        }
    }

    #[inline(always)]
    fn get(&self, kmer: &Kmer<K>) -> Option<&VertexState> {
        match self {
            Self::OneWord(map) => map.get(kmer.as_u128() as u64),
            Self::TwoWord(map) => map.get(kmer),
        }
    }

    #[inline(always)]
    fn get_mut(&mut self, kmer: &Kmer<K>) -> Option<&mut VertexState> {
        match self {
            Self::OneWord(map) => map.get_mut(kmer.as_u128() as u64),
            Self::TwoWord(map) => map.get_mut(kmer),
        }
    }

    fn keys_vec(&self) -> Vec<Kmer<K>> {
        match self {
            Self::OneWord(map) => map
                .entries
                .iter()
                .map(|&(key, _)| Kmer::<K>::from_bits(key as u128))
                .collect(),
            Self::TwoWord(map) => map.keys().copied().collect(),
        }
    }

    fn dense_key_state(&self, index: usize) -> Option<(Kmer<K>, VertexState)> {
        match self {
            Self::OneWord(map) => map
                .entries
                .get(index)
                .map(|&(key, state)| (Kmer::<K>::from_bits(key as u128), state)),
            Self::TwoWord(_) => None,
        }
    }

    fn dense_len(&self) -> Option<usize> {
        match self {
            Self::OneWord(map) => Some(map.entries.len()),
            Self::TwoWord(_) => None,
        }
    }

    #[inline(always)]
    fn state_or_default(&mut self, kmer: Kmer<K>) -> &mut VertexState {
        match self {
            Self::OneWord(map) => {
                let key = kmer.as_u128() as u64;
                map.state_or_default(key)
            }
            Self::TwoWord(map) => {
                let hash = local_vertex_hash(kmer);
                match map
                    .raw_entry_mut()
                    .from_hash(hash, |stored| *stored == kmer)
                {
                    RawEntryMut::Occupied(entry) => entry.into_mut(),
                    RawEntryMut::Vacant(entry) => {
                        entry
                            .insert_with_hasher(hash, kmer, VertexState::default(), |stored| {
                                local_vertex_hash(*stored)
                            })
                            .1
                    }
                }
            }
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
        .map(|&b| Base::from_ascii(b).complement().to_ascii())
        .collect()
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

fn minimal_rotation(label: &[u8]) -> Vec<u8> {
    let start = least_rotation_start(label);
    label[start..]
        .iter()
        .chain(label[..start].iter())
        .copied()
        .collect()
}

fn least_rotation_start(s: &[u8]) -> usize {
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

impl<const K: usize> LocalSubgraph<K> {
    pub fn from_bucket_path(
        path: impl AsRef<Path>,
        cutoff: u32,
    ) -> Result<Self, LocalSubgraphError> {
        Self::from_bucket_paths_with_capacity(&[path.as_ref().to_path_buf()], cutoff, 0, None)
    }

    pub fn from_manifest_entries(
        entries: &[BucketManifestEntry],
        cutoff: u32,
    ) -> Result<Self, LocalSubgraphError> {
        let paths = entries
            .iter()
            .map(|entry| entry.path.clone())
            .collect::<Vec<_>>();
        // Weak super-kmers overlap heavily. HumGut's C++ construction observes
        // about 2.15 unique vertices per record (2.4 in the largest bucket), so
        // reserving K/4 slots wastes cache and scan bandwidth on sparse tables.
        let vertex_capacity = entries
            .iter()
            .map(|entry| entry.records as usize)
            .sum::<usize>()
            .saturating_mul(5)
            / 2;
        Self::from_bucket_paths_with_capacity(&paths, cutoff, vertex_capacity, None)
    }

    pub(crate) fn from_manifest_entries_reusing(
        entries: &[BucketManifestEntry],
        cutoff: u32,
        vertices: Option<LocalVertexMap<K>>,
    ) -> Result<Self, LocalSubgraphError> {
        let paths = entries
            .iter()
            .map(|entry| entry.path.clone())
            .collect::<Vec<_>>();
        let vertex_capacity = entries
            .iter()
            .map(|entry| entry.records as usize)
            .sum::<usize>()
            .saturating_mul(5)
            / 2;
        Self::from_bucket_paths_with_capacity(&paths, cutoff, vertex_capacity, vertices)
    }

    fn from_bucket_paths_with_capacity(
        paths: &[std::path::PathBuf],
        cutoff: u32,
        vertex_capacity: usize,
        reusable_vertices: Option<LocalVertexMap<K>>,
    ) -> Result<Self, LocalSubgraphError> {
        if cutoff == 0 {
            return Err(LocalSubgraphError::InvalidCutoff);
        }
        let Some(first_path) = paths.first() else {
            return Err(LocalSubgraphError::EmptyBucketGroup);
        };
        let mut reader = BucketReader::open(first_path)?;
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

        reader.try_for_each_borrowed_packed_record(|record| {
            subgraph.add_borrowed_packed_record(record)
        })?;
        for path in &paths[1..] {
            let mut reader = BucketReader::open(path)?;
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

    pub fn contract_compact(&mut self) -> Result<Vec<LocalUnitig<K>>, LocalSubgraphError> {
        self.contract_internal(false)
    }

    pub fn contract_colored(
        &mut self,
        entries: &[BucketManifestEntry],
    ) -> Result<Vec<ColoredLocalUnitig<K>>, LocalSubgraphError> {
        self.contract_colored_with_known(entries, |_| None)
    }

    pub fn contract_colored_with_known<F>(
        &mut self,
        entries: &[BucketManifestEntry],
        color_is_known: F,
    ) -> Result<Vec<ColoredLocalUnitig<K>>, LocalSubgraphError>
    where
        F: Fn(u64) -> Option<ColorCoordinate>,
    {
        self.contract_colored_impl(entries, color_is_known, None)
    }

    pub(crate) fn contract_colored_resolved(
        &mut self,
        entries: &[BucketManifestEntry],
        repository: &ConcurrentColorRepository,
        worker: usize,
    ) -> Result<Vec<ColoredLocalUnitig<K>>, LocalSubgraphError> {
        let mut resolve = |color_hash: u64, sources: &[u32]| {
            repository
                .resolve_or_insert(color_hash, sources, worker)
                .map_err(LocalSubgraphError::Color)
        };
        self.contract_colored_impl(
            entries,
            |color_hash| repository.get(color_hash),
            Some(&mut resolve),
        )
    }

    fn contract_colored_impl<F>(
        &mut self,
        entries: &[BucketManifestEntry],
        color_is_known: F,
        mut resolve_color: Option<
            &mut dyn FnMut(u64, &[u32]) -> Result<ColorCoordinate, LocalSubgraphError>,
        >,
    ) -> Result<Vec<ColoredLocalUnitig<K>>, LocalSubgraphError>
    where
        F: Fn(u64) -> Option<ColorCoordinate>,
    {
        if !self.colored {
            return Err(LocalSubgraphError::MalformedRecord);
        }
        let mut representative_indices = HashMap::<u64, usize, FastBuildHasher>::default();
        let mut representatives = Vec::<Kmer<K>>::new();
        let mut unitigs = Vec::new();
        let mut unitig_hash_runs = Vec::new();
        let mut back_walk = UnitigWalk::default();
        let mut front_walk = UnitigWalk::default();
        if let Some(vertex_count) = self.vertices.dense_len()
            && std::env::var_os("CF3_RS_SORT_LOCAL_VERTICES").is_none()
        {
            for index in 0..vertex_count {
                let (v_hat, state) = self
                    .vertices
                    .dense_key_state(index)
                    .expect("dense vertex index is in bounds");
                self.contract_colored_vertex(
                    v_hat,
                    state,
                    &color_is_known,
                    &mut unitigs,
                    &mut unitig_hash_runs,
                    &mut representative_indices,
                    &mut representatives,
                    &mut back_walk,
                    &mut front_walk,
                )?;
            }
        } else {
            let mut vertices = self.vertices.keys_vec();
            if std::env::var_os("CF3_RS_SORT_LOCAL_VERTICES").is_some() {
                vertices.sort_unstable();
            }
            for v_hat in vertices {
                let Some(state) = self.vertices.get(&v_hat).copied() else {
                    continue;
                };
                self.contract_colored_vertex(
                    v_hat,
                    state,
                    &color_is_known,
                    &mut unitigs,
                    &mut unitig_hash_runs,
                    &mut representative_indices,
                    &mut representatives,
                    &mut back_walk,
                    &mut front_walk,
                )?;
            }
        }

        let mut wanted = WantedColorMap::<K>::with_capacity(representatives.len());
        for (index, &vertex) in representatives.iter().enumerate() {
            wanted.insert(vertex, index);
        }
        let mut source_sets = vec![Vec::<u32>::new(); representatives.len()];
        for entry in entries {
            let mut reader = BucketReader::open(&entry.path)?;
            reader.try_for_each_borrowed_packed_record(|record| {
                collect_wanted_color_relations::<K>(record, &wanted, &mut source_sets)
            })?;
        }

        unitigs
            .into_iter()
            .zip(unitig_hash_runs)
            .map(|(unitig, runs)| {
                let colors = runs
                    .into_iter()
                    .map(|(offset, color_hash, coordinate)| {
                        let (coordinate, sources) = match coordinate {
                            Some(coordinate) => (Some(coordinate), Vec::new()),
                            None => {
                                let &index = representative_indices
                                    .get(&color_hash)
                                    .ok_or(LocalSubgraphError::MissingVertex)?;
                                match resolve_color.as_deref_mut() {
                                    Some(resolve) => (
                                        Some(resolve(color_hash, &source_sets[index])?),
                                        Vec::new(),
                                    ),
                                    None => (None, source_sets[index].clone()),
                                }
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

    #[allow(clippy::too_many_arguments)]
    fn contract_colored_vertex<F>(
        &mut self,
        v_hat: Kmer<K>,
        state: VertexState,
        color_is_known: &F,
        unitigs: &mut Vec<LocalUnitig<K>>,
        unitig_hash_runs: &mut Vec<Vec<(u32, u64, Option<ColorCoordinate>)>>,
        representative_indices: &mut HashMap<u64, usize, FastBuildHasher>,
        representatives: &mut Vec<Kmer<K>>,
        back_walk: &mut UnitigWalk<K>,
        front_walk: &mut UnitigWalk<K>,
    ) -> Result<(), LocalSubgraphError>
    where
        F: Fn(u64) -> Option<ColorCoordinate>,
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

        let unitig = self.extract_maximal_unitig_compact(v_hat, back_walk, front_walk)?;
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

        let mut runs = Vec::<(u32, u64, Option<ColorCoordinate>)>::new();
        let mut record_hash = |offset: usize, vertex: Kmer<K>, color_hash: u64| {
            if runs
                .last()
                .is_none_or(|&(_, previous, _)| previous != color_hash)
            {
                let coordinate = color_is_known(color_hash);
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
        if unitig.is_cycle {
            let mut vertex = Kmer::<K>::from_ascii(&unitig.label[..K])?;
            for offset in 0..=unitig.label.len() - K {
                let canonical = vertex.canonical();
                let color_hash = self
                    .vertices
                    .get(&canonical)
                    .ok_or(LocalSubgraphError::MissingVertex)?
                    .color_hash();
                record_hash(offset, canonical, color_hash);
                if offset < unitig.label.len() - K {
                    vertex = vertex.roll_forward(Base::from_ascii(unitig.label[offset + K]));
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
        unitigs.push(unitig);
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
        if let Some(vertex_count) = self.vertices.dense_len()
            && std::env::var_os("CF3_RS_SORT_LOCAL_VERTICES").is_none()
        {
            for index in 0..vertex_count {
                let (v_hat, state) = self
                    .vertices
                    .dense_key_state(index)
                    .expect("dense vertex index is in bounds");
                self.contract_vertex(
                    v_hat,
                    state,
                    collect_vertices,
                    &mut unitigs,
                    &mut back_walk,
                    &mut front_walk,
                )?;
            }
        } else {
            let mut vertices = self.vertices.keys_vec();
            if std::env::var_os("CF3_RS_SORT_LOCAL_VERTICES").is_some() {
                vertices.sort_unstable();
            }
            for v_hat in vertices {
                let Some(state) = self.vertices.get(&v_hat).copied() else {
                    continue;
                };
                self.contract_vertex(
                    v_hat,
                    state,
                    collect_vertices,
                    &mut unitigs,
                    &mut back_walk,
                    &mut front_walk,
                )?;
            }
        }

        Ok(unitigs)
    }

    fn contract_vertex(
        &mut self,
        v_hat: Kmer<K>,
        state: VertexState,
        collect_vertices: bool,
        unitigs: &mut Vec<LocalUnitig<K>>,
        back_walk: &mut UnitigWalk<K>,
        front_walk: &mut UnitigWalk<K>,
    ) -> Result<(), LocalSubgraphError> {
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
        unitigs.push(unitig);
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
                .map(|&base| Base::from_ascii(base).complement().to_ascii()),
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

    fn extract_maximal_unitig_compact(
        &mut self,
        v_hat: Kmer<K>,
        back_walk: &mut UnitigWalk<K>,
        front_walk: &mut UnitigWalk<K>,
    ) -> Result<LocalUnitig<K>, LocalSubgraphError> {
        let back_term = self.walk_unitig(v_hat, Side::Back, true, back_walk)?;

        if back_walk.is_cycle {
            return Ok(LocalUnitig {
                label: canonical_cycle_label::<K>(back_walk.label.clone()),
                vertices: Vec::new(),
                color_hashes: Vec::new(),
                left_exit: None,
                right_exit: None,
                is_cycle: true,
            });
        }

        let front_term = self.walk_unitig(v_hat, Side::Front, true, front_walk)?;
        let mut label = Vec::with_capacity(front_walk.label.len() + back_walk.label.len() - K);
        label.extend(
            front_walk
                .label
                .iter()
                .rev()
                .map(|&base| Base::from_ascii(base).complement().to_ascii()),
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
            vertices: Vec::new(),
            color_hashes: Vec::new(),
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

fn canonical_vertices_from_label<const K: usize>(
    label: &[u8],
) -> Result<Vec<Kmer<K>>, LocalSubgraphError> {
    if label.len() < K {
        return Err(LocalSubgraphError::MalformedRecord);
    }
    let mut vertices = Vec::with_capacity(label.len() - K + 1);
    let mut vertex = Kmer::<K>::from_ascii(&label[..K])?;
    vertices.push(vertex.canonical());
    for &base in &label[K..] {
        vertex = vertex.roll_forward(Base::from_ascii(base));
        vertices.push(vertex.canonical());
    }
    Ok(vertices)
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
fn local_vertex_hash<const K: usize>(kmer: Kmer<K>) -> u64 {
    let words = kmer.words();
    hash_two_u64(words[0], words[1])
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
    use crate::buckets::{BucketEmitter, read_manifest};
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
        for source_id in [1, 2] {
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
        let entries = read_manifest(&dir).unwrap();
        let mut subgraph = LocalSubgraph::<3>::from_manifest_entries(&entries, 1).unwrap();
        let unitigs = subgraph.contract_colored(&entries).unwrap();
        assert_eq!(unitigs.len(), 1);
        assert_eq!(unitigs[0].colors.len(), 1);
        assert_eq!(unitigs[0].colors[0].offset, 0);
        assert_eq!(unitigs[0].colors[0].sources, [1, 2]);
        fs::remove_dir_all(dir).unwrap();
    }
}

use crate::Side;
#[cfg(test)]
use crate::buckets::BucketRecord;
use crate::buckets::{BucketError, BucketManifestEntry, BucketPackedRecord, BucketReader};
use crate::dna::Base;
use crate::hash::{FastBuildHasher, hash_two_u64};
use crate::kmer::{Kmer, KmerError};
use crate::state::VertexState;
use hashbrown::HashMap;
use hashbrown::hash_map::RawEntryMut;
use std::collections::HashSet;
use std::path::Path;

pub type LocalVertexMap<const K: usize> = HashMap<Kmer<K>, VertexState, FastBuildHasher>;
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
    pub left_exit: Option<(Kmer<K>, Side)>,
    pub right_exit: Option<(Kmer<K>, Side)>,
    pub is_cycle: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct UnitigWalk<const K: usize> {
    label: Vec<u8>,
    vertices: Vec<Kmer<K>>,
    is_cycle: bool,
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
    fn new(observed: Kmer<K>) -> Self {
        Self {
            observed,
            reverse: observed.reverse_complement(),
        }
    }

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
    fn init(v: DirectedKmer<K>, collect_vertices: bool) -> Self {
        Self {
            label: v.observed.to_ascii_string().into_bytes(),
            vertices: if collect_vertices {
                vec![v.canonical()]
            } else {
                Vec::new()
            },
            is_cycle: false,
        }
    }

    fn extend(
        &mut self,
        v: DirectedKmer<K>,
        base: Base,
        anchor: Kmer<K>,
        collect_vertices: bool,
    ) -> bool {
        if v.canonical() == anchor {
            self.is_cycle = true;
            return false;
        }

        self.label.push(base.to_ascii());
        if collect_vertices {
            self.vertices.push(v.canonical());
        }
        true
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
        Self::from_bucket_paths_with_capacity(&[path.as_ref().to_path_buf()], cutoff, 0)
    }

    pub fn from_manifest_entries(
        entries: &[BucketManifestEntry],
        cutoff: u32,
    ) -> Result<Self, LocalSubgraphError> {
        let paths = entries
            .iter()
            .map(|entry| entry.path.clone())
            .collect::<Vec<_>>();
        let vertices_per_record_hint = (K / 4).max(4);
        let vertex_capacity = entries
            .iter()
            .map(|entry| entry.records as usize)
            .sum::<usize>()
            .saturating_mul(vertices_per_record_hint);
        Self::from_bucket_paths_with_capacity(&paths, cutoff, vertex_capacity)
    }

    fn from_bucket_paths_with_capacity(
        paths: &[std::path::PathBuf],
        cutoff: u32,
        vertex_capacity: usize,
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
        let mut subgraph = Self {
            graph_id,
            colored,
            cutoff,
            vertices: HashMap::with_capacity_and_hasher(
                vertex_capacity,
                FastBuildHasher::default(),
            ),
            edges: HashSet::with_hasher(FastBuildHasher::default()),
            stats: LocalSubgraphStats::default(),
        };

        let mut record = BucketPackedRecord::default();
        reader
            .try_for_each_packed_record(&mut record, |record| subgraph.add_packed_record(record))?;
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
            reader.try_for_each_packed_record(&mut record, |record| {
                subgraph.add_packed_record(record)
            })?;
        }

        subgraph.stats.unique_vertices = subgraph.vertices.len() as u64;
        subgraph.stats.unique_edges = subgraph.edges.len() as u64;

        Ok(subgraph)
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

    fn contract_internal(
        &mut self,
        collect_vertices: bool,
    ) -> Result<Vec<LocalUnitig<K>>, LocalSubgraphError> {
        let mut unitigs = Vec::new();
        let mut vertices = self.vertices.keys().copied().collect::<Vec<_>>();
        if std::env::var_os("CF3_RS_SORT_LOCAL_VERTICES").is_some() {
            vertices.sort_unstable();
        }

        for v_hat in vertices {
            let Some(state) = self.vertices.get(&v_hat).copied() else {
                continue;
            };
            if state.is_visited() {
                continue;
            }
            if state.is_isolated(self.cutoff) {
                if let Some(st) = self.vertices.get_mut(&v_hat) {
                    st.mark_visited();
                }
                self.stats.isolated_vertices += 1;
                continue;
            }

            let unitig = self.extract_maximal_unitig(v_hat, collect_vertices)?;
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
        }

        Ok(unitigs)
    }

    fn extract_maximal_unitig(
        &mut self,
        v_hat: Kmer<K>,
        collect_vertices: bool,
    ) -> Result<LocalUnitig<K>, LocalSubgraphError> {
        let (back_walk, back_term) = self.walk_unitig(v_hat, Side::Back, collect_vertices)?;

        if back_walk.is_cycle {
            return Ok(LocalUnitig {
                label: canonical_cycle_label::<K>(back_walk.label.clone()),
                vertices: back_walk.vertices,
                left_exit: None,
                right_exit: None,
                is_cycle: true,
            });
        }

        let (front_walk, front_term) = self.walk_unitig(v_hat, Side::Front, collect_vertices)?;
        let vertices = if collect_vertices {
            let mut vertices = front_walk.vertices;
            vertices.reverse();
            vertices.extend(back_walk.vertices.into_iter().skip(1));
            vertices
        } else {
            Vec::new()
        };

        let mut label = reverse_complement_label(&front_walk.label);
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
    ) -> Result<(UnitigWalk<K>, WalkTermination<K>), LocalSubgraphError> {
        let icc_return_side = start_side.inverse();
        let mut v = if start_side == Side::Back {
            DirectedKmer::new(v_hat)
        } else {
            DirectedKmer::new(v_hat.reverse_complement())
        };
        let mut side = start_side;
        let mut walk = UnitigWalk::init(v, collect_vertices);

        loop {
            let canonical = v.canonical();
            let state = {
                let state = self
                    .vertices
                    .get_mut(&canonical)
                    .ok_or(LocalSubgraphError::MissingVertex)?;
                let copied = *state;
                state.mark_visited();
                copied
            };

            let mut edge = state.edge_at(side, self.cutoff);
            if edge == Base::N {
                return Ok((walk, WalkTermination::Branched));
            }
            if edge == Base::E {
                if !state.is_discontinuous(side) {
                    return Ok((walk, WalkTermination::DeadEnded));
                }
                return Ok((walk, WalkTermination::Exited(v)));
            }

            if side == Side::Front {
                edge = edge.complement();
            }
            v = v.roll_forward(edge);

            let next_state = *self
                .vertices
                .get(&v.canonical())
                .ok_or(LocalSubgraphError::MissingVertex)?;
            side = v.entrance_side();
            if next_state.is_branching_side(side, self.cutoff) {
                return Ok((walk, WalkTermination::Crossed));
            }
            if next_state.is_visited() {
                if v.canonical() == v_hat && side == icc_return_side {
                    walk.is_cycle = true;
                }
                return Ok((walk, WalkTermination::Crossed));
            }

            if !walk.extend(v, edge, v_hat, collect_vertices) {
                return Ok((walk, WalkTermination::Crossed));
            }
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

    fn add_packed_record(&mut self, record: &BucketPackedRecord) -> Result<(), LocalSubgraphError> {
        if record.graph_id != self.graph_id {
            return Err(LocalSubgraphError::GraphMismatch {
                expected: self.graph_id,
                got: record.graph_id,
            });
        }
        if record.len < K || record.len > record.words.len() * 32 {
            return Err(LocalSubgraphError::MalformedRecord);
        }

        self.stats.weak_superkmers += 1;
        self.stats.weak_superkmer_bases += record.len as u64;

        let last_vertex_offset = record.len - K;
        let mut observed_bits = 0u128;
        for idx in 0..K {
            observed_bits = (observed_bits << 2) | packed_base(&record.words, idx).bits() as u128;
        }
        let mut observed = Kmer::<K>::from_bits(observed_bits);
        let mut reverse = observed.reverse_complement();
        let mut prev = None;
        for offset in 0..=last_vertex_offset {
            let in_canonical_form = observed <= reverse;
            let canonical = if in_canonical_form { observed } else { reverse };
            let pred_base = if offset == 0 {
                Base::E
            } else {
                packed_base(&record.words, offset - 1)
            };
            let succ_base = if offset == last_vertex_offset {
                Base::E
            } else {
                packed_base(&record.words, offset + K)
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

    #[inline]
    fn vertex_state_or_default(&mut self, kmer: Kmer<K>) -> &mut VertexState {
        let hash = local_vertex_hash(kmer);
        match self
            .vertices
            .raw_entry_mut()
            .from_hash(hash, |key| *key == kmer)
        {
            RawEntryMut::Occupied(entry) => entry.into_mut(),
            RawEntryMut::Vacant(entry) => {
                entry
                    .insert_with_hasher(hash, kmer, VertexState::default(), |key| {
                        local_vertex_hash(*key)
                    })
                    .1
            }
        }
    }
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
fn local_vertex_hash<const K: usize>(kmer: Kmer<K>) -> u64 {
    let words = kmer.words();
    hash_two_u64(words[0], words[1])
}

#[derive(Debug)]
pub enum LocalSubgraphError {
    Bucket(BucketError),
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
            vertices: HashMap::with_hasher(FastBuildHasher::default()),
            edges: HashSet::with_hasher(FastBuildHasher::default()),
            stats: LocalSubgraphStats::default(),
        };

        subgraph.add_record(&record).unwrap();
        subgraph.stats.unique_vertices = subgraph.vertices.len() as u64;
        subgraph.stats.unique_edges = subgraph.edges.len() as u64;

        assert_eq!(subgraph.stats.observed_vertices, 3);
        assert_eq!(subgraph.stats.observed_edges, 2);
        assert_eq!(
            subgraph.stats.discontinuity_fronts + subgraph.stats.discontinuity_backs,
            2
        );
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
}

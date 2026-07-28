//! End-to-end uncolored graph construction.
//!
//! The production entry point streams local contraction into the external
//! discontinuity pipeline. A global in-memory contractor remains available
//! only as a debug/reference implementation.

use crate::Side;
use crate::buckets::{BucketError, BucketStore};
use crate::discontinuity::{
    spawn_background_dir_removal,
    SerialUncoloredCollator, emit_uncolored_external_discontinuity_inputs_with_threads_in_dir,
    report_process_memory, trim_process_allocations,
};
use crate::dna::{Base, complement_ascii};
use crate::hash::FastBuildHasher;
use crate::kmer::{Kmer, KmerError};
use crate::params::BuildParams;
use crate::state::VertexState;
use crate::subgraph::LocalSubgraphError;
use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

type FastHashMap<K, V> = HashMap<K, V, FastBuildHasher>;

/// Summary of a completed uncolored graph build.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UncoloredBuildStats {
    pub input_buckets: usize,
    pub bucket_records: u64,
    pub observed_edges: u64,
    pub retained_edges: u64,
    pub unitigs: u64,
    pub unitig_bases: u64,
    pub output_path: PathBuf,
}

/// Builds an uncolored compacted de Bruijn graph from partition buckets.
///
/// The external-memory path is used unless the debug contractor environment
/// switch is explicitly enabled.
pub fn build_uncolored_from_buckets<const K: usize>(
    params: &BuildParams,
    bucket_dir: impl AsRef<Path>,
) -> Result<UncoloredBuildStats, UncoloredBuildError> {
    if std::env::var_os("CF3_RS_DEBUG_GLOBAL_CONTRACTOR").is_some() {
        build_uncolored_with_debug_global_contractor::<K>(params, bucket_dir)
    } else {
        build_uncolored_with_serial_discontinuity_pipeline::<K>(params, bucket_dir)
    }
}

pub fn build_uncolored_with_debug_global_contractor<const K: usize>(
    params: &BuildParams,
    bucket_dir: impl AsRef<Path>,
) -> Result<UncoloredBuildStats, UncoloredBuildError> {
    if params.color {
        return Err(UncoloredBuildError::ColoredUnsupported);
    }

    let bucket_dir = bucket_dir.as_ref();
    let cutoff = params.cutoff();
    let mut bucket_records = 0u64;
    let mut observed_edges = 0u64;
    let mut graph = DebugGlobalCanonicalGraph::<K>::new(cutoff);

    let (store, entries) = BucketStore::open_dir(bucket_dir)?;
    for entry in &entries {
        let mut reader = store.reader(entry)?;
        let header = reader.header();
        if header.k != params.k || header.minimizer_len != params.minimizer_len || header.colored {
            return Err(UncoloredBuildError::BucketParamsMismatch {
                path: bucket_dir.to_path_buf(),
            });
        }

        while let Some(record) = reader.next_record()? {
            bucket_records += 1;
            observed_edges += graph.add_label(&record.label)?;
        }
    }

    let retained_edges = graph.unique_edges.len() as u64;
    let mut unitigs = graph.contract()?;
    normalize_unitigs(&mut unitigs);

    let output_path = PathBuf::from(format!("{}.fa", params.output_prefix));
    write_fasta(&output_path, &unitigs)?;

    Ok(UncoloredBuildStats {
        input_buckets: entries.len(),
        bucket_records,
        observed_edges,
        retained_edges,
        unitigs: unitigs.len() as u64,
        unitig_bases: unitigs.iter().map(|u| u.len() as u64).sum(),
        output_path,
    })
}

pub fn build_uncolored_with_serial_discontinuity_pipeline<const K: usize>(
    params: &BuildParams,
    bucket_dir: impl AsRef<Path>,
) -> Result<UncoloredBuildStats, UncoloredBuildError> {
    if params.color {
        return Err(UncoloredBuildError::ColoredUnsupported);
    }

    let local_start = Instant::now();
    report_process_memory("before local contraction");
    let output_name = Path::new(&params.output_prefix)
        .file_name()
        .and_then(|s| s.to_str())
        .filter(|s| !s.is_empty())
        .unwrap_or("cuttlefish3");
    let label_path =
        PathBuf::from(&params.work_dir).join(format!("{output_name}.cf3rs.lmtig-labels"));
    let bucket_dir = bucket_dir.as_ref().to_path_buf();
    let local_threads = params.local_workers();
    eprintln!("cuttlefish3-rs: local contraction using {local_threads} worker(s)");
    let output_path = PathBuf::from(format!("{}.fa", params.output_prefix));
    let mut inputs = emit_uncolored_external_discontinuity_inputs_with_threads_in_dir::<K>(
        &bucket_dir,
        params.cutoff(),
        local_threads,
        &label_path,
        Some(&output_path),
    )?;
    let local_elapsed = local_start.elapsed();
    report_process_memory("after local contraction before trim");
    eprintln!(
        "cuttlefish3-rs: local contraction emitted {} unitig(s), {} discontinuity exit(s)",
        inputs.stats.local_unitigs, inputs.stats.discontinuity_exits
    );
    eprintln!(
        "cuttlefish3-rs: local contraction phase completed in {:.3}s",
        local_elapsed.as_secs_f64()
    );
    // With buckets in containers this is a bulk delete of 129 files holding
    // hundreds of gigabytes, not the single leftover manifest the per-file
    // layout left behind -- local contraction unlinked those as it consumed
    // them. So it goes to the background unlinkers, like the edge-matrix and
    // expansion directories, and is joined before the build returns.
    let bucket_reclaim = (std::env::var_os("CF3_RS_KEEP_INTERMEDIATES").is_none())
        .then(|| spawn_background_dir_removal(bucket_dir.clone()));
    trim_process_allocations();
    report_process_memory("after local contraction trim");
    eprintln!("cuttlefish3-rs: collating final unitigs");
    let collation_start = Instant::now();
    report_process_memory("before collation");
    let coord_dir =
        PathBuf::from(&params.work_dir).join(format!("{output_name}.cf3rs.stitch-coords"));
    let final_dir =
        PathBuf::from(&params.work_dir).join(format!("{output_name}.cf3rs.final-unitigs"));
    eprintln!("cuttlefish3-rs: writing FASTA to {}", output_path.display());
    let post_local_threads = params.post_local_workers();
    eprintln!("cuttlefish3-rs: collation using {post_local_threads} worker(s)");
    let stats = SerialUncoloredCollator::collate_external_stitched_to_fasta_with_threads_in_dir(
        &mut inputs,
        post_local_threads,
        &coord_dir,
        &final_dir,
        &output_path,
    )?;
    let collation_elapsed = collation_start.elapsed();
    report_process_memory("after collation");
    eprintln!(
        "cuttlefish3-rs: collation and FASTA write completed in {:.3}s",
        collation_elapsed.as_secs_f64()
    );

    if let Some(handle) = bucket_reclaim {
        let _ = handle.join();
    }

    Ok(UncoloredBuildStats {
        input_buckets: inputs.stats.input_buckets,
        bucket_records: inputs.stats.weak_superkmers,
        observed_edges: inputs.stats.discontinuity_exits,
        retained_edges: inputs.stats.discontinuity_exits,
        unitigs: stats.emitted_unitigs,
        unitig_bases: stats.emitted_bases,
        output_path,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct CanonicalEdge<const K: usize> {
    from: Kmer<K>,
    to: Kmer<K>,
}

struct DebugGlobalCanonicalGraph<const K: usize> {
    cutoff: u32,
    vertices: FastHashMap<Kmer<K>, VertexState>,
    unique_edges: BTreeSet<CanonicalEdge<K>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct DirectedKmer<const K: usize> {
    observed: Kmer<K>,
}

impl<const K: usize> DirectedKmer<K> {
    fn new(observed: Kmer<K>) -> Self {
        Self { observed }
    }

    #[inline]
    fn canonical(self) -> Kmer<K> {
        self.observed.canonical()
    }

    #[inline]
    fn in_canonical_form(self) -> bool {
        self.observed.is_canonical()
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
        Self::new(self.observed.roll_forward(base))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct UnitigWalk<const K: usize> {
    label: Vec<u8>,
    anchor: Kmer<K>,
}

impl<const K: usize> UnitigWalk<K> {
    fn init(v: DirectedKmer<K>) -> Self {
        Self {
            label: v.observed.to_ascii_string().into_bytes(),
            anchor: v.canonical(),
        }
    }

    fn extend(&mut self, v: DirectedKmer<K>, base: Base) -> bool {
        if v.canonical() == self.anchor {
            return false;
        }

        self.label.push(base.to_ascii());
        true
    }
}

impl<const K: usize> DebugGlobalCanonicalGraph<K> {
    fn new(cutoff: u32) -> Self {
        Self {
            cutoff,
            vertices: FastHashMap::default(),
            unique_edges: BTreeSet::new(),
        }
    }

    fn add_label(&mut self, label: &[u8]) -> Result<u64, KmerError> {
        if label.len() < K {
            return Ok(0);
        }

        let last_vertex_offset = label.len() - K;
        let mut prev = None;
        let mut observed_edges = 0u64;
        for offset in 0..=last_vertex_offset {
            let directed = DirectedKmer::new(Kmer::<K>::from_ascii(&label[offset..offset + K])?);
            let canonical = directed.canonical();
            let pred_base = if offset == 0 {
                Base::E
            } else {
                Base::from_ascii(label[offset - 1])
            };
            let succ_base = if offset == last_vertex_offset {
                Base::E
            } else {
                Base::from_ascii(label[offset + K])
            };
            let mut front = if directed.in_canonical_form() {
                pred_base
            } else {
                succ_base.complement()
            };
            let mut back = if directed.in_canonical_form() {
                succ_base
            } else {
                pred_base.complement()
            };

            if offset > 0 && Some(canonical) == prev {
                if directed.in_canonical_form() {
                    front = Base::E;
                } else {
                    back = Base::E;
                }
            }

            self.vertices
                .entry(canonical)
                .or_default()
                .update_edges(front, back);

            if let Some(from) = prev {
                self.unique_edges.insert(CanonicalEdge {
                    from,
                    to: canonical,
                });
                observed_edges += 1;
            }
            prev = Some(canonical);
        }

        Ok(observed_edges)
    }

    fn contract(&mut self) -> Result<Vec<Vec<u8>>, UncoloredBuildError> {
        let mut unitigs = Vec::new();
        let mut vertices = self.vertices.keys().copied().collect::<Vec<_>>();
        vertices.sort_unstable();

        for v_hat in vertices {
            let Some(state) = self.vertices.get(&v_hat).copied() else {
                continue;
            };
            if state.is_visited() || state.is_isolated(self.cutoff) {
                continue;
            }

            unitigs.push(self.extract_maximal_unitig(v_hat)?);
        }

        Ok(unitigs)
    }

    fn extract_maximal_unitig(&mut self, v_hat: Kmer<K>) -> Result<Vec<u8>, UncoloredBuildError> {
        let (back_walk, back_is_cycle) = self.walk_unitig(v_hat, Side::Back)?;
        if back_is_cycle {
            return Ok(canonical_label(back_walk.label));
        }

        let (front_walk, _) = self.walk_unitig(v_hat, Side::Front)?;
        let mut label = reverse_complement_label(&front_walk.label);
        label.extend_from_slice(&back_walk.label[K..]);
        Ok(canonical_label(label))
    }

    fn walk_unitig(
        &mut self,
        v_hat: Kmer<K>,
        start_side: Side,
    ) -> Result<(UnitigWalk<K>, bool), UncoloredBuildError> {
        let icc_return_side = start_side.inverse();
        let mut v = if start_side == Side::Back {
            DirectedKmer::new(v_hat)
        } else {
            DirectedKmer::new(v_hat.reverse_complement())
        };
        let mut side = start_side;
        let mut walk = UnitigWalk::init(v);

        loop {
            let canonical = v.canonical();
            let state = *self
                .vertices
                .get(&canonical)
                .ok_or(UncoloredBuildError::MissingVertex)?;
            self.vertices
                .get_mut(&canonical)
                .ok_or(UncoloredBuildError::MissingVertex)?
                .mark_visited();

            let mut edge = state.edge_at(side, self.cutoff);
            if edge == Base::N || edge == Base::E {
                return Ok((walk, false));
            }

            if side == Side::Front {
                edge = edge.complement();
            }
            v = v.roll_forward(edge);

            let next_state = *self
                .vertices
                .get(&v.canonical())
                .ok_or(UncoloredBuildError::MissingVertex)?;
            side = v.entrance_side();
            if next_state.is_branching_side(side, self.cutoff) {
                return Ok((walk, false));
            }
            if next_state.is_visited() {
                return Ok((walk, v.canonical() == v_hat && side == icc_return_side));
            }

            if !walk.extend(v, edge) {
                return Ok((walk, false));
            }
            side = side.inverse();
        }
    }
}

fn normalize_unitigs(unitigs: &mut Vec<Vec<u8>>) {
    for label in unitigs.iter_mut() {
        let rc = reverse_complement_label(label);
        if rc < *label {
            *label = rc;
        }
    }

    unitigs.sort_unstable();
    unitigs.dedup();
}

fn reverse_complement_label(label: &[u8]) -> Vec<u8> {
    label
        .iter()
        .rev()
        .map(|&base| complement_ascii(base))
        .collect()
}

fn canonical_label(label: Vec<u8>) -> Vec<u8> {
    let rc = reverse_complement_label(&label);
    if rc < label { rc } else { label }
}

fn write_fasta(path: &Path, unitigs: &[Vec<u8>]) -> Result<(), UncoloredBuildError> {
    let file = File::create(path).map_err(|source| UncoloredBuildError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut out = BufWriter::new(file);

    for label in unitigs {
        writeln!(out, ">0").map_err(|source| UncoloredBuildError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        out.write_all(label)
            .and_then(|_| out.write_all(b"\n"))
            .map_err(|source| UncoloredBuildError::Io {
                path: path.to_path_buf(),
                source,
            })?;
    }

    out.flush().map_err(|source| UncoloredBuildError::Io {
        path: path.to_path_buf(),
        source,
    })
}

#[derive(Debug)]
pub enum UncoloredBuildError {
    Bucket(BucketError),
    DiscontinuityInput(crate::discontinuity::DiscontinuityInputError),
    SerialCollation(crate::discontinuity::SerialCollationError),
    SerialEdgeMatrix(crate::discontinuity::SerialEdgeMatrixError),
    LocalSubgraph(LocalSubgraphError),
    Kmer(KmerError),
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    ColoredUnsupported,
    BucketParamsMismatch {
        path: PathBuf,
    },
    MissingVertex,
}

impl From<BucketError> for UncoloredBuildError {
    fn from(value: BucketError) -> Self {
        Self::Bucket(value)
    }
}

impl From<crate::discontinuity::DiscontinuityInputError> for UncoloredBuildError {
    fn from(value: crate::discontinuity::DiscontinuityInputError) -> Self {
        Self::DiscontinuityInput(value)
    }
}

impl From<crate::discontinuity::SerialCollationError> for UncoloredBuildError {
    fn from(value: crate::discontinuity::SerialCollationError) -> Self {
        Self::SerialCollation(value)
    }
}

impl From<crate::discontinuity::SerialEdgeMatrixError> for UncoloredBuildError {
    fn from(value: crate::discontinuity::SerialEdgeMatrixError) -> Self {
        Self::SerialEdgeMatrix(value)
    }
}

impl From<LocalSubgraphError> for UncoloredBuildError {
    fn from(value: LocalSubgraphError) -> Self {
        Self::LocalSubgraph(value)
    }
}

impl From<KmerError> for UncoloredBuildError {
    fn from(value: KmerError) -> Self {
        Self::Kmer(value)
    }
}

impl std::fmt::Display for UncoloredBuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Bucket(err) => write!(f, "{err}"),
            Self::DiscontinuityInput(err) => write!(f, "{err}"),
            Self::SerialCollation(err) => write!(f, "{err}"),
            Self::SerialEdgeMatrix(err) => write!(f, "{err}"),
            Self::LocalSubgraph(err) => write!(f, "{err}"),
            Self::Kmer(err) => write!(f, "{err}"),
            Self::Io { path, source } => write!(f, "{}: {source}", path.display()),
            Self::ColoredUnsupported => {
                write!(
                    f,
                    "colored graph output is not implemented in the Rust path yet"
                )
            }
            Self::BucketParamsMismatch { path } => write!(
                f,
                "weak-superkmer bucket parameters do not match the build request: {}",
                path.display()
            ),
            Self::MissingVertex => write!(f, "canonical graph edge references a missing vertex"),
        }
    }
}

impl std::error::Error for UncoloredBuildError {}

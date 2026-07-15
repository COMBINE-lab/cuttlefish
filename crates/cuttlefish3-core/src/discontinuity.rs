use crate::DEFAULT_VERTEX_PARTITIONS;
use crate::Side;
use crate::buckets::{BucketError, BucketManifestEntry, read_manifest};
use crate::color::{
    ColorError, ColorRepositoryManifest, ColorRunSidecar, ColorRunSidecarWriter,
    ConcurrentColorRepository, ConcurrentColorRunSidecarWriter, append_color_runs,
    reverse_color_runs, reverse_color_runs_in_place,
};
use crate::dna::{Base, complement_ascii};
use crate::hash::{FastBuildHasher, hash_bytes, hash_u64, wyhash_u64};
use crate::kmer::Kmer;
use crate::state::{UnitigColor, VertexState};
use crate::subgraph::{LocalSubgraph, LocalSubgraphError, LocalVertexMap};
use leapfrog::LeapMap;
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};
use scc::{HashMap as SccHashMap, hash_map::Entry as SccEntry};
use std::cell::UnsafeCell;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::marker::PhantomData;
use std::mem::MaybeUninit;
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::{
    Arc, Mutex,
    atomic::{AtomicU8, AtomicU64, AtomicUsize, Ordering},
    mpsc,
};
use std::time::{Duration, Instant};
use xxhash_rust::xxh3::xxh3_64;

type FastHashMap<K, V> = HashMap<K, V, FastBuildHasher>;
type FastHashSet<T> = HashSet<T, FastBuildHasher>;

fn keep_intermediates() -> bool {
    std::env::var_os("CF3_RS_KEEP_INTERMEDIATES").is_some()
}

fn remove_serial_file(path: &Path) -> Result<(), SerialCollationError> {
    if keep_intermediates() {
        return Ok(());
    }
    match fs::remove_file(path) {
        Ok(()) => Ok(()),
        Err(source) if source.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(source) => Err(SerialCollationError::Io {
            path: path.to_path_buf(),
            source,
        }),
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiscontinuityInputs<const K: usize> {
    pub unitigs: Vec<DiscontinuityUnitig<K>>,
    labels: Vec<u8>,
    pub stats: DiscontinuityInputStats,
}

#[derive(Debug)]
pub struct ExternalDiscontinuityInputs<const K: usize> {
    unitig_path: PathBuf,
    label_path: PathBuf,
    unitigs: usize,
    compact_unitigs: bool,
    ranges: Vec<ExternalLocalUnitigRange>,
    edge_matrix: Option<BlockedEdgeMatrix<K>>,
    color_runs: Option<ColorRunSidecar>,
    local_unitig_buckets: Option<Vec<LocalUnitigBucketEntry>>,
    trivial_fasta: Option<(PathBuf, u64, u64)>,
    color_repository: Option<ColorRepositoryManifest>,
    pub stats: DiscontinuityInputStats,
}

impl<const K: usize> ExternalDiscontinuityInputs<K> {
    pub fn unitig_count(&self) -> usize {
        self.unitigs
    }

    pub fn color_runs(&self) -> Option<&ColorRunSidecar> {
        self.color_runs.as_ref()
    }

    pub fn color_repository(&self) -> Option<&ColorRepositoryManifest> {
        self.color_repository.as_ref()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ExternalLocalUnitigRange {
    start_unitig: usize,
    unitigs: usize,
    label_start: u64,
    label_len: u64,
    color_start: u64,
}

#[derive(Debug)]
struct LocalUnitigBucketEntry {
    bucket_id: u16,
    unitig_path: PathBuf,
    label_path: PathBuf,
    unitigs: usize,
    color_runs: Option<ColorRunSidecar>,
}

struct LocalUnitigBucketWriter {
    bucket_id: u16,
    unitig_path: PathBuf,
    label_path: PathBuf,
    unitigs: BufWriter<File>,
    labels: BufWriter<File>,
    color_runs: Option<ColorRunSidecarWriter>,
    unitig_count: usize,
}

impl LocalUnitigBucketWriter {
    fn create(dir: &Path, bucket_id: u16, colored: bool) -> Result<Self, DiscontinuityInputError> {
        let unitig_path = dir.join(format!("{bucket_id:03}.unitigs"));
        let label_path = dir.join(format!("{bucket_id:03}.labels"));
        let unitig_file =
            File::create(&unitig_path).map_err(|source| DiscontinuityInputError::Io {
                path: unitig_path.clone(),
                source,
            })?;
        let label_file =
            File::create(&label_path).map_err(|source| DiscontinuityInputError::Io {
                path: label_path.clone(),
                source,
            })?;
        let color_runs = colored
            .then(|| ColorRunSidecarWriter::create(dir.join(format!("{bucket_id:03}.colors"))))
            .transpose()?;
        Ok(Self {
            bucket_id,
            unitig_path,
            label_path,
            unitigs: BufWriter::with_capacity(1024 * 1024, unitig_file),
            labels: BufWriter::with_capacity(4 * 1024 * 1024, label_file),
            color_runs,
            unitig_count: 0,
        })
    }

    fn write<const K: usize>(
        &mut self,
        labels: &[u8],
        unitigs: &[DiscontinuityUnitig<K>],
        colors: Option<&[Vec<UnitigColor>]>,
    ) -> Result<usize, DiscontinuityInputError> {
        if colors.is_some_and(|runs| runs.len() != unitigs.len()) {
            return Err(DiscontinuityInputError::MissingColorRuns);
        }
        let base = self.unitig_count;
        self.labels
            .write_all(labels)
            .map_err(|source| DiscontinuityInputError::Io {
                path: self.label_path.clone(),
                source,
            })?;
        for (index, unitig) in unitigs.iter().enumerate() {
            write_discontinuity_unitig_record(&mut self.unitigs, &self.unitig_path, unitig, true)?;
            if let Some(writer) = self.color_runs.as_mut() {
                writer.write_unitig(&colors.expect("colored bucket has colors")[index])?;
            }
        }
        self.unitig_count += unitigs.len();
        Ok(base)
    }

    fn finish(mut self) -> Result<LocalUnitigBucketEntry, DiscontinuityInputError> {
        self.unitigs
            .flush()
            .map_err(|source| DiscontinuityInputError::Io {
                path: self.unitig_path.clone(),
                source,
            })?;
        self.labels
            .flush()
            .map_err(|source| DiscontinuityInputError::Io {
                path: self.label_path.clone(),
                source,
            })?;
        Ok(LocalUnitigBucketEntry {
            bucket_id: self.bucket_id,
            unitig_path: self.unitig_path,
            label_path: self.label_path,
            unitigs: self.unitig_count,
            color_runs: self
                .color_runs
                .map(ColorRunSidecarWriter::finish)
                .transpose()?,
        })
    }
}

fn finish_local_unitig_writers(
    writers: Vec<Mutex<LocalUnitigBucketWriter>>,
    threads: usize,
) -> Result<Vec<LocalUnitigBucketEntry>, DiscontinuityInputError> {
    let worker_count = threads.max(1).min(writers.len().max(1));
    let mut work = (0..worker_count).map(|_| Vec::new()).collect::<Vec<_>>();
    for (index, writer) in writers.into_iter().enumerate() {
        work[index % worker_count].push(writer);
    }
    let mut entries = std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(worker_count);
        for worker_writers in work {
            handles.push(scope.spawn(move || {
                let mut entries = Vec::with_capacity(worker_writers.len());
                for writer in worker_writers {
                    entries.push(
                        writer
                            .into_inner()
                            .map_err(|_| DiscontinuityInputError::WorkerPanic)?
                            .finish()?,
                    );
                }
                Ok::<_, DiscontinuityInputError>(entries)
            }));
        }
        let mut entries = Vec::new();
        for handle in handles {
            entries.extend(
                handle
                    .join()
                    .map_err(|_| DiscontinuityInputError::WorkerPanic)??,
            );
        }
        Ok::<_, DiscontinuityInputError>(entries)
    })?;
    entries.sort_unstable_by_key(|entry| entry.bucket_id);
    Ok(entries)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ExternalLabelRef {
    label_start: u64,
    label_len: u32,
}

#[cfg(all(target_os = "linux", target_env = "gnu"))]
unsafe extern "C" {
    fn malloc_trim(pad: usize) -> i32;
    fn getrlimit(resource: i32, limits: *mut RLimit) -> i32;
    fn setrlimit(resource: i32, limits: *const RLimit) -> i32;
}

#[cfg(all(target_os = "linux", target_env = "gnu"))]
#[repr(C)]
struct RLimit {
    current: u64,
    maximum: u64,
}

#[cfg(all(target_os = "linux", target_env = "gnu"))]
fn raise_open_file_limit() -> bool {
    const RLIMIT_NOFILE: i32 = 7;
    unsafe {
        let mut limits = RLimit {
            current: 0,
            maximum: 0,
        };
        if getrlimit(RLIMIT_NOFILE, &mut limits) != 0 {
            return false;
        }
        limits.current = limits.maximum;
        setrlimit(RLIMIT_NOFILE, &limits) == 0
    }
}

#[cfg(not(all(target_os = "linux", target_env = "gnu")))]
fn raise_open_file_limit() -> bool {
    false
}

pub fn trim_process_allocations() {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    unsafe {
        let _ = malloc_trim(0);
    }
}

pub fn report_process_memory(label: &str) {
    if std::env::var_os("CF3_RS_PROFILE_RSS").is_none() {
        return;
    }
    let Ok(status) = std::fs::read_to_string("/proc/self/status") else {
        return;
    };
    let mut rss_kib = None;
    let mut hwm_kib = None;
    for line in status.lines() {
        if let Some(value) = line.strip_prefix("VmRSS:") {
            rss_kib = parse_status_kib(value);
        } else if let Some(value) = line.strip_prefix("VmHWM:") {
            hwm_kib = parse_status_kib(value);
        }
    }
    if let (Some(rss), Some(hwm)) = (rss_kib, hwm_kib) {
        eprintln!("cuttlefish3-rs: rss {label}: current {rss} KiB, peak {hwm} KiB");
    }
}

fn parse_status_kib(value: &str) -> Option<u64> {
    value
        .split_whitespace()
        .next()
        .and_then(|value| value.parse::<u64>().ok())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct DiscontinuityInputStats {
    pub input_buckets: usize,
    pub weak_superkmers: u64,
    pub local_unitigs: u64,
    pub discontinuity_exits: u64,
    pub unitig_bases: u64,
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DiscontinuityUnitig<const K: usize> {
    pub label_start: u64,
    left_vertex: Kmer<K>,
    right_vertex: Kmer<K>,
    pub label_len: u32,
    flags: u8,
}

const DISCONTINUITY_LEFT_EXIT: u8 = 1 << 0;
const DISCONTINUITY_LEFT_BACK: u8 = 1 << 1;
const DISCONTINUITY_RIGHT_EXIT: u8 = 1 << 2;
const DISCONTINUITY_RIGHT_BACK: u8 = 1 << 3;
const DISCONTINUITY_CYCLE: u8 = 1 << 4;

impl<const K: usize> DiscontinuityInputs<K> {
    pub fn empty(stats: DiscontinuityInputStats) -> Self {
        Self {
            unitigs: Vec::new(),
            labels: Vec::new(),
            stats,
        }
    }

    pub fn from_unitigs(
        unitigs: impl IntoIterator<Item = OwnedDiscontinuityUnitig<K>>,
        stats: DiscontinuityInputStats,
    ) -> Self {
        let mut inputs = Self::empty(stats);
        for unitig in unitigs {
            inputs.push_unitig(unitig);
        }
        inputs
    }

    pub fn push_unitig(&mut self, unitig: OwnedDiscontinuityUnitig<K>) -> usize {
        let index = self.unitigs.len();
        let label_start = self.labels.len();
        let label_len = unitig.label.len();
        assert!(u64::try_from(label_start).is_ok());
        assert!(u32::try_from(label_len).is_ok());
        self.labels.extend_from_slice(&unitig.label);
        self.unitigs.push(DiscontinuityUnitig {
            label_start: label_start as u64,
            left_vertex: unitig
                .left_exit
                .map(|endpoint| endpoint.vertex)
                .unwrap_or_else(Kmer::zero),
            right_vertex: unitig
                .right_exit
                .map(|endpoint| endpoint.vertex)
                .unwrap_or_else(Kmer::zero),
            label_len: label_len as u32,
            flags: discontinuity_unitig_flags(unitig.left_exit, unitig.right_exit, unitig.is_cycle),
        });
        index
    }

    #[inline]
    pub fn label(&self, unitig_index: usize) -> &[u8] {
        self.unitigs[unitig_index].label(self)
    }

    #[inline]
    pub fn try_label(&self, unitig_index: usize) -> Option<&[u8]> {
        self.unitigs
            .get(unitig_index)
            .map(|unitig| unitig.label(self))
    }
}

impl<const K: usize> DiscontinuityUnitig<K> {
    #[inline]
    pub fn label<'a>(&self, inputs: &'a DiscontinuityInputs<K>) -> &'a [u8] {
        let start = self.label_start as usize;
        let end = start + self.label_len as usize;
        &inputs.labels[start..end]
    }

    #[inline]
    pub fn left_exit(&self) -> Option<DiscontinuityEndpoint<K>> {
        (self.flags & DISCONTINUITY_LEFT_EXIT != 0).then(|| DiscontinuityEndpoint {
            vertex: self.left_vertex,
            side: if self.flags & DISCONTINUITY_LEFT_BACK != 0 {
                Side::Back
            } else {
                Side::Front
            },
        })
    }

    #[inline]
    pub fn right_exit(&self) -> Option<DiscontinuityEndpoint<K>> {
        (self.flags & DISCONTINUITY_RIGHT_EXIT != 0).then(|| DiscontinuityEndpoint {
            vertex: self.right_vertex,
            side: if self.flags & DISCONTINUITY_RIGHT_BACK != 0 {
                Side::Back
            } else {
                Side::Front
            },
        })
    }

    #[inline]
    pub fn is_cycle(&self) -> bool {
        self.flags & DISCONTINUITY_CYCLE != 0
    }
}

struct ExternalDiscontinuityReader<const K: usize> {
    label_file: File,
    unitig_path: PathBuf,
    label_path: PathBuf,
    unitigs: usize,
    compact_unitigs: bool,
}

impl<const K: usize> ExternalDiscontinuityReader<K> {
    fn open(inputs: &ExternalDiscontinuityInputs<K>) -> Result<Self, SerialCollationError> {
        let label_file =
            File::open(&inputs.label_path).map_err(|source| SerialCollationError::Io {
                path: inputs.label_path.clone(),
                source,
            })?;
        Ok(Self {
            label_file,
            unitig_path: inputs.unitig_path.clone(),
            label_path: inputs.label_path.clone(),
            unitigs: inputs.unitigs,
            compact_unitigs: inputs.compact_unitigs,
        })
    }

    fn read_label(
        &self,
        unitig: &DiscontinuityUnitig<K>,
        scratch: &mut Vec<u8>,
    ) -> Result<(), SerialCollationError> {
        scratch.resize(unitig.label_len as usize, 0);
        self.label_file
            .read_exact_at(scratch, unitig.label_start)
            .map_err(|source| SerialCollationError::Io {
                path: self.label_path.clone(),
                source,
            })
    }

    fn iter(&self) -> Result<ExternalDiscontinuityIter<K>, SerialCollationError> {
        let file = File::open(&self.unitig_path).map_err(|source| SerialCollationError::Io {
            path: self.unitig_path.clone(),
            source,
        })?;
        Ok(ExternalDiscontinuityIter {
            input: BufReader::with_capacity(1024 * 1024, file),
            path: self.unitig_path.clone(),
            remaining: self.unitigs,
            compact_unitigs: self.compact_unitigs,
        })
    }
}

struct ExternalDiscontinuityIter<const K: usize> {
    input: BufReader<File>,
    path: PathBuf,
    remaining: usize,
    compact_unitigs: bool,
}

impl<const K: usize> ExternalDiscontinuityIter<K> {
    fn next_unitig(&mut self) -> Result<DiscontinuityUnitig<K>, SerialCollationError> {
        if self.remaining == 0 {
            return Err(SerialCollationError::MalformedCoordBucket(
                self.path.clone(),
            ));
        }
        self.remaining -= 1;
        read_discontinuity_unitig_from_reader(&mut self.input, &self.path, self.compact_unitigs)
    }
}

impl<const K: usize> Iterator for ExternalDiscontinuityIter<K> {
    type Item = Result<DiscontinuityUnitig<K>, SerialCollationError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            return None;
        }
        Some(self.next_unitig())
    }
}

fn read_discontinuity_unitig_from_reader<const K: usize>(
    input: &mut BufReader<File>,
    path: &Path,
    compact: bool,
) -> Result<DiscontinuityUnitig<K>, SerialCollationError> {
    if compact {
        let mut bytes = [0u8; 8];
        input
            .read_exact(&mut bytes)
            .map_err(|source| SerialCollationError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        return Ok(DiscontinuityUnitig {
            label_start: 0,
            left_vertex: Kmer::zero(),
            right_vertex: Kmer::zero(),
            label_len: u32::from_le_bytes(bytes[..4].try_into().expect("label length")),
            flags: bytes[4],
        });
    }
    let mut out = MaybeUninit::<DiscontinuityUnitig<K>>::zeroed();
    let bytes = unsafe {
        std::slice::from_raw_parts_mut(
            out.as_mut_ptr().cast::<u8>(),
            std::mem::size_of::<DiscontinuityUnitig<K>>(),
        )
    };
    input
        .read_exact(bytes)
        .map_err(|source| SerialCollationError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(unsafe { out.assume_init() })
}

const fn external_unitig_record_len<const K: usize>(compact: bool) -> usize {
    if compact {
        8
    } else {
        std::mem::size_of::<DiscontinuityUnitig<K>>()
    }
}

fn read_discontinuity_unitig_from_reader_for_input<const K: usize>(
    input: &mut BufReader<File>,
    path: &Path,
) -> Result<DiscontinuityUnitig<K>, DiscontinuityInputError> {
    let mut out = MaybeUninit::<DiscontinuityUnitig<K>>::zeroed();
    let bytes = unsafe {
        std::slice::from_raw_parts_mut(
            out.as_mut_ptr().cast::<u8>(),
            std::mem::size_of::<DiscontinuityUnitig<K>>(),
        )
    };
    input
        .read_exact(bytes)
        .map_err(|source| DiscontinuityInputError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(unsafe { out.assume_init() })
}

#[inline]
fn discontinuity_unitig_flags<const K: usize>(
    left_exit: Option<DiscontinuityEndpoint<K>>,
    right_exit: Option<DiscontinuityEndpoint<K>>,
    is_cycle: bool,
) -> u8 {
    let mut flags = 0u8;
    if let Some(endpoint) = left_exit {
        flags |= DISCONTINUITY_LEFT_EXIT;
        if endpoint.side == Side::Back {
            flags |= DISCONTINUITY_LEFT_BACK;
        }
    }
    if let Some(endpoint) = right_exit {
        flags |= DISCONTINUITY_RIGHT_EXIT;
        if endpoint.side == Side::Back {
            flags |= DISCONTINUITY_RIGHT_BACK;
        }
    }
    if is_cycle {
        flags |= DISCONTINUITY_CYCLE;
    }
    flags
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OwnedDiscontinuityUnitig<const K: usize> {
    pub graph_id: usize,
    pub label: Vec<u8>,
    pub left_exit: Option<DiscontinuityEndpoint<K>>,
    pub right_exit: Option<DiscontinuityEndpoint<K>>,
    pub is_cycle: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DiscontinuityEndpoint<const K: usize> {
    pub vertex: Kmer<K>,
    pub side: Side,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SerialEdgeMatrix<const K: usize> {
    vertex_partitions: usize,
    blocks: Vec<Vec<Vec<DiscontinuityEdge<K>>>>,
    stats: SerialEdgeMatrixStats,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SerialEdgeMatrixStats {
    pub edges: u64,
    pub phi_edges: u64,
    pub diagonal_edges: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiscontinuityEdge<const K: usize> {
    pub first: MatrixEndpoint<K>,
    pub second: MatrixEndpoint<K>,
    pub weight: u64,
    pub unitig_bucket: u16,
    pub unitig_index: usize,
    pub unitig_exit_side: Side,
    pub phantom_unitig: Option<DiscontinuityEndpoint<K>>,
    pub swapped: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MatrixEndpoint<const K: usize> {
    Phi,
    Vertex(DiscontinuityEndpoint<K>),
}

impl<const K: usize> SerialEdgeMatrix<K> {
    pub fn new(vertex_partitions: usize) -> Result<Self, SerialEdgeMatrixError> {
        if vertex_partitions == 0 || !vertex_partitions.is_power_of_two() {
            return Err(SerialEdgeMatrixError::InvalidPartitionCount(
                vertex_partitions,
            ));
        }

        let partition_count = vertex_partitions + 1;
        Ok(Self {
            vertex_partitions,
            blocks: vec![vec![Vec::new(); partition_count]; partition_count],
            stats: SerialEdgeMatrixStats::default(),
        })
    }

    pub fn from_inputs(
        inputs: &DiscontinuityInputs<K>,
        vertex_partitions: usize,
    ) -> Result<Self, SerialEdgeMatrixError> {
        let mut matrix = Self::new(vertex_partitions)?;
        for (unitig_index, unitig) in inputs.unitigs.iter().enumerate() {
            match (unitig.left_exit(), unitig.right_exit()) {
                (Some(left), Some(right)) => matrix.add_edge(
                    MatrixEndpoint::Vertex(left),
                    MatrixEndpoint::Vertex(right),
                    1,
                    unitig_index,
                ),
                (Some(endpoint), None) => matrix.add_edge_with_orientation(
                    MatrixEndpoint::Phi,
                    MatrixEndpoint::Vertex(endpoint),
                    1,
                    unitig_index,
                    Side::Front,
                ),
                (None, Some(endpoint)) => matrix.add_edge_with_orientation(
                    MatrixEndpoint::Phi,
                    MatrixEndpoint::Vertex(endpoint),
                    1,
                    unitig_index,
                    Side::Back,
                ),
                (None, None) => {}
            }
        }
        Ok(matrix)
    }

    pub fn add_edge(
        &mut self,
        first: MatrixEndpoint<K>,
        second: MatrixEndpoint<K>,
        weight: u64,
        unitig_index: usize,
    ) {
        self.add_edge_with_orientation(first, second, weight, unitig_index, Side::Back);
    }

    pub fn add_edge_with_orientation(
        &mut self,
        first: MatrixEndpoint<K>,
        second: MatrixEndpoint<K>,
        weight: u64,
        unitig_index: usize,
        unitig_exit_side: Side,
    ) {
        self.add_edge_with_orientation_and_phantom(
            first,
            second,
            weight,
            unitig_index,
            unitig_exit_side,
            None,
        );
    }

    pub fn add_edge_with_orientation_and_phantom(
        &mut self,
        first: MatrixEndpoint<K>,
        second: MatrixEndpoint<K>,
        weight: u64,
        unitig_index: usize,
        unitig_exit_side: Side,
        phantom_unitig: Option<DiscontinuityEndpoint<K>>,
    ) {
        let first_partition = self.partition(first);
        let second_partition = self.partition(second);
        let swapped = first_partition > second_partition;
        let unitig_exit_side = if swapped {
            unitig_exit_side.inverse()
        } else {
            unitig_exit_side
        };
        let (row, col, first, second) = if swapped {
            (second_partition, first_partition, second, first)
        } else {
            (first_partition, second_partition, first, second)
        };

        self.stats.edges += 1;
        if first.is_phi() || second.is_phi() {
            self.stats.phi_edges += 1;
        }
        if row == col {
            self.stats.diagonal_edges += 1;
        }

        self.blocks[row][col].push(DiscontinuityEdge {
            first,
            second,
            weight,
            unitig_bucket: 0,
            unitig_index,
            unitig_exit_side,
            phantom_unitig,
            swapped,
        });
    }

    #[inline]
    pub const fn vertex_partitions(&self) -> usize {
        self.vertex_partitions
    }

    #[inline]
    pub const fn partition_count(&self) -> usize {
        self.vertex_partitions + 1
    }

    #[inline]
    pub const fn stats(&self) -> SerialEdgeMatrixStats {
        self.stats
    }

    #[inline]
    pub fn partition(&self, endpoint: MatrixEndpoint<K>) -> usize {
        match endpoint {
            MatrixEndpoint::Phi => 0,
            MatrixEndpoint::Vertex(endpoint) => {
                ((endpoint.vertex.hash64(0) as usize) & (self.vertex_partitions - 1)) + 1
            }
        }
    }

    #[inline]
    pub fn block(&self, row: usize, col: usize) -> &[DiscontinuityEdge<K>] {
        assert!(row <= col);
        &self.blocks[row][col]
    }

    pub fn edges(&self) -> impl Iterator<Item = &DiscontinuityEdge<K>> {
        self.blocks
            .iter()
            .enumerate()
            .flat_map(|(row, blocks)| blocks[row..].iter())
            .flat_map(|block| block.iter())
    }
}

impl<const K: usize> MatrixEndpoint<K> {
    #[inline]
    pub const fn is_phi(self) -> bool {
        matches!(self, Self::Phi)
    }
}

const BLOCKED_EDGE_WRITE_BUFFER_BYTES: usize = 256 * 1024;

#[derive(Debug)]
struct BlockedEdgeBlock {
    path: PathBuf,
    edges: usize,
    record_len: usize,
    buffer: Vec<u8>,
}

#[derive(Debug)]
struct BlockedEdgeMatrix<const K: usize> {
    dir: PathBuf,
    vertex_partitions: usize,
    blocks: Vec<BlockedEdgeBlock>,
    stats: SerialEdgeMatrixStats,
    phantom: PhantomData<[(); K]>,
}

struct PreparedBlockedEdge {
    block: usize,
    bytes: [u8; 72],
    phi: bool,
    diagonal: bool,
}

enum BlockedReadTask<'a> {
    File {
        file: &'a File,
        path: &'a Path,
        offset: u64,
        len: usize,
    },
    Memory(&'a [u8]),
}

struct ConcurrentBlockedAppend {
    path: PathBuf,
    record_len: usize,
    buffer: Mutex<Vec<u8>>,
    edges: AtomicUsize,
}

struct ConcurrentBlockedEdgeWriters {
    partition_count: usize,
    blocks: Vec<ConcurrentBlockedAppend>,
    phi_edges: AtomicU64,
    diagonal_edges: AtomicU64,
}

impl ConcurrentBlockedEdgeWriters {
    fn new<const K: usize>(matrix: &BlockedEdgeMatrix<K>) -> Self {
        Self {
            partition_count: matrix.partition_count(),
            blocks: matrix
                .blocks
                .iter()
                .map(|block| ConcurrentBlockedAppend {
                    path: block.path.clone(),
                    record_len: block.record_len,
                    buffer: Mutex::new(Vec::new()),
                    edges: AtomicUsize::new(0),
                })
                .collect(),
            phi_edges: AtomicU64::new(0),
            diagonal_edges: AtomicU64::new(0),
        }
    }

    fn add(&self, edge: &PreparedBlockedEdge) -> Result<(), SerialCollationError> {
        let block = &self.blocks[edge.block];
        let mut buffer = block
            .buffer
            .lock()
            .unwrap_or_else(|poison| poison.into_inner());
        if buffer.len() + block.record_len > BLOCKED_EDGE_WRITE_BUFFER_BYTES {
            append_blocked_bytes(&block.path, &buffer)?;
            buffer.clear();
        }
        buffer.extend_from_slice(&edge.bytes[..block.record_len]);
        block.edges.fetch_add(1, Ordering::Relaxed);
        Ok(())
    }

    fn add_batch(&self, edges: &[PreparedBlockedEdge]) -> Result<(), SerialCollationError> {
        if edges.is_empty() {
            return Ok(());
        }
        let block = &self.blocks[edges[0].block];
        let mut buffer = block
            .buffer
            .lock()
            .unwrap_or_else(|poison| poison.into_inner());
        for edge in edges {
            debug_assert_eq!(edge.block, edges[0].block);
            if buffer.len() + block.record_len > BLOCKED_EDGE_WRITE_BUFFER_BYTES {
                append_blocked_bytes(&block.path, &buffer)?;
                buffer.clear();
            }
            buffer.extend_from_slice(&edge.bytes[..block.record_len]);
        }
        block.edges.fetch_add(edges.len(), Ordering::Relaxed);
        Ok(())
    }

    fn add_prepared_edges<const K: usize>(
        &self,
        edges: &mut [PreparedBlockedEdge],
        unitig_base: usize,
        unitig_bucket: u16,
    ) -> Result<(), SerialCollationError> {
        if edges.is_empty() {
            return Ok(());
        }
        edges.sort_unstable_by_key(|edge| edge.block);
        let unitig_off = blocked_edge_unitig_offset::<K>();
        let mut start = 0;
        while start < edges.len() {
            let block_id = edges[start].block;
            let mut end = start + 1;
            while end < edges.len() && edges[end].block == block_id {
                end += 1;
            }
            let block = &self.blocks[block_id];
            let mut buffer = block
                .buffer
                .lock()
                .unwrap_or_else(|poison| poison.into_inner());
            for edge in &edges[start..end] {
                let mut bytes = edge.bytes;
                add_unitig_base_to_encoded_edge::<K>(&mut bytes, unitig_off, unitig_base);
                set_encoded_edge_unitig_bucket::<K>(&mut bytes, unitig_bucket);
                if buffer.len() + block.record_len > BLOCKED_EDGE_WRITE_BUFFER_BYTES {
                    append_blocked_bytes(&block.path, &buffer)?;
                    buffer.clear();
                }
                buffer.extend_from_slice(&bytes[..block.record_len]);
            }
            block.edges.fetch_add(end - start, Ordering::Relaxed);
            start = end;
        }
        self.phi_edges.fetch_add(
            edges.iter().filter(|edge| edge.phi).count() as u64,
            Ordering::Relaxed,
        );
        self.diagonal_edges.fetch_add(
            edges.iter().filter(|edge| edge.diagonal).count() as u64,
            Ordering::Relaxed,
        );
        Ok(())
    }

    fn finish_into<const K: usize>(
        &self,
        matrix: &mut BlockedEdgeMatrix<K>,
    ) -> Result<(), SerialCollationError> {
        let mut edges = 0u64;
        for (block, target) in self.blocks.iter().zip(&mut matrix.blocks) {
            let mut buffer = block
                .buffer
                .lock()
                .unwrap_or_else(|poison| poison.into_inner());
            if !buffer.is_empty() {
                append_blocked_bytes(&block.path, &buffer)?;
                buffer.clear();
            }
            target.edges = block.edges.load(Ordering::Relaxed);
            edges += target.edges as u64;
        }
        matrix.stats.edges = edges;
        matrix.stats.phi_edges = self.phi_edges.load(Ordering::Relaxed);
        matrix.stats.diagonal_edges = self.diagonal_edges.load(Ordering::Relaxed);
        Ok(())
    }

    fn flush_column(&self, col: usize) -> Result<(), SerialCollationError> {
        for row in 0..=col {
            self.flush_block(row * self.partition_count + col)?;
        }
        Ok(())
    }

    fn flush_block(&self, block_id: usize) -> Result<(), SerialCollationError> {
        let block = &self.blocks[block_id];
        let mut buffer = block
            .buffer
            .lock()
            .unwrap_or_else(|poison| poison.into_inner());
        if !buffer.is_empty() {
            append_blocked_bytes(&block.path, &buffer)?;
            buffer.clear();
        }
        Ok(())
    }

    fn take_added_edges(&self, row: usize, col: usize) -> usize {
        self.blocks[row * self.partition_count + col]
            .edges
            .swap(0, Ordering::Relaxed)
    }
}

fn append_blocked_bytes(path: &Path, bytes: &[u8]) -> Result<(), SerialCollationError> {
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .map_err(|source| SerialCollationError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    file.write_all(bytes)
        .map_err(|source| SerialCollationError::Io {
            path: path.to_path_buf(),
            source,
        })
}

impl<const K: usize> BlockedEdgeMatrix<K> {
    fn create(dir: &Path, vertex_partitions: usize) -> Result<Self, SerialCollationError> {
        if dir.exists() {
            fs::remove_dir_all(dir).map_err(|source| SerialCollationError::Io {
                path: dir.to_path_buf(),
                source,
            })?;
        }
        fs::create_dir_all(dir).map_err(|source| SerialCollationError::Io {
            path: dir.to_path_buf(),
            source,
        })?;
        let partition_count = vertex_partitions + 1;
        let record_len = discontinuity_edge_record_len::<K>();
        let mut blocks = Vec::with_capacity(partition_count * partition_count);
        for row in 0..partition_count {
            for col in 0..partition_count {
                blocks.push(BlockedEdgeBlock {
                    path: dir.join(format!("{row}.{col}.edge")),
                    edges: 0,
                    record_len,
                    buffer: Vec::new(),
                });
            }
        }
        Ok(Self {
            dir: dir.to_path_buf(),
            vertex_partitions,
            blocks,
            stats: SerialEdgeMatrixStats::default(),
            phantom: PhantomData,
        })
    }

    #[inline]
    fn partition_count(&self) -> usize {
        self.vertex_partitions + 1
    }

    #[inline]
    fn partition(&self, endpoint: MatrixEndpoint<K>) -> usize {
        edge_matrix_partition(self.vertex_partitions, endpoint)
    }

    #[inline]
    fn block_index(&self, row: usize, col: usize) -> usize {
        row * self.partition_count() + col
    }

    fn add_edge_record(
        &mut self,
        mut edge: DiscontinuityEdge<K>,
    ) -> Result<(), SerialCollationError> {
        let first_partition = self.partition(edge.first);
        let second_partition = self.partition(edge.second);
        if first_partition > second_partition {
            std::mem::swap(&mut edge.first, &mut edge.second);
            edge.unitig_exit_side = edge.unitig_exit_side.inverse();
            edge.swapped = !edge.swapped;
        }
        let row = first_partition.min(second_partition);
        let col = first_partition.max(second_partition);
        self.stats.edges += 1;
        self.stats.phi_edges += u64::from(edge.first.is_phi() || edge.second.is_phi());
        self.stats.diagonal_edges += u64::from(row == col);
        let idx = self.block_index(row, col);
        let record = encode_discontinuity_edge(&edge);
        let block = &mut self.blocks[idx];
        if block.buffer.len() + block.record_len > BLOCKED_EDGE_WRITE_BUFFER_BYTES {
            flush_blocked_edge_block(block)?;
        }
        block.buffer.extend_from_slice(&record[..block.record_len]);
        block.edges += 1;
        if block.buffer.len() + block.record_len > BLOCKED_EDGE_WRITE_BUFFER_BYTES {
            flush_blocked_edge_block(block)?;
        }
        Ok(())
    }

    fn add_prepared_edges(
        &mut self,
        edges: &[PreparedBlockedEdge],
        unitig_base: usize,
    ) -> Result<(), SerialCollationError> {
        let unitig_off = blocked_edge_unitig_offset::<K>();
        for edge in edges {
            let mut bytes = edge.bytes;
            add_unitig_base_to_encoded_edge::<K>(&mut bytes, unitig_off, unitig_base);
            let block = &mut self.blocks[edge.block];
            if block.buffer.len() + block.record_len > BLOCKED_EDGE_WRITE_BUFFER_BYTES {
                flush_blocked_edge_block(block)?;
            }
            block.buffer.extend_from_slice(&bytes[..block.record_len]);
            block.edges += 1;
            self.stats.edges += 1;
            self.stats.phi_edges += u64::from(edge.phi);
            self.stats.diagonal_edges += u64::from(edge.diagonal);
        }
        Ok(())
    }

    fn add_prepared_edges_absolute(
        &mut self,
        edges: &[PreparedBlockedEdge],
    ) -> Result<(), SerialCollationError> {
        for edge in edges {
            let block = &mut self.blocks[edge.block];
            if block.buffer.len() + block.record_len > BLOCKED_EDGE_WRITE_BUFFER_BYTES {
                flush_blocked_edge_block(block)?;
            }
            block
                .buffer
                .extend_from_slice(&edge.bytes[..block.record_len]);
            block.edges += 1;
            self.stats.edges += 1;
            self.stats.phi_edges += u64::from(edge.phi);
            self.stats.diagonal_edges += u64::from(edge.diagonal);
        }
        Ok(())
    }

    fn flush_block(&mut self, row: usize, col: usize) -> Result<(), SerialCollationError> {
        let idx = self.block_index(row, col);
        flush_blocked_edge_block(&mut self.blocks[idx])
    }

    fn flush_all(&mut self) -> Result<(), SerialCollationError> {
        for block in &mut self.blocks {
            flush_blocked_edge_block(block)?;
        }
        Ok(())
    }

    fn flush_all_with_threads(&mut self, threads: usize) -> Result<(), SerialCollationError> {
        let workers = threads.max(1).min(self.blocks.len());
        if workers == 1 {
            return self.flush_all();
        }
        let chunk = self.blocks.len().div_ceil(workers);
        std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(workers);
            for blocks in self.blocks.chunks_mut(chunk) {
                handles.push(scope.spawn(move || {
                    for block in blocks {
                        flush_blocked_edge_block(block)?;
                    }
                    Ok::<_, SerialCollationError>(())
                }));
            }
            for handle in handles {
                handle
                    .join()
                    .map_err(|_| SerialCollationError::WorkerPanic)??;
            }
            Ok(())
        })
    }

    fn read_flushed_block(
        &self,
        row: usize,
        col: usize,
    ) -> Result<Vec<DiscontinuityEdge<K>>, SerialCollationError> {
        let idx = self.block_index(row, col);
        let block = &self.blocks[idx];
        if block.edges == 0 {
            return Ok(Vec::new());
        }
        let bytes = match fs::read(&block.path) {
            Ok(bytes) => bytes,
            Err(source) if source.kind() == std::io::ErrorKind::NotFound => Vec::new(),
            Err(source) => {
                return Err(SerialCollationError::Io {
                    path: block.path.clone(),
                    source,
                });
            }
        };
        if bytes.len() + block.buffer.len() != block.edges * block.record_len {
            return Err(SerialCollationError::MalformedCoordBucket(
                block.path.clone(),
            ));
        }
        let mut records = Vec::with_capacity(block.edges);
        for encoded in bytes.chunks_exact(block.record_len) {
            records.push(decode_discontinuity_edge::<K>(encoded));
        }
        for encoded in block.buffer.chunks_exact(block.record_len) {
            records.push(decode_discontinuity_edge::<K>(encoded));
        }
        Ok(records)
    }

    fn load_column(
        &mut self,
        col: usize,
        threads: usize,
    ) -> Result<SerialEdgeMatrix<K>, SerialCollationError> {
        let mut column = SerialEdgeMatrix::new(self.vertex_partitions)
            .map_err(|_| SerialCollationError::MalformedCoordBucket(self.dir.clone()))?;
        self.load_column_into(col, threads, &mut column)?;
        Ok(column)
    }

    fn load_column_into(
        &mut self,
        col: usize,
        threads: usize,
        column: &mut SerialEdgeMatrix<K>,
    ) -> Result<(), SerialCollationError> {
        for row in 0..=col {
            self.flush_block(row, col)?;
        }
        let workers = threads.max(1).min(col + 1);
        let chunk = (col + 1).div_ceil(workers);
        let block_groups = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(workers);
            let matrix = &*self;
            for row_start in (0..=col).step_by(chunk) {
                let row_end = (row_start + chunk).min(col + 1);
                handles.push(scope.spawn(move || {
                    let mut blocks = Vec::with_capacity(row_end - row_start);
                    for row in row_start..row_end {
                        blocks.push((row, matrix.read_flushed_block(row, col)?));
                    }
                    Ok::<_, SerialCollationError>(blocks)
                }));
            }
            let mut groups = Vec::with_capacity(handles.len());
            for handle in handles {
                groups.push(
                    handle
                        .join()
                        .map_err(|_| SerialCollationError::WorkerPanic)??,
                );
            }
            Ok::<_, SerialCollationError>(groups)
        })?;
        column.stats = SerialEdgeMatrixStats::default();
        for group in block_groups {
            for (row, edges) in group {
                column.stats.edges += edges.len() as u64;
                column.stats.phi_edges += edges
                    .iter()
                    .filter(|edge| edge.first.is_phi() || edge.second.is_phi())
                    .count() as u64;
                column.stats.diagonal_edges += u64::from(row == col) * edges.len() as u64;
                let old = std::mem::replace(&mut column.blocks[row][col], edges);
                drop(old);
            }
        }
        Ok(())
    }

    fn read_flushed_row(
        &self,
        row: usize,
        threads: usize,
    ) -> Result<Vec<Vec<DiscontinuityEdge<K>>>, SerialCollationError> {
        let partition_count = self.partition_count();
        let block_count = partition_count - row;
        let workers = threads.max(1).min(block_count);
        let chunk = block_count.div_ceil(workers);
        let groups = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(workers);
            for col_start in (row..partition_count).step_by(chunk) {
                let col_end = (col_start + chunk).min(partition_count);
                handles.push(scope.spawn(move || {
                    let mut blocks = Vec::with_capacity(col_end - col_start);
                    for col in col_start..col_end {
                        blocks.push((col, self.read_flushed_block(row, col)?));
                    }
                    Ok::<_, SerialCollationError>(blocks)
                }));
            }
            let mut groups = Vec::with_capacity(handles.len());
            for handle in handles {
                groups.push(
                    handle
                        .join()
                        .map_err(|_| SerialCollationError::WorkerPanic)??,
                );
            }
            Ok::<_, SerialCollationError>(groups)
        })?;
        let mut row_blocks = (0..block_count).map(|_| Vec::new()).collect::<Vec<_>>();
        for group in groups {
            for (col, edges) in group {
                row_blocks[col - row] = edges;
            }
        }
        Ok(row_blocks)
    }
}

#[derive(Debug)]
struct ExternalBlockedContraction<const K: usize> {
    vertex_partitions: usize,
    expansion_matrix: BlockedEdgeMatrix<K>,
    pub compressed_diagonal_edges: Vec<Vec<DiscontinuityEdge<K>>>,
    meta_vertex_dir: PathBuf,
    meta_vertex_count: u64,
    meta_vertices_per_partition: Vec<usize>,
    stats: FullSerialContractionStats,
}

#[derive(Default)]
struct BlockedContractTimings {
    flush: Duration,
    clear: Duration,
    diagonal: Duration,
    read: Duration,
    scan: Duration,
    finish: Duration,
}

fn flush_blocked_edge_block(block: &mut BlockedEdgeBlock) -> Result<(), SerialCollationError> {
    if block.buffer.is_empty() {
        return Ok(());
    }
    let file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&block.path)
        .map_err(|source| SerialCollationError::Io {
            path: block.path.clone(),
            source,
        })?;
    let mut writer = BufWriter::with_capacity(BLOCKED_EDGE_WRITE_BUFFER_BYTES, file);
    writer
        .write_all(&block.buffer)
        .map_err(|source| SerialCollationError::Io {
            path: block.path.clone(),
            source,
        })?;
    writer.flush().map_err(|source| SerialCollationError::Io {
        path: block.path.clone(),
        source,
    })?;
    block.buffer.clear();
    Ok(())
}

fn encode_discontinuity_edge<const K: usize>(edge: &DiscontinuityEdge<K>) -> [u8; 72] {
    let mut bytes = [0u8; 72];
    if K <= 31 {
        let endpoint_bits = |endpoint: MatrixEndpoint<K>| match endpoint {
            MatrixEndpoint::Phi => u64::MAX,
            MatrixEndpoint::Vertex(endpoint) => endpoint.vertex.as_u128() as u64,
        };
        bytes[..8].copy_from_slice(&endpoint_bits(edge.first).to_le_bytes());
        bytes[8..16].copy_from_slice(&endpoint_bits(edge.second).to_le_bytes());
        let weight = u16::try_from(edge.weight).expect("discontinuity edge weight fits u16");
        bytes[16..18].copy_from_slice(&weight.to_le_bytes());
        let unitig_index = if edge.unitig_index == usize::MAX {
            u32::MAX
        } else {
            u32::try_from(edge.unitig_index).expect("local unitig index fits compact edge")
        };
        bytes[18..22].copy_from_slice(&unitig_index.to_le_bytes());
        bytes[22] = u8::from(
            matches!(edge.first, MatrixEndpoint::Vertex(endpoint) if endpoint.side == Side::Back),
        ) | (u8::from(
            matches!(edge.second, MatrixEndpoint::Vertex(endpoint) if endpoint.side == Side::Back),
        ) << 1)
            | (u8::from(edge.unitig_exit_side == Side::Back) << 2)
            | (u8::from(edge.swapped) << 3)
            | (u8::from(edge.phantom_unitig.is_some()) << 4);
        bytes[23..25].copy_from_slice(&edge.unitig_bucket.to_le_bytes());
        return bytes;
    }
    let kmer_bytes = discontinuity_edge_kmer_bytes::<K>();
    let endpoint_len = kmer_bytes + 2;
    encode_matrix_endpoint(edge.first, &mut bytes[0..endpoint_len], kmer_bytes);
    encode_matrix_endpoint(
        edge.second,
        &mut bytes[endpoint_len..2 * endpoint_len],
        kmer_bytes,
    );
    let weight_off = 2 * endpoint_len;
    let unitig_off = weight_off + 8;
    let flags_off = unitig_off + 8;
    bytes[weight_off..unitig_off].copy_from_slice(&edge.weight.to_le_bytes());
    bytes[unitig_off..flags_off].copy_from_slice(&(edge.unitig_index as u64).to_le_bytes());
    bytes[flags_off] = edge.unitig_exit_side as u8;
    bytes[flags_off + 1] = u8::from(edge.swapped);
    if let Some(endpoint) = edge.phantom_unitig {
        bytes[flags_off + 2] = 1;
        let phantom_off = flags_off + 3;
        bytes[phantom_off..phantom_off + kmer_bytes]
            .copy_from_slice(&endpoint.vertex.as_u128().to_le_bytes()[..kmer_bytes]);
        bytes[phantom_off + kmer_bytes] = endpoint.side as u8;
    }
    let bucket_off = discontinuity_edge_record_len::<K>() - 2;
    bytes[bucket_off..bucket_off + 2].copy_from_slice(&edge.unitig_bucket.to_le_bytes());
    bytes
}

fn encode_matrix_endpoint<const K: usize>(
    endpoint: MatrixEndpoint<K>,
    dst: &mut [u8],
    kmer_bytes: usize,
) {
    if let MatrixEndpoint::Vertex(endpoint) = endpoint {
        dst[0] = 1;
        dst[1..1 + kmer_bytes]
            .copy_from_slice(&endpoint.vertex.as_u128().to_le_bytes()[..kmer_bytes]);
        dst[1 + kmer_bytes] = endpoint.side as u8;
    }
}

fn decode_discontinuity_edge<const K: usize>(bytes: &[u8]) -> DiscontinuityEdge<K> {
    if K <= 31 {
        let read_u64 = |range: std::ops::Range<usize>| {
            u64::from_le_bytes(bytes[range].try_into().expect("eight-byte edge field"))
        };
        let flags = bytes[22];
        let endpoint = |bits: u64, side_bit: u8| {
            if bits == u64::MAX {
                MatrixEndpoint::Phi
            } else {
                MatrixEndpoint::Vertex(DiscontinuityEndpoint {
                    vertex: Kmer::from_bits(bits as u128),
                    side: if flags & side_bit == 0 {
                        Side::Front
                    } else {
                        Side::Back
                    },
                })
            }
        };
        let first = endpoint(read_u64(0..8), 1);
        let second = endpoint(read_u64(8..16), 2);
        let raw_unitig = u32::from_le_bytes(bytes[18..22].try_into().unwrap());
        let phantom_unitig = if flags & (1 << 4) == 0 {
            None
        } else {
            match (first, second) {
                (MatrixEndpoint::Vertex(endpoint), MatrixEndpoint::Phi)
                | (MatrixEndpoint::Phi, MatrixEndpoint::Vertex(endpoint)) => Some(endpoint),
                _ => None,
            }
        };
        return DiscontinuityEdge {
            first,
            second,
            weight: u64::from(u16::from_le_bytes(bytes[16..18].try_into().unwrap())),
            unitig_bucket: u16::from_le_bytes(bytes[23..25].try_into().unwrap()),
            unitig_index: if raw_unitig == u32::MAX {
                usize::MAX
            } else {
                raw_unitig as usize
            },
            unitig_exit_side: if flags & (1 << 2) == 0 {
                Side::Front
            } else {
                Side::Back
            },
            phantom_unitig,
            swapped: flags & (1 << 3) != 0,
        };
    }
    let kmer_bytes = discontinuity_edge_kmer_bytes::<K>();
    let endpoint_len = kmer_bytes + 2;
    let endpoint = |src: &[u8]| {
        if src[0] == 0 {
            MatrixEndpoint::Phi
        } else {
            let mut bits = [0u8; 16];
            bits[..kmer_bytes].copy_from_slice(&src[1..1 + kmer_bytes]);
            MatrixEndpoint::Vertex(DiscontinuityEndpoint {
                vertex: Kmer::from_bits(u128::from_le_bytes(bits)),
                side: decode_side(src[1 + kmer_bytes]),
            })
        }
    };
    let weight_off = 2 * endpoint_len;
    let unitig_off = weight_off + 8;
    let flags_off = unitig_off + 8;
    let mut weight = [0u8; 8];
    weight.copy_from_slice(&bytes[weight_off..unitig_off]);
    let mut unitig = [0u8; 8];
    unitig.copy_from_slice(&bytes[unitig_off..flags_off]);
    let phantom_unitig = (bytes[flags_off + 2] != 0).then(|| {
        let phantom_off = flags_off + 3;
        let mut bits = [0u8; 16];
        bits[..kmer_bytes].copy_from_slice(&bytes[phantom_off..phantom_off + kmer_bytes]);
        DiscontinuityEndpoint {
            vertex: Kmer::from_bits(u128::from_le_bytes(bits)),
            side: decode_side(bytes[phantom_off + kmer_bytes]),
        }
    });
    DiscontinuityEdge {
        first: endpoint(&bytes[0..endpoint_len]),
        second: endpoint(&bytes[endpoint_len..2 * endpoint_len]),
        weight: u64::from_le_bytes(weight),
        unitig_bucket: u16::from_le_bytes(
            bytes[discontinuity_edge_record_len::<K>() - 2..discontinuity_edge_record_len::<K>()]
                .try_into()
                .unwrap(),
        ),
        unitig_index: u64::from_le_bytes(unitig) as usize,
        unitig_exit_side: decode_side(bytes[flags_off]),
        phantom_unitig,
        swapped: bytes[flags_off + 1] != 0,
    }
}

#[inline]
fn decode_partition_incoming<const K: usize>(
    bytes: &[u8],
) -> Option<(Kmer<K>, PartitionOtherEnd<K>)> {
    if K <= 31 {
        let first = u64::from_le_bytes(bytes[..8].try_into().unwrap());
        let second = u64::from_le_bytes(bytes[8..16].try_into().unwrap());
        if second == u64::MAX {
            return None;
        }
        let flags = bytes[22];
        let endpoint = if first == u64::MAX {
            MatrixEndpoint::Phi
        } else {
            MatrixEndpoint::Vertex(DiscontinuityEndpoint {
                vertex: Kmer::from_bits(first as u128),
                side: if flags & 1 == 0 {
                    Side::Front
                } else {
                    Side::Back
                },
            })
        };
        return Some((
            Kmer::from_bits(second as u128),
            PartitionOtherEnd {
                endpoint,
                side_at_current: if flags & 2 == 0 {
                    Side::Front
                } else {
                    Side::Back
                },
                weight: u64::from(u16::from_le_bytes(bytes[16..18].try_into().unwrap())),
                in_same_part: false,
                processed: false,
            },
        ));
    }
    let kmer_bytes = discontinuity_edge_kmer_bytes::<K>();
    let endpoint_len = kmer_bytes + 2;
    let decode_endpoint = |src: &[u8]| {
        if src[0] == 0 {
            MatrixEndpoint::Phi
        } else {
            let mut bits = [0u8; 16];
            bits[..kmer_bytes].copy_from_slice(&src[1..1 + kmer_bytes]);
            MatrixEndpoint::Vertex(DiscontinuityEndpoint {
                vertex: Kmer::from_bits(u128::from_le_bytes(bits)),
                side: decode_side(src[1 + kmer_bytes]),
            })
        }
    };
    let lower = decode_endpoint(&bytes[..endpoint_len]);
    let MatrixEndpoint::Vertex(current) = decode_endpoint(&bytes[endpoint_len..2 * endpoint_len])
    else {
        return None;
    };
    let weight_off = 2 * endpoint_len;
    let mut weight = [0u8; 8];
    weight.copy_from_slice(&bytes[weight_off..weight_off + 8]);
    Some((
        current.vertex,
        PartitionOtherEnd {
            endpoint: lower,
            side_at_current: current.side,
            weight: u64::from_le_bytes(weight),
            in_same_part: false,
            processed: false,
        },
    ))
}

#[inline]
const fn discontinuity_edge_kmer_bytes<const K: usize>() -> usize {
    (2 * K).div_ceil(8)
}

#[inline]
const fn discontinuity_edge_record_len<const K: usize>() -> usize {
    if K <= 31 {
        25
    } else {
        3 * discontinuity_edge_kmer_bytes::<K>() + 24
    }
}

#[inline]
const fn blocked_edge_unitig_offset<const K: usize>() -> usize {
    if K <= 31 {
        18
    } else {
        2 * (discontinuity_edge_kmer_bytes::<K>() + 2) + 8
    }
}

fn add_unitig_base_to_encoded_edge<const K: usize>(
    bytes: &mut [u8; 72],
    unitig_off: usize,
    unitig_base: usize,
) {
    if K <= 31 {
        let local = u32::from_le_bytes(bytes[unitig_off..unitig_off + 4].try_into().unwrap());
        if local != u32::MAX {
            let index = unitig_base
                .checked_add(local as usize)
                .and_then(|index| u32::try_from(index).ok())
                .expect("global local-unitig index fits compact edge");
            bytes[unitig_off..unitig_off + 4].copy_from_slice(&index.to_le_bytes());
        }
    } else {
        let local = u64::from_le_bytes(bytes[unitig_off..unitig_off + 8].try_into().unwrap());
        if local != u64::MAX {
            let index = unitig_base + local as usize;
            bytes[unitig_off..unitig_off + 8].copy_from_slice(&(index as u64).to_le_bytes());
        }
    }
}

#[inline]
fn set_encoded_edge_unitig_bucket<const K: usize>(bytes: &mut [u8; 72], bucket: u16) {
    let offset = discontinuity_edge_record_len::<K>() - 2;
    bytes[offset..offset + 2].copy_from_slice(&bucket.to_le_bytes());
}

#[inline]
fn decode_side(value: u8) -> Side {
    if value == Side::Front as u8 {
        Side::Front
    } else {
        Side::Back
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SerialContraction {
    pub components: Vec<SerialComponent>,
    pub stats: SerialContractionStats,
}

pub const DISCONTINUITY_PARALLELIZATION_OPPORTUNITIES: &[&str] = &[
    "diagonal block compression per partition can run concurrently with non-diagonal column scans",
    "non-diagonal blocks in a partition column can be scanned independently with a shared vertex table",
    "compressed diagonal-chain edges can be processed independently after the column scan",
    "false-phantom filtering is a parallel scan over unprocessed partition-table entries",
    "different output edge/meta-vertex buckets can use batched thread-local buffers before external-memory flush",
];

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiagonalCompression<const K: usize> {
    pub partition: usize,
    pub edges: Vec<DiscontinuityEdge<K>>,
    pub expansion_edges: Vec<DiscontinuityEdge<K>>,
    pub meta_vertices: Vec<SerialMetaVertex<K>>,
    pub stats: DiagonalCompressionStats,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SerialMetaVertex<const K: usize> {
    pub vertex: Kmer<K>,
    pub partition: usize,
    pub entry_side: Side,
    pub weight: u64,
    pub is_cycle: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct DiagonalCompressionStats {
    pub input_edges: u64,
    pub compressed_edges: u64,
    pub meta_vertices: u64,
    pub isolated_cordless_cycles: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PartitionContraction<const K: usize> {
    pub partition: usize,
    pub edges: Vec<DiscontinuityEdge<K>>,
    pub expansion_edges: Vec<DiscontinuityEdge<K>>,
    pub compressed_diagonal_edges: Vec<DiscontinuityEdge<K>>,
    pub meta_vertices: Vec<SerialMetaVertex<K>>,
    pub stats: PartitionContractionStats,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PartitionContractionStats {
    pub input_non_diagonal_edges: u64,
    pub compressed_diagonal_edges: u64,
    pub output_edges: u64,
    pub meta_vertices: u64,
    pub phantom_edges: u64,
    pub isolated_cordless_cycles: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FullSerialDiscontinuityContraction<const K: usize> {
    pub vertex_partitions: usize,
    pub final_edges: Vec<DiscontinuityEdge<K>>,
    pub expansion_edges: Vec<DiscontinuityEdge<K>>,
    pub expansion_matrix: SerialEdgeMatrix<K>,
    pub compressed_diagonal_edges: Vec<Vec<DiscontinuityEdge<K>>>,
    pub meta_vertices: Vec<SerialMetaVertex<K>>,
    pub partitions: Vec<PartitionContraction<K>>,
    pub stats: FullSerialContractionStats,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct FullSerialContractionStats {
    pub partitions: u64,
    pub input_edges: u64,
    pub partition_output_edges: u64,
    pub reinserted_edges: u64,
    pub final_edges: u64,
    pub meta_vertices: u64,
    pub phantom_edges: u64,
    pub isolated_cordless_cycles: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PathInfo<const K: usize> {
    pub path_id: Kmer<K>,
    pub rank: u64,
    pub exit_side: Side,
    pub is_cycle: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VertexPathInfo<const K: usize> {
    pub vertex: Kmer<K>,
    pub info: PathInfo<K>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EdgePathInfo<const K: usize> {
    pub unitig_index: usize,
    pub phantom_unitig: Option<DiscontinuityEndpoint<K>>,
    pub info: PathInfo<K>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SerialExpansion<const K: usize> {
    pub vertices: Vec<VertexPathInfo<K>>,
    pub edges: Vec<EdgePathInfo<K>>,
    pub stats: SerialExpansionStats,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SerialExpansionStats {
    pub seed_vertices: u64,
    pub inferred_vertices: u64,
    pub edge_path_infos: u64,
    pub unresolved_edges: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SerialCollation {
    pub unitigs: Vec<Vec<u8>>,
    pub stats: SerialCollationStats,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SerialCollationStats {
    pub input_path_infos: u64,
    pub emitted_unitigs: u64,
    pub emitted_bases: u64,
    pub missing_unitig_labels: u64,
    pub direct_local_unitigs: u64,
    pub stitched_discontinuity_unitigs: u64,
}

#[derive(Debug)]
pub enum SerialCollationError {
    Io {
        path: std::path::PathBuf,
        source: std::io::Error,
    },
    MalformedCoordBucket(std::path::PathBuf),
    WorkerPanic,
    Color(ColorError),
}

impl std::fmt::Display for SerialCollationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io { path, source } => write!(f, "{}: {source}", path.display()),
            Self::MalformedCoordBucket(path) => {
                write!(
                    f,
                    "malformed stitched coordinate bucket: {}",
                    path.display()
                )
            }
            Self::WorkerPanic => write!(f, "stitched coordinate worker thread panicked"),
            Self::Color(err) => write!(f, "{err}"),
        }
    }
}

impl std::error::Error for SerialCollationError {}

impl From<ColorError> for SerialCollationError {
    fn from(value: ColorError) -> Self {
        Self::Color(value)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SerialComponent {
    pub endpoints: usize,
    pub edges: u64,
    pub weight: u64,
    pub phi_edges: u64,
    pub cyclic: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SerialContractionStats {
    pub input_edges: u64,
    pub components: u64,
    pub phi_edges: u64,
    pub cyclic_components: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct EndpointKey<const K: usize>(MatrixEndpoint<K>);

pub struct SerialDiscontinuityContractor;

struct PartitionColumnScan<const K: usize> {
    table: FastHashMap<Kmer<K>, PartitionOtherEnd<K>>,
    edges: Vec<DiscontinuityEdge<K>>,
    meta_vertices: Vec<SerialMetaVertex<K>>,
    input_non_diagonal_edges: u64,
}

struct PartitionColumnScanOutput<const K: usize> {
    edges: Vec<DiscontinuityEdge<K>>,
    meta_vertices: Vec<SerialMetaVertex<K>>,
    input_non_diagonal_edges: u64,
}

#[derive(Debug, Clone, Copy)]
struct PartitionColumnIncoming<const K: usize> {
    vertex: Kmer<K>,
    end: PartitionOtherEnd<K>,
}

fn partition_column_vertex_shard<const K: usize>(vertex: Kmer<K>, shard_mask: usize) -> usize {
    let bits = vertex.as_u128();
    let mixed = hash_u64((bits as u64) ^ ((bits >> 64) as u64), 0);
    (mixed as usize) & shard_mask
}

impl SerialDiscontinuityContractor {
    pub fn compress_diagonal_block<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
    ) -> DiagonalCompression<K> {
        let mut ends = FastHashMap::with_hasher(FastBuildHasher::default());
        Self::compress_diagonal_block_with_ends(matrix, partition, &mut ends)
    }

    fn compress_diagonal_block_with_ends<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
        ends: &mut FastHashMap<Kmer<K>, DiagonalOtherEnd<K>>,
    ) -> DiagonalCompression<K> {
        assert!(partition > 0 && partition < matrix.partition_count());

        let diagonal_edges = matrix.block(partition, partition);
        ends.clear();
        ends.reserve(diagonal_edges.len().saturating_mul(2));
        let mut expansion_edges = Vec::with_capacity(diagonal_edges.len());
        let mut meta_vertices = Vec::new();
        let mut isolated_cordless_cycles = 0u64;

        for edge in diagonal_edges {
            let (x, y) = match (edge.first, edge.second) {
                (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) => (x, y),
                _ => continue,
            };

            let end_x = ends.get(&x.vertex).copied();
            let end_y = ends.get(&y.vertex).copied();

            let u = end_x
                .map(|end| DiscontinuityEndpoint {
                    vertex: end.vertex,
                    side: end.side_at_vertex,
                })
                .unwrap_or(x);
            let v = end_y
                .map(|end| DiscontinuityEndpoint {
                    vertex: end.vertex,
                    side: end.side_at_vertex,
                })
                .unwrap_or(y);
            let weight =
                end_x.map_or(0, |end| end.weight) + edge.weight + end_y.map_or(0, |end| end.weight);

            assert_eq!(matrix.partition(MatrixEndpoint::Vertex(u)), partition);
            assert_eq!(matrix.partition(MatrixEndpoint::Vertex(v)), partition);

            if u.vertex == v.vertex {
                ends.insert(
                    u.vertex,
                    DiagonalOtherEnd {
                        vertex: Kmer::zero(),
                        side_at_vertex: Side::Back,
                        side_at_current: Side::Front,
                        weight: 1,
                        unitig_index: 0,
                        unitig_exit_side: Side::Back,
                        is_phi: true,
                    },
                );
                meta_vertices.push(SerialMetaVertex {
                    vertex: u.vertex,
                    partition,
                    entry_side: Side::Front.inverse(),
                    weight: 1,
                    is_cycle: true,
                });
                isolated_cordless_cycles += 1;
            } else if u.vertex == y.vertex && v.vertex == x.vertex {
                ends.insert(
                    u.vertex,
                    DiagonalOtherEnd {
                        vertex: Kmer::zero(),
                        side_at_vertex: Side::Back,
                        side_at_current: u.side.inverse(),
                        weight: 1,
                        unitig_index: 0,
                        unitig_exit_side: Side::Back,
                        is_phi: true,
                    },
                );
                ends.insert(
                    v.vertex,
                    DiagonalOtherEnd {
                        vertex: Kmer::zero(),
                        side_at_vertex: Side::Back,
                        side_at_current: v.side.inverse(),
                        weight: 1,
                        unitig_index: 0,
                        unitig_exit_side: Side::Back,
                        is_phi: true,
                    },
                );
                meta_vertices.push(SerialMetaVertex {
                    vertex: u.vertex,
                    partition,
                    entry_side: u.side,
                    weight: 1,
                    is_cycle: true,
                });
                isolated_cordless_cycles += 1;
            } else {
                expansion_edges.push(DiscontinuityEdge {
                    first: MatrixEndpoint::Vertex(DiscontinuityEndpoint {
                        vertex: u.vertex,
                        side: u.side,
                    }),
                    second: MatrixEndpoint::Vertex(DiscontinuityEndpoint {
                        vertex: v.vertex,
                        side: v.side,
                    }),
                    weight,
                    unitig_bucket: if weight == 1 { edge.unitig_bucket } else { 0 },
                    unitig_index: if weight == 1 { edge.unitig_index } else { 0 },
                    unitig_exit_side: if weight == 1 {
                        edge.unitig_exit_side
                    } else {
                        Side::Back
                    },
                    phantom_unitig: None,
                    swapped: false,
                });
                ends.insert(
                    u.vertex,
                    DiagonalOtherEnd {
                        vertex: v.vertex,
                        side_at_vertex: v.side,
                        side_at_current: u.side,
                        weight,
                        unitig_index: 0,
                        unitig_exit_side: Side::Back,
                        is_phi: false,
                    },
                );
                ends.insert(
                    v.vertex,
                    DiagonalOtherEnd {
                        vertex: u.vertex,
                        side_at_vertex: u.side,
                        side_at_current: v.side,
                        weight,
                        unitig_index: 0,
                        unitig_exit_side: Side::Back,
                        is_phi: false,
                    },
                );
            }
        }

        let mut compressed_edges = Vec::new();
        for (&vertex, &end) in ends.iter() {
            if end.is_phi {
                continue;
            }

            if vertex < end.vertex {
                let Some(reverse) = ends.get(&end.vertex) else {
                    continue;
                };
                if reverse.vertex == vertex {
                    compressed_edges.push(DiscontinuityEdge {
                        first: MatrixEndpoint::Vertex(DiscontinuityEndpoint {
                            vertex,
                            side: end.side_at_current,
                        }),
                        second: MatrixEndpoint::Vertex(DiscontinuityEndpoint {
                            vertex: end.vertex,
                            side: end.side_at_vertex,
                        }),
                        weight: end.weight,
                        unitig_bucket: 0,
                        unitig_index: end.unitig_index,
                        unitig_exit_side: end.unitig_exit_side,
                        phantom_unitig: None,
                        swapped: false,
                    });
                }
            }
        }
        DiagonalCompression {
            partition,
            stats: DiagonalCompressionStats {
                input_edges: diagonal_edges.len() as u64,
                compressed_edges: compressed_edges.len() as u64,
                meta_vertices: meta_vertices.len() as u64,
                isolated_cordless_cycles,
            },
            edges: compressed_edges,
            expansion_edges,
            meta_vertices,
        }
    }

    pub fn contract_partition<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
    ) -> PartitionContraction<K> {
        Self::contract_partition_with_threads(matrix, partition, 1)
    }

    pub fn contract_partition_with_threads<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
        threads: usize,
    ) -> PartitionContraction<K> {
        assert!(partition > 0 && partition < matrix.partition_count());

        let (diagonal, scan) = if threads > 1 {
            std::thread::scope(|scope| {
                let diagonal = scope.spawn(|| Self::compress_diagonal_block(matrix, partition));
                let scan = Self::scan_partition_column_with_threads(matrix, partition, threads);
                (
                    diagonal
                        .join()
                        .expect("diagonal compression worker panicked"),
                    scan,
                )
            })
        } else {
            (
                Self::compress_diagonal_block(matrix, partition),
                Self::scan_partition_column_with_threads(matrix, partition, threads),
            )
        };

        let mut table = scan.table;
        let mut edges = scan.edges;
        let mut expansion_edges = Vec::new();
        let mut meta_vertices = scan.meta_vertices;
        meta_vertices.extend(diagonal.meta_vertices);
        let mut phantom_edges = 0u64;

        // Parallelization note: compressed diagonal-chain edges are independent
        // after the non-diagonal table is built. The serial loop makes false
        // phantom creation explicit for tests.
        for edge in &diagonal.edges {
            let (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) = (edge.first, edge.second)
            else {
                continue;
            };

            if !table.contains_key(&x.vertex) {
                phantom_edges += 1;
                let phantom = DiscontinuityEndpoint {
                    vertex: x.vertex,
                    side: x.side.inverse(),
                };
                expansion_edges.push(join_other_ends_with_phantom(
                    MatrixEndpoint::Phi,
                    MatrixEndpoint::Vertex(phantom),
                    1,
                    Some(phantom),
                ));
                table.insert(
                    x.vertex,
                    PartitionOtherEnd {
                        endpoint: MatrixEndpoint::Phi,
                        side_at_current: x.side.inverse(),
                        weight: 1,
                        in_same_part: false,
                        processed: false,
                    },
                );
            }

            if !table.contains_key(&y.vertex) {
                phantom_edges += 1;
                let phantom = DiscontinuityEndpoint {
                    vertex: y.vertex,
                    side: y.side.inverse(),
                };
                expansion_edges.push(join_other_ends_with_phantom(
                    MatrixEndpoint::Phi,
                    MatrixEndpoint::Vertex(phantom),
                    1,
                    Some(phantom),
                ));
                table.insert(
                    y.vertex,
                    PartitionOtherEnd {
                        endpoint: MatrixEndpoint::Phi,
                        side_at_current: y.side.inverse(),
                        weight: 1,
                        in_same_part: false,
                        processed: false,
                    },
                );
            }

            let x_end = table.get(&x.vertex).copied().unwrap();
            let y_end = table.get(&y.vertex).copied().unwrap();
            if x_end.endpoint.is_phi() && y_end.endpoint.is_phi() {
                meta_vertices.push(two_weight_meta_vertex(
                    x.vertex,
                    partition,
                    x.side.inverse(),
                    x_end.weight,
                    edge.weight + y_end.weight,
                    false,
                ));
            } else {
                edges.push(join_other_ends(
                    x_end.endpoint,
                    y_end.endpoint,
                    x_end.weight + edge.weight + y_end.weight,
                ));
            }

            if let Some(end) = table.get_mut(&x.vertex) {
                end.processed = true;
            }
            if let Some(end) = table.get_mut(&y.vertex) {
                end.processed = true;
            }
        }

        // Parallelization note: this is a scan over table entries in C++.
        let unprocessed = table
            .iter()
            .filter_map(|(&vertex, &end)| (!end.processed).then_some((vertex, end)))
            .collect::<Vec<_>>();
        for (vertex, end) in unprocessed {
            phantom_edges += 1;
            let phantom = DiscontinuityEndpoint {
                vertex,
                side: end.side_at_current.inverse(),
            };
            expansion_edges.push(join_other_ends_with_phantom(
                MatrixEndpoint::Phi,
                MatrixEndpoint::Vertex(phantom),
                1,
                Some(phantom),
            ));

            if end.endpoint.is_phi() {
                meta_vertices.push(two_weight_meta_vertex(
                    vertex,
                    partition,
                    end.side_at_current,
                    end.weight,
                    1,
                    false,
                ));
            } else {
                edges.push(join_other_ends(
                    MatrixEndpoint::Phi,
                    end.endpoint,
                    1 + end.weight,
                ));
            }
        }

        edges.sort_by_key(|edge| {
            (
                endpoint_sort_key(edge.first),
                endpoint_sort_key(edge.second),
                edge.weight,
            )
        });
        meta_vertices.sort_by_key(|meta| {
            (
                meta.partition,
                meta.vertex.as_u128(),
                meta.entry_side as u8,
                meta.weight,
                meta.is_cycle,
            )
        });

        PartitionContraction {
            partition,
            compressed_diagonal_edges: diagonal.expansion_edges,
            expansion_edges,
            stats: PartitionContractionStats {
                input_non_diagonal_edges: scan.input_non_diagonal_edges,
                compressed_diagonal_edges: diagonal.stats.compressed_edges,
                output_edges: edges.len() as u64,
                meta_vertices: meta_vertices.len() as u64,
                phantom_edges,
                isolated_cordless_cycles: diagonal.stats.isolated_cordless_cycles,
            },
            edges,
            meta_vertices,
        }
    }

    fn scan_partition_column<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
    ) -> PartitionColumnScan<K> {
        let column_edges = (0..partition)
            .map(|row| matrix.block(row, partition).len())
            .sum::<usize>();
        let mut table = FastHashMap::<Kmer<K>, PartitionOtherEnd<K>>::with_capacity_and_hasher(
            column_edges,
            FastBuildHasher::default(),
        );
        let mut edges = Vec::with_capacity(column_edges / 2);
        let mut meta_vertices = Vec::new();
        let mut input_non_diagonal_edges = 0u64;

        // Parallelization note: the C++ implementation scans blocks in this
        // column concurrently against a shared table. This scan is now isolated
        // from diagonal compression and can be split by row range next.
        for row in 0..partition {
            for edge in matrix.block(row, partition) {
                let (lower, current) = match endpoint_in_partition(matrix, edge, partition) {
                    Some(pair) => pair,
                    None => continue,
                };
                input_non_diagonal_edges += 1;

                let incoming = PartitionOtherEnd {
                    endpoint: lower,
                    side_at_current: current.side,
                    weight: edge.weight,
                    in_same_part: false,
                    processed: false,
                };

                absorb_partition_other_end(
                    current.vertex,
                    incoming,
                    partition,
                    &mut table,
                    &mut edges,
                    &mut meta_vertices,
                );
            }
        }

        PartitionColumnScan {
            table,
            edges,
            meta_vertices,
            input_non_diagonal_edges,
        }
    }

    fn scan_partition_column_into<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
        table: &mut FastHashMap<Kmer<K>, PartitionOtherEnd<K>>,
    ) -> PartitionColumnScanOutput<K> {
        table.clear();
        let column_edges = (0..partition)
            .map(|row| matrix.block(row, partition).len())
            .sum::<usize>();
        table.reserve(column_edges.saturating_sub(table.capacity()));
        let mut edges = Vec::with_capacity(column_edges / 2);
        let mut meta_vertices = Vec::new();
        let mut input_non_diagonal_edges = 0u64;

        for row in 0..partition {
            for edge in matrix.block(row, partition) {
                let (lower, current) = match endpoint_in_partition(matrix, edge, partition) {
                    Some(pair) => pair,
                    None => continue,
                };
                input_non_diagonal_edges += 1;

                let incoming = PartitionOtherEnd {
                    endpoint: lower,
                    side_at_current: current.side,
                    weight: edge.weight,
                    in_same_part: false,
                    processed: false,
                };

                absorb_partition_other_end(
                    current.vertex,
                    incoming,
                    partition,
                    table,
                    &mut edges,
                    &mut meta_vertices,
                );
            }
        }

        PartitionColumnScanOutput {
            edges,
            meta_vertices,
            input_non_diagonal_edges,
        }
    }

    fn contract_partition_with_reusable_table<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
        threads: usize,
        table: &mut FastHashMap<Kmer<K>, PartitionOtherEnd<K>>,
    ) -> PartitionContraction<K> {
        assert!(partition > 0 && partition < matrix.partition_count());

        let (diagonal, scan) = if threads > 1 {
            std::thread::scope(|scope| {
                let diagonal = scope.spawn(|| Self::compress_diagonal_block(matrix, partition));
                let scan = Self::scan_partition_column_into(matrix, partition, table);
                (
                    diagonal
                        .join()
                        .expect("diagonal compression worker panicked"),
                    scan,
                )
            })
        } else {
            (
                Self::compress_diagonal_block(matrix, partition),
                Self::scan_partition_column_into(matrix, partition, table),
            )
        };

        let mut edges = scan.edges;
        let mut expansion_edges = Vec::new();
        let mut meta_vertices = scan.meta_vertices;
        meta_vertices.extend(diagonal.meta_vertices);
        let mut phantom_edges = 0u64;

        for edge in &diagonal.edges {
            let (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) = (edge.first, edge.second)
            else {
                continue;
            };

            if !table.contains_key(&x.vertex) {
                phantom_edges += 1;
                let phantom = DiscontinuityEndpoint {
                    vertex: x.vertex,
                    side: x.side.inverse(),
                };
                expansion_edges.push(join_other_ends_with_phantom(
                    MatrixEndpoint::Phi,
                    MatrixEndpoint::Vertex(phantom),
                    1,
                    Some(phantom),
                ));
                table.insert(
                    x.vertex,
                    PartitionOtherEnd {
                        endpoint: MatrixEndpoint::Phi,
                        side_at_current: x.side.inverse(),
                        weight: 1,
                        in_same_part: false,
                        processed: false,
                    },
                );
            }

            if !table.contains_key(&y.vertex) {
                phantom_edges += 1;
                let phantom = DiscontinuityEndpoint {
                    vertex: y.vertex,
                    side: y.side.inverse(),
                };
                expansion_edges.push(join_other_ends_with_phantom(
                    MatrixEndpoint::Phi,
                    MatrixEndpoint::Vertex(phantom),
                    1,
                    Some(phantom),
                ));
                table.insert(
                    y.vertex,
                    PartitionOtherEnd {
                        endpoint: MatrixEndpoint::Phi,
                        side_at_current: y.side.inverse(),
                        weight: 1,
                        in_same_part: false,
                        processed: false,
                    },
                );
            }

            let x_end = table.get(&x.vertex).copied().unwrap();
            let y_end = table.get(&y.vertex).copied().unwrap();
            if x_end.endpoint.is_phi() && y_end.endpoint.is_phi() {
                meta_vertices.push(two_weight_meta_vertex(
                    x.vertex,
                    partition,
                    x.side.inverse(),
                    x_end.weight,
                    edge.weight + y_end.weight,
                    false,
                ));
            } else {
                edges.push(join_other_ends(
                    x_end.endpoint,
                    y_end.endpoint,
                    x_end.weight + edge.weight + y_end.weight,
                ));
            }

            if let Some(end) = table.get_mut(&x.vertex) {
                end.processed = true;
            }
            if let Some(end) = table.get_mut(&y.vertex) {
                end.processed = true;
            }
        }

        for (&vertex, &end) in table.iter() {
            if end.processed {
                continue;
            }
            phantom_edges += 1;
            let phantom = DiscontinuityEndpoint {
                vertex,
                side: end.side_at_current.inverse(),
            };
            expansion_edges.push(join_other_ends_with_phantom(
                MatrixEndpoint::Phi,
                MatrixEndpoint::Vertex(phantom),
                1,
                Some(phantom),
            ));

            if end.endpoint.is_phi() {
                meta_vertices.push(two_weight_meta_vertex(
                    vertex,
                    partition,
                    end.side_at_current,
                    end.weight,
                    1,
                    false,
                ));
            } else {
                edges.push(join_other_ends(
                    MatrixEndpoint::Phi,
                    end.endpoint,
                    1 + end.weight,
                ));
            }
        }

        PartitionContraction {
            partition,
            compressed_diagonal_edges: diagonal.expansion_edges,
            expansion_edges,
            stats: PartitionContractionStats {
                input_non_diagonal_edges: scan.input_non_diagonal_edges,
                compressed_diagonal_edges: diagonal.stats.compressed_edges,
                output_edges: edges.len() as u64,
                meta_vertices: meta_vertices.len() as u64,
                phantom_edges,
                isolated_cordless_cycles: diagonal.stats.isolated_cordless_cycles,
            },
            edges,
            meta_vertices,
        }
    }

    fn scan_partition_column_atomic<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
        threads: usize,
        atomic_table: &AtomicPartitionTable<K>,
        pool: &ThreadPool,
    ) -> PartitionColumnScanOutput<K> {
        atomic_table.clear(threads, pool);
        const EDGES_PER_TASK: usize = 16 * 1024;
        let tasks = (0..partition)
            .flat_map(|row| matrix.block(row, partition).chunks(EDGES_PER_TASK))
            .collect::<Vec<_>>();
        let outputs = pool.install(|| {
            tasks
                .into_par_iter()
                .map(|input| {
                    let mut edges = Vec::with_capacity(input.len() / 2);
                    let mut meta_vertices = Vec::new();
                    let mut input_edges = 0u64;
                    for edge in input {
                        let Some((lower, current)) = endpoint_in_partition(matrix, edge, partition)
                        else {
                            continue;
                        };
                        input_edges += 1;
                        atomic_table.absorb(
                            current.vertex,
                            PartitionOtherEnd {
                                endpoint: lower,
                                side_at_current: current.side,
                                weight: edge.weight,
                                in_same_part: false,
                                processed: false,
                            },
                            partition,
                            &mut edges,
                            &mut meta_vertices,
                        );
                    }
                    (edges, meta_vertices, input_edges)
                })
                .collect::<Vec<_>>()
        });
        let mut edges = Vec::new();
        let mut meta_vertices = Vec::new();
        let mut input_non_diagonal_edges = 0;
        for (worker_edges, worker_meta, worker_input) in outputs {
            edges.extend(worker_edges);
            meta_vertices.extend(worker_meta);
            input_non_diagonal_edges += worker_input;
        }
        PartitionColumnScanOutput {
            edges,
            meta_vertices,
            input_non_diagonal_edges,
        }
    }

    fn contract_blocked_partition_atomic<const K: usize>(
        matrix: &mut BlockedEdgeMatrix<K>,
        appenders: &ConcurrentBlockedEdgeWriters,
        diagonal_matrix: &mut SerialEdgeMatrix<K>,
        partition: usize,
        threads: usize,
        table: &AtomicPartitionTable<K>,
        diagonal_ends: &mut FastHashMap<Kmer<K>, DiagonalOtherEnd<K>>,
        pool: &ThreadPool,
        timings: &mut BlockedContractTimings,
    ) -> Result<(PartitionContraction<K>, u64), SerialCollationError> {
        let phase = Instant::now();
        appenders.flush_column(partition)?;
        for row in 0..=partition {
            matrix.flush_block(row, partition)?;
            let added = appenders.take_added_edges(row, partition);
            if added != 0 {
                let block_id = matrix.block_index(row, partition);
                matrix.blocks[block_id].edges += added;
                matrix.stats.edges += added as u64;
                matrix.stats.phi_edges += u64::from(row == 0) * added as u64;
                matrix.stats.diagonal_edges += u64::from(row == partition) * added as u64;
            }
        }
        timings.flush += phase.elapsed();
        let phase = Instant::now();
        table.clear(threads, pool);
        timings.clear += phase.elapsed();

        let diagonal_read_started = Instant::now();
        let diagonal_edges = matrix.read_flushed_block(partition, partition)?;
        diagonal_matrix.blocks[partition][partition] = diagonal_edges;
        let diagonal_read_elapsed = diagonal_read_started.elapsed();

        let record_len = discontinuity_edge_record_len::<K>();
        let block_ids = (0..partition)
            .filter_map(|row| {
                let block_id = matrix.block_index(row, partition);
                (matrix.blocks[block_id].edges != 0).then_some(block_id)
            })
            .collect::<Vec<_>>();
        let records_per_task = (1024 * 1024 / record_len).max(1);
        let bytes_per_task = records_per_task * record_len;
        let mut block_files = Vec::with_capacity(block_ids.len());
        for &block_id in &block_ids {
            let block = &matrix.blocks[block_id];
            let expected_file_bytes = block
                .edges
                .checked_mul(record_len)
                .and_then(|bytes| bytes.checked_sub(block.buffer.len()))
                .ok_or_else(|| SerialCollationError::MalformedCoordBucket(block.path.clone()))?;
            if expected_file_bytes == 0 {
                block_files.push((block_id, None, 0));
                continue;
            }
            let file = File::open(&block.path).map_err(|source| SerialCollationError::Io {
                path: block.path.clone(),
                source,
            })?;
            let file_bytes = file
                .metadata()
                .map_err(|source| SerialCollationError::Io {
                    path: block.path.clone(),
                    source,
                })?
                .len() as usize;
            if file_bytes != expected_file_bytes || file_bytes % record_len != 0 {
                return Err(SerialCollationError::MalformedCoordBucket(
                    block.path.clone(),
                ));
            }
            block_files.push((block_id, Some(file), file_bytes));
        }
        let mut tasks = Vec::new();
        for (block_id, file, file_bytes) in &block_files {
            if let Some(file) = file {
                let path = matrix.blocks[*block_id].path.as_path();
                for offset in (0..*file_bytes).step_by(bytes_per_task) {
                    tasks.push(BlockedReadTask::File {
                        file,
                        path,
                        offset: offset as u64,
                        len: bytes_per_task.min(*file_bytes - offset),
                    });
                }
            }
            tasks.extend(
                matrix.blocks[*block_id]
                    .buffer
                    .chunks(bytes_per_task)
                    .map(BlockedReadTask::Memory),
            );
        }
        let ((outputs, scan_elapsed), (diagonal, diagonal_elapsed)) = pool.install(|| {
            rayon::join(
                || {
                    let scan_started = Instant::now();
                    let outputs = tasks
                        .into_par_iter()
                        .map_init(
                            || vec![0u8; bytes_per_task],
                            |scratch, task| {
                                let bytes = match task {
                                    BlockedReadTask::File {
                                        file,
                                        path,
                                        offset,
                                        len,
                                    } => {
                                        file.read_exact_at(&mut scratch[..len], offset).map_err(
                                            |source| SerialCollationError::Io {
                                                path: path.to_path_buf(),
                                                source,
                                            },
                                        )?;
                                        &scratch[..len]
                                    }
                                    BlockedReadTask::Memory(bytes) => bytes,
                                };
                                let mut prepared = Vec::with_capacity(bytes.len() / record_len / 2);
                                let mut meta_vertices = Vec::new();
                                let mut input_edges = 0u64;
                                for encoded in bytes.chunks_exact(record_len) {
                                    let Some((vertex, incoming)) =
                                        decode_partition_incoming::<K>(encoded)
                                    else {
                                        continue;
                                    };
                                    input_edges += 1;
                                    table.absorb_prepared(
                                        vertex,
                                        incoming,
                                        partition,
                                        matrix.vertex_partitions,
                                        &mut prepared,
                                        &mut meta_vertices,
                                    );
                                }
                                let emitted = emit_prepared_edge_batch(&mut prepared, appenders)?;
                                Ok::<_, SerialCollationError>((meta_vertices, input_edges, emitted))
                            },
                        )
                        .collect::<Result<Vec<_>, SerialCollationError>>();
                    (outputs, scan_started.elapsed())
                },
                || {
                    let diagonal_started = Instant::now();
                    let diagonal = Self::compress_diagonal_block_with_ends(
                        diagonal_matrix,
                        partition,
                        diagonal_ends,
                    );
                    (diagonal, diagonal_started.elapsed())
                },
            )
        });
        let outputs = outputs?;
        diagonal_matrix.blocks[partition][partition].clear();
        timings.scan += scan_elapsed;
        timings.diagonal += diagonal_read_elapsed + diagonal_elapsed;
        let mut meta_vertices = Vec::new();
        let mut input_non_diagonal_edges = 0;
        let mut pre_reinserted_edges = 0u64;
        for (local_meta, local_input, local_emitted) in outputs {
            meta_vertices.extend(local_meta);
            input_non_diagonal_edges += local_input;
            pre_reinserted_edges += local_emitted;
        }
        let scan = PartitionColumnScanOutput {
            edges: Vec::new(),
            meta_vertices,
            input_non_diagonal_edges,
        };
        let phase = Instant::now();
        let contracted =
            Self::finish_atomic_partition(partition, threads, table, pool, diagonal, scan);
        timings.finish += phase.elapsed();
        let mut contracted = contracted;
        contracted.stats.output_edges += pre_reinserted_edges;
        Ok((contracted, pre_reinserted_edges))
    }

    fn contract_partition_with_owned_tables<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
        threads: usize,
        tables: &mut OwnedPartitionTables<K>,
    ) -> PartitionContraction<K> {
        let (diagonal, scan) = if threads > 1 {
            std::thread::scope(|scope| {
                let diagonal = scope.spawn(|| Self::compress_diagonal_block(matrix, partition));
                let scan = Self::scan_partition_column_owned(matrix, partition, threads, tables);
                (
                    diagonal
                        .join()
                        .expect("diagonal compression worker panicked"),
                    scan,
                )
            })
        } else {
            let diagonal = Self::compress_diagonal_block(matrix, partition);
            let scan = Self::scan_partition_column_owned(matrix, partition, 1, tables);
            (diagonal, scan)
        };

        let mut edges = scan.edges;
        let mut expansion_edges = Vec::new();
        let mut meta_vertices = scan.meta_vertices;
        meta_vertices.extend(diagonal.meta_vertices);
        let mut phantom_edges = 0u64;

        for edge in &diagonal.edges {
            let (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) = (edge.first, edge.second)
            else {
                continue;
            };
            for endpoint in [x, y] {
                if !tables.contains_key(&endpoint.vertex) {
                    phantom_edges += 1;
                    let phantom = DiscontinuityEndpoint {
                        vertex: endpoint.vertex,
                        side: endpoint.side.inverse(),
                    };
                    expansion_edges.push(join_other_ends_with_phantom(
                        MatrixEndpoint::Phi,
                        MatrixEndpoint::Vertex(phantom),
                        1,
                        Some(phantom),
                    ));
                    tables.insert(
                        endpoint.vertex,
                        PartitionOtherEnd {
                            endpoint: MatrixEndpoint::Phi,
                            side_at_current: endpoint.side.inverse(),
                            weight: 1,
                            in_same_part: false,
                            processed: false,
                        },
                    );
                }
            }

            let x_end = *tables.get(&x.vertex).expect("inserted missing x endpoint");
            let y_end = *tables.get(&y.vertex).expect("inserted missing y endpoint");
            if x_end.endpoint.is_phi() && y_end.endpoint.is_phi() {
                meta_vertices.push(two_weight_meta_vertex(
                    x.vertex,
                    partition,
                    x.side.inverse(),
                    x_end.weight,
                    edge.weight + y_end.weight,
                    false,
                ));
            } else {
                edges.push(join_other_ends(
                    x_end.endpoint,
                    y_end.endpoint,
                    x_end.weight + edge.weight + y_end.weight,
                ));
            }
            tables
                .get_mut(&x.vertex)
                .expect("x endpoint exists")
                .processed = true;
            tables
                .get_mut(&y.vertex)
                .expect("y endpoint exists")
                .processed = true;
        }

        for map in &tables.maps {
            for (&vertex, &end) in map {
                if end.processed {
                    continue;
                }
                phantom_edges += 1;
                let phantom = DiscontinuityEndpoint {
                    vertex,
                    side: end.side_at_current.inverse(),
                };
                expansion_edges.push(join_other_ends_with_phantom(
                    MatrixEndpoint::Phi,
                    MatrixEndpoint::Vertex(phantom),
                    1,
                    Some(phantom),
                ));
                if end.endpoint.is_phi() {
                    meta_vertices.push(two_weight_meta_vertex(
                        vertex,
                        partition,
                        end.side_at_current,
                        end.weight,
                        1,
                        false,
                    ));
                } else {
                    edges.push(join_other_ends(
                        MatrixEndpoint::Phi,
                        end.endpoint,
                        1 + end.weight,
                    ));
                }
            }
        }

        PartitionContraction {
            partition,
            compressed_diagonal_edges: diagonal.expansion_edges,
            expansion_edges,
            stats: PartitionContractionStats {
                input_non_diagonal_edges: scan.input_non_diagonal_edges,
                compressed_diagonal_edges: diagonal.stats.compressed_edges,
                output_edges: edges.len() as u64,
                meta_vertices: meta_vertices.len() as u64,
                phantom_edges,
                isolated_cordless_cycles: diagonal.stats.isolated_cordless_cycles,
            },
            edges,
            meta_vertices,
        }
    }

    fn contract_partition_with_atomic_table<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
        threads: usize,
        table: &AtomicPartitionTable<K>,
        pool: &ThreadPool,
    ) -> PartitionContraction<K> {
        let (diagonal, scan) = std::thread::scope(|scope| {
            let diagonal = scope.spawn(|| Self::compress_diagonal_block(matrix, partition));
            let scan = Self::scan_partition_column_atomic(matrix, partition, threads, table, pool);
            (
                diagonal
                    .join()
                    .expect("diagonal compression worker panicked"),
                scan,
            )
        });
        Self::finish_atomic_partition(partition, threads, table, pool, diagonal, scan)
    }

    fn finish_atomic_partition<const K: usize>(
        partition: usize,
        threads: usize,
        table: &AtomicPartitionTable<K>,
        pool: &ThreadPool,
        diagonal: DiagonalCompression<K>,
        scan: PartitionColumnScanOutput<K>,
    ) -> PartitionContraction<K> {
        let mut edges = scan.edges;
        let mut expansion_edges = Vec::new();
        let mut meta_vertices = scan.meta_vertices;
        meta_vertices.extend(diagonal.meta_vertices);
        let mut phantom_edges = 0u64;

        let diagonal_chunk = diagonal.edges.len().div_ceil(threads.max(1)).max(1);
        let diagonal_outputs = pool.install(|| {
            diagonal
                .edges
                .par_chunks(diagonal_chunk)
                .map(|diagonal_edges| {
                    let mut local_edges = Vec::new();
                    let mut local_expansion = Vec::new();
                    let mut local_meta = Vec::new();
                    let mut local_phantoms = 0u64;
                    for edge in diagonal_edges {
                        let (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) =
                            (edge.first, edge.second)
                        else {
                            continue;
                        };
                        let x_stored = table.get(x.vertex);
                        let y_stored = table.get(y.vertex);
                        let endpoint = |vertex: DiscontinuityEndpoint<K>| PartitionOtherEnd {
                            endpoint: MatrixEndpoint::Phi,
                            side_at_current: vertex.side.inverse(),
                            weight: 1,
                            in_same_part: false,
                            processed: true,
                        };
                        let x_end = x_stored.unwrap_or_else(|| endpoint(x));
                        let y_end = y_stored.unwrap_or_else(|| endpoint(y));
                        for (stored, endpoint) in [(x_stored, x), (y_stored, y)] {
                            if stored.is_none() {
                                local_phantoms += 1;
                                let phantom = DiscontinuityEndpoint {
                                    vertex: endpoint.vertex,
                                    side: endpoint.side.inverse(),
                                };
                                local_expansion.push(join_other_ends_with_phantom(
                                    MatrixEndpoint::Phi,
                                    MatrixEndpoint::Vertex(phantom),
                                    1,
                                    Some(phantom),
                                ));
                            }
                        }
                        if x_end.endpoint.is_phi() && y_end.endpoint.is_phi() {
                            local_meta.push(two_weight_meta_vertex(
                                x.vertex,
                                partition,
                                x.side.inverse(),
                                x_end.weight,
                                edge.weight + y_end.weight,
                                false,
                            ));
                        } else {
                            local_edges.push(join_other_ends(
                                x_end.endpoint,
                                y_end.endpoint,
                                x_end.weight + edge.weight + y_end.weight,
                            ));
                        }
                        if x_stored.is_some() {
                            table.mark_processed(x.vertex);
                        }
                        if y_stored.is_some() {
                            table.mark_processed(y.vertex);
                        }
                    }
                    (local_edges, local_expansion, local_meta, local_phantoms)
                })
                .collect::<Vec<_>>()
        });
        for (local_edges, local_expansion, local_meta, local_phantoms) in diagonal_outputs {
            edges.extend(local_edges);
            expansion_edges.extend(local_expansion);
            meta_vertices.extend(local_meta);
            phantom_edges += local_phantoms;
        }

        let scan_chunk = table.slots.len().div_ceil(threads.max(1)).max(1);
        let phantom_outputs = pool.install(|| {
            table
                .slots
                .par_chunks(scan_chunk)
                .map(|slots| {
                    let mut local_edges = Vec::new();
                    let mut local_expansion = Vec::new();
                    let mut local_meta = Vec::new();
                    let mut local_phantoms = 0u64;
                    for slot in slots {
                        let key = slot.key.load(Ordering::Relaxed);
                        if key == AtomicPartitionTable::<K>::EMPTY {
                            continue;
                        }
                        let end = unsafe { (*slot.value.get()).assume_init() }.unpack::<K>();
                        if end.processed {
                            continue;
                        }
                        let vertex = Kmer::from_bits(key as u128);
                        local_phantoms += 1;
                        let phantom = DiscontinuityEndpoint {
                            vertex,
                            side: end.side_at_current.inverse(),
                        };
                        local_expansion.push(join_other_ends_with_phantom(
                            MatrixEndpoint::Phi,
                            MatrixEndpoint::Vertex(phantom),
                            1,
                            Some(phantom),
                        ));
                        if end.endpoint.is_phi() {
                            local_meta.push(two_weight_meta_vertex(
                                vertex,
                                partition,
                                end.side_at_current,
                                end.weight,
                                1,
                                false,
                            ));
                        } else {
                            local_edges.push(join_other_ends(
                                MatrixEndpoint::Phi,
                                end.endpoint,
                                1 + end.weight,
                            ));
                        }
                    }
                    (local_edges, local_expansion, local_meta, local_phantoms)
                })
                .collect::<Vec<_>>()
        });
        for (local_edges, local_expansion, local_meta, local_phantoms) in phantom_outputs {
            edges.extend(local_edges);
            expansion_edges.extend(local_expansion);
            meta_vertices.extend(local_meta);
            phantom_edges += local_phantoms;
        }
        PartitionContraction {
            partition,
            compressed_diagonal_edges: diagonal.expansion_edges,
            expansion_edges,
            stats: PartitionContractionStats {
                input_non_diagonal_edges: scan.input_non_diagonal_edges,
                compressed_diagonal_edges: diagonal.stats.compressed_edges,
                output_edges: edges.len() as u64,
                meta_vertices: meta_vertices.len() as u64,
                phantom_edges,
                isolated_cordless_cycles: diagonal.stats.isolated_cordless_cycles,
            },
            edges,
            meta_vertices,
        }
    }

    fn contract_partition_with_flat_table<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
        table: &mut FlatPartitionTable<K>,
    ) -> PartitionContraction<K> {
        let (diagonal, scan) = std::thread::scope(|scope| {
            let diagonal = scope.spawn(|| Self::compress_diagonal_block(matrix, partition));
            let scan = Self::scan_partition_column_flat(matrix, partition, table);
            (
                diagonal
                    .join()
                    .expect("diagonal compression worker panicked"),
                scan,
            )
        });
        let mut edges = scan.edges;
        let mut expansion_edges = Vec::new();
        let mut meta_vertices = scan.meta_vertices;
        meta_vertices.extend(diagonal.meta_vertices);
        let mut phantom_edges = 0u64;
        for edge in &diagonal.edges {
            let (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) = (edge.first, edge.second)
            else {
                continue;
            };
            for endpoint in [x, y] {
                if table.get(endpoint.vertex).is_none() {
                    phantom_edges += 1;
                    let phantom = DiscontinuityEndpoint {
                        vertex: endpoint.vertex,
                        side: endpoint.side.inverse(),
                    };
                    expansion_edges.push(join_other_ends_with_phantom(
                        MatrixEndpoint::Phi,
                        MatrixEndpoint::Vertex(phantom),
                        1,
                        Some(phantom),
                    ));
                    table.insert(
                        endpoint.vertex,
                        PartitionOtherEnd {
                            endpoint: MatrixEndpoint::Phi,
                            side_at_current: endpoint.side.inverse(),
                            weight: 1,
                            in_same_part: false,
                            processed: false,
                        },
                    );
                }
            }
            let x_end = table.get(x.vertex).expect("inserted x endpoint");
            let y_end = table.get(y.vertex).expect("inserted y endpoint");
            if x_end.endpoint.is_phi() && y_end.endpoint.is_phi() {
                meta_vertices.push(two_weight_meta_vertex(
                    x.vertex,
                    partition,
                    x.side.inverse(),
                    x_end.weight,
                    edge.weight + y_end.weight,
                    false,
                ));
            } else {
                edges.push(join_other_ends(
                    x_end.endpoint,
                    y_end.endpoint,
                    x_end.weight + edge.weight + y_end.weight,
                ));
            }
            table.mark_processed(x.vertex);
            table.mark_processed(y.vertex);
        }
        for &idx in &table.occupied {
            let end = unsafe { table.values[idx].assume_init() };
            if end.processed {
                continue;
            }
            let vertex = Kmer::from_bits(table.keys[idx] as u128);
            phantom_edges += 1;
            let phantom = DiscontinuityEndpoint {
                vertex,
                side: end.side_at_current.inverse(),
            };
            expansion_edges.push(join_other_ends_with_phantom(
                MatrixEndpoint::Phi,
                MatrixEndpoint::Vertex(phantom),
                1,
                Some(phantom),
            ));
            if end.endpoint.is_phi() {
                meta_vertices.push(two_weight_meta_vertex(
                    vertex,
                    partition,
                    end.side_at_current,
                    end.weight,
                    1,
                    false,
                ));
            } else {
                edges.push(join_other_ends(
                    MatrixEndpoint::Phi,
                    end.endpoint,
                    1 + end.weight,
                ));
            }
        }
        PartitionContraction {
            partition,
            compressed_diagonal_edges: diagonal.expansion_edges,
            expansion_edges,
            stats: PartitionContractionStats {
                input_non_diagonal_edges: scan.input_non_diagonal_edges,
                compressed_diagonal_edges: diagonal.stats.compressed_edges,
                output_edges: edges.len() as u64,
                meta_vertices: meta_vertices.len() as u64,
                phantom_edges,
                isolated_cordless_cycles: diagonal.stats.isolated_cordless_cycles,
            },
            edges,
            meta_vertices,
        }
    }

    fn scan_partition_column_flat<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
        table: &mut FlatPartitionTable<K>,
    ) -> PartitionColumnScanOutput<K> {
        table.clear();
        let column_edges = (0..partition)
            .map(|row| matrix.block(row, partition).len())
            .sum::<usize>();
        let mut edges = Vec::with_capacity(column_edges / 2);
        let mut meta_vertices = Vec::new();
        let mut input_non_diagonal_edges = 0;
        for row in 0..partition {
            for edge in matrix.block(row, partition) {
                let Some((lower, current)) = endpoint_in_partition(matrix, edge, partition) else {
                    continue;
                };
                input_non_diagonal_edges += 1;
                table.absorb(
                    current.vertex,
                    PartitionOtherEnd {
                        endpoint: lower,
                        side_at_current: current.side,
                        weight: edge.weight,
                        in_same_part: false,
                        processed: false,
                    },
                    partition,
                    &mut edges,
                    &mut meta_vertices,
                );
            }
        }
        PartitionColumnScanOutput {
            edges,
            meta_vertices,
            input_non_diagonal_edges,
        }
    }

    fn scan_partition_column_owned<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
        threads: usize,
        tables: &mut OwnedPartitionTables<K>,
    ) -> PartitionColumnScanOutput<K> {
        tables.clear();
        let routing_workers = threads.max(1).min(partition);
        let rows_per_worker = partition.div_ceil(routing_workers);
        let owner_count = tables.maps.len();
        let owner_mask = tables.mask;
        let routed = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(routing_workers);
            for row_start in (0..partition).step_by(rows_per_worker) {
                let row_end = (row_start + rows_per_worker).min(partition);
                handles.push(scope.spawn(move || {
                    let mut owners = (0..owner_count)
                        .map(|_| Vec::<PartitionColumnIncoming<K>>::new())
                        .collect::<Vec<_>>();
                    let mut input_edges = 0u64;
                    for row in row_start..row_end {
                        for edge in matrix.block(row, partition) {
                            let Some((lower, current)) =
                                endpoint_in_partition(matrix, edge, partition)
                            else {
                                continue;
                            };
                            input_edges += 1;
                            let owner = partition_column_vertex_shard(current.vertex, owner_mask);
                            owners[owner].push(PartitionColumnIncoming {
                                vertex: current.vertex,
                                end: PartitionOtherEnd {
                                    endpoint: lower,
                                    side_at_current: current.side,
                                    weight: edge.weight,
                                    in_same_part: false,
                                    processed: false,
                                },
                            });
                        }
                    }
                    (owners, input_edges)
                }));
            }
            handles
                .into_iter()
                .map(|handle| handle.join().expect("partition routing worker panicked"))
                .collect::<Vec<_>>()
        });

        let input_non_diagonal_edges = routed.iter().map(|(_, count)| *count).sum();
        let mut owner_inputs = (0..owner_count)
            .map(|_| Vec::<PartitionColumnIncoming<K>>::new())
            .collect::<Vec<_>>();
        for (mut worker_owners, _) in routed {
            for (owner, incoming) in worker_owners.iter_mut().enumerate() {
                owner_inputs[owner].append(incoming);
            }
        }

        let outputs = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(owner_count);
            for (map, incoming) in tables.maps.iter_mut().zip(owner_inputs) {
                handles.push(scope.spawn(move || {
                    map.reserve(incoming.len().saturating_sub(map.capacity()));
                    let mut edges = Vec::with_capacity(incoming.len() / 2);
                    let mut meta_vertices = Vec::new();
                    for incoming in incoming {
                        absorb_partition_other_end(
                            incoming.vertex,
                            incoming.end,
                            partition,
                            map,
                            &mut edges,
                            &mut meta_vertices,
                        );
                    }
                    (edges, meta_vertices)
                }));
            }
            handles
                .into_iter()
                .map(|handle| handle.join().expect("partition owner worker panicked"))
                .collect::<Vec<_>>()
        });
        let mut edges = Vec::new();
        let mut meta_vertices = Vec::new();
        for (worker_edges, worker_meta) in outputs {
            edges.extend(worker_edges);
            meta_vertices.extend(worker_meta);
        }
        PartitionColumnScanOutput {
            edges,
            meta_vertices,
            input_non_diagonal_edges,
        }
    }

    fn scan_partition_column_with_threads<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
        threads: usize,
    ) -> PartitionColumnScan<K> {
        let column_edges = (0..partition)
            .map(|row| matrix.block(row, partition).len())
            .sum::<usize>();
        if threads <= 1 || partition < 2 || column_edges < 64 * 1024 {
            return Self::scan_partition_column(matrix, partition);
        }
        if std::env::var_os("CF3_RS_LEAPMAP_CONTRACTOR").is_some() {
            return Self::scan_partition_column_leapmap_with_threads(
                matrix,
                partition,
                threads,
                column_edges,
            );
        }
        if std::env::var_os("CF3_RS_SCC_CONTRACTOR").is_some() {
            return Self::scan_partition_column_scc_with_threads(
                matrix,
                partition,
                threads,
                column_edges,
            );
        }
        if std::env::var_os("CF3_RS_PARALLEL_CONTRACTOR").is_none() {
            return Self::scan_partition_column(matrix, partition);
        }
        let workers = threads.max(1).min(partition);
        let shard_count = (workers * 8).next_power_of_two().max(8);
        let shard_mask = shard_count - 1;
        let mut shards = (0..shard_count)
            .map(|_| Vec::<PartitionColumnIncoming<K>>::new())
            .collect::<Vec<_>>();
        let mut input_non_diagonal_edges = 0u64;
        for row in 0..partition {
            for edge in matrix.block(row, partition) {
                let (lower, current) = match endpoint_in_partition(matrix, edge, partition) {
                    Some(pair) => pair,
                    None => continue,
                };
                input_non_diagonal_edges += 1;
                let incoming = PartitionOtherEnd {
                    endpoint: lower,
                    side_at_current: current.side,
                    weight: edge.weight,
                    in_same_part: false,
                    processed: false,
                };
                let shard_id = partition_column_vertex_shard(current.vertex, shard_mask);
                shards[shard_id].push(PartitionColumnIncoming {
                    vertex: current.vertex,
                    end: incoming,
                });
            }
        }

        let shards_per_worker = shard_count.div_ceil(workers);
        let worker_outputs = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(workers);
            for shard_start in (0..shard_count).step_by(shards_per_worker) {
                let shard_end = (shard_start + shards_per_worker).min(shard_count);
                let shard_slice = &shards[shard_start..shard_end];
                handles.push(scope.spawn(move || {
                    let shard_edges = shard_slice.iter().map(Vec::len).sum::<usize>();
                    let mut table =
                        FastHashMap::<Kmer<K>, PartitionOtherEnd<K>>::with_capacity_and_hasher(
                            shard_edges,
                            FastBuildHasher::default(),
                        );
                    let mut edges = Vec::with_capacity(shard_edges / 2);
                    let mut meta_vertices = Vec::new();
                    for incoming in shard_slice.iter().flatten() {
                        absorb_partition_other_end(
                            incoming.vertex,
                            incoming.end,
                            partition,
                            &mut table,
                            &mut edges,
                            &mut meta_vertices,
                        );
                    }
                    (table, edges, meta_vertices)
                }));
            }

            let mut outputs = Vec::with_capacity(handles.len());
            for handle in handles {
                outputs.push(
                    handle
                        .join()
                        .expect("partition column scan worker panicked"),
                );
            }
            outputs
        });

        let mut table = FastHashMap::<Kmer<K>, PartitionOtherEnd<K>>::with_capacity_and_hasher(
            column_edges,
            FastBuildHasher::default(),
        );
        let mut edges = Vec::with_capacity(column_edges / 2);
        let mut meta_vertices = Vec::new();
        for (worker_table, worker_edges, worker_meta_vertices) in worker_outputs {
            edges.extend(worker_edges);
            meta_vertices.extend(worker_meta_vertices);
            table.extend(worker_table);
        }
        edges.shrink_to_fit();
        meta_vertices.shrink_to_fit();

        PartitionColumnScan {
            table,
            edges,
            meta_vertices,
            input_non_diagonal_edges,
        }
    }

    fn scan_partition_column_leapmap_with_threads<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
        threads: usize,
        column_edges: usize,
    ) -> PartitionColumnScan<K> {
        let workers = threads.max(1).min(partition);
        let shard_count = (workers * 8).next_power_of_two().max(8);
        let shard_mask = shard_count - 1;
        let mut shards = (0..shard_count)
            .map(|_| Vec::<PartitionColumnIncoming<K>>::new())
            .collect::<Vec<_>>();
        let mut input_non_diagonal_edges = 0u64;
        for row in 0..partition {
            for edge in matrix.block(row, partition) {
                let (lower, current) = match endpoint_in_partition(matrix, edge, partition) {
                    Some(pair) => pair,
                    None => continue,
                };
                input_non_diagonal_edges += 1;
                let incoming = PartitionOtherEnd {
                    endpoint: lower,
                    side_at_current: current.side,
                    weight: edge.weight,
                    in_same_part: false,
                    processed: false,
                };
                let shard_id = partition_column_vertex_shard(current.vertex, shard_mask);
                shards[shard_id].push(PartitionColumnIncoming {
                    vertex: current.vertex,
                    end: incoming,
                });
            }
        }

        let shards_per_worker = shard_count.div_ceil(workers);
        let worker_outputs = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(workers);
            for shard_start in (0..shard_count).step_by(shards_per_worker) {
                let shard_end = (shard_start + shards_per_worker).min(shard_count);
                let shard_slice = &shards[shard_start..shard_end];
                handles.push(scope.spawn(move || {
                    let shard_edges = shard_slice.iter().map(Vec::len).sum::<usize>();
                    let table =
                        LeapMap::<Kmer<K>, usize, FastBuildHasher>::with_capacity_and_hasher(
                            shard_edges.saturating_mul(2).max(8),
                            FastBuildHasher::default(),
                        );
                    let mut payloads = Vec::<PartitionOtherEnd<K>>::with_capacity(shard_edges);
                    let mut edges = Vec::with_capacity(shard_edges / 2);
                    let mut meta_vertices = Vec::new();
                    for incoming in shard_slice.iter().flatten() {
                        match table.get(&incoming.vertex) {
                            None => {
                                let idx = payloads.len();
                                payloads.push(incoming.end);
                                table.insert(incoming.vertex, idx);
                            }
                            Some(mut existing_ref) => {
                                let idx = existing_ref
                                    .value()
                                    .expect("leapmap partition column entry disappeared");
                                let existing = payloads[idx];
                                let mut updated = existing;
                                if updated.in_same_part {
                                    updated = incoming.end;
                                }
                                updated.processed = true;
                                payloads[idx] = updated;
                                absorb_partition_collision_outputs(
                                    incoming.vertex,
                                    incoming.end,
                                    existing,
                                    partition,
                                    &mut edges,
                                    &mut meta_vertices,
                                );
                            }
                        }
                    }
                    let mut compact_table =
                        FastHashMap::<Kmer<K>, PartitionOtherEnd<K>>::with_capacity_and_hasher(
                            table.len(),
                            FastBuildHasher::default(),
                        );
                    for (vertex, idx) in table {
                        compact_table.insert(vertex, payloads[idx]);
                    }
                    (compact_table, edges, meta_vertices)
                }));
            }

            let mut outputs = Vec::with_capacity(handles.len());
            for handle in handles {
                outputs.push(
                    handle
                        .join()
                        .expect("leapmap partition column scan worker panicked"),
                );
            }
            outputs
        });

        let mut table = FastHashMap::<Kmer<K>, PartitionOtherEnd<K>>::with_capacity_and_hasher(
            column_edges,
            FastBuildHasher::default(),
        );
        let mut edges = Vec::with_capacity(column_edges / 2);
        let mut meta_vertices = Vec::new();
        for (worker_table, worker_edges, worker_meta_vertices) in worker_outputs {
            table.extend(worker_table);
            edges.extend(worker_edges);
            meta_vertices.extend(worker_meta_vertices);
        }
        edges.shrink_to_fit();
        meta_vertices.shrink_to_fit();

        PartitionColumnScan {
            table,
            edges,
            meta_vertices,
            input_non_diagonal_edges,
        }
    }

    fn scan_partition_column_scc_with_threads<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        partition: usize,
        threads: usize,
        column_edges: usize,
    ) -> PartitionColumnScan<K> {
        let workers = threads.max(1).min(partition);
        let shard_count = (workers * 8).next_power_of_two().max(8);
        let shard_mask = shard_count - 1;
        let mut shards = (0..shard_count)
            .map(|_| Vec::<PartitionColumnIncoming<K>>::new())
            .collect::<Vec<_>>();
        let mut input_non_diagonal_edges = 0u64;
        for row in 0..partition {
            for edge in matrix.block(row, partition) {
                let (lower, current) = match endpoint_in_partition(matrix, edge, partition) {
                    Some(pair) => pair,
                    None => continue,
                };
                input_non_diagonal_edges += 1;
                let incoming = PartitionOtherEnd {
                    endpoint: lower,
                    side_at_current: current.side,
                    weight: edge.weight,
                    in_same_part: false,
                    processed: false,
                };
                let shard_id = partition_column_vertex_shard(current.vertex, shard_mask);
                shards[shard_id].push(PartitionColumnIncoming {
                    vertex: current.vertex,
                    end: incoming,
                });
            }
        }

        let shards_per_worker = shard_count.div_ceil(workers);
        let worker_outputs = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(workers);
            for shard_start in (0..shard_count).step_by(shards_per_worker) {
                let shard_end = (shard_start + shards_per_worker).min(shard_count);
                let shard_slice = &shards[shard_start..shard_end];
                handles.push(scope.spawn(move || {
                    let shard_edges = shard_slice.iter().map(Vec::len).sum::<usize>();
                    let table =
                        SccHashMap::<Kmer<K>, PartitionOtherEnd<K>, FastBuildHasher>::with_capacity_and_hasher(
                            shard_edges.max(8),
                            FastBuildHasher::default(),
                        );
                    let mut edges = Vec::with_capacity(shard_edges / 2);
                    let mut meta_vertices = Vec::new();
                    for incoming in shard_slice.iter().flatten() {
                        match table.entry_sync(incoming.vertex) {
                            SccEntry::Vacant(vacant) => {
                                vacant.insert_entry(incoming.end);
                            }
                            SccEntry::Occupied(mut occupied) => {
                                let existing = *occupied.get();
                                let current_end = occupied.get_mut();
                                if current_end.in_same_part {
                                    *current_end = incoming.end;
                                }
                                current_end.processed = true;
                                absorb_partition_collision_outputs(
                                    incoming.vertex,
                                    incoming.end,
                                    existing,
                                    partition,
                                    &mut edges,
                                    &mut meta_vertices,
                                );
                            }
                        }
                    }
                    let mut compact_table =
                        FastHashMap::<Kmer<K>, PartitionOtherEnd<K>>::with_capacity_and_hasher(
                            shard_edges,
                            FastBuildHasher::default(),
                        );
                    table.iter_sync(|vertex, end| {
                        compact_table.insert(*vertex, *end);
                        true
                    });
                    (compact_table, edges, meta_vertices)
                }));
            }

            let mut outputs = Vec::with_capacity(handles.len());
            for handle in handles {
                outputs.push(
                    handle
                        .join()
                        .expect("scc partition column scan worker panicked"),
                );
            }
            outputs
        });

        let mut table = FastHashMap::<Kmer<K>, PartitionOtherEnd<K>>::with_capacity_and_hasher(
            column_edges,
            FastBuildHasher::default(),
        );
        let mut edges = Vec::with_capacity(column_edges / 2);
        let mut meta_vertices = Vec::new();
        for (worker_table, worker_edges, worker_meta_vertices) in worker_outputs {
            table.extend(worker_table);
            edges.extend(worker_edges);
            meta_vertices.extend(worker_meta_vertices);
        }
        edges.shrink_to_fit();
        meta_vertices.shrink_to_fit();

        PartitionColumnScan {
            table,
            edges,
            meta_vertices,
            input_non_diagonal_edges,
        }
    }

    pub fn contract_all_partitions<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
    ) -> FullSerialDiscontinuityContraction<K> {
        Self::contract_all_partitions_with_threads(matrix, 1)
    }

    pub fn contract_all_partitions_with_threads<const K: usize>(
        matrix: &SerialEdgeMatrix<K>,
        threads: usize,
    ) -> FullSerialDiscontinuityContraction<K> {
        Self::contract_all_partitions_owned(matrix.clone(), threads)
    }

    fn contract_all_partitions_owned<const K: usize>(
        working: SerialEdgeMatrix<K>,
        threads: usize,
    ) -> FullSerialDiscontinuityContraction<K> {
        Self::contract_all_partitions_owned_impl::<K, true>(working, threads)
    }

    fn contract_blocked_external<const K: usize>(
        mut working: BlockedEdgeMatrix<K>,
        threads: usize,
    ) -> Result<ExternalBlockedContraction<K>, SerialCollationError> {
        let setup_started = Instant::now();
        let partition_count = working.partition_count();
        let total_partitions = partition_count - 1;
        let mut compressed_diagonal_edges =
            (0..partition_count).map(|_| Vec::new()).collect::<Vec<_>>();
        let meta_vertex_dir = working.dir.join("contracted-meta-path-info");
        if meta_vertex_dir.exists() {
            fs::remove_dir_all(&meta_vertex_dir).map_err(|source| SerialCollationError::Io {
                path: meta_vertex_dir.clone(),
                source,
            })?;
        }
        fs::create_dir_all(&meta_vertex_dir).map_err(|source| SerialCollationError::Io {
            path: meta_vertex_dir.clone(),
            source,
        })?;
        let mut meta_vertex_count = 0u64;
        let mut meta_vertices_per_partition = vec![0usize; partition_count];
        let mut reusable_table =
            FastHashMap::<Kmer<K>, PartitionOtherEnd<K>>::with_hasher(FastBuildHasher::default());
        let mut diagonal_ends =
            FastHashMap::<Kmer<K>, DiagonalOtherEnd<K>>::with_hasher(FastBuildHasher::default());
        let mut column = SerialEdgeMatrix::new(working.vertex_partitions)
            .map_err(|_| SerialCollationError::MalformedCoordBucket(working.dir.clone()))?;
        let max_partition_vertices = (1..partition_count)
            .map(|col| {
                let column_edges = (0..=col)
                    .map(|row| working.blocks[working.block_index(row, col)].edges)
                    .sum::<usize>();
                column_edges + working.blocks[working.block_index(col, col)].edges
            })
            .max()
            .unwrap_or(1);
        let atomic_table =
            (K <= 31).then(|| AtomicPartitionTable::<K>::with_max_entries(max_partition_vertices));
        let contraction_pool = ThreadPoolBuilder::new()
            .num_threads(threads.max(1))
            .build()
            .map_err(|_| SerialCollationError::WorkerPanic)?;
        let appenders = ConcurrentBlockedEdgeWriters::new(&working);
        let setup_elapsed = setup_started.elapsed();
        if let Some(table) = &atomic_table {
            eprintln!(
                "cuttlefish3-rs: contraction table setup {:.3}s; max entries {}, capacity {}, slot {} byte(s)",
                setup_elapsed.as_secs_f64(),
                max_partition_vertices,
                table.slots.len(),
                std::mem::size_of::<AtomicPartitionSlot>(),
            );
        }
        let mut stats = FullSerialContractionStats {
            input_edges: working.stats.edges,
            ..FullSerialContractionStats::default()
        };
        let started = Instant::now();
        let mut block_load_elapsed = Duration::default();
        let mut partition_contract_elapsed = Duration::default();
        let mut reinsert_elapsed = Duration::default();
        let mut direct_timings = BlockedContractTimings::default();
        let mut meta_write_elapsed = Duration::default();
        eprintln!(
            "cuttlefish3-rs: contracting {} blocked discontinuity partition(s) with {} worker(s)",
            total_partitions, threads
        );

        for (completed, partition) in (1..partition_count).rev().enumerate() {
            let phase_started = Instant::now();
            let (mut contracted, pre_reinserted_edges) = if let Some(table) = &atomic_table {
                Self::contract_blocked_partition_atomic(
                    &mut working,
                    &appenders,
                    &mut column,
                    partition,
                    threads,
                    table,
                    &mut diagonal_ends,
                    &contraction_pool,
                    &mut direct_timings,
                )?
            } else {
                let load_started = Instant::now();
                working.load_column_into(partition, threads, &mut column)?;
                block_load_elapsed += load_started.elapsed();
                (
                    Self::contract_partition_with_reusable_table(
                        &column,
                        partition,
                        threads,
                        &mut reusable_table,
                    ),
                    0,
                )
            };
            partition_contract_elapsed += phase_started.elapsed();
            stats.partition_output_edges += contracted.stats.output_edges;
            stats.phantom_edges += contracted.stats.phantom_edges;
            stats.isolated_cordless_cycles += contracted.stats.isolated_cordless_cycles;
            stats.reinserted_edges += pre_reinserted_edges;
            compressed_diagonal_edges[partition] =
                std::mem::take(&mut contracted.compressed_diagonal_edges);

            let phase_started = Instant::now();
            let vertex_partitions = working.vertex_partitions;
            let reinserted = contraction_pool.install(|| {
                let (expansion_result, reinsert_result) = rayon::join(
                    || {
                        emit_contracted_edge_chunks(
                            contracted.expansion_edges,
                            vertex_partitions,
                            None,
                            &appenders,
                        )
                    },
                    || {
                        emit_contracted_edge_chunks(
                            contracted.edges,
                            vertex_partitions,
                            Some(partition),
                            &appenders,
                        )
                    },
                );
                expansion_result?;
                reinsert_result
            })?;
            stats.reinserted_edges += reinserted;
            reinsert_elapsed += phase_started.elapsed();
            let meta_write_started = Instant::now();
            write_meta_vertex_bucket_parallel(
                &meta_vertex_dir,
                partition,
                &contracted.meta_vertices,
                &contraction_pool,
            )?;
            meta_write_elapsed += meta_write_started.elapsed();
            meta_vertex_count += contracted.meta_vertices.len() as u64;
            meta_vertices_per_partition[partition] += contracted.meta_vertices.len();
            for row in 0..=partition {
                column.blocks[row][partition].clear();
            }
            column.stats = SerialEdgeMatrixStats::default();
            report_discontinuity_contraction_progress(completed + 1, total_partitions, started);
        }
        let flush_started = Instant::now();
        for row in 0..partition_count {
            for col in row..partition_count {
                let block_id = working.block_index(row, col);
                appenders.flush_block(block_id)?;
                let added = appenders.take_added_edges(row, col);
                if added != 0 {
                    working.blocks[block_id].edges += added;
                    working.stats.edges += added as u64;
                    working.stats.phi_edges += u64::from(row == 0) * added as u64;
                    working.stats.diagonal_edges += u64::from(row == col) * added as u64;
                }
            }
        }
        working.flush_all_with_threads(threads)?;
        let block_flush_elapsed = flush_started.elapsed();
        eprintln!(
            "cuttlefish3-rs: blocked contraction detail: setup {:.3}s, block load/decode {:.3}s, partition contraction {:.3}s, reinsertion {:.3}s, meta write {:.3}s, block flush {:.3}s",
            setup_elapsed.as_secs_f64(),
            block_load_elapsed.as_secs_f64(),
            partition_contract_elapsed.as_secs_f64(),
            reinsert_elapsed.as_secs_f64(),
            meta_write_elapsed.as_secs_f64(),
            block_flush_elapsed.as_secs_f64(),
        );
        if atomic_table.is_some() {
            eprintln!(
                "cuttlefish3-rs: fused contraction phases: column flush {:.3}s, table clear {:.3}s, diagonal {:.3}s, raw read {:.3}s, table scan {:.3}s, finalize {:.3}s",
                direct_timings.flush.as_secs_f64(),
                direct_timings.clear.as_secs_f64(),
                direct_timings.diagonal.as_secs_f64(),
                direct_timings.read.as_secs_f64(),
                direct_timings.scan.as_secs_f64(),
                direct_timings.finish.as_secs_f64(),
            );
        }
        stats.partitions = total_partitions as u64;
        stats.meta_vertices = meta_vertex_count;
        Ok(ExternalBlockedContraction {
            vertex_partitions: working.vertex_partitions,
            expansion_matrix: working,
            compressed_diagonal_edges,
            meta_vertex_dir,
            meta_vertex_count,
            meta_vertices_per_partition,
            stats,
        })
    }

    fn contract_all_partitions_owned_impl<const K: usize, const SORT_EXPANSION_EDGES: bool>(
        mut working: SerialEdgeMatrix<K>,
        threads: usize,
    ) -> FullSerialDiscontinuityContraction<K> {
        let mut partitions = Vec::new();
        let mut final_edges = Vec::new();
        let mut meta_vertices = Vec::new();
        let mut compressed_diagonal_edges = (0..working.partition_count())
            .map(|_| Vec::new())
            .collect::<Vec<_>>();
        let mut reusable_table =
            FastHashMap::<Kmer<K>, PartitionOtherEnd<K>>::with_hasher(FastBuildHasher::default());
        let mut stats = FullSerialContractionStats {
            input_edges: working.stats().edges,
            ..FullSerialContractionStats::default()
        };
        let total_partitions = working.partition_count().saturating_sub(1);
        let started = Instant::now();
        if total_partitions >= 16 {
            eprintln!(
                "cuttlefish3-rs: contracting {} discontinuity partition(s) with {} worker(s)",
                total_partitions, threads
            );
        }

        for (completed, partition) in (1..working.partition_count()).rev().enumerate() {
            let mut contracted = if SORT_EXPANSION_EDGES {
                Self::contract_partition_with_threads(&working, partition, threads)
            } else {
                Self::contract_partition_with_reusable_table(
                    &working,
                    partition,
                    threads,
                    &mut reusable_table,
                )
            };
            stats.partition_output_edges += contracted.stats.output_edges;
            stats.phantom_edges += contracted.stats.phantom_edges;
            stats.isolated_cordless_cycles += contracted.stats.isolated_cordless_cycles;
            if !SORT_EXPANSION_EDGES && contracted.partition < compressed_diagonal_edges.len() {
                compressed_diagonal_edges[contracted.partition] =
                    std::mem::take(&mut contracted.compressed_diagonal_edges);
            }
            for edge in &contracted.expansion_edges {
                working.add_edge_with_orientation_and_phantom(
                    edge.first,
                    edge.second,
                    edge.weight,
                    edge.unitig_index,
                    edge.unitig_exit_side,
                    edge.phantom_unitig,
                );
            }

            for edge in &contracted.edges {
                if max_endpoint_partition(&working, edge.first, edge.second) < partition {
                    working.add_edge_with_orientation_and_phantom(
                        edge.first,
                        edge.second,
                        edge.weight,
                        edge.unitig_index,
                        edge.unitig_exit_side,
                        edge.phantom_unitig,
                    );
                    stats.reinserted_edges += 1;
                } else if SORT_EXPANSION_EDGES {
                    final_edges.push(edge.clone());
                }
            }

            meta_vertices.extend(contracted.meta_vertices.iter().cloned());
            if SORT_EXPANSION_EDGES {
                if contracted.partition < compressed_diagonal_edges.len() {
                    compressed_diagonal_edges[contracted.partition] =
                        contracted.compressed_diagonal_edges.clone();
                }
                partitions.push(contracted);
            }
            report_discontinuity_contraction_progress(completed + 1, total_partitions, started);
        }

        if SORT_EXPANSION_EDGES {
            final_edges.sort_by_key(|edge| {
                (
                    endpoint_sort_key(edge.first),
                    endpoint_sort_key(edge.second),
                    edge.weight,
                )
            });
            meta_vertices.sort_by_key(|meta| {
                (
                    meta.partition,
                    meta.vertex.as_u128(),
                    meta.entry_side as u8,
                    meta.weight,
                    meta.is_cycle,
                )
            });
        }

        stats.partitions = total_partitions as u64;
        stats.final_edges = final_edges.len() as u64;
        stats.meta_vertices = meta_vertices.len() as u64;

        let expansion_edges = if SORT_EXPANSION_EDGES {
            let mut expansion_edges = working.edges().cloned().collect::<Vec<_>>();
            expansion_edges.sort_by_key(|edge| {
                (
                    endpoint_sort_key(edge.first),
                    endpoint_sort_key(edge.second),
                    edge.weight,
                    edge.unitig_index,
                )
            });
            expansion_edges
        } else {
            Vec::new()
        };

        FullSerialDiscontinuityContraction {
            vertex_partitions: working.vertex_partitions(),
            final_edges,
            expansion_edges,
            expansion_matrix: working,
            compressed_diagonal_edges,
            meta_vertices,
            partitions,
            stats,
        }
    }

    pub fn contract<const K: usize>(matrix: &SerialEdgeMatrix<K>) -> SerialContraction {
        let mut endpoint_ids = FastHashMap::<EndpointKey<K>, usize>::with_capacity_and_hasher(
            matrix.stats.edges as usize * 2,
            FastBuildHasher::default(),
        );
        let mut adjacency = Vec::<Vec<(usize, u64, bool)>>::new();

        for edge in matrix.edges() {
            let first = endpoint_id(&mut endpoint_ids, &mut adjacency, edge.first);
            let second = endpoint_id(&mut endpoint_ids, &mut adjacency, edge.second);
            let phi_edge = edge.first.is_phi() || edge.second.is_phi();
            adjacency[first].push((second, edge.weight, phi_edge));
            adjacency[second].push((first, edge.weight, phi_edge));
        }

        let mut seen = vec![false; adjacency.len()];
        let mut components = Vec::new();
        for start in 0..adjacency.len() {
            if seen[start] {
                continue;
            }

            let mut stack = vec![start];
            let mut vertices = FastHashSet::default();
            let mut edge_visits = 0u64;
            let mut weight_visits = 0u64;
            let mut phi_visits = 0u64;
            seen[start] = true;

            while let Some(v) = stack.pop() {
                vertices.insert(v);
                for &(next, weight, phi_edge) in &adjacency[v] {
                    edge_visits += 1;
                    weight_visits += weight;
                    phi_visits += u64::from(phi_edge);
                    if !seen[next] {
                        seen[next] = true;
                        stack.push(next);
                    }
                }
            }

            let edges = edge_visits / 2;
            let phi_edges = phi_visits / 2;
            let cyclic = edges > 0 && vertices.iter().all(|&v| adjacency[v].len() == 2);
            components.push(SerialComponent {
                endpoints: vertices.len(),
                edges,
                weight: weight_visits / 2,
                phi_edges,
                cyclic,
            });
        }

        components.sort_by_key(|component| {
            (
                std::cmp::Reverse(component.edges),
                std::cmp::Reverse(component.weight),
                component.endpoints,
            )
        });

        let stats = SerialContractionStats {
            input_edges: matrix.stats.edges,
            components: components.len() as u64,
            phi_edges: matrix.stats.phi_edges,
            cyclic_components: components
                .iter()
                .filter(|component| component.cyclic)
                .count() as u64,
        };

        SerialContraction { components, stats }
    }
}

pub struct SerialDiscontinuityExpander;

impl SerialDiscontinuityExpander {
    pub fn infer<const K: usize>(
        source: PathInfo<K>,
        source_side: Side,
        target_side: Side,
        weight: u64,
    ) -> PathInfo<K> {
        let rank = if source_side == source.exit_side {
            source.rank + weight
        } else {
            source.rank.saturating_sub(weight)
        };
        let exit_side = if source_side == source.exit_side {
            target_side.inverse()
        } else {
            target_side
        };

        PathInfo {
            path_id: source.path_id,
            rank,
            exit_side,
            is_cycle: source.is_cycle,
        }
    }

    #[inline(always)]
    fn infer_compact(
        source: CompactExpansionPathInfo,
        source_side: Side,
        target_side: Side,
        weight: u64,
    ) -> CompactExpansionPathInfo {
        let rank = if source_side == source.exit_side() {
            source.rank() + weight
        } else {
            source.rank().saturating_sub(weight)
        };
        let exit_side = if source_side == source.exit_side() {
            target_side.inverse()
        } else {
            target_side
        };
        CompactExpansionPathInfo::new(source.path_id, rank, exit_side, source.is_cycle())
    }

    fn expand_non_diagonal_raw<const K: usize>(
        matrix: &BlockedEdgeMatrix<K>,
        partition: usize,
        map: &ExpansionPathInfoTable<K>,
        vertex_path_info_dir: &Path,
        edge_writers: &ConcurrentStitchedCoordWriters<'_>,
        ranges: &[ExternalLocalUnitigRange],
        range_index: &ExternalRangeIndex,
        ranges_per_bucket: usize,
        error_path: &Path,
        pool: &ThreadPool,
    ) -> Result<
        (
            Vec<(StitchedCoordRecord, DiscontinuityEndpoint<K>)>,
            u64,
            u64,
            u64,
        ),
        SerialCollationError,
    > {
        let blocks = (partition + 1..matrix.partition_count())
            .filter_map(|col| {
                let block_id = matrix.block_index(partition, col);
                let block = &matrix.blocks[block_id];
                (block.edges != 0).then_some((col, block_id))
            })
            .collect::<Vec<_>>();
        let block_bytes = pool.install(|| {
            blocks
                .par_iter()
                .copied()
                .map(|(col, block_id)| {
                    let block = &matrix.blocks[block_id];
                    let bytes = match fs::read(&block.path) {
                        Ok(bytes) => bytes,
                        Err(source) if source.kind() == std::io::ErrorKind::NotFound => Vec::new(),
                        Err(source) => {
                            return Err(SerialCollationError::Io {
                                path: block.path.clone(),
                                source,
                            });
                        }
                    };
                    if bytes.len() % block.record_len != 0 {
                        return Err(SerialCollationError::MalformedCoordBucket(
                            block.path.clone(),
                        ));
                    }
                    Ok::<_, SerialCollationError>((col, block_id, bytes))
                })
                .collect::<Result<Vec<_>, SerialCollationError>>()
        })?;
        let mut tasks = Vec::new();
        for (col, block_id, bytes) in &block_bytes {
            let block = &matrix.blocks[*block_id];
            let record_len = block.record_len;
            let bytes_per_task = (1024 * 1024 / record_len).max(1) * record_len;
            for chunk in bytes.chunks(bytes_per_task) {
                tasks.push((*col, record_len, chunk));
            }
            for chunk in block.buffer.chunks(bytes_per_task) {
                tasks.push((*col, record_len, chunk));
            }
        }
        let outputs = pool.install(|| {
            tasks
                .into_par_iter()
                .map_init(
                    || {
                        (
                            Vec::<u8>::new(),
                            (0..edge_writers.writers.len())
                                .map(|_| Vec::<StitchedCoordRecord>::new())
                                .collect::<Vec<_>>(),
                            Vec::<usize>::new(),
                        )
                    },
                    |(vertex_records, path_records, used_path_buckets),
                     (col, record_len, bytes)| {
                        vertex_records.clear();
                        debug_assert!(used_path_buckets.is_empty());
                        let mut phantoms = Vec::new();
                        let mut unresolved = 0u64;
                        let mut inferred = 0u64;
                        let mut emitted = 0u64;
                        for encoded in bytes.chunks_exact(record_len) {
                            let edge = decode_discontinuity_edge::<K>(encoded);
                            let (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) =
                                (edge.first, edge.second)
                            else {
                                continue;
                            };
                            let Some(x_info) = map.get_compact(x.vertex) else {
                                unresolved += 1;
                                continue;
                            };
                            let y_info = Self::infer_compact(x_info, x.side, y.side, edge.weight);
                            if y_info.rank() > 0 {
                                append_encoded_compact_vertex_path_info::<K>(
                                    vertex_records,
                                    y.vertex,
                                    y_info,
                                );
                                inferred += 1;
                            }
                            if edge.weight == 1 {
                                let record = stitched_record_from_compact_edge_path_info(
                                    &edge,
                                    compact_edge_path_info(&edge, x_info, y_info),
                                    error_path,
                                )?;
                                if let Some(phantom) = edge.phantom_unitig {
                                    phantoms.push((record, phantom));
                                } else {
                                    let bucket_id = edge_path_info_bucket(
                                        &edge,
                                        ranges,
                                        range_index,
                                        ranges_per_bucket,
                                    )
                                    .ok_or_else(|| {
                                        SerialCollationError::MalformedCoordBucket(
                                            error_path.to_path_buf(),
                                        )
                                    })?;
                                    let bucket = &mut path_records[bucket_id];
                                    if bucket.is_empty() {
                                        used_path_buckets.push(bucket_id);
                                    }
                                    bucket.push(record);
                                }
                                emitted += 1;
                            }
                        }
                        if !vertex_records.is_empty() {
                            let path = vertex_path_info_bucket_path(vertex_path_info_dir, col);
                            let mut file = OpenOptions::new()
                                .create(true)
                                .append(true)
                                .open(&path)
                                .map_err(|source| SerialCollationError::Io {
                                    path: path.clone(),
                                    source,
                                })?;
                            file.write_all(&vertex_records)
                                .map_err(|source| SerialCollationError::Io { path, source })?;
                        }
                        for bucket_id in used_path_buckets.drain(..) {
                            edge_writers.write_path_records(bucket_id, &path_records[bucket_id])?;
                            path_records[bucket_id].clear();
                        }
                        Ok::<_, SerialCollationError>((phantoms, unresolved, inferred, emitted))
                    },
                )
                .collect::<Result<Vec<_>, SerialCollationError>>()
        })?;
        let mut phantoms = Vec::new();
        let mut unresolved = 0;
        let mut inferred = 0;
        let mut emitted = 0;
        for (local_phantoms, local_unresolved, local_inferred, local_emitted) in outputs {
            phantoms.extend(local_phantoms);
            unresolved += local_unresolved;
            inferred += local_inferred;
            emitted += local_emitted;
        }
        Ok((phantoms, unresolved, inferred, emitted))
    }

    pub fn expand<const K: usize>(
        contraction: &FullSerialDiscontinuityContraction<K>,
    ) -> SerialExpansion<K> {
        Self::expand_impl(contraction, false)
    }

    pub fn expand_with_original_edges<const K: usize>(
        contraction: &FullSerialDiscontinuityContraction<K>,
    ) -> SerialExpansion<K> {
        Self::expand_impl(contraction, true)
    }

    pub fn expand_cpp_ordered<const K: usize>(
        contraction: &FullSerialDiscontinuityContraction<K>,
    ) -> SerialExpansion<K> {
        Self::expand_cpp_ordered_impl::<K, true, true>(contraction)
    }

    #[allow(dead_code)]
    fn expand_cpp_ordered_external<const K: usize>(
        contraction: &FullSerialDiscontinuityContraction<K>,
    ) -> SerialExpansion<K> {
        Self::expand_cpp_ordered_impl::<K, false, false>(contraction)
    }

    fn expand_cpp_ordered_external_to_range_buckets<const K: usize>(
        contraction: &ExternalBlockedContraction<K>,
        ranges: &[ExternalLocalUnitigRange],
        ranges_per_bucket: usize,
        bucket_count: usize,
        error_path: &Path,
        expansion_dir: &Path,
        threads: usize,
    ) -> Result<RangeBucketedExpansion<K>, SerialCollationError> {
        if expansion_dir.exists() {
            fs::remove_dir_all(expansion_dir).map_err(|source| SerialCollationError::Io {
                path: expansion_dir.to_path_buf(),
                source,
            })?;
        }
        fs::create_dir_all(expansion_dir).map_err(|source| SerialCollationError::Io {
            path: expansion_dir.to_path_buf(),
            source,
        })?;
        let expansion_pool = ThreadPoolBuilder::new()
            .num_threads(threads.max(1))
            .build()
            .map_err(|_| SerialCollationError::WorkerPanic)?;

        let partition_count = contraction.vertex_partitions + 1;
        let vertex_path_info_dir = contraction.meta_vertex_dir.clone();
        let mut vertex_writers =
            VertexPathInfoBucketWriters::<K>::open_existing(&vertex_path_info_dir, partition_count);
        let meta_vertices_per_partition = &contraction.meta_vertices_per_partition;
        let seed_vertices = contraction.meta_vertex_count;

        let mut diagonal_by_partition = contraction.compressed_diagonal_edges.clone();
        diagonal_by_partition.resize_with(partition_count, Vec::new);
        let matrix = &contraction.expansion_matrix;
        let edge_path_info_dir = expansion_dir.join("P_e");
        fs::create_dir_all(&edge_path_info_dir).map_err(|source| SerialCollationError::Io {
            path: edge_path_info_dir.clone(),
            source,
        })?;
        let edge_writers =
            ConcurrentStitchedCoordWriters::new(&edge_path_info_dir, bucket_count, threads);
        let range_index = ExternalRangeIndex::new(ranges);
        let mut phantom_records = Vec::new();
        let mut unresolved_edges = 0u64;
        let mut inferred_vertices = 0u64;
        let mut edge_path_infos = 0u64;
        let max_partition_entries = (1..partition_count)
            .map(|partition| {
                let incoming = (1..partition)
                    .map(|row| matrix.blocks[matrix.block_index(row, partition)].edges)
                    .sum::<usize>();
                let diagonal = matrix.blocks[matrix.block_index(partition, partition)].edges;
                incoming + 2 * diagonal + meta_vertices_per_partition[partition]
            })
            .max()
            .unwrap_or(1);
        let map = ExpansionPathInfoTable::<K>::with_max_entries(max_partition_entries);
        eprintln!(
            "cuttlefish3-rs: expansion table max entries {}, capacity {}, slot {} byte(s)",
            max_partition_entries,
            map.capacity(),
            map.slot_size(),
        );
        let mut row_load_elapsed = Duration::default();
        let mut path_info_load_elapsed = Duration::default();
        let mut path_info_read_elapsed = Duration::default();
        let mut path_info_insert_elapsed = Duration::default();
        let mut path_info_clear_elapsed = Duration::default();
        let mut compressed_diagonal_elapsed = Duration::default();
        let mut non_diagonal_elapsed = Duration::default();
        let mut original_edge_elapsed = Duration::default();
        let expansion_work_started = Instant::now();

        for partition in 1..partition_count {
            let row_load_started = Instant::now();
            let row_blocks = if K <= 31 {
                vec![matrix.read_flushed_block(partition, partition)?]
            } else {
                matrix.read_flushed_row(partition, threads)?
            };
            row_load_elapsed += row_load_started.elapsed();
            let diagonal_block = &row_blocks[0];
            let phase_started = Instant::now();
            map.clear();
            path_info_clear_elapsed += phase_started.elapsed();
            vertex_writers.flush_bucket(partition)?;
            read_vertex_path_info_bucket_into::<K>(
                &vertex_path_info_dir,
                partition,
                error_path,
                &map,
                &expansion_pool,
                &mut path_info_read_elapsed,
                &mut path_info_insert_elapsed,
            )?;
            path_info_load_elapsed += phase_started.elapsed();

            let phase_started = Instant::now();
            for edge in diagonal_by_partition[partition].iter().rev() {
                let (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) =
                    (edge.first, edge.second)
                else {
                    continue;
                };

                if K <= 31 {
                    if let Some(y_info) = map.get_compact(y.vertex) {
                        let x_info = Self::infer_compact(y_info, y.side, x.side, edge.weight);
                        if x_info.rank() > 0 && map.insert_compact(x.vertex, x_info) {
                            inferred_vertices += 1;
                        }
                    } else if let Some(x_info) = map.get_compact(x.vertex) {
                        let y_info = Self::infer_compact(x_info, x.side, y.side, edge.weight);
                        if y_info.rank() > 0 && map.insert_compact(y.vertex, y_info) {
                            inferred_vertices += 1;
                        }
                    } else {
                        unresolved_edges += 1;
                    }
                } else if let Some(y_info) = map.get(y.vertex) {
                    let x_info = Self::infer(y_info, y.side, x.side, edge.weight);
                    if x_info.rank > 0 && map.insert(x.vertex, x_info) {
                        inferred_vertices += 1;
                    }
                } else if let Some(x_info) = map.get(x.vertex) {
                    let y_info = Self::infer(x_info, x.side, y.side, edge.weight);
                    if y_info.rank > 0 && map.insert(y.vertex, y_info) {
                        inferred_vertices += 1;
                    }
                } else {
                    unresolved_edges += 1;
                }
            }
            compressed_diagonal_elapsed += phase_started.elapsed();

            let phase_started = Instant::now();
            let non_diagonal_cols = partition_count.saturating_sub(partition + 1);
            let non_diagonal_edges = row_blocks[1..].iter().map(Vec::len).sum::<usize>();
            let workers = threads.max(1).min(non_diagonal_cols.max(1));
            if K <= 31 {
                let (local_phantoms, local_unresolved, local_inferred, local_edge_infos) =
                    Self::expand_non_diagonal_raw(
                        matrix,
                        partition,
                        &map,
                        &vertex_path_info_dir,
                        &edge_writers,
                        ranges,
                        &range_index,
                        ranges_per_bucket,
                        error_path,
                        &expansion_pool,
                    )?;
                phantom_records.extend(local_phantoms);
                unresolved_edges += local_unresolved;
                inferred_vertices += local_inferred;
                edge_path_infos += local_edge_infos;
            } else if workers == 1 || non_diagonal_cols < 2 || non_diagonal_edges < 64 * 1024 {
                for col in partition + 1..partition_count {
                    for edge in &row_blocks[col - partition] {
                        let (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) =
                            (edge.first, edge.second)
                        else {
                            continue;
                        };
                        let Some(x_info) = map.get(x.vertex) else {
                            unresolved_edges += 1;
                            continue;
                        };
                        let y_info = Self::infer(x_info, x.side, y.side, edge.weight);
                        if y_info.rank > 0 {
                            vertex_writers.write_record(
                                col,
                                &VertexPathInfo {
                                    vertex: y.vertex,
                                    info: y_info,
                                },
                            )?;
                            inferred_vertices += 1;
                        }
                        if edge.weight == 1 {
                            push_edge_path_record_to_range_bucket_writer(
                                edge,
                                edge_path_info(edge, x_info, y_info),
                                ranges,
                                &range_index,
                                ranges_per_bucket,
                                &edge_writers,
                                &mut phantom_records,
                                error_path,
                            )?;
                            edge_path_infos += 1;
                        }
                    }
                }
            } else {
                for col in partition + 1..partition_count {
                    vertex_writers.flush_bucket(col)?;
                }
                let chunk = non_diagonal_cols.div_ceil(workers);
                let worker_outputs = expansion_pool.install(|| {
                    (0..workers)
                        .into_par_iter()
                        .map(|worker_id| {
                            let col_start = partition + 1 + worker_id * chunk;
                            if col_start >= partition_count {
                                return Ok(None);
                            }
                            let col_end = (col_start + chunk).min(partition_count);
                            let map = &map;
                            let row_blocks = &row_blocks;
                            let edge_writers = &edge_writers;
                            let vertex_path_info_dir = &vertex_path_info_dir;
                            let mut range_records = (0..bucket_count)
                                .map(|_| Vec::<StitchedCoordRecord>::new())
                                .collect::<Vec<_>>();
                            let mut local_phantoms =
                                Vec::<(StitchedCoordRecord, DiscontinuityEndpoint<K>)>::new();
                            let mut local_unresolved = 0u64;
                            let mut local_inferred = 0u64;
                            let mut local_edge_infos = 0u64;

                            for col in col_start..col_end {
                                let vertex_path =
                                    vertex_path_info_bucket_path(vertex_path_info_dir, col);
                                let mut vertex_out: Option<BufWriter<File>> = None;
                                for edge in &row_blocks[col - partition] {
                                    let (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) =
                                        (edge.first, edge.second)
                                    else {
                                        continue;
                                    };
                                    let Some(x_info) = map.get(x.vertex) else {
                                        local_unresolved += 1;
                                        continue;
                                    };
                                    let y_info = Self::infer(x_info, x.side, y.side, edge.weight);
                                    if y_info.rank > 0 {
                                        if vertex_out.is_none() {
                                            let file = OpenOptions::new()
                                                .create(true)
                                                .append(true)
                                                .open(&vertex_path)
                                                .map_err(|source| SerialCollationError::Io {
                                                    path: vertex_path.clone(),
                                                    source,
                                                })?;
                                            vertex_out = Some(BufWriter::with_capacity(
                                                VERTEX_PATH_INFO_WRITE_BUFFER,
                                                file,
                                            ));
                                        }
                                        vertex_out
                                            .as_mut()
                                            .expect("vertex writer was just created")
                                            .write_all(
                                                &encoded_vertex_path_info_record(&VertexPathInfo {
                                                    vertex: y.vertex,
                                                    info: y_info,
                                                })[..vertex_path_info_record_len::<K>()],
                                            )
                                            .map_err(|source| SerialCollationError::Io {
                                                path: vertex_path.clone(),
                                                source,
                                            })?;
                                        local_inferred += 1;
                                    }
                                    if edge.weight == 1 {
                                        let info = edge_path_info(edge, x_info, y_info);
                                        let record = stitched_record_from_edge_path_info(
                                            edge, info, error_path,
                                        )?;
                                        if let Some(phantom) = edge.phantom_unitig {
                                            local_phantoms.push((record, phantom));
                                        } else {
                                            let bucket_id = edge_path_info_bucket(
                                                edge,
                                                ranges,
                                                &range_index,
                                                ranges_per_bucket,
                                            )
                                            .ok_or_else(|| {
                                                SerialCollationError::MalformedCoordBucket(
                                                    error_path.to_path_buf(),
                                                )
                                            })?;
                                            range_records
                                                .get_mut(bucket_id)
                                                .ok_or_else(|| {
                                                    SerialCollationError::MalformedCoordBucket(
                                                        error_path.to_path_buf(),
                                                    )
                                                })?
                                                .push(record);
                                        }
                                        local_edge_infos += 1;
                                    }
                                }
                                if let Some(mut out) = vertex_out {
                                    out.flush().map_err(|source| SerialCollationError::Io {
                                        path: vertex_path,
                                        source,
                                    })?;
                                }
                            }

                            for (bucket_id, records) in range_records.iter().enumerate() {
                                edge_writers.write_path_records(bucket_id, records)?;
                            }

                            Ok::<_, SerialCollationError>(Some((
                                local_phantoms,
                                local_unresolved,
                                local_inferred,
                                local_edge_infos,
                            )))
                        })
                        .collect::<Result<Vec<_>, SerialCollationError>>()
                })?;

                for (local_phantoms, local_unresolved, local_inferred, local_edge_infos) in
                    worker_outputs.into_iter().flatten()
                {
                    unresolved_edges += local_unresolved;
                    inferred_vertices += local_inferred;
                    edge_path_infos += local_edge_infos;
                    phantom_records.extend(local_phantoms);
                }
            }
            non_diagonal_elapsed += phase_started.elapsed();

            let phase_started = Instant::now();
            let diagonal_records_per_chunk =
                (1024 * 1024 / std::mem::size_of::<DiscontinuityEdge<K>>()).max(1);
            let diagonal_outputs = expansion_pool.install(|| {
                diagonal_block
                    .par_chunks(diagonal_records_per_chunk)
                    .map(|edges| {
                        let mut records = (0..bucket_count)
                            .map(|_| None::<Vec<StitchedCoordRecord>>)
                            .collect::<Vec<_>>();
                        let mut used_buckets = Vec::new();
                        let mut phantoms = Vec::new();
                        let mut local_unresolved = 0u64;
                        let mut local_emitted = 0u64;
                        for edge in edges {
                            if edge.weight != 1 {
                                continue;
                            }
                            let (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) =
                                (edge.first, edge.second)
                            else {
                                continue;
                            };
                            let record = if K <= 31 {
                                let Some(x_info) = map.get_compact(x.vertex) else {
                                    local_unresolved += 1;
                                    continue;
                                };
                                let y_info =
                                    Self::infer_compact(x_info, x.side, y.side, edge.weight);
                                stitched_record_from_compact_edge_path_info(
                                    edge,
                                    compact_edge_path_info(edge, x_info, y_info),
                                    error_path,
                                )?
                            } else {
                                let Some(x_info) = map.get(x.vertex) else {
                                    local_unresolved += 1;
                                    continue;
                                };
                                let y_info = Self::infer(x_info, x.side, y.side, edge.weight);
                                stitched_record_from_edge_path_info(
                                    edge,
                                    edge_path_info(edge, x_info, y_info),
                                    error_path,
                                )?
                            };
                            if let Some(phantom) = edge.phantom_unitig {
                                phantoms.push((record, phantom));
                            } else {
                                let bucket_id = edge_path_info_bucket(
                                    edge,
                                    ranges,
                                    &range_index,
                                    ranges_per_bucket,
                                )
                                .ok_or_else(|| {
                                    SerialCollationError::MalformedCoordBucket(
                                        error_path.to_path_buf(),
                                    )
                                })?;
                                let bucket = &mut records[bucket_id];
                                if bucket.is_none() {
                                    *bucket = Some(Vec::new());
                                    used_buckets.push(bucket_id);
                                }
                                bucket
                                    .as_mut()
                                    .expect("edge bucket was just initialized")
                                    .push(record);
                            }
                            local_emitted += 1;
                        }
                        for bucket_id in used_buckets {
                            edge_writers.write_path_records(
                                bucket_id,
                                records[bucket_id]
                                    .as_deref()
                                    .expect("used edge bucket is initialized"),
                            )?;
                        }
                        Ok::<_, SerialCollationError>((phantoms, local_unresolved, local_emitted))
                    })
                    .collect::<Result<Vec<_>, SerialCollationError>>()
            })?;
            for (phantoms, local_unresolved, local_emitted) in diagonal_outputs {
                phantom_records.extend(phantoms);
                unresolved_edges += local_unresolved;
                edge_path_infos += local_emitted;
            }

            let phi_block = &matrix.blocks[matrix.block_index(0, partition)];
            if phi_block.edges != 0 {
                let bytes = match fs::read(&phi_block.path) {
                    Ok(bytes) => bytes,
                    Err(source) if source.kind() == std::io::ErrorKind::NotFound => Vec::new(),
                    Err(source) => {
                        return Err(SerialCollationError::Io {
                            path: phi_block.path.clone(),
                            source,
                        });
                    }
                };
                if bytes.len() + phi_block.buffer.len() != phi_block.edges * phi_block.record_len {
                    return Err(SerialCollationError::MalformedCoordBucket(
                        phi_block.path.clone(),
                    ));
                }
                let records_per_chunk = (1024 * 1024 / phi_block.record_len).max(1);
                let bytes_per_chunk = records_per_chunk * phi_block.record_len;
                let phi_outputs = expansion_pool.install(|| {
                    bytes
                        .par_chunks(bytes_per_chunk)
                        .chain(phi_block.buffer.par_chunks(bytes_per_chunk))
                        .map(|chunk| {
                            let mut records = (0..bucket_count)
                                .map(|_| None::<Vec<StitchedCoordRecord>>)
                                .collect::<Vec<_>>();
                            let mut used_buckets = Vec::new();
                            let mut phantoms = Vec::new();
                            let mut local_unresolved = 0u64;
                            let mut local_emitted = 0u64;
                            for encoded in chunk.chunks_exact(phi_block.record_len) {
                                let edge = decode_discontinuity_edge::<K>(encoded);
                                if edge.weight != 1 {
                                    continue;
                                }
                                let (MatrixEndpoint::Phi, MatrixEndpoint::Vertex(v)) =
                                    (edge.first, edge.second)
                                else {
                                    continue;
                                };
                                let record = if K <= 31 {
                                    let Some(v_info) = map.get_compact(v.vertex) else {
                                        local_unresolved += 1;
                                        continue;
                                    };
                                    stitched_record_from_compact_edge_path_info(
                                        &edge,
                                        compact_phi_edge_path_info(&edge, v_info),
                                        error_path,
                                    )?
                                } else {
                                    let Some(v_info) = map.get(v.vertex) else {
                                        local_unresolved += 1;
                                        continue;
                                    };
                                    stitched_record_from_edge_path_info(
                                        &edge,
                                        phi_edge_path_info(&edge, v_info),
                                        error_path,
                                    )?
                                };
                                if let Some(phantom) = edge.phantom_unitig {
                                    phantoms.push((record, phantom));
                                } else {
                                    let bucket_id = edge_path_info_bucket(
                                        &edge,
                                        ranges,
                                        &range_index,
                                        ranges_per_bucket,
                                    )
                                    .ok_or_else(|| {
                                        SerialCollationError::MalformedCoordBucket(
                                            error_path.to_path_buf(),
                                        )
                                    })?;
                                    let bucket = &mut records[bucket_id];
                                    if bucket.is_none() {
                                        *bucket = Some(Vec::new());
                                        used_buckets.push(bucket_id);
                                    }
                                    bucket
                                        .as_mut()
                                        .expect("edge bucket was just initialized")
                                        .push(record);
                                }
                                local_emitted += 1;
                            }
                            for bucket_id in used_buckets {
                                edge_writers.write_path_records(
                                    bucket_id,
                                    records[bucket_id]
                                        .as_deref()
                                        .expect("used edge bucket is initialized"),
                                )?;
                            }
                            Ok::<_, SerialCollationError>((
                                phantoms,
                                local_unresolved,
                                local_emitted,
                            ))
                        })
                        .collect::<Result<Vec<_>, SerialCollationError>>()
                })?;
                for (phantoms, local_unresolved, local_emitted) in phi_outputs {
                    phantom_records.extend(phantoms);
                    unresolved_edges += local_unresolved;
                    edge_path_infos += local_emitted;
                }
            }
            original_edge_elapsed += phase_started.elapsed();
        }

        drop(vertex_writers);
        eprintln!(
            "cuttlefish3-rs: blocked expansion detail: row load/decode {:.3}s, propagation/emission {:.3}s",
            row_load_elapsed.as_secs_f64(),
            expansion_work_started
                .elapsed()
                .saturating_sub(row_load_elapsed)
                .as_secs_f64()
        );
        eprintln!(
            "cuttlefish3-rs: path-info load phases: clear {:.3}s, read {:.3}s, decode/insert {:.3}s",
            path_info_clear_elapsed.as_secs_f64(),
            path_info_read_elapsed.as_secs_f64(),
            path_info_insert_elapsed.as_secs_f64(),
        );
        eprintln!(
            "cuttlefish3-rs: blocked expansion phases: path-info load {:.3}s, compressed diagonal {:.3}s, non-diagonal {:.3}s, original diagonal/phi {:.3}s",
            path_info_load_elapsed.as_secs_f64(),
            compressed_diagonal_elapsed.as_secs_f64(),
            non_diagonal_elapsed.as_secs_f64(),
            original_edge_elapsed.as_secs_f64(),
        );
        let mut record_manifest = edge_writers.finish(&expansion_pool)?;
        record_manifest.sort_by(|left, right| {
            left.bucket_id
                .cmp(&right.bucket_id)
                .then_with(|| left.path.cmp(&right.path))
        });

        Ok(RangeBucketedExpansion {
            stats: SerialExpansionStats {
                seed_vertices,
                inferred_vertices,
                edge_path_infos,
                unresolved_edges,
            },
            records: Vec::new(),
            record_manifest,
            phantom_records,
        })
    }

    fn expand_cpp_ordered_impl<
        const K: usize,
        const COLLECT_VERTICES: bool,
        const SORT_DEDUP_EDGES: bool,
    >(
        contraction: &FullSerialDiscontinuityContraction<K>,
    ) -> SerialExpansion<K> {
        let partition_count = contraction.vertex_partitions + 1;
        let mut vertex_buckets = vec![Vec::<VertexPathInfo<K>>::new(); partition_count];
        for meta in &contraction.meta_vertices {
            if meta.partition < partition_count {
                vertex_buckets[meta.partition].push(VertexPathInfo {
                    vertex: meta.vertex,
                    info: PathInfo {
                        path_id: meta.vertex,
                        rank: meta.weight,
                        exit_side: meta.entry_side,
                        is_cycle: meta.is_cycle,
                    },
                });
            }
        }
        let seed_vertices = vertex_buckets.iter().map(Vec::len).sum::<usize>() as u64;

        let mut diagonal_by_partition = contraction.compressed_diagonal_edges.clone();
        diagonal_by_partition.resize_with(partition_count, Vec::new);
        if diagonal_by_partition.iter().all(Vec::is_empty) {
            for partition in &contraction.partitions {
                if partition.partition < partition_count {
                    diagonal_by_partition[partition.partition]
                        .extend(partition.compressed_diagonal_edges.iter().cloned());
                }
            }
        }

        let mut non_diagonal_by_row =
            vec![Vec::<(usize, DiscontinuityEdge<K>)>::new(); partition_count];
        let mut diagonal_by_row = vec![Vec::<DiscontinuityEdge<K>>::new(); partition_count];
        let mut phi_by_partition = vec![Vec::<DiscontinuityEdge<K>>::new(); partition_count];
        for edge in &contraction.expansion_edges {
            let (row, col) =
                edge_matrix_row_col(contraction.vertex_partitions, edge.first, edge.second);
            if row == 0 {
                if col < partition_count {
                    phi_by_partition[col].push(edge.clone());
                }
            } else if row == col {
                if row < partition_count {
                    diagonal_by_row[row].push(edge.clone());
                }
            } else if row < partition_count {
                non_diagonal_by_row[row].push((col, edge.clone()));
            }
        }
        for row in 1..partition_count {
            non_diagonal_by_row[row].sort_by_key(|(col, edge)| {
                (
                    *col,
                    endpoint_sort_key(edge.first),
                    endpoint_sort_key(edge.second),
                    edge.weight,
                    edge.unitig_index,
                )
            });
            diagonal_by_row[row].sort_by_key(|edge| {
                (
                    endpoint_sort_key(edge.first),
                    endpoint_sort_key(edge.second),
                    edge.weight,
                    edge.unitig_index,
                )
            });
            phi_by_partition[row].sort_by_key(|edge| {
                (
                    endpoint_sort_key(edge.first),
                    endpoint_sort_key(edge.second),
                    edge.weight,
                    edge.unitig_index,
                )
            });
        }
        let mut edges = Vec::new();
        let mut unresolved_edges = 0u64;
        let mut inferred_vertices = 0u64;
        let mut map = FastHashMap::<Kmer<K>, PathInfo<K>>::default();

        for partition in 1..partition_count {
            map.clear();
            for vertex_info in &vertex_buckets[partition] {
                map.insert(vertex_info.vertex, vertex_info.info);
            }

            for edge in diagonal_by_partition[partition].iter().rev() {
                let (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) =
                    (edge.first, edge.second)
                else {
                    continue;
                };

                if let Some(&y_info) = map.get(&y.vertex) {
                    let x_info = Self::infer(y_info, y.side, x.side, edge.weight);
                    if x_info.rank > 0 && !map.contains_key(&x.vertex) {
                        map.insert(x.vertex, x_info);
                        inferred_vertices += 1;
                    }
                } else if let Some(&x_info) = map.get(&x.vertex) {
                    let y_info = Self::infer(x_info, x.side, y.side, edge.weight);
                    if y_info.rank > 0 && !map.contains_key(&y.vertex) {
                        map.insert(y.vertex, y_info);
                        inferred_vertices += 1;
                    }
                } else {
                    unresolved_edges += 1;
                }
            }

            for (col, edge) in &non_diagonal_by_row[partition] {
                let (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) =
                    (edge.first, edge.second)
                else {
                    continue;
                };
                let Some(&x_info) = map.get(&x.vertex) else {
                    unresolved_edges += 1;
                    continue;
                };
                let y_info = Self::infer(x_info, x.side, y.side, edge.weight);
                if y_info.rank > 0 {
                    vertex_buckets[*col].push(VertexPathInfo {
                        vertex: y.vertex,
                        info: y_info,
                    });
                    inferred_vertices += 1;
                }
                if edge.weight == 1 {
                    edges.push(EdgePathInfo {
                        unitig_index: edge.unitig_index,
                        phantom_unitig: edge.phantom_unitig,
                        info: edge_path_info(edge, x_info, y_info),
                    });
                }
            }

            for edge in &diagonal_by_row[partition] {
                if edge.weight != 1 {
                    continue;
                }
                let (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) =
                    (edge.first, edge.second)
                else {
                    continue;
                };
                let Some(&x_info) = map.get(&x.vertex) else {
                    unresolved_edges += 1;
                    continue;
                };
                let y_info = Self::infer(x_info, x.side, y.side, edge.weight);
                edges.push(EdgePathInfo {
                    unitig_index: edge.unitig_index,
                    phantom_unitig: edge.phantom_unitig,
                    info: edge_path_info(edge, x_info, y_info),
                });
            }

            for edge in &phi_by_partition[partition] {
                if edge.weight != 1 {
                    continue;
                }
                let (MatrixEndpoint::Phi, MatrixEndpoint::Vertex(v)) = (edge.first, edge.second)
                else {
                    continue;
                };
                if let Some(&v_info) = map.get(&v.vertex) {
                    edges.push(EdgePathInfo {
                        unitig_index: edge.unitig_index,
                        phantom_unitig: edge.phantom_unitig,
                        info: phi_edge_path_info(edge, v_info),
                    });
                } else {
                    unresolved_edges += 1;
                }
            }
        }

        if SORT_DEDUP_EDGES {
            edges.sort_by_key(|edge| {
                (
                    edge.unitig_index,
                    edge.info.path_id.as_u128(),
                    edge.info.rank,
                    edge.info.exit_side as u8,
                )
            });
            edges.dedup();
        }

        let vertices = if COLLECT_VERTICES {
            let mut vertices = vertex_buckets
                .into_iter()
                .flatten()
                .collect::<Vec<VertexPathInfo<K>>>();
            vertices.sort_by_key(|entry| {
                (
                    entry.vertex.as_u128(),
                    entry.info.path_id.as_u128(),
                    entry.info.rank,
                    entry.info.exit_side as u8,
                )
            });
            vertices.dedup();
            vertices
        } else {
            Vec::new()
        };

        SerialExpansion {
            stats: SerialExpansionStats {
                seed_vertices,
                inferred_vertices,
                edge_path_infos: edges.len() as u64,
                unresolved_edges,
            },
            vertices,
            edges,
        }
    }

    fn expand_impl<const K: usize>(
        contraction: &FullSerialDiscontinuityContraction<K>,
        include_original_edges: bool,
    ) -> SerialExpansion<K> {
        let mut vertex_info = FastHashMap::<Kmer<K>, PathInfo<K>>::with_capacity_and_hasher(
            contraction.meta_vertices.len(),
            FastBuildHasher::default(),
        );
        let mut edges = Vec::new();
        let mut unresolved_edges = 0u64;

        for meta in &contraction.meta_vertices {
            vertex_info.insert(
                meta.vertex,
                PathInfo {
                    path_id: meta.vertex,
                    rank: meta.weight,
                    exit_side: meta.entry_side,
                    is_cycle: meta.is_cycle,
                },
            );
        }
        let seed_vertices = vertex_info.len() as u64;

        let mut expansion_edges =
            if include_original_edges && !contraction.expansion_edges.is_empty() {
                contraction.expansion_edges.clone()
            } else {
                contraction.final_edges.clone()
            };
        expansion_edges.sort_by_key(|edge| {
            (
                endpoint_sort_key(edge.first),
                endpoint_sort_key(edge.second),
                edge.weight,
                edge.unitig_index,
            )
        });

        let mut changed = true;
        while changed {
            changed = false;
            for edge in &expansion_edges {
                match (edge.first, edge.second) {
                    (MatrixEndpoint::Phi, MatrixEndpoint::Vertex(v))
                    | (MatrixEndpoint::Vertex(v), MatrixEndpoint::Phi) => {
                        if edge.weight == 1 {
                            if let Some(&info) = vertex_info.get(&v.vertex) {
                                edges.push(EdgePathInfo {
                                    unitig_index: edge.unitig_index,
                                    phantom_unitig: edge.phantom_unitig,
                                    info: phi_edge_path_info(edge, info),
                                });
                            }
                        }
                    }
                    (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) => {
                        let x_info = vertex_info.get(&x.vertex).copied();
                        let y_info = vertex_info.get(&y.vertex).copied();
                        match (x_info, y_info) {
                            (Some(x_info), None) => {
                                let inferred = Self::infer(x_info, x.side, y.side, edge.weight);
                                if inferred.rank > 0 {
                                    vertex_info.insert(y.vertex, inferred);
                                    changed = true;
                                }
                            }
                            (None, Some(y_info)) => {
                                let inferred = Self::infer(y_info, y.side, x.side, edge.weight);
                                if inferred.rank > 0 {
                                    vertex_info.insert(x.vertex, inferred);
                                    changed = true;
                                }
                            }
                            (Some(x_info), Some(y_info)) if edge.weight == 1 => {
                                if x_info.path_id == y_info.path_id {
                                    edges.push(EdgePathInfo {
                                        unitig_index: edge.unitig_index,
                                        phantom_unitig: edge.phantom_unitig,
                                        info: edge_path_info(edge, x_info, y_info),
                                    });
                                }
                            }
                            _ => {}
                        }
                    }
                    _ => {}
                }
            }
        }

        for edge in &expansion_edges {
            match (edge.first, edge.second) {
                (MatrixEndpoint::Phi, MatrixEndpoint::Vertex(v))
                | (MatrixEndpoint::Vertex(v), MatrixEndpoint::Phi) => {
                    if !vertex_info.contains_key(&v.vertex) {
                        unresolved_edges += 1;
                    }
                }
                (MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y)) => {
                    if !vertex_info.contains_key(&x.vertex) || !vertex_info.contains_key(&y.vertex)
                    {
                        unresolved_edges += 1;
                    }
                }
                _ => {}
            }
        }

        edges.sort_by_key(|edge| {
            (
                edge.unitig_index,
                edge.info.path_id.as_u128(),
                edge.info.rank,
                edge.info.exit_side as u8,
            )
        });
        edges.dedup();

        let mut vertices = vertex_info
            .into_iter()
            .map(|(vertex, info)| VertexPathInfo { vertex, info })
            .collect::<Vec<_>>();
        vertices.sort_by_key(|entry| {
            (
                entry.vertex.as_u128(),
                entry.info.path_id.as_u128(),
                entry.info.rank,
                entry.info.exit_side as u8,
            )
        });

        SerialExpansion {
            stats: SerialExpansionStats {
                seed_vertices,
                inferred_vertices: (vertices.len() as u64).saturating_sub(seed_vertices),
                edge_path_infos: edges.len() as u64,
                unresolved_edges,
            },
            vertices,
            edges,
        }
    }
}

pub struct SerialUncoloredCollator;

impl SerialUncoloredCollator {
    pub fn collate<const K: usize>(
        inputs: &DiscontinuityInputs<K>,
        expansion: &SerialExpansion<K>,
    ) -> SerialCollation {
        Self::collate_with_threads(inputs, expansion, 1)
    }

    pub fn collate_with_threads<const K: usize>(
        inputs: &DiscontinuityInputs<K>,
        expansion: &SerialExpansion<K>,
        threads: usize,
    ) -> SerialCollation {
        eprintln!(
            "cuttlefish3-rs: collating {} path edge info record(s) against {} local unitig(s)",
            expansion.edges.len(),
            inputs.unitigs.len()
        );
        let mut records = Vec::<CollationRecord<K>>::new();
        let mut missing_unitig_labels = 0u64;
        let mut path_unitigs = vec![false; inputs.unitigs.len()];

        let collect_started = Instant::now();
        for edge in &expansion.edges {
            let Some(label) = inputs.try_label(edge.unitig_index) else {
                missing_unitig_labels += 1;
                continue;
            };
            path_unitigs[edge.unitig_index] = true;
            records.push(CollationRecord {
                info: edge.info,
                label: label.to_vec(),
            });
        }
        eprintln!(
            "cuttlefish3-rs: collation collected path records in {:.3}s",
            collect_started.elapsed().as_secs_f64()
        );

        let record_sort_started = Instant::now();
        records.sort_by_key(|record| {
            (
                record.info.path_id.as_u128(),
                record.info.rank,
                record.info.exit_side as u8,
            )
        });
        eprintln!(
            "cuttlefish3-rs: collation sorted path records in {:.3}s",
            record_sort_started.elapsed().as_secs_f64()
        );

        let expanded_started = Instant::now();
        let mut expanded_unitigs = Vec::new();
        let mut start = 0;
        while start < records.len() {
            let path_id = records[start].info.path_id;
            let is_cycle = records[start].info.is_cycle;
            let mut end = start + 1;
            while end < records.len() && records[end].info.path_id == path_id {
                end += 1;
            }

            let mut label = Vec::new();
            if end - start == 2
                && !is_cycle
                && records[start].info.rank == 0
                && records[start + 1].info.rank == 0
            {
                append_or_init::<K>(
                    &mut label,
                    oriented_label(
                        &records[start].label,
                        records[start].info.exit_side == Side::Front,
                    ),
                );
                append_or_init::<K>(
                    &mut label,
                    oriented_label(
                        &records[start + 1].label,
                        records[start + 1].info.exit_side != Side::Front,
                    ),
                );
            } else {
                for record in &records[start..end] {
                    append_or_init::<K>(
                        &mut label,
                        oriented_label(&record.label, record.info.exit_side == Side::Front),
                    );
                }
            }

            if is_cycle && label.len() >= K {
                label.truncate(label.len() - (K - 1));
            }

            expanded_unitigs.push(canonical_label(label));
            start = end;
        }
        expanded_unitigs.sort_unstable();
        expanded_unitigs.dedup();
        eprintln!(
            "cuttlefish3-rs: collation built expanded path unitigs in {:.3}s",
            expanded_started.elapsed().as_secs_f64()
        );
        eprintln!(
            "cuttlefish3-rs: collated {} expanded path unitig(s)",
            expanded_unitigs.len()
        );

        let stitch_started = Instant::now();
        let mut stitched_unitigs = stitch_discontinuity_paths::<K>(inputs, &[], threads);
        stitched_unitigs.sort_unstable();
        stitched_unitigs.dedup();
        let stitched_discontinuity_unitigs = stitched_unitigs.len() as u64;
        eprintln!(
            "cuttlefish3-rs: collation stitched discontinuity paths in {:.3}s",
            stitch_started.elapsed().as_secs_f64()
        );
        eprintln!(
            "cuttlefish3-rs: added {} stitched discontinuity unitig(s)",
            stitched_discontinuity_unitigs
        );

        let suppress_started = Instant::now();
        let suppressed_expanded =
            suppress_labels_contained_in_sources(&mut expanded_unitigs, &stitched_unitigs, threads);
        eprintln!(
            "cuttlefish3-rs: collation suppressed contained expanded labels in {:.3}s",
            suppress_started.elapsed().as_secs_f64()
        );
        if suppressed_expanded > 0 {
            eprintln!(
                "cuttlefish3-rs: suppressed {} expanded path unitig(s) contained in stitched paths",
                suppressed_expanded
            );
        }

        let mut unitigs = Vec::with_capacity(
            expanded_unitigs.len() + stitched_unitigs.len() + inputs.unitigs.len(),
        );
        let mut unitig_origins = Vec::new();
        let merge_started = Instant::now();
        for unitig in expanded_unitigs {
            unitigs.push(unitig);
            unitig_origins.push("expanded");
        }
        for unitig in stitched_unitigs {
            unitigs.push(unitig);
            unitig_origins.push("stitched");
        }

        let mut direct_local_unitigs = 0u64;
        for (index, unitig) in inputs.unitigs.iter().enumerate() {
            if !path_unitigs[index] && unitig.left_exit().is_none() && unitig.right_exit().is_none()
            {
                unitigs.push(canonical_label(unitig.label(inputs).to_vec()));
                unitig_origins.push("direct");
                direct_local_unitigs += 1;
            }
        }
        eprintln!(
            "cuttlefish3-rs: collation merged direct local unitigs in {:.3}s",
            merge_started.elapsed().as_secs_f64()
        );
        eprintln!(
            "cuttlefish3-rs: added {} direct local unitig(s)",
            direct_local_unitigs
        );

        if let Ok(path) = std::env::var("CF3_RS_DUMP_COLLATION_CANDIDATES") {
            if !path.is_empty() {
                if let Err(err) = dump_collation_candidates(
                    std::path::Path::new(&path),
                    &unitigs,
                    &unitig_origins,
                ) {
                    eprintln!("cuttlefish3-rs: failed to dump collation candidates: {err}");
                }
            }
        }

        let final_sort_started = Instant::now();
        unitigs = bucketed_maximal_unitig_reduce(unitigs, threads);
        eprintln!(
            "cuttlefish3-rs: collation bucket reduce completed in {:.3}s",
            final_sort_started.elapsed().as_secs_f64()
        );
        if std::env::var_os("CF3_RS_RETAIN_MAXIMAL_LABELS").is_some() {
            retain_maximal_labels_indexed_with_threads(&mut unitigs, threads);
        }

        SerialCollation {
            stats: SerialCollationStats {
                input_path_infos: expansion.edges.len() as u64,
                emitted_unitigs: unitigs.len() as u64,
                emitted_bases: unitigs.iter().map(|unitig| unitig.len() as u64).sum(),
                missing_unitig_labels,
                direct_local_unitigs,
                stitched_discontinuity_unitigs,
            },
            unitigs,
        }
    }

    pub fn collate_path_info_only_with_threads<const K: usize>(
        inputs: &DiscontinuityInputs<K>,
        expansion: &SerialExpansion<K>,
        threads: usize,
    ) -> SerialCollation {
        eprintln!(
            "cuttlefish3-rs: collating {} path edge info record(s) against {} local unitig(s)",
            expansion.edges.len(),
            inputs.unitigs.len()
        );
        let mut records = Vec::<CollationRecord<K>>::with_capacity(expansion.edges.len());
        let mut missing_unitig_labels = 0u64;
        let mut path_unitigs = vec![false; inputs.unitigs.len()];

        let collect_started = Instant::now();
        for edge in &expansion.edges {
            let Some(label) = inputs.try_label(edge.unitig_index) else {
                missing_unitig_labels += 1;
                continue;
            };
            path_unitigs[edge.unitig_index] = true;
            records.push(CollationRecord {
                info: edge.info,
                label: label.to_vec(),
            });
        }
        eprintln!(
            "cuttlefish3-rs: path-info collation collected records in {:.3}s",
            collect_started.elapsed().as_secs_f64()
        );

        let record_sort_started = Instant::now();
        records.sort_by_key(|record| {
            (
                record.info.path_id.as_u128(),
                record.info.rank,
                record.info.exit_side as u8,
            )
        });
        eprintln!(
            "cuttlefish3-rs: path-info collation sorted records in {:.3}s",
            record_sort_started.elapsed().as_secs_f64()
        );

        let build_started = Instant::now();
        let mut unitigs = Vec::new();
        let mut start = 0;
        while start < records.len() {
            let path_id = records[start].info.path_id;
            let is_cycle = records[start].info.is_cycle;
            let mut end = start + 1;
            while end < records.len() && records[end].info.path_id == path_id {
                end += 1;
            }

            let mut label = Vec::new();
            if end - start == 2
                && !is_cycle
                && records[start].info.rank == 0
                && records[start + 1].info.rank == 0
            {
                append_or_init::<K>(
                    &mut label,
                    oriented_label(
                        &records[start].label,
                        records[start].info.exit_side == Side::Front,
                    ),
                );
                append_or_init::<K>(
                    &mut label,
                    oriented_label(
                        &records[start + 1].label,
                        records[start + 1].info.exit_side != Side::Front,
                    ),
                );
            } else {
                for record in &records[start..end] {
                    append_or_init::<K>(
                        &mut label,
                        oriented_label(&record.label, record.info.exit_side == Side::Front),
                    );
                }
            }

            if label.len() >= K {
                if is_cycle {
                    label.truncate(label.len().saturating_sub(K - 1));
                    unitigs.push(normalize_stitched_cycle::<K>(&label));
                } else {
                    unitigs.push(canonical_label(label));
                }
            }
            start = end;
        }
        let stitched_discontinuity_unitigs = unitigs.len() as u64;
        eprintln!(
            "cuttlefish3-rs: path-info collation built path unitigs in {:.3}s",
            build_started.elapsed().as_secs_f64()
        );
        eprintln!(
            "cuttlefish3-rs: added {} path-info discontinuity unitig(s)",
            stitched_discontinuity_unitigs
        );

        let merge_started = Instant::now();
        let mut direct_local_unitigs = 0u64;
        for (index, unitig) in inputs.unitigs.iter().enumerate() {
            if !path_unitigs[index] && unitig.left_exit().is_none() && unitig.right_exit().is_none()
            {
                unitigs.push(canonical_label(unitig.label(inputs).to_vec()));
                direct_local_unitigs += 1;
            }
        }
        eprintln!(
            "cuttlefish3-rs: path-info collation merged direct local unitigs in {:.3}s",
            merge_started.elapsed().as_secs_f64()
        );
        eprintln!(
            "cuttlefish3-rs: added {} direct local unitig(s)",
            direct_local_unitigs
        );

        let final_sort_started = Instant::now();
        unitigs = bucketed_maximal_unitig_reduce(unitigs, threads);
        eprintln!(
            "cuttlefish3-rs: path-info collation bucket reduce completed in {:.3}s",
            final_sort_started.elapsed().as_secs_f64()
        );

        SerialCollation {
            stats: SerialCollationStats {
                input_path_infos: expansion.edges.len() as u64,
                emitted_unitigs: unitigs.len() as u64,
                emitted_bases: unitigs.iter().map(|unitig| unitig.len() as u64).sum(),
                missing_unitig_labels,
                direct_local_unitigs,
                stitched_discontinuity_unitigs,
            },
            unitigs,
        }
    }

    pub fn collate_stitched_with_threads<const K: usize>(
        inputs: &DiscontinuityInputs<K>,
        threads: usize,
    ) -> SerialCollation {
        Self::collate_stitched_with_threads_impl(inputs, threads, None)
            .expect("in-memory stitched collation cannot fail")
    }

    pub fn collate_stitched_with_threads_in_dir<const K: usize>(
        inputs: &DiscontinuityInputs<K>,
        threads: usize,
        coord_dir: &Path,
    ) -> Result<SerialCollation, SerialCollationError> {
        Self::collate_stitched_with_threads_impl(inputs, threads, Some(coord_dir))
    }

    pub fn collate_stitched_to_fasta_with_threads_in_dir<const K: usize>(
        inputs: &DiscontinuityInputs<K>,
        threads: usize,
        coord_dir: &Path,
        final_dir: &Path,
        output_path: &Path,
    ) -> Result<SerialCollationStats, SerialCollationError> {
        eprintln!(
            "cuttlefish3-rs: collating {} local unitig(s) by discontinuity-end stitching",
            inputs.unitigs.len()
        );

        let spill_started = Instant::now();
        let mut final_buckets = FinalUnitigBucketWriters::create(final_dir, threads)?;
        let stitched_discontinuity_unitigs = stitch_discontinuity_paths_to_final_buckets::<K>(
            inputs,
            &[],
            threads,
            coord_dir,
            &mut final_buckets,
        )?;

        let mut direct_local_unitigs = 0u64;
        for unitig in &inputs.unitigs {
            if unitig.left_exit().is_none() && unitig.right_exit().is_none() {
                let label = canonical_label(unitig.label(inputs).to_vec());
                final_buckets.write_label(&label)?;
                direct_local_unitigs += 1;
            }
        }
        let manifest = final_buckets.finish()?;
        eprintln!(
            "cuttlefish3-rs: collation stitched and spilled final-label candidates in {:.3}s",
            spill_started.elapsed().as_secs_f64()
        );
        eprintln!(
            "cuttlefish3-rs: added {} stitched discontinuity unitig(s)",
            stitched_discontinuity_unitigs
        );
        eprintln!(
            "cuttlefish3-rs: added {} direct local unitig(s)",
            direct_local_unitigs
        );

        let reduce_started = Instant::now();
        let (emitted_unitigs, emitted_bases) =
            reduce_final_unitig_buckets_to_fasta(&manifest, output_path)?;
        eprintln!(
            "cuttlefish3-rs: collation external bucket reduce and FASTA write completed in {:.3}s",
            reduce_started.elapsed().as_secs_f64()
        );

        Ok(SerialCollationStats {
            input_path_infos: 0,
            emitted_unitigs,
            emitted_bases,
            missing_unitig_labels: 0,
            direct_local_unitigs,
            stitched_discontinuity_unitigs,
        })
    }

    pub fn collate_external_stitched_to_fasta_with_threads_in_dir<const K: usize>(
        inputs: &mut ExternalDiscontinuityInputs<K>,
        threads: usize,
        coord_dir: &Path,
        final_dir: &Path,
        output_path: &Path,
    ) -> Result<SerialCollationStats, SerialCollationError> {
        eprintln!(
            "cuttlefish3-rs: collating {} local unitig(s) by external discontinuity-end stitching",
            inputs.unitig_count()
        );

        let spill_started = Instant::now();
        let colored = inputs.color_runs.is_some()
            || inputs
                .local_unitig_buckets
                .as_ref()
                .is_some_and(|buckets| buckets.iter().any(|bucket| bucket.color_runs.is_some()));
        let mut final_buckets = FinalUnitigBucketWriters::create_direct(output_path, colored)?;
        let mut direct_local_unitigs = 0u64;
        if let Some((path, records, bases)) = inputs.trivial_fasta.take() {
            final_buckets.append_direct_fasta_file(&path, records, bases)?;
            remove_serial_file(&path)?;
            direct_local_unitigs = records;
        }
        let use_cpp_path_info = std::env::var_os("CF3_RS_ENDPOINT_STITCH").is_none();
        let mut direct_local_complete = false;
        let stitched_discontinuity_unitigs = if use_cpp_path_info {
            let result = collate_external_cpp_path_info_to_final_buckets::<K>(
                inputs,
                threads,
                coord_dir,
                &mut final_buckets,
            )?;
            direct_local_unitigs += result.direct_local_unitigs;
            direct_local_complete = result.direct_local_unitigs_complete;
            result.stitched_unitigs
        } else {
            stitch_external_discontinuity_paths_to_final_buckets::<K>(
                inputs,
                threads,
                coord_dir,
                &mut final_buckets,
            )?
        };

        if !direct_local_complete {
            let reader = ExternalDiscontinuityReader::open(inputs)?;
            let mut color_reader = None;
            let mut range_index = 0usize;
            let mut scratch = Vec::new();
            for (unitig_index, unitig) in reader.iter()?.enumerate() {
                let unitig = unitig?;
                if inputs
                    .ranges
                    .get(range_index)
                    .is_some_and(|range| range.start_unitig == unitig_index)
                {
                    color_reader = inputs
                        .color_runs
                        .as_ref()
                        .map(|sidecar| sidecar.reader_at(inputs.ranges[range_index].color_start))
                        .transpose()?;
                    range_index += 1;
                }
                let colors = color_reader
                    .as_mut()
                    .map(|reader| reader.read_next())
                    .transpose()?;
                if unitig.left_exit().is_none() && unitig.right_exit().is_none() {
                    reader.read_label(&unitig, &mut scratch)?;
                    let reverse = reverse_complement_is_less(&scratch);
                    let label = canonical_label(scratch.clone());
                    if let Some(mut colors) = colors {
                        if reverse {
                            colors = reverse_color_runs(
                                &colors,
                                (unitig.label_len as usize - K + 1) as u32,
                            );
                        }
                        final_buckets.write_colored_label(&label, &colors)?;
                    } else {
                        final_buckets.write_label(&label)?;
                    }
                    direct_local_unitigs += 1;
                }
            }
        }
        if !keep_intermediates() {
            remove_serial_file(&inputs.unitig_path)?;
            remove_serial_file(&inputs.label_path)?;
            if let Some(sidecar) = inputs.color_runs.as_ref() {
                remove_serial_file(&sidecar.run_path)?;
            }
        }
        let manifest = final_buckets.finish()?;
        eprintln!(
            "cuttlefish3-rs: external collation stitched and spilled final-label candidates in {:.3}s",
            spill_started.elapsed().as_secs_f64()
        );
        eprintln!(
            "cuttlefish3-rs: added {} stitched discontinuity unitig(s)",
            stitched_discontinuity_unitigs
        );
        eprintln!(
            "cuttlefish3-rs: added {} direct local unitig(s)",
            direct_local_unitigs
        );

        let reduce_started = Instant::now();
        let (emitted_unitigs, emitted_bases) =
            reduce_final_unitig_buckets_to_fasta(&manifest, output_path)?;
        if !keep_intermediates() {
            let _ = fs::remove_dir(coord_dir);
            let _ = fs::remove_dir(final_dir);
        }
        eprintln!(
            "cuttlefish3-rs: collation external bucket reduce and FASTA write completed in {:.3}s",
            reduce_started.elapsed().as_secs_f64()
        );

        Ok(SerialCollationStats {
            input_path_infos: 0,
            emitted_unitigs,
            emitted_bases,
            missing_unitig_labels: 0,
            direct_local_unitigs,
            stitched_discontinuity_unitigs,
        })
    }

    fn collate_stitched_with_threads_impl<const K: usize>(
        inputs: &DiscontinuityInputs<K>,
        threads: usize,
        coord_dir: Option<&Path>,
    ) -> Result<SerialCollation, SerialCollationError> {
        eprintln!(
            "cuttlefish3-rs: collating {} local unitig(s) by discontinuity-end stitching",
            inputs.unitigs.len()
        );

        let stitch_started = Instant::now();
        let mut unitigs = if let Some(coord_dir) = coord_dir {
            stitch_discontinuity_paths_with_coord_dir::<K>(inputs, &[], threads, coord_dir)?
        } else {
            stitch_discontinuity_paths::<K>(inputs, &[], threads)
        };
        let stitched_discontinuity_unitigs = unitigs.len() as u64;
        eprintln!(
            "cuttlefish3-rs: collation stitched discontinuity paths in {:.3}s",
            stitch_started.elapsed().as_secs_f64()
        );
        eprintln!(
            "cuttlefish3-rs: added {} stitched discontinuity unitig(s)",
            stitched_discontinuity_unitigs
        );

        let merge_started = Instant::now();
        let mut direct_local_unitigs = 0u64;
        for unitig in &inputs.unitigs {
            if unitig.left_exit().is_none() && unitig.right_exit().is_none() {
                unitigs.push(canonical_label(unitig.label(inputs).to_vec()));
                direct_local_unitigs += 1;
            }
        }
        eprintln!(
            "cuttlefish3-rs: collation merged direct local unitigs in {:.3}s",
            merge_started.elapsed().as_secs_f64()
        );
        eprintln!(
            "cuttlefish3-rs: added {} direct local unitig(s)",
            direct_local_unitigs
        );

        let final_sort_started = Instant::now();
        unitigs = bucketed_maximal_unitig_reduce(unitigs, threads);
        eprintln!(
            "cuttlefish3-rs: collation bucket reduce completed in {:.3}s",
            final_sort_started.elapsed().as_secs_f64()
        );
        if std::env::var_os("CF3_RS_RETAIN_MAXIMAL_LABELS").is_some() {
            retain_maximal_labels_indexed_with_threads(&mut unitigs, threads);
        }

        Ok(SerialCollation {
            stats: SerialCollationStats {
                input_path_infos: 0,
                emitted_unitigs: unitigs.len() as u64,
                emitted_bases: unitigs.iter().map(|unitig| unitig.len() as u64).sum(),
                missing_unitig_labels: 0,
                direct_local_unitigs,
                stitched_discontinuity_unitigs,
            },
            unitigs,
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CollationRecord<const K: usize> {
    info: PathInfo<K>,
    label: Vec<u8>,
}

fn append_or_init<const K: usize>(label: &mut Vec<u8>, next: Vec<u8>) {
    if label.is_empty() {
        label.extend_from_slice(&next);
    } else {
        label.extend_from_slice(&next[K..]);
    }
}

fn oriented_label(label: &[u8], reverse_complement: bool) -> Vec<u8> {
    if reverse_complement {
        reverse_complement_label(label)
    } else {
        label.to_vec()
    }
}

fn reverse_complement_label(label: &[u8]) -> Vec<u8> {
    label
        .iter()
        .rev()
        .map(|&base| complement_ascii(base))
        .collect()
}

fn canonical_label(label: Vec<u8>) -> Vec<u8> {
    if reverse_complement_is_less(&label) {
        reverse_complement_label(&label)
    } else {
        label
    }
}

fn bucketed_maximal_unitig_reduce(mut unitigs: Vec<Vec<u8>>, threads: usize) -> Vec<Vec<u8>> {
    if unitigs.len() <= 2048 {
        unitigs.sort_unstable();
        unitigs.dedup();
        return unitigs;
    }

    let bucket_count = (threads.max(1) * 1024)
        .next_power_of_two()
        .clamp(1024, 16_384);
    let bucket_mask = bucket_count - 1;
    let mut buckets = (0..bucket_count).map(|_| Vec::new()).collect::<Vec<_>>();
    for unitig in unitigs {
        let bucket = (hash_bytes(&unitig, 0) as usize) & bucket_mask;
        buckets[bucket].push(unitig);
    }

    let workers = threads.max(1).min(bucket_count);
    if workers == 1 {
        for bucket in &mut buckets {
            bucket.sort_unstable();
            bucket.dedup();
        }
    } else {
        let chunk_size = bucket_count.div_ceil(workers);
        std::thread::scope(|scope| {
            for chunk in buckets.chunks_mut(chunk_size) {
                scope.spawn(move || {
                    for bucket in chunk {
                        bucket.sort_unstable();
                        bucket.dedup();
                    }
                });
            }
        });
    }

    let total = buckets.iter().map(Vec::len).sum();
    let mut reduced = Vec::with_capacity(total);
    for mut bucket in buckets {
        reduced.append(&mut bucket);
    }
    reduced
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FinalUnitigBucketEntry {
    bucket_id: usize,
    path: PathBuf,
    records: u64,
    bases: u64,
    colored: bool,
    direct_output: bool,
}

struct FinalUnitigBucketWriters {
    dir: PathBuf,
    bucket_mask: usize,
    next_bucket: usize,
    writers: Vec<Option<FinalUnitigBucketWriter>>,
    records: Vec<u64>,
    bases: Vec<u64>,
    open_writers: usize,
    colored: bool,
    direct_output: Option<BufWriter<File>>,
    direct_output_path: Option<PathBuf>,
    direct_records: u64,
    direct_record_id_highwater: u64,
    direct_bases: u64,
}

impl FinalUnitigBucketWriters {
    fn create(dir: &Path, threads: usize) -> Result<Self, SerialCollationError> {
        Self::create_with_color(dir, threads, false)
    }

    fn create_with_color(
        dir: &Path,
        threads: usize,
        colored: bool,
    ) -> Result<Self, SerialCollationError> {
        if dir.exists() {
            fs::remove_dir_all(dir).map_err(|source| SerialCollationError::Io {
                path: dir.to_path_buf(),
                source,
            })?;
        }
        fs::create_dir_all(dir).map_err(|source| SerialCollationError::Io {
            path: dir.to_path_buf(),
            source,
        })?;
        let bucket_count = final_unitig_bucket_count(threads);
        let mut writers = Vec::with_capacity(bucket_count);
        writers.resize_with(bucket_count, || None);
        Ok(Self {
            dir: dir.to_path_buf(),
            bucket_mask: bucket_count - 1,
            next_bucket: 0,
            writers,
            records: vec![0; bucket_count],
            bases: vec![0; bucket_count],
            open_writers: 0,
            colored,
            direct_output: None,
            direct_output_path: None,
            direct_records: 0,
            direct_record_id_highwater: 0,
            direct_bases: 0,
        })
    }

    fn create_direct(output_path: &Path, colored: bool) -> Result<Self, SerialCollationError> {
        let file = File::create(output_path).map_err(|source| SerialCollationError::Io {
            path: output_path.to_path_buf(),
            source,
        })?;
        Ok(Self {
            dir: output_path.parent().unwrap_or(Path::new(".")).to_path_buf(),
            bucket_mask: 0,
            next_bucket: 0,
            writers: Vec::new(),
            records: Vec::new(),
            bases: Vec::new(),
            open_writers: 0,
            colored,
            direct_output: Some(BufWriter::with_capacity(8 * 1024 * 1024, file)),
            direct_output_path: Some(output_path.to_path_buf()),
            direct_records: 0,
            direct_record_id_highwater: 0,
            direct_bases: 0,
        })
    }

    fn write_label(&mut self, label: &[u8]) -> Result<(), SerialCollationError> {
        if self.direct_output.is_some() {
            return self.write_direct_record(label, &[]);
        }
        let bucket_id = self.next_bucket;
        self.next_bucket = (self.next_bucket + 1) & self.bucket_mask;
        if self.writers[bucket_id].is_none() {
            self.evict_writer_if_needed(bucket_id)?;
            self.writers[bucket_id] = Some(FinalUnitigBucketWriter::open(
                &self.dir,
                bucket_id,
                self.colored,
            )?);
            self.open_writers += 1;
        }
        self.writers[bucket_id]
            .as_mut()
            .expect("final unitig bucket writer was just created")
            .write_record(label, &[])?;
        self.records[bucket_id] += 1;
        self.bases[bucket_id] += label.len() as u64;
        Ok(())
    }

    fn write_colored_label(
        &mut self,
        label: &[u8],
        colors: &[UnitigColor],
    ) -> Result<(), SerialCollationError> {
        if !self.colored {
            return Err(SerialCollationError::MalformedCoordBucket(self.dir.clone()));
        }
        if self.direct_output.is_some() {
            return self.write_direct_record(label, colors);
        }
        let bucket_id = self.next_bucket;
        self.next_bucket = (self.next_bucket + 1) & self.bucket_mask;
        if self.writers[bucket_id].is_none() {
            self.evict_writer_if_needed(bucket_id)?;
            self.writers[bucket_id] =
                Some(FinalUnitigBucketWriter::open(&self.dir, bucket_id, true)?);
            self.open_writers += 1;
        }
        self.writers[bucket_id]
            .as_mut()
            .expect("final unitig bucket writer was just created")
            .write_record(label, colors)?;
        self.records[bucket_id] += 1;
        self.bases[bucket_id] += label.len() as u64;
        Ok(())
    }

    fn write_direct_record(
        &mut self,
        label: &[u8],
        colors: &[UnitigColor],
    ) -> Result<(), SerialCollationError> {
        self.direct_records += 1;
        self.direct_record_id_highwater += 1;
        self.direct_bases += label.len() as u64;
        let output_path = self
            .direct_output_path
            .as_ref()
            .expect("direct output path exists")
            .clone();
        let out = self.direct_output.as_mut().expect("direct output exists");
        let mut header = Vec::with_capacity(4 + colors.len() * 12);
        header.extend_from_slice(b">0");
        for color in colors {
            header.push(b' ');
            append_decimal_u64(&mut header, color.raw());
        }
        header.push(b'\n');
        out.write_all(&header)
            .and_then(|_| out.write_all(label))
            .and_then(|_| out.write_all(b"\n"))
            .map_err(|source| SerialCollationError::Io {
                path: output_path,
                source,
            })
    }

    fn write_direct_batch(
        &mut self,
        bytes: &[u8],
        records: u64,
        bases: u64,
    ) -> Result<(), SerialCollationError> {
        let output_path = self
            .direct_output_path
            .as_ref()
            .expect("direct output path exists")
            .clone();
        self.direct_output
            .as_mut()
            .expect("direct output exists")
            .write_all(bytes)
            .map_err(|source| SerialCollationError::Io {
                path: output_path,
                source,
            })?;
        self.direct_records += records;
        self.direct_bases += bases;
        Ok(())
    }

    fn append_direct_fasta_file(
        &mut self,
        path: &Path,
        records: u64,
        bases: u64,
    ) -> Result<(), SerialCollationError> {
        let file = File::open(path).map_err(|source| SerialCollationError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        std::io::copy(
            &mut BufReader::with_capacity(8 * 1024 * 1024, file),
            self.direct_output.as_mut().expect("direct output exists"),
        )
        .map_err(|source| SerialCollationError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        self.direct_records += records;
        self.direct_record_id_highwater += records;
        self.direct_bases += bases;
        Ok(())
    }

    fn prepare_parallel_direct_output(
        &mut self,
    ) -> Result<(File, PathBuf, u64), SerialCollationError> {
        let path = self
            .direct_output_path
            .clone()
            .expect("direct output path exists");
        let output = self.direct_output.as_mut().expect("direct output exists");
        output.flush().map_err(|source| SerialCollationError::Io {
            path: path.clone(),
            source,
        })?;
        let file = output
            .get_ref()
            .try_clone()
            .map_err(|source| SerialCollationError::Io {
                path: path.clone(),
                source,
            })?;
        let offset = file
            .metadata()
            .map_err(|source| SerialCollationError::Io {
                path: path.clone(),
                source,
            })?
            .len();
        Ok((file, path, offset))
    }

    fn evict_writer_if_needed(
        &mut self,
        requested_bucket_id: usize,
    ) -> Result<(), SerialCollationError> {
        if self.open_writers < MAX_OPEN_STITCH_ENDPOINT_WRITERS {
            return Ok(());
        }
        let evict_bucket_id = self
            .writers
            .iter()
            .enumerate()
            .find_map(|(bucket_id, writer)| {
                (bucket_id != requested_bucket_id && writer.is_some()).then_some(bucket_id)
            })
            .unwrap_or(requested_bucket_id);
        if let Some(mut writer) = self.writers[evict_bucket_id].take() {
            writer.flush()?;
            self.open_writers -= 1;
        }
        Ok(())
    }

    fn finish(mut self) -> Result<Vec<FinalUnitigBucketEntry>, SerialCollationError> {
        if let Some(mut output) = self.direct_output.take() {
            output.flush().map_err(|source| SerialCollationError::Io {
                path: self
                    .direct_output_path
                    .clone()
                    .expect("direct output path exists"),
                source,
            })?;
            return Ok(vec![FinalUnitigBucketEntry {
                bucket_id: 0,
                path: self.direct_output_path.expect("direct output path exists"),
                records: self.direct_records,
                bases: self.direct_bases,
                colored: self.colored,
                direct_output: true,
            }]);
        }
        let mut manifest = Vec::new();
        for writer in self.writers.iter_mut().flatten() {
            writer.flush()?;
        }
        for (bucket_id, &records) in self.records.iter().enumerate() {
            if records != 0 {
                manifest.push(FinalUnitigBucketEntry {
                    bucket_id,
                    path: self.dir.join(format!("{bucket_id:05}.fub")),
                    records,
                    bases: self.bases[bucket_id],
                    colored: self.colored,
                    direct_output: false,
                });
            }
        }
        manifest.sort_by_key(|entry| entry.bucket_id);
        Ok(manifest)
    }
}

struct FinalUnitigBucketWriter {
    path: PathBuf,
    out: BufWriter<File>,
    colored: bool,
}

impl FinalUnitigBucketWriter {
    fn open(dir: &Path, bucket_id: usize, colored: bool) -> Result<Self, SerialCollationError> {
        let path = dir.join(format!("{bucket_id:05}.fub"));
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)
            .map_err(|source| SerialCollationError::Io {
                path: path.clone(),
                source,
            })?;
        Ok(Self {
            path,
            out: BufWriter::with_capacity(FINAL_UNITIG_BUCKET_WRITE_BUFFER, file),
            colored,
        })
    }

    fn write_record(
        &mut self,
        label: &[u8],
        colors: &[UnitigColor],
    ) -> Result<(), SerialCollationError> {
        let len = u32::try_from(label.len())
            .map_err(|_| SerialCollationError::MalformedCoordBucket(self.path.clone()))?;
        self.out
            .write_all(&len.to_le_bytes())
            .and_then(|_| {
                if self.colored {
                    let count = u32::try_from(colors.len())
                        .map_err(|_| std::io::Error::from(std::io::ErrorKind::InvalidData))?;
                    self.out.write_all(&count.to_le_bytes())?;
                }
                self.out.write_all(label)
            })
            .and_then(|_| {
                if self.colored {
                    for color in colors {
                        self.out.write_all(&color.raw().to_le_bytes())?;
                    }
                }
                Ok(())
            })
            .map_err(|source| SerialCollationError::Io {
                path: self.path.clone(),
                source,
            })?;
        Ok(())
    }

    fn flush(&mut self) -> Result<(), SerialCollationError> {
        self.out.flush().map_err(|source| SerialCollationError::Io {
            path: self.path.clone(),
            source,
        })
    }
}

fn reduce_final_unitig_buckets_to_fasta(
    manifest: &[FinalUnitigBucketEntry],
    output_path: &Path,
) -> Result<(u64, u64), SerialCollationError> {
    if let [entry] = manifest
        && entry.direct_output
    {
        if entry.path != output_path {
            return Err(SerialCollationError::MalformedCoordBucket(
                entry.path.clone(),
            ));
        }
        return Ok((entry.records, entry.bases));
    }
    let file = File::create(output_path).map_err(|source| SerialCollationError::Io {
        path: output_path.to_path_buf(),
        source,
    })?;
    let mut out = file;
    let mut fasta_buffer = Vec::with_capacity(8 * 1024 * 1024);
    let mut emitted = 0u64;
    let mut bases = 0u64;

    for entry in manifest {
        let bucket = read_final_unitig_bucket(entry)?;
        remove_serial_file(&entry.path)?;
        let bytes = &bucket.bytes;
        for (label_start, label_len, color_start, color_count) in bucket.labels {
            let label = &bytes[label_start..label_start + label_len];
            emitted += 1;
            bases += label.len() as u64;
            fasta_buffer.extend_from_slice(b">0");
            for color in &bucket.colors[color_start..color_start + color_count] {
                fasta_buffer.push(b' ');
                append_decimal_u64(&mut fasta_buffer, color.raw());
            }
            fasta_buffer.push(b'\n');
            fasta_buffer.extend_from_slice(&label);
            fasta_buffer.push(b'\n');
            if fasta_buffer.len() >= 8 * 1024 * 1024 {
                out.write_all(&fasta_buffer)
                    .map_err(|source| SerialCollationError::Io {
                        path: output_path.to_path_buf(),
                        source,
                    })?;
                fasta_buffer.clear();
            }
        }
    }

    out.write_all(&fasta_buffer)
        .map_err(|source| SerialCollationError::Io {
            path: output_path.to_path_buf(),
            source,
        })?;

    out.flush().map_err(|source| SerialCollationError::Io {
        path: output_path.to_path_buf(),
        source,
    })?;
    Ok((emitted, bases))
}

#[inline]
fn append_decimal_u64(output: &mut Vec<u8>, mut value: u64) {
    let mut digits = [0u8; 20];
    let mut start = digits.len();
    loop {
        start -= 1;
        digits[start] = b'0' + (value % 10) as u8;
        value /= 10;
        if value == 0 {
            break;
        }
    }
    output.extend_from_slice(&digits[start..]);
}

struct FinalUnitigBucketData {
    bytes: Vec<u8>,
    labels: Vec<(usize, usize, usize, usize)>,
    colors: Vec<UnitigColor>,
}

fn read_final_unitig_bucket(
    entry: &FinalUnitigBucketEntry,
) -> Result<FinalUnitigBucketData, SerialCollationError> {
    let bytes = fs::read(&entry.path).map_err(|source| SerialCollationError::Io {
        path: entry.path.clone(),
        source,
    })?;
    let mut labels = Vec::with_capacity(entry.records as usize);
    let mut colors = Vec::new();
    let mut cursor = 0usize;
    for _ in 0..entry.records {
        let len_end = cursor
            .checked_add(4)
            .ok_or_else(|| SerialCollationError::MalformedCoordBucket(entry.path.clone()))?;
        let len_bytes = bytes
            .get(cursor..len_end)
            .ok_or_else(|| SerialCollationError::MalformedCoordBucket(entry.path.clone()))?;
        let len =
            u32::from_le_bytes(len_bytes.try_into().expect("four-byte label length")) as usize;
        let color_count = if entry.colored {
            let count_end = len_end + 4;
            let count_bytes = bytes
                .get(len_end..count_end)
                .ok_or_else(|| SerialCollationError::MalformedCoordBucket(entry.path.clone()))?;
            cursor = count_end;
            u32::from_le_bytes(count_bytes.try_into().expect("four-byte color count")) as usize
        } else {
            cursor = len_end;
            0
        };
        let label_end = cursor
            .checked_add(len)
            .ok_or_else(|| SerialCollationError::MalformedCoordBucket(entry.path.clone()))?;
        bytes
            .get(cursor..label_end)
            .ok_or_else(|| SerialCollationError::MalformedCoordBucket(entry.path.clone()))?;
        let label_start = cursor;
        cursor = label_end;
        let color_start = colors.len();
        for _ in 0..color_count {
            let color_end = cursor + 8;
            let raw = bytes
                .get(cursor..color_end)
                .ok_or_else(|| SerialCollationError::MalformedCoordBucket(entry.path.clone()))?;
            let raw = u64::from_le_bytes(raw.try_into().expect("eight-byte unitig color"));
            colors.push(UnitigColor::new(
                (raw & 0xff_ffff) as u32,
                crate::state::ColorCoordinate::from_u40(raw >> 24),
            ));
            cursor = color_end;
        }
        labels.push((label_start, len, color_start, color_count));
    }
    if cursor != bytes.len() {
        return Err(SerialCollationError::MalformedCoordBucket(
            entry.path.clone(),
        ));
    }
    Ok(FinalUnitigBucketData {
        bytes,
        labels,
        colors,
    })
}

fn reverse_complement_is_less(label: &[u8]) -> bool {
    for (idx, &base) in label.iter().enumerate() {
        let rc_base = complement_ascii(label[label.len() - 1 - idx]);
        if rc_base != base {
            return rc_base < base;
        }
    }
    false
}

fn retain_maximal_labels(unitigs: &mut Vec<Vec<u8>>) {
    let mut keep = vec![true; unitigs.len()];
    for i in 0..unitigs.len() {
        for j in 0..unitigs.len() {
            if i == j || unitigs[i].len() >= unitigs[j].len() {
                continue;
            }
            if contains_label(&unitigs[j], &unitigs[i])
                || contains_label(&reverse_complement_label(&unitigs[j]), &unitigs[i])
            {
                keep[i] = false;
                break;
            }
        }
    }

    let mut idx = 0usize;
    unitigs.retain(|_| {
        let retain = keep[idx];
        idx += 1;
        retain
    });
}

fn retain_maximal_labels_indexed_with_threads(unitigs: &mut Vec<Vec<u8>>, threads: usize) {
    let started = Instant::now();
    let workers = threads.max(1).min(unitigs.len().max(1));
    eprintln!(
        "cuttlefish3-rs: retaining maximal labels among {} unitig(s) with indexed containment and {} worker(s)",
        unitigs.len(),
        workers
    );
    if unitigs.len() <= 2048 {
        retain_maximal_labels(unitigs);
        return;
    }

    let mut lengths = unitigs.iter().map(Vec::len).collect::<Vec<_>>();
    lengths.sort_unstable();
    lengths.dedup();

    let mut keep = vec![true; unitigs.len()];
    let mut patterns_by_len = FastHashMap::<usize, FastHashMap<u128, Vec<usize>>>::default();
    for &len in &lengths {
        let mut patterns = FastHashMap::<u128, Vec<usize>>::default();
        for index in unitigs
            .iter()
            .enumerate()
            .filter_map(|(index, unitig)| (unitig.len() == len).then_some(index))
        {
            patterns
                .entry(hash_label(&unitigs[index]))
                .or_default()
                .push(index);
        }
        patterns_by_len.insert(len, patterns);
    }

    let powers = lengths
        .iter()
        .map(|&len| {
            (
                len,
                (hash_pow(HASH_BASE_1, len), hash_pow(HASH_BASE_2, len)),
            )
        })
        .collect::<FastHashMap<_, _>>();

    if workers == 1 {
        for haystack in unitigs.iter() {
            mark_contained_candidates_for_haystack(
                haystack,
                &lengths,
                &patterns_by_len,
                &powers,
                unitigs,
                &mut keep,
            );
            let reverse = reverse_complement_label(haystack);
            mark_contained_candidates_for_haystack(
                &reverse,
                &lengths,
                &patterns_by_len,
                &powers,
                unitigs,
                &mut keep,
            );
        }
    } else {
        let chunk_len = unitigs.len().div_ceil(workers);
        let local_keeps = std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for chunk in unitigs.chunks(chunk_len) {
                let lengths = &lengths;
                let patterns_by_len = &patterns_by_len;
                let powers = &powers;
                let unitigs = &*unitigs;
                handles.push(scope.spawn(move || {
                    let mut local_keep = vec![true; unitigs.len()];
                    for haystack in chunk {
                        mark_contained_candidates_for_haystack(
                            haystack,
                            &lengths,
                            &patterns_by_len,
                            &powers,
                            unitigs,
                            &mut local_keep,
                        );
                        let reverse = reverse_complement_label(haystack);
                        mark_contained_candidates_for_haystack(
                            &reverse,
                            &lengths,
                            &patterns_by_len,
                            &powers,
                            unitigs,
                            &mut local_keep,
                        );
                    }
                    local_keep
                }));
            }

            handles
                .into_iter()
                .map(|handle| handle.join().expect("maximal-label filter worker panicked"))
                .collect::<Vec<_>>()
        });
        for local_keep in local_keeps {
            for (dst, src) in keep.iter_mut().zip(local_keep) {
                *dst &= src;
            }
        }
    }

    if let Ok(path) = std::env::var("CF3_RS_DUMP_MAXIMAL_FILTER_REMOVALS") {
        if !path.is_empty() {
            if let Err(err) = dump_removed_labels(std::path::Path::new(&path), unitigs, &keep) {
                eprintln!("cuttlefish3-rs: failed to dump maximal-filter removals: {err}");
            }
        }
    }

    let mut idx = 0usize;
    unitigs.retain(|_| {
        let retain = keep[idx];
        idx += 1;
        retain
    });
    eprintln!(
        "cuttlefish3-rs: maximal-label filter completed in {:.1}s",
        started.elapsed().as_secs_f64()
    );
}

fn dump_removed_labels(
    path: &std::path::Path,
    unitigs: &[Vec<u8>],
    keep: &[bool],
) -> std::io::Result<()> {
    use std::io::Write;

    let mut out = std::fs::File::create(path)?;
    for (index, (unitig, &retain)) in unitigs.iter().zip(keep).enumerate() {
        if !retain {
            writeln!(out, ">{index}_len_{}", unitig.len())?;
            out.write_all(unitig)?;
            writeln!(out)?;
        }
    }
    Ok(())
}

fn dump_collation_candidates(
    path: &std::path::Path,
    unitigs: &[Vec<u8>],
    origins: &[&'static str],
) -> std::io::Result<()> {
    use std::io::Write;

    let mut out = std::fs::File::create(path)?;
    for (index, (unitig, origin)) in unitigs.iter().zip(origins).enumerate() {
        writeln!(out, ">{index}_{origin}_len_{}", unitig.len())?;
        out.write_all(unitig)?;
        writeln!(out)?;
    }
    Ok(())
}

fn suppress_labels_contained_in_sources(
    labels: &mut Vec<Vec<u8>>,
    sources: &[Vec<u8>],
    threads: usize,
) -> usize {
    if labels.is_empty() || sources.is_empty() {
        return 0;
    }

    let mut lengths = labels.iter().map(Vec::len).collect::<Vec<_>>();
    lengths.sort_unstable();
    lengths.dedup();

    let mut patterns_by_len = FastHashMap::<usize, FastHashMap<u128, Vec<usize>>>::default();
    for &len in &lengths {
        let mut patterns = FastHashMap::<u128, Vec<usize>>::default();
        for (index, label) in labels
            .iter()
            .enumerate()
            .filter(|(_, label)| label.len() == len)
        {
            patterns.entry(hash_label(label)).or_default().push(index);
        }
        patterns_by_len.insert(len, patterns);
    }

    let powers = lengths
        .iter()
        .map(|&len| {
            (
                len,
                (hash_pow(HASH_BASE_1, len), hash_pow(HASH_BASE_2, len)),
            )
        })
        .collect::<FastHashMap<_, _>>();

    let workers = threads.max(1).min(sources.len().max(1));
    let mut keep = vec![true; labels.len()];
    if workers == 1 {
        for source in sources {
            mark_contained_candidates_for_haystack(
                source,
                &lengths,
                &patterns_by_len,
                &powers,
                labels,
                &mut keep,
            );
            let reverse = reverse_complement_label(source);
            mark_contained_candidates_for_haystack(
                &reverse,
                &lengths,
                &patterns_by_len,
                &powers,
                labels,
                &mut keep,
            );
        }
    } else {
        let chunk_len = sources.len().div_ceil(workers);
        let local_keeps = std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for chunk in sources.chunks(chunk_len) {
                let lengths = &lengths;
                let patterns_by_len = &patterns_by_len;
                let powers = &powers;
                let labels = &*labels;
                handles.push(scope.spawn(move || {
                    let mut local_keep = vec![true; labels.len()];
                    for source in chunk {
                        mark_contained_candidates_for_haystack(
                            source,
                            lengths,
                            patterns_by_len,
                            powers,
                            labels,
                            &mut local_keep,
                        );
                        let reverse = reverse_complement_label(source);
                        mark_contained_candidates_for_haystack(
                            &reverse,
                            lengths,
                            patterns_by_len,
                            powers,
                            labels,
                            &mut local_keep,
                        );
                    }
                    local_keep
                }));
            }

            handles
                .into_iter()
                .map(|handle| handle.join().expect("stitched-containment worker panicked"))
                .collect::<Vec<_>>()
        });
        for local_keep in local_keeps {
            for (dst, src) in keep.iter_mut().zip(local_keep) {
                *dst &= src;
            }
        }
    }

    let before = labels.len();
    let mut idx = 0usize;
    labels.retain(|_| {
        let retain = keep[idx];
        idx += 1;
        retain
    });
    before - labels.len()
}

fn contains_label(haystack: &[u8], needle: &[u8]) -> bool {
    needle.len() <= haystack.len()
        && haystack
            .windows(needle.len())
            .any(|window| window == needle)
}

const HASH_BASE_1: u64 = 1_099_511_628_211;
const HASH_BASE_2: u64 = 1_000_000_007;

fn mark_contained_candidates_for_haystack(
    haystack: &[u8],
    lengths: &[usize],
    patterns_by_len: &FastHashMap<usize, FastHashMap<u128, Vec<usize>>>,
    powers: &FastHashMap<usize, (u64, u64)>,
    unitigs: &[Vec<u8>],
    keep: &mut [bool],
) {
    if haystack.is_empty() {
        return;
    }
    let mut prefix1 = Vec::with_capacity(haystack.len() + 1);
    let mut prefix2 = Vec::with_capacity(haystack.len() + 1);
    prefix1.push(0u64);
    prefix2.push(0u64);
    for &base in haystack {
        let code = hash_base_code(base);
        prefix1.push(
            prefix1
                .last()
                .copied()
                .unwrap()
                .wrapping_mul(HASH_BASE_1)
                .wrapping_add(code),
        );
        prefix2.push(
            prefix2
                .last()
                .copied()
                .unwrap()
                .wrapping_mul(HASH_BASE_2)
                .wrapping_add(code),
        );
    }

    for &len in lengths {
        if len >= haystack.len() {
            continue;
        }
        let Some(patterns) = patterns_by_len.get(&len) else {
            continue;
        };
        let Some(&(pow1, pow2)) = powers.get(&len) else {
            continue;
        };
        for start in 0..=haystack.len() - len {
            let end = start + len;
            let h1 = prefix1[end].wrapping_sub(prefix1[start].wrapping_mul(pow1));
            let h2 = prefix2[end].wrapping_sub(prefix2[start].wrapping_mul(pow2));
            let hash = ((h1 as u128) << 64) | h2 as u128;
            let Some(indices) = patterns.get(&hash) else {
                continue;
            };
            let window = &haystack[start..end];
            for &index in indices {
                if keep[index] && unitigs[index].as_slice() == window {
                    keep[index] = false;
                }
            }
        }
    }
}

fn hash_label(label: &[u8]) -> u128 {
    let mut h1 = 0u64;
    let mut h2 = 0u64;
    for &base in label {
        let code = hash_base_code(base);
        h1 = h1.wrapping_mul(HASH_BASE_1).wrapping_add(code);
        h2 = h2.wrapping_mul(HASH_BASE_2).wrapping_add(code);
    }
    ((h1 as u128) << 64) | h2 as u128
}

fn hash_pow(base: u64, len: usize) -> u64 {
    let mut pow = 1u64;
    for _ in 0..len {
        pow = pow.wrapping_mul(base);
    }
    pow
}

fn hash_base_code(base: u8) -> u64 {
    match base {
        b'A' => 1,
        b'C' => 2,
        b'G' => 3,
        b'T' => 4,
        _ => 5,
    }
}

fn phi_edge_path_info<const K: usize>(
    edge: &DiscontinuityEdge<K>,
    vertex_info: PathInfo<K>,
) -> PathInfo<K> {
    let rank = if vertex_info.rank == 1 {
        0
    } else {
        vertex_info.rank
    };
    let exit_side = if rank == 0 {
        edge.unitig_exit_side
    } else {
        edge.unitig_exit_side.inverse()
    };

    PathInfo {
        path_id: vertex_info.path_id,
        rank,
        exit_side,
        is_cycle: vertex_info.is_cycle,
    }
}

#[inline(always)]
fn compact_phi_edge_path_info<const K: usize>(
    edge: &DiscontinuityEdge<K>,
    vertex_info: CompactExpansionPathInfo,
) -> CompactExpansionPathInfo {
    compact_phi_edge_path_info_with_exit(edge.unitig_exit_side, vertex_info)
}

#[inline(always)]
fn compact_phi_edge_path_info_with_exit(
    unitig_exit_side: Side,
    vertex_info: CompactExpansionPathInfo,
) -> CompactExpansionPathInfo {
    let rank = if vertex_info.rank() == 1 {
        0
    } else {
        vertex_info.rank()
    };
    let exit_side = if rank == 0 {
        unitig_exit_side
    } else {
        unitig_exit_side.inverse()
    };
    CompactExpansionPathInfo::new(vertex_info.path_id, rank, exit_side, vertex_info.is_cycle())
}

fn edge_path_info<const K: usize>(
    edge: &DiscontinuityEdge<K>,
    first_info: PathInfo<K>,
    second_info: PathInfo<K>,
) -> PathInfo<K> {
    let rank = first_info.rank.min(second_info.rank);
    let exit_side = if rank == first_info.rank {
        edge.unitig_exit_side
    } else {
        edge.unitig_exit_side.inverse()
    };

    PathInfo {
        path_id: first_info.path_id,
        rank,
        exit_side,
        is_cycle: first_info.is_cycle,
    }
}

#[inline(always)]
fn compact_edge_path_info<const K: usize>(
    edge: &DiscontinuityEdge<K>,
    first_info: CompactExpansionPathInfo,
    second_info: CompactExpansionPathInfo,
) -> CompactExpansionPathInfo {
    compact_edge_path_info_with_exit(edge.unitig_exit_side, first_info, second_info)
}

#[inline(always)]
fn compact_edge_path_info_with_exit(
    unitig_exit_side: Side,
    first_info: CompactExpansionPathInfo,
    second_info: CompactExpansionPathInfo,
) -> CompactExpansionPathInfo {
    let rank = first_info.rank().min(second_info.rank());
    let exit_side = if rank == first_info.rank() {
        unitig_exit_side
    } else {
        unitig_exit_side.inverse()
    };
    CompactExpansionPathInfo::new(first_info.path_id, rank, exit_side, first_info.is_cycle())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LabelEnd {
    Left,
    Right,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct HalfEnd {
    unitig_index: u32,
}

impl HalfEnd {
    #[inline]
    fn unitig_index(self) -> usize {
        self.unitig_index as usize
    }
}

#[inline]
fn label_end_for_node(node: usize) -> LabelEnd {
    if node & 1 == 0 {
        LabelEnd::Left
    } else {
        LabelEnd::Right
    }
}

#[inline]
fn reverse_for_stitch_node(node: usize) -> bool {
    node & 1 == 1
}

const STITCH_NO_NODE: u32 = u32::MAX;

#[inline]
fn stitch_node(node: usize) -> u32 {
    u32::try_from(node).expect("stitch node index exceeds u32")
}

#[inline]
fn stitch_node_index(node: u32) -> usize {
    node as usize
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
struct StitchAdjacency {
    to: usize,
    unitig_index: Option<usize>,
}

#[derive(Debug, Clone, Copy, Default)]
struct StitchAdjacencyList {
    edges: [StitchAdjacency; 3],
    len: u8,
}

impl StitchAdjacencyList {
    fn push(&mut self, edge: StitchAdjacency) {
        let len = self.len as usize;
        debug_assert!(len < self.edges.len());
        if len < self.edges.len() {
            self.edges[len] = edge;
            self.len += 1;
        }
    }

    fn len(&self) -> usize {
        self.len as usize
    }

    fn iter(&self) -> impl Iterator<Item = &StitchAdjacency> {
        self.edges[..self.len()].iter()
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct StitchVertexEnds {
    fronts: [u32; 2],
    backs: [u32; 2],
    nodes: [u32; 2],
    front_len: u8,
    back_len: u8,
    total_len: u8,
}

impl StitchVertexEnds {
    fn push(&mut self, side: Side, node: u32) {
        if self.total_len < 2 {
            self.nodes[self.total_len as usize] = node;
        }
        self.total_len = self.total_len.saturating_add(1);

        match side {
            Side::Front => {
                if self.front_len < 2 {
                    self.fronts[self.front_len as usize] = node;
                }
                self.front_len = self.front_len.saturating_add(1);
            }
            Side::Back => {
                if self.back_len < 2 {
                    self.backs[self.back_len as usize] = node;
                }
                self.back_len = self.back_len.saturating_add(1);
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct StitchEndpointRecord<const K: usize> {
    vertex: Kmer<K>,
    side: Side,
    node: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StitchEndpointBucketEntry {
    path: PathBuf,
    records: u64,
}

const MAX_OPEN_STITCH_ENDPOINT_WRITERS: usize = 512;

struct StitchEndpointBucketWriters {
    dir: PathBuf,
    bucket_mask: usize,
    writers: Vec<Option<StitchEndpointBucketWriter>>,
    records: Vec<u64>,
    open_writers: usize,
}

impl StitchEndpointBucketWriters {
    fn create(dir: &Path, threads: usize) -> Result<Self, SerialCollationError> {
        if dir.exists() {
            fs::remove_dir_all(dir).map_err(|source| SerialCollationError::Io {
                path: dir.to_path_buf(),
                source,
            })?;
        }
        fs::create_dir_all(dir).map_err(|source| SerialCollationError::Io {
            path: dir.to_path_buf(),
            source,
        })?;
        let bucket_count = stitch_endpoint_bucket_count(threads);
        let mut writers = Vec::with_capacity(bucket_count);
        writers.resize_with(bucket_count, || None);
        Ok(Self {
            dir: dir.to_path_buf(),
            bucket_mask: bucket_count - 1,
            writers,
            records: vec![0; bucket_count],
            open_writers: 0,
        })
    }

    fn write_record<const K: usize>(
        &mut self,
        record: &StitchEndpointRecord<K>,
    ) -> Result<(), SerialCollationError> {
        let bits = record.vertex.as_u128();
        let bucket_id = (hash_bytes(&bits.to_le_bytes(), 0) as usize) & self.bucket_mask;
        if self.writers[bucket_id].is_none() {
            self.evict_writer_if_needed(bucket_id)?;
            self.writers[bucket_id] = Some(StitchEndpointBucketWriter::open(
                &self.dir,
                bucket_id,
                self.records[bucket_id],
            )?);
            self.open_writers += 1;
        }
        let writer = self.writers[bucket_id]
            .as_mut()
            .expect("stitch endpoint bucket writer was just created");
        writer.write_record(bits, record.side, record.node)?;
        self.records[bucket_id] += 1;
        Ok(())
    }

    fn evict_writer_if_needed(
        &mut self,
        requested_bucket_id: usize,
    ) -> Result<(), SerialCollationError> {
        if self.open_writers < MAX_OPEN_STITCH_ENDPOINT_WRITERS {
            return Ok(());
        }
        let evict_bucket_id = self
            .writers
            .iter()
            .enumerate()
            .find_map(|(bucket_id, writer)| {
                (bucket_id != requested_bucket_id && writer.is_some()).then_some(bucket_id)
            })
            .unwrap_or(requested_bucket_id);
        if let Some(mut writer) = self.writers[evict_bucket_id].take() {
            writer.flush()?;
            self.open_writers -= 1;
        }
        Ok(())
    }

    fn finish(mut self) -> Result<Vec<StitchEndpointBucketEntry>, SerialCollationError> {
        let mut manifest = Vec::new();
        for writer in self.writers.iter_mut().flatten() {
            writer.flush()?;
        }
        for (bucket_id, &records) in self.records.iter().enumerate() {
            if records != 0 {
                manifest.push(StitchEndpointBucketEntry {
                    path: self.dir.join(format!("{bucket_id:05}.seb")),
                    records,
                });
            }
        }
        Ok(manifest)
    }
}

struct StitchEndpointBucketWriter {
    path: PathBuf,
    out: BufWriter<File>,
    records: u64,
}

impl StitchEndpointBucketWriter {
    fn open(dir: &Path, bucket_id: usize, records: u64) -> Result<Self, SerialCollationError> {
        let path = dir.join(format!("{bucket_id:05}.seb"));
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)
            .map_err(|source| SerialCollationError::Io {
                path: path.clone(),
                source,
            })?;
        Ok(Self {
            path,
            out: BufWriter::with_capacity(STITCH_COORD_SHARD_WRITE_BUFFER, file),
            records,
        })
    }

    fn write_record(
        &mut self,
        vertex_bits: u128,
        side: Side,
        node: u32,
    ) -> Result<(), SerialCollationError> {
        let side = match side {
            Side::Front => 0u8,
            Side::Back => 1u8,
        };
        let mut bytes = [0u8; 21];
        bytes[..16].copy_from_slice(&vertex_bits.to_le_bytes());
        bytes[16] = side;
        bytes[17..21].copy_from_slice(&node.to_le_bytes());
        self.out
            .write_all(&bytes)
            .map_err(|source| SerialCollationError::Io {
                path: self.path.clone(),
                source,
            })?;
        self.records += 1;
        Ok(())
    }

    fn flush(&mut self) -> Result<(), SerialCollationError> {
        self.out.flush().map_err(|source| SerialCollationError::Io {
            path: self.path.clone(),
            source,
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StitchedPath {
    label: Vec<u8>,
    is_cycle: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct StitchedCoordRecord {
    path_id: u64,
    rank: u64,
    unitig_index: u32,
    reverse: bool,
    is_cycle: bool,
}

#[derive(Clone, Copy)]
struct DenseLocalPathInfo {
    path_id: u64,
    rank_and_flags: u64,
}

impl DenseLocalPathInfo {
    const EMPTY: Self = Self {
        path_id: 0,
        rank_and_flags: u64::MAX,
    };

    fn from_record(record: StitchedCoordRecord) -> Self {
        Self {
            path_id: record.path_id,
            rank_and_flags: (record.rank << 2)
                | u64::from(record.reverse)
                | (u64::from(record.is_cycle) << 1),
        }
    }

    fn to_record(self, unitig_index: usize) -> Option<StitchedCoordRecord> {
        (self.rank_and_flags != u64::MAX).then(|| StitchedCoordRecord {
            path_id: self.path_id,
            rank: self.rank_and_flags >> 2,
            unitig_index: unitig_index as u32,
            reverse: self.rank_and_flags & 1 != 0,
            is_cycle: self.rank_and_flags & 2 != 0,
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StitchedCoordBucketEntry {
    bucket_id: usize,
    records: u64,
    path: PathBuf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct MaterializedStitchedCoordRecord {
    path_id: u64,
    rank: u64,
    label_offset: u64,
    label_len: u32,
    reverse: bool,
    is_cycle: bool,
    color_index: u32,
    color_count: u32,
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct LoadedMaterializedStitchedCoordRecord {
    path_id: u64,
    label_offset: u32,
    color_start: u32,
    rank: u16,
    label_len: u16,
    color_count: u16,
    flags: u16,
}

const _: () = assert!(std::mem::size_of::<LoadedMaterializedStitchedCoordRecord>() == 24);

impl LoadedMaterializedStitchedCoordRecord {
    const REVERSE_FLAG: u16 = 1;
    const CYCLE_FLAG: u16 = 1 << 1;

    #[inline(always)]
    fn reverse(self) -> bool {
        self.flags & Self::REVERSE_FLAG != 0
    }

    #[inline(always)]
    fn is_cycle(self) -> bool {
        self.flags & Self::CYCLE_FLAG != 0
    }

    #[inline(always)]
    fn color_count(self) -> u32 {
        u32::from(self.color_count)
    }

    fn new(
        path_id: u64,
        rank: u64,
        label_offset: u32,
        label_len: u32,
        reverse: bool,
        is_cycle: bool,
        color_start: u32,
        color_count: u32,
    ) -> Self {
        Self {
            path_id,
            label_offset,
            color_start,
            rank: u16::try_from(rank).expect("materialized rank fits C++ weight_t"),
            label_len: u16::try_from(label_len)
                .expect("materialized label length fits C++ uni_len_t"),
            color_count: u16::try_from(color_count)
                .expect("materialized color count fits C++ uni_len_t"),
            flags: u16::from(reverse) * Self::REVERSE_FLAG | u16::from(is_cycle) * Self::CYCLE_FLAG,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct MaterializedStitchedCoordBucketEntry {
    bucket_id: usize,
    records: u64,
    label_bytes: u64,
    coord_path: PathBuf,
    label_path: PathBuf,
    color_path: Option<PathBuf>,
    color_runs: u64,
}

const STITCH_COORD_MAGIC: &[u8; 8] = b"CF3SCB2\0";
const MATERIALIZED_STITCH_COORD_MAGIC: &[u8; 8] = b"CF3MCB2\0";
const STITCH_COORD_HEADER_LEN: u64 = 32;
const STITCH_PATH_INFO_RECORD_LEN: u64 = 24;
const STITCH_COORD_RECORD_LEN: u64 = 24;
const MATERIALIZED_STITCH_COORD_SHARD_WRITE_BUFFER: usize = 16 * 1024;
const FINAL_UNITIG_BUCKET_WRITE_BUFFER: usize = 128 * 1024;
const STITCH_COORD_SHARD_WRITE_BUFFER: usize = 1024 * 1024;
const STITCH_COORD_RECORD_WRITE_BUFFER: usize = 1024 * 1024;
const EDGE_PATH_INFO_WORKER_BUFFER: usize = 128 * 1024;
const STITCH_COORD_REVERSE_FLAG: u8 = 1;
const STITCH_COORD_CYCLE_FLAG: u8 = 2;
const MAX_OPEN_STITCH_PATH_INFO_WRITERS: usize = 768;
const MAX_OPEN_MATERIALIZED_STITCH_WRITERS: usize = 768;
const MAX_OPEN_MATERIALIZED_STITCH_WRITERS_PER_SHARD: usize = 8;
const STITCH_PATH_INFO_BUCKET_TARGET: usize = 32;
const DEFAULT_INMEM_STITCH_PATH_INFO_LIMIT: usize = 768 * 1024 * 1024;
fn stitch_discontinuity_paths<const K: usize>(
    inputs: &DiscontinuityInputs<K>,
    skip_unitigs: &[bool],
    threads: usize,
) -> Vec<Vec<u8>> {
    stitch_discontinuity_paths_impl::<K>(inputs, skip_unitigs, threads, None)
        .expect("in-memory discontinuity stitching cannot fail")
}

fn stitch_discontinuity_paths_with_coord_dir<const K: usize>(
    inputs: &DiscontinuityInputs<K>,
    skip_unitigs: &[bool],
    threads: usize,
    coord_dir: &Path,
) -> Result<Vec<Vec<u8>>, SerialCollationError> {
    stitch_discontinuity_paths_impl::<K>(inputs, skip_unitigs, threads, Some(coord_dir))
}

fn stitch_discontinuity_paths_to_final_buckets<const K: usize>(
    inputs: &DiscontinuityInputs<K>,
    skip_unitigs: &[bool],
    threads: usize,
    coord_dir: &Path,
    final_buckets: &mut FinalUnitigBucketWriters,
) -> Result<u64, SerialCollationError> {
    let debug_stitch = std::env::var_os("CF3_RS_DEBUG_STITCH").is_some();
    let endpoint_dir = coord_dir.join("endpoints");
    let mut endpoint_writers = StitchEndpointBucketWriters::create(&endpoint_dir, threads)?;
    let build_started = Instant::now();
    let mut half_ends = Vec::<HalfEnd>::new();

    for (unitig_index, unitig) in inputs.unitigs.iter().enumerate() {
        if skip_unitigs.get(unitig_index).copied().unwrap_or(false) {
            continue;
        }
        if unitig.left_exit().is_none() && unitig.right_exit().is_none() {
            continue;
        }

        let (left_endpoint, right_endpoint) = endpoints_by_label_end(unitig);
        let left_node = half_ends.len();
        half_ends.push(HalfEnd {
            unitig_index: u32::try_from(unitig_index).expect("unitig index exceeds u32"),
        });
        if let Some(endpoint) = left_endpoint {
            endpoint_writers.write_record(&StitchEndpointRecord {
                vertex: endpoint.vertex,
                side: endpoint.side,
                node: stitch_node(left_node),
            })?;
        }

        let right_node = half_ends.len();
        half_ends.push(HalfEnd {
            unitig_index: u32::try_from(unitig_index).expect("unitig index exceeds u32"),
        });
        if let Some(endpoint) = right_endpoint {
            endpoint_writers.write_record(&StitchEndpointRecord {
                vertex: endpoint.vertex,
                side: endpoint.side,
                node: stitch_node(right_node),
            })?;
        }
    }
    let half_end_elapsed = build_started.elapsed();

    let mut join_neighbor = vec![STITCH_NO_NODE; half_ends.len()];
    let endpoint_sort_started = Instant::now();
    let endpoint_manifest = endpoint_writers.finish()?;
    let endpoint_sort_elapsed = endpoint_sort_started.elapsed();
    let endpoint_join_started = Instant::now();
    join_neighbors_from_endpoint_buckets::<K>(&endpoint_manifest, &mut join_neighbor, threads)?;
    let endpoint_join_elapsed = endpoint_join_started.elapsed();

    if debug_stitch {
        eprintln!("stitch inputs: {} unitigs", inputs.unitigs.len());
        for (node, half_end) in half_ends.iter().enumerate() {
            let (left_endpoint, right_endpoint) =
                endpoints_by_label_end(&inputs.unitigs[half_end.unitig_index()]);
            let label_end = label_end_for_node(node);
            let endpoint = match label_end {
                LabelEnd::Left => left_endpoint,
                LabelEnd::Right => right_endpoint,
            };
            eprintln!(
                "node {node}: unitig={} end={:?} endpoint={:?} join={}",
                half_end.unitig_index(),
                label_end,
                endpoint.map(debug_endpoint),
                join_neighbor[node],
            );
        }
    }

    let component_started = Instant::now();
    let component_starts = stitch_component_starts_from_neighbors(&half_ends, &join_neighbor);
    let component_elapsed = component_started.elapsed();
    let walk_started = Instant::now();
    let manifest = write_materialized_stitched_coord_buckets(
        inputs,
        coord_dir,
        threads,
        &half_ends,
        &join_neighbor,
        &component_starts,
    )?;
    let emitted = reduce_materialized_stitched_coord_bucket_files_to_final::<K>(
        &manifest,
        &[],
        threads,
        final_buckets,
    )?;
    let walk_elapsed = walk_started.elapsed();

    eprintln!(
        "cuttlefish3-rs: stitch detail: half-ends {:.3}s, endpoint sort {:.3}s, endpoint join {:.3}s, components {:.3}s, materialized walk+reduce {:.3}s",
        half_end_elapsed.as_secs_f64(),
        endpoint_sort_elapsed.as_secs_f64(),
        endpoint_join_elapsed.as_secs_f64(),
        component_elapsed.as_secs_f64(),
        walk_elapsed.as_secs_f64(),
    );
    Ok(emitted)
}

fn stitch_external_discontinuity_paths_to_final_buckets<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    threads: usize,
    coord_dir: &Path,
    final_buckets: &mut FinalUnitigBucketWriters,
) -> Result<u64, SerialCollationError> {
    report_process_memory("external stitch start");
    let endpoint_dir = coord_dir.join("endpoints");
    let mut endpoint_writers = StitchEndpointBucketWriters::create(&endpoint_dir, threads)?;
    let reader = ExternalDiscontinuityReader::open(inputs)?;
    let build_started = Instant::now();
    let mut half_ends = Vec::<HalfEnd>::new();
    let use_direct_label_reads = std::env::var_os("CF3_RS_DIRECT_STITCH_LABEL_READS").is_some();
    let use_ordered_labels = std::env::var_os("CF3_RS_ORDERED_STITCH_LABELS").is_some();
    let collect_label_refs = use_direct_label_reads || use_ordered_labels;
    let mut label_refs = if collect_label_refs {
        Vec::<ExternalLabelRef>::with_capacity(inputs.unitig_count())
    } else {
        Vec::<ExternalLabelRef>::new()
    };

    for (unitig_index, unitig) in reader.iter()?.enumerate() {
        let unitig = unitig?;
        if collect_label_refs {
            label_refs.push(ExternalLabelRef {
                label_start: unitig.label_start,
                label_len: unitig.label_len,
            });
        }
        if unitig.left_exit().is_none() && unitig.right_exit().is_none() {
            continue;
        }

        let (left_endpoint, right_endpoint) = endpoints_by_label_end(&unitig);
        let left_node = half_ends.len();
        half_ends.push(HalfEnd {
            unitig_index: u32::try_from(unitig_index).expect("unitig index exceeds u32"),
        });
        if let Some(endpoint) = left_endpoint {
            endpoint_writers.write_record(&StitchEndpointRecord {
                vertex: endpoint.vertex,
                side: endpoint.side,
                node: stitch_node(left_node),
            })?;
        }

        let right_node = half_ends.len();
        half_ends.push(HalfEnd {
            unitig_index: u32::try_from(unitig_index).expect("unitig index exceeds u32"),
        });
        if let Some(endpoint) = right_endpoint {
            endpoint_writers.write_record(&StitchEndpointRecord {
                vertex: endpoint.vertex,
                side: endpoint.side,
                node: stitch_node(right_node),
            })?;
        }
    }
    let half_end_elapsed = build_started.elapsed();
    report_process_memory("external stitch after half-ends");

    let mut join_neighbor = vec![STITCH_NO_NODE; half_ends.len()];
    let endpoint_sort_started = Instant::now();
    let endpoint_manifest = endpoint_writers.finish()?;
    let endpoint_sort_elapsed = endpoint_sort_started.elapsed();
    let endpoint_join_started = Instant::now();
    join_neighbors_from_endpoint_buckets::<K>(&endpoint_manifest, &mut join_neighbor, threads)?;
    let endpoint_join_elapsed = endpoint_join_started.elapsed();
    report_process_memory("external stitch after endpoint join");

    let walk_started = Instant::now();
    let component_elapsed;
    let emitted = if use_direct_label_reads {
        let component_started = Instant::now();
        let component_starts = stitch_component_starts_from_neighbors(&half_ends, &join_neighbor);
        component_elapsed = component_started.elapsed();
        emit_external_stitched_labels_to_final_buckets::<K>(
            inputs,
            threads,
            &label_refs,
            &half_ends,
            &join_neighbor,
            &component_starts,
            final_buckets,
        )?
    } else if use_ordered_labels {
        let component_started = Instant::now();
        let component_starts = stitch_component_starts_from_neighbors(&half_ends, &join_neighbor);
        component_elapsed = component_started.elapsed();
        emit_external_ordered_stitched_labels_to_final_buckets::<K>(
            inputs,
            coord_dir,
            &label_refs,
            &half_ends,
            &join_neighbor,
            &component_starts,
            final_buckets,
        )?
    } else {
        let manifest =
            if std::env::var_os("CF3_RS_NEIGHBOR_STITCH_PATH_INFO").is_some() {
                let component_started = Instant::now();
                let manifest = write_external_materialized_stitched_coord_buckets_from_neighbors::<
                    K,
                >(
                    inputs, coord_dir, threads, &half_ends, &join_neighbor
                )?;
                component_elapsed = component_started.elapsed();
                report_process_memory("external stitch after materialized coord write");
                manifest
            } else {
                let component_started = Instant::now();
                let component_starts =
                    stitch_component_starts_from_neighbors(&half_ends, &join_neighbor);
                component_elapsed = component_started.elapsed();
                write_external_materialized_stitched_coord_buckets::<K>(
                    inputs,
                    coord_dir,
                    threads,
                    &half_ends,
                    &join_neighbor,
                    &component_starts,
                )?
            };
        report_process_memory("external stitch before materialized coord reduce");
        let reduce_started = Instant::now();
        let emitted = reduce_materialized_stitched_coord_bucket_files_to_final::<K>(
            &manifest,
            &[],
            threads,
            final_buckets,
        )?;
        report_process_memory("external stitch after materialized coord reduce");
        eprintln!(
            "cuttlefish3-rs: external materialized stitch reduce {:.3}s",
            reduce_started.elapsed().as_secs_f64()
        );
        emitted
    };
    let walk_elapsed = walk_started.elapsed();

    eprintln!(
        "cuttlefish3-rs: external stitch detail: half-ends {:.3}s, endpoint sort {:.3}s, endpoint join {:.3}s, components {:.3}s, stitch label materialization/reduce {:.3}s",
        half_end_elapsed.as_secs_f64(),
        endpoint_sort_elapsed.as_secs_f64(),
        endpoint_join_elapsed.as_secs_f64(),
        component_elapsed.as_secs_f64(),
        walk_elapsed.as_secs_f64(),
    );
    Ok(emitted)
}

struct ExternalCppPathInfoCollation {
    stitched_unitigs: u64,
    direct_local_unitigs: u64,
    direct_local_unitigs_complete: bool,
}

fn collate_external_cpp_path_info_to_final_buckets<const K: usize>(
    inputs: &mut ExternalDiscontinuityInputs<K>,
    threads: usize,
    coord_dir: &Path,
    final_buckets: &mut FinalUnitigBucketWriters,
) -> Result<ExternalCppPathInfoCollation, SerialCollationError> {
    eprintln!(
        "cuttlefish3-rs: deriving discontinuity path-info by C++-style contraction/expansion"
    );

    let matrix_started = Instant::now();
    let matrix = inputs
        .edge_matrix
        .take()
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;
    let matrix_dir = matrix.dir.clone();
    eprintln!(
        "cuttlefish3-rs: discontinuity edge matrix handed off in {:.3}s; {} edge(s), {} phi edge(s), {} diagonal edge(s)",
        matrix_started.elapsed().as_secs_f64(),
        matrix.stats.edges,
        matrix.stats.phi_edges,
        matrix.stats.diagonal_edges
    );

    let contract_started = Instant::now();
    report_process_memory("before C++-style discontinuity contraction");
    let contraction = SerialDiscontinuityContractor::contract_blocked_external(matrix, threads)?;
    report_process_memory("after C++-style discontinuity contraction");
    eprintln!(
        "cuttlefish3-rs: discontinuity graph contracted in {:.3}s; {} meta-vertices, {} final edge(s), {} reinserted edge(s)",
        contract_started.elapsed().as_secs_f64(),
        contraction.stats.meta_vertices,
        contraction.stats.final_edges,
        contraction.stats.reinserted_edges
    );
    let (ranges_per_bucket, path_info_bucket_count) = if let Some(buckets) =
        inputs.local_unitig_buckets.as_ref()
    {
        (1, buckets.len().max(1))
    } else {
        let ranges_per_bucket = ranges_per_path_info_bucket(inputs.ranges.len(), threads.max(1));
        (
            ranges_per_bucket,
            inputs.ranges.len().div_ceil(ranges_per_bucket).max(1),
        )
    };
    let expansion_dir = coord_dir.with_file_name(format!(
        "{}.cpp-expansion",
        coord_dir
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("coord")
    ));
    let expand_started = Instant::now();
    report_process_memory("before C++-style discontinuity expansion");
    let expansion = SerialDiscontinuityExpander::expand_cpp_ordered_external_to_range_buckets(
        &contraction,
        &inputs.ranges,
        ranges_per_bucket,
        path_info_bucket_count,
        &inputs.unitig_path,
        &expansion_dir,
        threads,
    )?;
    report_process_memory("after C++-style discontinuity expansion");
    eprintln!(
        "cuttlefish3-rs: contracted graph expanded in {:.3}s; {} path-info edge record(s), {} unresolved edge(s)",
        expand_started.elapsed().as_secs_f64(),
        expansion.stats.edge_path_infos,
        expansion.stats.unresolved_edges
    );
    drop(contraction);
    let _ = fs::remove_dir_all(&matrix_dir);
    trim_process_allocations();
    report_process_memory("after discontinuity contraction release");

    let map_started = Instant::now();
    report_process_memory("before C++-style path-info map");
    let materialized = map_external_cpp_path_info_buckets_to_max_unitig_buckets::<K>(
        inputs,
        coord_dir,
        threads,
        ranges_per_bucket,
        expansion,
        final_buckets,
    )?;
    report_process_memory("after C++-style path-info map");
    let _ = fs::remove_dir_all(&expansion_dir);
    eprintln!(
        "cuttlefish3-rs: C++-style path-info map completed in {:.3}s",
        map_started.elapsed().as_secs_f64()
    );
    let reduce_started = Instant::now();
    report_process_memory("before C++-style path-info reduce");
    let emitted = reduce_materialized_stitched_coord_bucket_files_to_final::<K>(
        &materialized.manifest,
        &materialized.retained,
        threads,
        final_buckets,
    )?;
    report_process_memory("after C++-style path-info reduce");
    eprintln!(
        "cuttlefish3-rs: C++-style path-info reduce {:.3}s",
        reduce_started.elapsed().as_secs_f64()
    );
    Ok(ExternalCppPathInfoCollation {
        stitched_unitigs: emitted,
        direct_local_unitigs: materialized.direct_local_unitigs,
        direct_local_unitigs_complete: materialized.direct_local_unitigs_complete,
    })
}

fn prepare_unitig_blocked_edge<const K: usize>(
    unitig: &DiscontinuityUnitig<K>,
    unitig_index: usize,
    vertex_partitions: usize,
) -> Option<PreparedBlockedEdge> {
    let mut edge = match (unitig.left_exit(), unitig.right_exit()) {
        (Some(left), Some(right)) => DiscontinuityEdge {
            first: MatrixEndpoint::Vertex(left),
            second: MatrixEndpoint::Vertex(right),
            weight: 1,
            unitig_bucket: 0,
            unitig_index,
            unitig_exit_side: Side::Back,
            phantom_unitig: None,
            swapped: false,
        },
        (Some(endpoint), None) => DiscontinuityEdge {
            first: MatrixEndpoint::Phi,
            second: MatrixEndpoint::Vertex(endpoint),
            weight: 1,
            unitig_bucket: 0,
            unitig_index,
            unitig_exit_side: Side::Front,
            phantom_unitig: None,
            swapped: false,
        },
        (None, Some(endpoint)) => DiscontinuityEdge {
            first: MatrixEndpoint::Phi,
            second: MatrixEndpoint::Vertex(endpoint),
            weight: 1,
            unitig_bucket: 0,
            unitig_index,
            unitig_exit_side: Side::Back,
            phantom_unitig: None,
            swapped: false,
        },
        (None, None) => return None,
    };
    let first_partition = edge_matrix_partition(vertex_partitions, edge.first);
    let second_partition = edge_matrix_partition(vertex_partitions, edge.second);
    if first_partition > second_partition {
        std::mem::swap(&mut edge.first, &mut edge.second);
        edge.unitig_exit_side = edge.unitig_exit_side.inverse();
        edge.swapped = true;
    }
    let row = first_partition.min(second_partition);
    let col = first_partition.max(second_partition);
    Some(PreparedBlockedEdge {
        block: row * (vertex_partitions + 1) + col,
        bytes: encode_discontinuity_edge(&edge),
        phi: row == 0,
        diagonal: row == col,
    })
}

fn prepare_existing_blocked_edge<const K: usize>(
    mut edge: DiscontinuityEdge<K>,
    vertex_partitions: usize,
) -> PreparedBlockedEdge {
    let first_partition = edge_matrix_partition(vertex_partitions, edge.first);
    let second_partition = edge_matrix_partition(vertex_partitions, edge.second);
    if first_partition > second_partition {
        std::mem::swap(&mut edge.first, &mut edge.second);
        edge.unitig_exit_side = edge.unitig_exit_side.inverse();
        edge.swapped = !edge.swapped;
    }
    let row = first_partition.min(second_partition);
    let col = first_partition.max(second_partition);
    PreparedBlockedEdge {
        block: row * (vertex_partitions + 1) + col,
        bytes: encode_discontinuity_edge(&edge),
        phi: row == 0,
        diagonal: row == col,
    }
}

fn emit_contracted_edge_chunks<const K: usize>(
    edges: Vec<DiscontinuityEdge<K>>,
    vertex_partitions: usize,
    max_partition_exclusive: Option<usize>,
    appenders: &ConcurrentBlockedEdgeWriters,
) -> Result<u64, SerialCollationError> {
    const EDGES_PER_BATCH: usize = 8 * 1024;
    edges
        .par_chunks(EDGES_PER_BATCH)
        .map(|chunk| {
            let mut prepared = Vec::with_capacity(chunk.len());
            for edge in chunk {
                if let Some(limit) = max_partition_exclusive {
                    let first = edge_matrix_partition(vertex_partitions, edge.first);
                    let second = edge_matrix_partition(vertex_partitions, edge.second);
                    if first.max(second) >= limit {
                        continue;
                    }
                }
                prepared.push(prepare_existing_blocked_edge(
                    edge.clone(),
                    vertex_partitions,
                ));
            }
            prepared.sort_unstable_by_key(|edge| edge.block);
            let mut start = 0;
            while start < prepared.len() {
                let block = prepared[start].block;
                let mut end = start + 1;
                while end < prepared.len() && prepared[end].block == block {
                    end += 1;
                }
                appenders.add_batch(&prepared[start..end])?;
                start = end;
            }
            Ok::<_, SerialCollationError>(prepared.len() as u64)
        })
        .try_reduce(|| 0, |left, right| Ok(left + right))
}

fn emit_contracted_edge_batch<const K: usize>(
    edges: Vec<DiscontinuityEdge<K>>,
    vertex_partitions: usize,
    appenders: &ConcurrentBlockedEdgeWriters,
) -> Result<u64, SerialCollationError> {
    let mut prepared = edges
        .into_iter()
        .map(|edge| prepare_existing_blocked_edge(edge, vertex_partitions))
        .collect::<Vec<_>>();
    prepared.sort_unstable_by_key(|edge| edge.block);
    let mut start = 0;
    while start < prepared.len() {
        let block = prepared[start].block;
        let mut end = start + 1;
        while end < prepared.len() && prepared[end].block == block {
            end += 1;
        }
        appenders.add_batch(&prepared[start..end])?;
        start = end;
    }
    Ok(prepared.len() as u64)
}

fn emit_prepared_edge_batch(
    prepared: &mut [PreparedBlockedEdge],
    appenders: &ConcurrentBlockedEdgeWriters,
) -> Result<u64, SerialCollationError> {
    prepared.sort_unstable_by_key(|edge| edge.block);
    let mut start = 0;
    while start < prepared.len() {
        let block = prepared[start].block;
        let mut end = start + 1;
        while end < prepared.len() && prepared[end].block == block {
            end += 1;
        }
        appenders.add_batch(&prepared[start..end])?;
        start = end;
    }
    Ok(prepared.len() as u64)
}

fn serial_collation_to_input_error(error: SerialCollationError) -> DiscontinuityInputError {
    match error {
        SerialCollationError::Io { path, source } => DiscontinuityInputError::Io { path, source },
        SerialCollationError::MalformedCoordBucket(path) => DiscontinuityInputError::Io {
            path,
            source: std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "malformed blocked edge matrix",
            ),
        },
        SerialCollationError::WorkerPanic => DiscontinuityInputError::WorkerPanic,
        SerialCollationError::Color(err) => DiscontinuityInputError::Color(err),
    }
}

struct ExternalCppPathInfoMaterialized {
    manifest: Vec<MaterializedStitchedCoordBucketEntry>,
    retained: Vec<Vec<PendingMaterializedBucket>>,
    direct_local_unitigs: u64,
    direct_local_unitigs_complete: bool,
}

struct RangeBucketedExpansion<const K: usize> {
    stats: SerialExpansionStats,
    records: Vec<Vec<StitchedCoordRecord>>,
    record_manifest: Vec<StitchedCoordBucketEntry>,
    phantom_records: Vec<(StitchedCoordRecord, DiscontinuityEndpoint<K>)>,
}

const VERTEX_PATH_INFO_WRITE_BUFFER: usize = 1024 * 1024;

#[inline]
const fn vertex_path_info_record_len<const K: usize>() -> usize {
    if K <= 31 {
        std::mem::size_of::<CompactVertexPathInfoRecord>()
    } else {
        2 * discontinuity_edge_kmer_bytes::<K>() + 9
    }
}

#[repr(C)]
#[derive(Clone, Copy)]
struct CompactVertexPathInfoRecord {
    vertex: u64,
    path_id: u64,
    rank_and_flags: u64,
}

const _: () = assert!(std::mem::size_of::<CompactVertexPathInfoRecord>() == 24);

struct VertexPathInfoBucketWriters<const K: usize> {
    dir: PathBuf,
    writers: Vec<Option<BufWriter<File>>>,
    buffers: Vec<Vec<u8>>,
    phantom: PhantomData<[(); K]>,
}

impl<const K: usize> VertexPathInfoBucketWriters<K> {
    fn create(dir: &Path, bucket_count: usize) -> Result<Self, SerialCollationError> {
        if dir.exists() {
            fs::remove_dir_all(dir).map_err(|source| SerialCollationError::Io {
                path: dir.to_path_buf(),
                source,
            })?;
        }
        fs::create_dir_all(dir).map_err(|source| SerialCollationError::Io {
            path: dir.to_path_buf(),
            source,
        })?;
        let mut writers = Vec::with_capacity(bucket_count);
        writers.resize_with(bucket_count, || None);
        let mut buffers = Vec::with_capacity(bucket_count);
        buffers.resize_with(bucket_count, Vec::new);
        Ok(Self {
            dir: dir.to_path_buf(),
            writers,
            buffers,
            phantom: PhantomData,
        })
    }

    fn open_existing(dir: &Path, bucket_count: usize) -> Self {
        let mut writers = Vec::with_capacity(bucket_count);
        writers.resize_with(bucket_count, || None);
        let mut buffers = Vec::with_capacity(bucket_count);
        buffers.resize_with(bucket_count, Vec::new);
        Self {
            dir: dir.to_path_buf(),
            writers,
            buffers,
            phantom: PhantomData,
        }
    }

    fn write_record(
        &mut self,
        bucket_id: usize,
        record: &VertexPathInfo<K>,
    ) -> Result<(), SerialCollationError> {
        if bucket_id >= self.writers.len() {
            return Err(SerialCollationError::MalformedCoordBucket(self.dir.clone()));
        }
        self.buffers[bucket_id].extend_from_slice(
            &encoded_vertex_path_info_record::<K>(record)[..vertex_path_info_record_len::<K>()],
        );
        if self.buffers[bucket_id].len() >= VERTEX_PATH_INFO_WRITE_BUFFER {
            self.write_buffer(bucket_id, false)?;
        }
        Ok(())
    }

    fn flush_bucket(&mut self, bucket_id: usize) -> Result<(), SerialCollationError> {
        self.write_buffer(bucket_id, true)
    }

    fn write_buffer(
        &mut self,
        bucket_id: usize,
        flush_writer: bool,
    ) -> Result<(), SerialCollationError> {
        if bucket_id >= self.writers.len() {
            return Err(SerialCollationError::MalformedCoordBucket(self.dir.clone()));
        }
        if self.buffers[bucket_id].is_empty() {
            if flush_writer && let Some(writer) = self.writers[bucket_id].as_mut() {
                writer.flush().map_err(|source| SerialCollationError::Io {
                    path: vertex_path_info_bucket_path(&self.dir, bucket_id),
                    source,
                })?;
            }
            return Ok(());
        }
        if self.writers[bucket_id].is_none() {
            let path = vertex_path_info_bucket_path(&self.dir, bucket_id);
            let file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(&path)
                .map_err(|source| SerialCollationError::Io {
                    path: path.clone(),
                    source,
                })?;
            self.writers[bucket_id] = Some(BufWriter::with_capacity(
                VERTEX_PATH_INFO_WRITE_BUFFER,
                file,
            ));
        }
        let writer = self.writers[bucket_id]
            .as_mut()
            .expect("vertex path-info writer was just created");
        writer
            .write_all(&self.buffers[bucket_id])
            .map_err(|source| SerialCollationError::Io {
                path: vertex_path_info_bucket_path(&self.dir, bucket_id),
                source,
            })?;
        if flush_writer {
            writer.flush().map_err(|source| SerialCollationError::Io {
                path: vertex_path_info_bucket_path(&self.dir, bucket_id),
                source,
            })?;
        }
        self.buffers[bucket_id].clear();
        Ok(())
    }
}

fn vertex_path_info_bucket_path(dir: &Path, bucket_id: usize) -> PathBuf {
    dir.join(format!("{bucket_id:05}.pv"))
}

fn encoded_vertex_path_info_record<const K: usize>(record: &VertexPathInfo<K>) -> [u8; 41] {
    if K <= 31 {
        let mut bytes = [0u8; 41];
        bytes[..8].copy_from_slice(&(record.vertex.as_u128() as u64).to_le_bytes());
        bytes[8..16].copy_from_slice(&(record.info.path_id.as_u128() as u64).to_le_bytes());
        let flags =
            u64::from(record.info.exit_side == Side::Back) | (u64::from(record.info.is_cycle) << 1);
        bytes[16..24].copy_from_slice(&((record.info.rank << 2) | flags).to_le_bytes());
        return bytes;
    }
    let kmer_bytes = discontinuity_edge_kmer_bytes::<K>();
    let path_off = kmer_bytes;
    let rank_off = 2 * kmer_bytes;
    let flags_off = rank_off + 8;
    let mut bytes = [0u8; 41];
    bytes[..kmer_bytes].copy_from_slice(&record.vertex.as_u128().to_le_bytes()[..kmer_bytes]);
    bytes[path_off..rank_off]
        .copy_from_slice(&record.info.path_id.as_u128().to_le_bytes()[..kmer_bytes]);
    bytes[rank_off..flags_off].copy_from_slice(&record.info.rank.to_le_bytes());
    let mut flags = 0u8;
    if record.info.exit_side == Side::Back {
        flags |= 1;
    }
    if record.info.is_cycle {
        flags |= 2;
    }
    bytes[flags_off] = flags;
    bytes
}

fn write_meta_vertex_bucket_parallel<const K: usize>(
    dir: &Path,
    bucket_id: usize,
    records: &[SerialMetaVertex<K>],
    pool: &ThreadPool,
) -> Result<(), SerialCollationError> {
    if records.is_empty() {
        return Ok(());
    }
    let path = vertex_path_info_bucket_path(dir, bucket_id);
    let file = File::create(&path).map_err(|source| SerialCollationError::Io {
        path: path.clone(),
        source,
    })?;
    let record_len = vertex_path_info_record_len::<K>();
    file.set_len((records.len() * record_len) as u64)
        .map_err(|source| SerialCollationError::Io {
            path: path.clone(),
            source,
        })?;
    let records_per_chunk = (1024 * 1024 / record_len).max(1);
    pool.install(|| {
        records
            .par_chunks(records_per_chunk)
            .enumerate()
            .try_for_each(|(chunk_id, chunk)| {
                let mut encoded = Vec::with_capacity(chunk.len() * record_len);
                for meta in chunk {
                    let record = VertexPathInfo {
                        vertex: meta.vertex,
                        info: PathInfo {
                            path_id: meta.vertex,
                            rank: meta.weight,
                            exit_side: meta.entry_side,
                            is_cycle: meta.is_cycle,
                        },
                    };
                    encoded
                        .extend_from_slice(&encoded_vertex_path_info_record(&record)[..record_len]);
                }
                let offset = (chunk_id * records_per_chunk * record_len) as u64;
                file.write_all_at(&encoded, offset)
                    .map_err(|source| SerialCollationError::Io {
                        path: path.clone(),
                        source,
                    })
            })
    })
}

#[inline(always)]
fn append_encoded_compact_vertex_path_info<const K: usize>(
    output: &mut Vec<u8>,
    vertex: Kmer<K>,
    info: CompactExpansionPathInfo,
) {
    debug_assert!(K <= 31);
    output.extend_from_slice(&(vertex.as_u128() as u64).to_le_bytes());
    output.extend_from_slice(&info.path_id.to_le_bytes());
    output.extend_from_slice(&info.rank_and_flags.to_le_bytes());
}

fn decoded_vertex_path_info_record<const K: usize>(bytes: &[u8]) -> VertexPathInfo<K> {
    if K <= 31 {
        let vertex = u64::from_le_bytes(bytes[..8].try_into().expect("u64 vertex field"));
        let path_id = u64::from_le_bytes(bytes[8..16].try_into().expect("u64 path field"));
        let rank_and_flags = u64::from_le_bytes(bytes[16..24].try_into().expect("u64 rank field"));
        return VertexPathInfo {
            vertex: Kmer::from_bits(vertex as u128),
            info: PathInfo {
                path_id: Kmer::from_bits(path_id as u128),
                rank: rank_and_flags >> 2,
                exit_side: if rank_and_flags & 1 == 0 {
                    Side::Front
                } else {
                    Side::Back
                },
                is_cycle: rank_and_flags & 2 != 0,
            },
        };
    }
    let kmer_bytes = discontinuity_edge_kmer_bytes::<K>();
    let path_off = kmer_bytes;
    let rank_off = 2 * kmer_bytes;
    let flags_off = rank_off + 8;
    let mut vertex = [0u8; 16];
    vertex[..kmer_bytes].copy_from_slice(&bytes[..kmer_bytes]);
    let mut path_id = [0u8; 16];
    path_id[..kmer_bytes].copy_from_slice(&bytes[path_off..rank_off]);
    let mut rank = [0u8; 8];
    rank.copy_from_slice(&bytes[rank_off..flags_off]);
    let flags = bytes[flags_off];
    VertexPathInfo {
        vertex: Kmer::from_bits(u128::from_le_bytes(vertex)),
        info: PathInfo {
            path_id: Kmer::from_bits(u128::from_le_bytes(path_id)),
            rank: u64::from_le_bytes(rank),
            exit_side: if flags & 1 == 0 {
                Side::Front
            } else {
                Side::Back
            },
            is_cycle: flags & 2 != 0,
        },
    }
}

fn read_vertex_path_info_bucket_into<const K: usize>(
    dir: &Path,
    bucket_id: usize,
    malformed_path: &Path,
    map: &ExpansionPathInfoTable<K>,
    pool: &ThreadPool,
    read_elapsed: &mut Duration,
    insert_elapsed: &mut Duration,
) -> Result<(), SerialCollationError> {
    let path = vertex_path_info_bucket_path(dir, bucket_id);
    if !path.exists() {
        return Ok(());
    }
    let record_len = vertex_path_info_record_len::<K>();
    let byte_len = fs::metadata(&path)
        .map_err(|source| SerialCollationError::Io {
            path: path.clone(),
            source,
        })?
        .len() as usize;
    if byte_len % record_len != 0 {
        return Err(SerialCollationError::MalformedCoordBucket(
            malformed_path.to_path_buf(),
        ));
    }
    if K <= 31 {
        let records_per_chunk = (1024 * 1024 / record_len).max(1);
        let record_count = byte_len / record_len;
        let next_record = AtomicUsize::new(0);
        let file = File::open(&path).map_err(|source| SerialCollationError::Io {
            path: path.clone(),
            source,
        })?;
        let worker_count = pool.current_num_threads().max(1).min(record_count.max(1));
        let worker_times = pool.install(|| {
            (0..worker_count)
                .into_par_iter()
                .map(|_| {
                    let empty = CompactVertexPathInfoRecord {
                        vertex: 0,
                        path_id: 0,
                        rank_and_flags: 0,
                    };
                    let mut records = vec![empty; records_per_chunk];
                    let mut read_time = Duration::default();
                    let mut insert_time = Duration::default();
                    loop {
                        let start = next_record.fetch_add(records_per_chunk, Ordering::Relaxed);
                        if start >= record_count {
                            break;
                        }
                        let count = records_per_chunk.min(record_count - start);
                        let bytes = unsafe {
                            std::slice::from_raw_parts_mut(
                                records.as_mut_ptr().cast::<u8>(),
                                count * record_len,
                            )
                        };
                        let started = Instant::now();
                        file.read_exact_at(bytes, (start * record_len) as u64)
                            .map_err(|source| SerialCollationError::Io {
                                path: path.clone(),
                                source,
                            })?;
                        read_time += started.elapsed();
                        let started = Instant::now();
                        for &record in &records[..count] {
                            map.insert_compact_record(record);
                        }
                        insert_time += started.elapsed();
                    }
                    Ok::<_, SerialCollationError>((read_time, insert_time))
                })
                .collect::<Result<Vec<_>, SerialCollationError>>()
        })?;
        let total_read = worker_times.iter().map(|times| times.0).sum::<Duration>();
        let total_insert = worker_times.iter().map(|times| times.1).sum::<Duration>();
        *read_elapsed += total_read / worker_count as u32;
        *insert_elapsed += total_insert / worker_count as u32;
        return Ok(());
    }
    let phase = Instant::now();
    let bytes = fs::read(&path).map_err(|source| SerialCollationError::Io {
        path: path.clone(),
        source,
    })?;
    *read_elapsed += phase.elapsed();
    let records_per_chunk = (1024 * 1024 / record_len).max(1);
    let phase = Instant::now();
    pool.install(|| {
        bytes
            .par_chunks(records_per_chunk * record_len)
            .for_each(|block| {
                for chunk in block.chunks_exact(record_len) {
                    map.insert_encoded(chunk);
                }
            });
    });
    *insert_elapsed += phase.elapsed();
    Ok(())
}

fn push_edge_path_record_to_range_bucket_writer<const K: usize>(
    edge: &DiscontinuityEdge<K>,
    info: PathInfo<K>,
    ranges: &[ExternalLocalUnitigRange],
    range_index: &ExternalRangeIndex,
    ranges_per_bucket: usize,
    writers: &ConcurrentStitchedCoordWriters<'_>,
    phantom_records: &mut Vec<(StitchedCoordRecord, DiscontinuityEndpoint<K>)>,
    error_path: &Path,
) -> Result<(), SerialCollationError> {
    let record = stitched_record_from_edge_path_info(edge, info, error_path)?;
    if let Some(phantom) = edge.phantom_unitig {
        phantom_records.push((record, phantom));
        return Ok(());
    }
    let bucket_id = edge_path_info_bucket(edge, ranges, range_index, ranges_per_bucket)
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(error_path.to_path_buf()))?;
    writers.write_record(bucket_id, record)
}

fn push_edge_path_record_to_range_buffers<const K: usize>(
    edge: &DiscontinuityEdge<K>,
    info: PathInfo<K>,
    ranges: &[ExternalLocalUnitigRange],
    range_index: &ExternalRangeIndex,
    ranges_per_bucket: usize,
    buffers: &mut [Vec<StitchedCoordRecord>],
    phantom_records: &mut Vec<(StitchedCoordRecord, DiscontinuityEndpoint<K>)>,
    error_path: &Path,
) -> Result<(), SerialCollationError> {
    let record = stitched_record_from_edge_path_info(edge, info, error_path)?;
    if let Some(phantom) = edge.phantom_unitig {
        phantom_records.push((record, phantom));
        return Ok(());
    }
    let bucket_id = edge_path_info_bucket(edge, ranges, range_index, ranges_per_bucket)
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(error_path.to_path_buf()))?;
    buffers[bucket_id].push(record);
    Ok(())
}

fn push_compact_edge_path_record_to_range_buffers<const K: usize>(
    edge: &DiscontinuityEdge<K>,
    info: CompactExpansionPathInfo,
    ranges: &[ExternalLocalUnitigRange],
    range_index: &ExternalRangeIndex,
    ranges_per_bucket: usize,
    buffers: &mut [Vec<StitchedCoordRecord>],
    phantom_records: &mut Vec<(StitchedCoordRecord, DiscontinuityEndpoint<K>)>,
    error_path: &Path,
) -> Result<(), SerialCollationError> {
    let record = stitched_record_from_compact_edge_path_info(edge, info, error_path)?;
    if let Some(phantom) = edge.phantom_unitig {
        phantom_records.push((record, phantom));
        return Ok(());
    }
    let bucket_id = edge_path_info_bucket(edge, ranges, range_index, ranges_per_bucket)
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(error_path.to_path_buf()))?;
    buffers[bucket_id].push(record);
    Ok(())
}

#[inline]
fn edge_path_info_bucket<const K: usize>(
    edge: &DiscontinuityEdge<K>,
    ranges: &[ExternalLocalUnitigRange],
    range_index: &ExternalRangeIndex,
    ranges_per_bucket: usize,
) -> Option<usize> {
    if edge.unitig_bucket != 0 {
        Some(edge.unitig_bucket as usize - 1)
    } else {
        range_index
            .find(ranges, edge.unitig_index)
            .map(|range_id| range_id / ranges_per_bucket)
    }
}

fn stitched_record_from_edge_path_info<const K: usize>(
    edge: &DiscontinuityEdge<K>,
    info: PathInfo<K>,
    error_path: &Path,
) -> Result<StitchedCoordRecord, SerialCollationError> {
    let path_id = u64::try_from(info.path_id.as_u128())
        .map_err(|_| SerialCollationError::MalformedCoordBucket(error_path.to_path_buf()))?;
    let unitig_index = u32::try_from(edge.unitig_index)
        .map_err(|_| SerialCollationError::MalformedCoordBucket(error_path.to_path_buf()))?;
    Ok(StitchedCoordRecord {
        path_id,
        rank: info.rank,
        unitig_index,
        reverse: info.exit_side == Side::Front,
        is_cycle: info.is_cycle,
    })
}

#[inline(always)]
fn stitched_record_from_compact_edge_path_info<const K: usize>(
    edge: &DiscontinuityEdge<K>,
    info: CompactExpansionPathInfo,
    error_path: &Path,
) -> Result<StitchedCoordRecord, SerialCollationError> {
    let unitig_index = u32::try_from(edge.unitig_index)
        .map_err(|_| SerialCollationError::MalformedCoordBucket(error_path.to_path_buf()))?;
    Ok(StitchedCoordRecord {
        path_id: info.path_id,
        rank: info.rank(),
        unitig_index,
        reverse: info.exit_side() == Side::Front,
        is_cycle: info.is_cycle(),
    })
}

#[allow(dead_code)]
fn write_external_cpp_path_info_materialized_coord_buckets<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    coord_dir: &Path,
    threads: usize,
    expansion: &SerialExpansion<K>,
    final_buckets: &mut FinalUnitigBucketWriters,
) -> Result<ExternalCppPathInfoMaterialized, SerialCollationError> {
    let _ = final_buckets;
    write_external_cpp_path_info_materialized_coord_buckets_bucketed::<K>(
        inputs, coord_dir, threads, expansion,
    )
}

#[allow(dead_code)]
fn write_external_cpp_path_info_materialized_coord_buckets_from_bucketed<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    coord_dir: &Path,
    threads: usize,
    ranges_per_bucket: usize,
    expansion: RangeBucketedExpansion<K>,
) -> Result<ExternalCppPathInfoMaterialized, SerialCollationError> {
    if coord_dir.exists() {
        fs::remove_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
            path: coord_dir.to_path_buf(),
            source,
        })?;
    }
    fs::create_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
        path: coord_dir.to_path_buf(),
        source,
    })?;

    let max_unitig_bucket_count = materialized_stitched_coord_bucket_count(threads);
    let max_unitig_bucket_mask = max_unitig_bucket_count - 1;
    let mut phantom_writers = MaterializedStitchedCoordShardWriters::new(
        coord_dir,
        threads.max(1),
        max_unitig_bucket_count,
    );
    for (record, phantom) in expansion.phantom_records {
        let label = phantom_unitig_label(phantom);
        let bucket_id = stitched_coord_bucket(record.path_id, max_unitig_bucket_mask);
        phantom_writers.write_materialized_record(bucket_id, &record, &label)?;
    }

    let path_info_records = expansion.records;
    let workers = threads.max(1).min(path_info_records.len().max(1));
    let mut manifest = if workers == 1 || path_info_records.len() < 2 {
        let mut writers =
            MaterializedStitchedCoordShardWriters::new(coord_dir, 0, max_unitig_bucket_count);
        for (bucket_id, records) in path_info_records.iter().enumerate() {
            materialize_external_stitched_coord_range_group_from_records::<K>(
                inputs,
                max_unitig_bucket_mask,
                ranges_per_bucket,
                bucket_id,
                records,
                &mut writers,
            )?;
        }
        writers.finish()?
    } else {
        let next_bucket = AtomicUsize::new(0);
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for worker_id in 0..workers {
                let next_bucket = &next_bucket;
                let path_info_records = &path_info_records;
                handles.push(scope.spawn(move || {
                    let mut writers = MaterializedStitchedCoordShardWriters::new(
                        coord_dir,
                        worker_id,
                        max_unitig_bucket_count,
                    );
                    loop {
                        let bucket_id = next_bucket.fetch_add(1, Ordering::Relaxed);
                        let Some(records) = path_info_records.get(bucket_id) else {
                            break;
                        };
                        materialize_external_stitched_coord_range_group_from_records::<K>(
                            inputs,
                            max_unitig_bucket_mask,
                            ranges_per_bucket,
                            bucket_id,
                            records,
                            &mut writers,
                        )?;
                    }
                    writers.finish()
                }));
            }

            let mut manifest = Vec::new();
            for handle in handles {
                manifest.extend(
                    handle
                        .join()
                        .map_err(|_| SerialCollationError::WorkerPanic)??,
                );
            }
            Ok::<_, SerialCollationError>(manifest)
        })?
    };
    manifest.extend(phantom_writers.finish()?);
    manifest.sort_by(|left, right| {
        left.bucket_id
            .cmp(&right.bucket_id)
            .then_with(|| left.coord_path.cmp(&right.coord_path))
    });
    Ok(ExternalCppPathInfoMaterialized {
        manifest,
        retained: Vec::new(),
        direct_local_unitigs: 0,
        direct_local_unitigs_complete: false,
    })
}

fn map_external_cpp_path_info_buckets_to_max_unitig_buckets<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    coord_dir: &Path,
    threads: usize,
    ranges_per_bucket: usize,
    expansion: RangeBucketedExpansion<K>,
    final_buckets: &mut FinalUnitigBucketWriters,
) -> Result<ExternalCppPathInfoMaterialized, SerialCollationError> {
    const DIRECT_FINAL_BATCH_RECORDS: usize = 4096;
    if coord_dir.exists() {
        fs::remove_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
            path: coord_dir.to_path_buf(),
            source,
        })?;
    }
    fs::create_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
        path: coord_dir.to_path_buf(),
        source,
    })?;

    let max_unitig_bucket_count = 1024;
    eprintln!(
        "cuttlefish3-rs: materializing final coordinates into {max_unitig_bucket_count} bucket(s)"
    );
    let max_unitig_bucket_mask = max_unitig_bucket_count - 1;
    let shared_writers = SharedMaterializedWriters::new(coord_dir, max_unitig_bucket_count);
    let mut phantom_writers =
        SharedMaterializedBatch::new(&shared_writers, max_unitig_bucket_count);
    for (record, phantom) in expansion.phantom_records {
        let label = phantom_unitig_label(phantom);
        let bucket_id = stitched_coord_bucket(record.path_id, max_unitig_bucket_mask);
        phantom_writers.write_materialized_record(bucket_id, &record, &label)?;
    }
    phantom_writers.finish()?;

    let mut manifest = Vec::new();
    let mut direct_local_unitigs = 0u64;
    let path_info_manifest = expansion.record_manifest;
    if let Some(local_buckets) = inputs.local_unitig_buckets.as_ref() {
        direct_local_unitigs = map_local_unitig_buckets_to_max_unitig_buckets::<K>(
            inputs,
            local_buckets,
            path_info_manifest,
            threads,
            max_unitig_bucket_mask,
            &shared_writers,
            final_buckets,
        )?;
        let materialized = shared_writers.finish()?;
        manifest = materialized.manifest;
        manifest.sort_by(|left, right| {
            left.bucket_id
                .cmp(&right.bucket_id)
                .then_with(|| left.coord_path.cmp(&right.coord_path))
        });
        return Ok(ExternalCppPathInfoMaterialized {
            manifest,
            retained: materialized.retained,
            direct_local_unitigs,
            direct_local_unitigs_complete: true,
        });
    }
    if !path_info_manifest.is_empty() {
        let path_info_bucket_count = inputs.ranges.len().div_ceil(ranges_per_bucket).max(1);
        let mut groups = vec![Vec::<StitchedCoordBucketEntry>::new(); path_info_bucket_count];
        for entry in path_info_manifest {
            if entry.bucket_id >= groups.len() {
                return Err(SerialCollationError::MalformedCoordBucket(entry.path));
            }
            groups[entry.bucket_id].push(entry);
        }

        let workers = threads.max(1).min(groups.len().max(1));
        if workers == 1 || groups.len() < 2 {
            let mut writers =
                SharedMaterializedBatch::new(&shared_writers, max_unitig_bucket_count);
            for (bucket_id, group) in groups.iter().enumerate() {
                let records = read_stitched_coord_bucket_group(group)?;
                let mut emit_direct =
                    |unitig: FinalUnitigRecord| -> Result<(), SerialCollationError> {
                        if unitig.colors.is_empty() {
                            final_buckets.write_label(&unitig.label)?;
                        } else {
                            final_buckets.write_colored_label(&unitig.label, &unitig.colors)?;
                        }
                        direct_local_unitigs += 1;
                        Ok(())
                    };
                map_external_cpp_path_info_range_bucket_owned::<K, _, _>(
                    inputs,
                    max_unitig_bucket_mask,
                    ranges_per_bucket,
                    bucket_id,
                    records,
                    &mut writers,
                    &mut emit_direct,
                )?;
            }
            writers.finish()?;
        } else {
            let next_bucket = AtomicUsize::new(0);
            let next_direct_record = AtomicU64::new(final_buckets.direct_record_id_highwater + 1);
            let (tx, rx) =
                mpsc::sync_channel::<Result<EncodedFinalBatch, SerialCollationError>>(workers * 2);
            let mut first_error = None;
            std::thread::scope(|scope| {
                let mut handles = Vec::new();
                for _worker_id in 0..workers {
                    let next_bucket = &next_bucket;
                    let next_direct_record = &next_direct_record;
                    let groups = &groups;
                    let tx = tx.clone();
                    let shared_writers = &shared_writers;
                    handles.push(scope.spawn(move || {
                        let mut writers =
                            SharedMaterializedBatch::new(&shared_writers, max_unitig_bucket_count);
                        let mut direct_batch = Vec::with_capacity(DIRECT_FINAL_BATCH_RECORDS);
                        let mut emit_direct =
                            |unitig: FinalUnitigRecord| -> Result<(), SerialCollationError> {
                                direct_batch.push(unitig);
                                if direct_batch.len() >= DIRECT_FINAL_BATCH_RECORDS {
                                    let batch = std::mem::take(&mut direct_batch);
                                    let first_record = next_direct_record
                                        .fetch_add(batch.len() as u64, Ordering::Relaxed);
                                    tx.send(Ok(encode_final_unitig_batch(batch, first_record)))
                                        .map_err(|_| SerialCollationError::WorkerPanic)?;
                                }
                                Ok(())
                            };

                        loop {
                            let bucket_id = next_bucket.fetch_add(1, Ordering::Relaxed);
                            let Some(group) = groups.get(bucket_id) else {
                                break;
                            };
                            let records = read_stitched_coord_bucket_group(group)?;
                            map_external_cpp_path_info_range_bucket_owned::<K, _, _>(
                                inputs,
                                max_unitig_bucket_mask,
                                ranges_per_bucket,
                                bucket_id,
                                records,
                                &mut writers,
                                &mut emit_direct,
                            )?;
                        }
                        if !direct_batch.is_empty() {
                            let first_record = next_direct_record
                                .fetch_add(direct_batch.len() as u64, Ordering::Relaxed);
                            tx.send(Ok(encode_final_unitig_batch(direct_batch, first_record)))
                                .map_err(|_| SerialCollationError::WorkerPanic)?;
                        }
                        writers.finish()
                    }));
                }
                drop(tx);

                for result in rx {
                    match result {
                        Ok(batch) if first_error.is_none() => {
                            if let Err(err) = final_buckets.write_direct_batch(
                                &batch.bytes,
                                batch.records,
                                batch.bases,
                            ) {
                                first_error = Some(err);
                            } else {
                                direct_local_unitigs += batch.records;
                            }
                        }
                        Ok(_) => {}
                        Err(err) if first_error.is_none() => first_error = Some(err),
                        Err(_) => {}
                    }
                }

                for handle in handles {
                    handle
                        .join()
                        .map_err(|_| SerialCollationError::WorkerPanic)??;
                }
                Ok::<_, SerialCollationError>(())
            })?;
            final_buckets.direct_record_id_highwater =
                next_direct_record.load(Ordering::Relaxed).saturating_sub(1);
            if let Some(err) = first_error {
                return Err(err);
            }
        }

        let materialized = shared_writers.finish()?;
        manifest = materialized.manifest;
        manifest.sort_by(|left, right| {
            left.bucket_id
                .cmp(&right.bucket_id)
                .then_with(|| left.coord_path.cmp(&right.coord_path))
        });
        return Ok(ExternalCppPathInfoMaterialized {
            manifest,
            retained: materialized.retained,
            direct_local_unitigs,
            direct_local_unitigs_complete: true,
        });
    }

    let path_info_records = expansion.records;
    let workers = threads.max(1).min(path_info_records.len().max(1));

    if workers == 1 || path_info_records.len() < 2 {
        let mut writers =
            MaterializedStitchedCoordShardWriters::new(coord_dir, 0, max_unitig_bucket_count);
        for (bucket_id, records) in path_info_records.iter().enumerate() {
            let mut emit_direct = |unitig: FinalUnitigRecord| -> Result<(), SerialCollationError> {
                if unitig.colors.is_empty() {
                    final_buckets.write_label(&unitig.label)?;
                } else {
                    final_buckets.write_colored_label(&unitig.label, &unitig.colors)?;
                }
                direct_local_unitigs += 1;
                Ok(())
            };
            map_external_cpp_path_info_range_bucket::<K, _, _>(
                inputs,
                max_unitig_bucket_mask,
                ranges_per_bucket,
                bucket_id,
                records,
                &mut writers,
                &mut emit_direct,
            )?;
        }
        manifest.extend(writers.finish()?);
    } else {
        let next_bucket = AtomicUsize::new(0);
        let (tx, rx) =
            mpsc::sync_channel::<Result<Vec<FinalUnitigRecord>, SerialCollationError>>(workers * 2);
        let mut first_error = None;
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for worker_id in 0..workers {
                let next_bucket = &next_bucket;
                let path_info_records = &path_info_records;
                let tx = tx.clone();
                handles.push(scope.spawn(move || {
                    let mut writers = MaterializedStitchedCoordShardWriters::new(
                        coord_dir,
                        worker_id,
                        max_unitig_bucket_count,
                    );
                    let mut direct_batch = Vec::with_capacity(256);
                    let mut emit_direct =
                        |unitig: FinalUnitigRecord| -> Result<(), SerialCollationError> {
                            direct_batch.push(unitig);
                            if direct_batch.len() >= 256 {
                                let batch = std::mem::take(&mut direct_batch);
                                tx.send(Ok(batch))
                                    .map_err(|_| SerialCollationError::WorkerPanic)?;
                            }
                            Ok(())
                        };

                    loop {
                        let bucket_id = next_bucket.fetch_add(1, Ordering::Relaxed);
                        let Some(records) = path_info_records.get(bucket_id) else {
                            break;
                        };
                        map_external_cpp_path_info_range_bucket::<K, _, _>(
                            inputs,
                            max_unitig_bucket_mask,
                            ranges_per_bucket,
                            bucket_id,
                            records,
                            &mut writers,
                            &mut emit_direct,
                        )?;
                    }
                    if !direct_batch.is_empty() {
                        tx.send(Ok(direct_batch))
                            .map_err(|_| SerialCollationError::WorkerPanic)?;
                    }
                    writers.finish()
                }));
            }
            drop(tx);

            for result in rx {
                match result {
                    Ok(unitigs) if first_error.is_none() => {
                        for unitig in unitigs {
                            let write = if unitig.colors.is_empty() {
                                final_buckets.write_label(&unitig.label)
                            } else {
                                final_buckets.write_colored_label(&unitig.label, &unitig.colors)
                            };
                            if let Err(err) = write {
                                first_error = Some(err);
                                break;
                            }
                            direct_local_unitigs += 1;
                        }
                    }
                    Ok(_) => {}
                    Err(err) if first_error.is_none() => first_error = Some(err),
                    Err(_) => {}
                }
            }

            for handle in handles {
                manifest.extend(
                    handle
                        .join()
                        .map_err(|_| SerialCollationError::WorkerPanic)??,
                );
            }
            Ok::<_, SerialCollationError>(())
        })?;
        if let Some(err) = first_error {
            return Err(err);
        }
    }

    let materialized = shared_writers.finish()?;
    manifest.extend(materialized.manifest);
    manifest.sort_by(|left, right| {
        left.bucket_id
            .cmp(&right.bucket_id)
            .then_with(|| left.coord_path.cmp(&right.coord_path))
    });
    Ok(ExternalCppPathInfoMaterialized {
        manifest,
        retained: materialized.retained,
        direct_local_unitigs,
        direct_local_unitigs_complete: true,
    })
}

fn map_local_unitig_buckets_to_max_unitig_buckets<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    local_buckets: &[LocalUnitigBucketEntry],
    path_info_manifest: Vec<StitchedCoordBucketEntry>,
    threads: usize,
    max_unitig_bucket_mask: usize,
    shared_writers: &SharedMaterializedWriters<'_>,
    final_buckets: &mut FinalUnitigBucketWriters,
) -> Result<u64, SerialCollationError> {
    const DIRECT_FINAL_BATCH_RECORDS: usize = 4096;
    let mut groups = (0..local_buckets.len())
        .map(|_| Vec::<StitchedCoordBucketEntry>::new())
        .collect::<Vec<_>>();
    for entry in path_info_manifest {
        let Some(group) = groups.get_mut(entry.bucket_id) else {
            return Err(SerialCollationError::MalformedCoordBucket(entry.path));
        };
        group.push(entry);
    }

    let workers = threads.max(1).min(local_buckets.len().max(1));
    let next_bucket = AtomicUsize::new(0);
    let next_direct_record = AtomicU64::new(final_buckets.direct_record_id_highwater + 1);
    let (tx, rx) =
        mpsc::sync_channel::<Result<EncodedFinalBatch, SerialCollationError>>(workers * 2);
    let mut direct_local_unitigs = 0u64;
    let mut first_error = None;
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for _ in 0..workers {
            let tx = tx.clone();
            let next_bucket = &next_bucket;
            let next_direct_record = &next_direct_record;
            let groups = &groups;
            handles.push(scope.spawn(move || {
                let mut writers =
                    SharedMaterializedBatch::new(shared_writers, max_unitig_bucket_mask + 1);
                let mut path_info = Vec::new();
                let mut direct_batch = None::<DirectFinalBatchBuilder>;
                let mut emit_direct =
                    |label: &[u8], colors: &[UnitigColor]| -> Result<(), SerialCollationError> {
                        let batch = direct_batch.get_or_insert_with(|| {
                            let first_record = next_direct_record
                                .fetch_add(DIRECT_FINAL_BATCH_RECORDS as u64, Ordering::Relaxed);
                            DirectFinalBatchBuilder::new(first_record)
                        });
                        batch.push(label, colors);
                        if batch.records as usize >= DIRECT_FINAL_BATCH_RECORDS {
                            tx.send(Ok(direct_batch.take().expect("batch exists").finish()))
                                .map_err(|_| SerialCollationError::WorkerPanic)?;
                        }
                        Ok(())
                    };
                loop {
                    let bucket_index = next_bucket.fetch_add(1, Ordering::Relaxed);
                    let Some(bucket) = local_buckets.get(bucket_index) else {
                        break;
                    };
                    read_stitched_coord_bucket_group_dense(
                        &groups[bucket_index],
                        bucket.unitigs,
                        &bucket.unitig_path,
                        &mut path_info,
                    )?;
                    map_local_unitig_bucket::<K, _, _>(
                        inputs,
                        bucket,
                        &path_info,
                        max_unitig_bucket_mask,
                        &mut writers,
                        &mut emit_direct,
                    )?;
                }
                if let Some(direct_batch) = direct_batch {
                    tx.send(Ok(direct_batch.finish()))
                        .map_err(|_| SerialCollationError::WorkerPanic)?;
                }
                writers.finish()
            }));
        }
        drop(tx);
        for result in rx {
            match result {
                Ok(batch) if first_error.is_none() => {
                    if let Err(err) =
                        final_buckets.write_direct_batch(&batch.bytes, batch.records, batch.bases)
                    {
                        first_error = Some(err);
                    } else {
                        direct_local_unitigs += batch.records;
                    }
                }
                Ok(_) => {}
                Err(err) if first_error.is_none() => first_error = Some(err),
                Err(_) => {}
            }
        }
        for handle in handles {
            if let Err(err) = handle
                .join()
                .map_err(|_| SerialCollationError::WorkerPanic)?
                && first_error.is_none()
            {
                first_error = Some(err);
            }
        }
        Ok::<_, SerialCollationError>(())
    })?;
    final_buckets.direct_record_id_highwater =
        next_direct_record.load(Ordering::Relaxed).saturating_sub(1);
    if let Some(err) = first_error {
        return Err(err);
    }
    Ok(direct_local_unitigs)
}

fn map_local_unitig_bucket<const K: usize, F, W>(
    inputs: &ExternalDiscontinuityInputs<K>,
    bucket: &LocalUnitigBucketEntry,
    path_info_by_unitig: &[DenseLocalPathInfo],
    max_unitig_bucket_mask: usize,
    writers: &mut W,
    emit_direct: &mut F,
) -> Result<(), SerialCollationError>
where
    F: FnMut(&[u8], &[UnitigColor]) -> Result<(), SerialCollationError>,
    W: MaterializedRecordSink,
{
    let unitig_file =
        File::open(&bucket.unitig_path).map_err(|source| SerialCollationError::Io {
            path: bucket.unitig_path.clone(),
            source,
        })?;
    let label_file = File::open(&bucket.label_path).map_err(|source| SerialCollationError::Io {
        path: bucket.label_path.clone(),
        source,
    })?;
    let mut unitig_input = BufReader::with_capacity(1024 * 1024, unitig_file);
    let mut label_input = BufReader::with_capacity(4 * 1024 * 1024, label_file);
    let mut color_input = bucket
        .color_runs
        .as_ref()
        .map(|sidecar| sidecar.reader_at(0))
        .transpose()?;
    let mut label = Vec::new();
    let mut colors = Vec::new();
    let mut discard = vec![0u8; 64 * 1024];

    for (unitig_index, &dense_record) in path_info_by_unitig.iter().enumerate() {
        let unitig = read_discontinuity_unitig_from_reader::<K>(
            &mut unitig_input,
            &bucket.unitig_path,
            true,
        )?;
        let has_colors = if let Some(reader) = color_input.as_mut() {
            reader.read_next_into(&mut colors)?;
            true
        } else {
            false
        };
        if let Some(record) = dense_record.to_record(unitig_index) {
            label.resize(unitig.label_len as usize, 0);
            label_input
                .read_exact(&mut label)
                .map_err(|source| SerialCollationError::Io {
                    path: bucket.label_path.clone(),
                    source,
                })?;
            let max_bucket = stitched_coord_bucket(record.path_id, max_unitig_bucket_mask);
            write_external_materialized_record(
                writers,
                max_bucket,
                &record,
                &label,
                has_colors.then_some(colors.as_slice()),
            )?;
        } else if unitig.left_exit().is_none() && unitig.right_exit().is_none() {
            label.resize(unitig.label_len as usize, 0);
            label_input
                .read_exact(&mut label)
                .map_err(|source| SerialCollationError::Io {
                    path: bucket.label_path.clone(),
                    source,
                })?;
            let reverse = reverse_complement_is_less(&label);
            if reverse {
                let reversed_label = reverse_complement_label(&label);
                let reversed_colors = if has_colors && !colors.is_empty() {
                    reverse_color_runs(&colors, (unitig.label_len as usize - K + 1) as u32)
                } else {
                    Vec::new()
                };
                emit_direct(&reversed_label, &reversed_colors)?;
            } else {
                let direct_colors = if has_colors { colors.as_slice() } else { &[] };
                emit_direct(&label, direct_colors)?;
            }
        } else {
            read_and_discard_exact(
                &mut label_input,
                &bucket.label_path,
                &mut discard,
                unitig.label_len as u64,
            )?;
        }
    }
    let _ = inputs;
    Ok(())
}

#[allow(dead_code)]
fn write_external_cpp_path_info_materialized_coord_buckets_bucketed<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    coord_dir: &Path,
    threads: usize,
    expansion: &SerialExpansion<K>,
) -> Result<ExternalCppPathInfoMaterialized, SerialCollationError> {
    if coord_dir.exists() {
        fs::remove_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
            path: coord_dir.to_path_buf(),
            source,
        })?;
    }
    fs::create_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
        path: coord_dir.to_path_buf(),
        source,
    })?;

    let max_unitig_bucket_count = materialized_stitched_coord_bucket_count(threads);
    let max_unitig_bucket_mask = max_unitig_bucket_count - 1;
    let ranges_per_bucket = ranges_per_path_info_bucket(inputs.ranges.len(), 1);
    let path_info_bucket_count = inputs.ranges.len().div_ceil(ranges_per_bucket).max(1);
    let mut path_info_records = (0..path_info_bucket_count)
        .map(|_| Vec::<StitchedCoordRecord>::new())
        .collect::<Vec<_>>();
    let mut phantom_writers = MaterializedStitchedCoordShardWriters::new(
        coord_dir,
        threads.max(1),
        max_unitig_bucket_count,
    );

    for edge in &expansion.edges {
        let path_id = u64::try_from(edge.info.path_id.as_u128())
            .map_err(|_| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;
        let unitig_index = u32::try_from(edge.unitig_index)
            .map_err(|_| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;
        let record = StitchedCoordRecord {
            path_id,
            rank: edge.info.rank,
            unitig_index,
            reverse: edge.info.exit_side == Side::Front,
            is_cycle: edge.info.is_cycle,
        };

        if let Some(phantom) = edge.phantom_unitig {
            let label = phantom_unitig_label(phantom);
            let bucket_id = stitched_coord_bucket(path_id, max_unitig_bucket_mask);
            phantom_writers.write_materialized_record(bucket_id, &record, &label)?;
            continue;
        }

        let range_id =
            external_range_id_for_unitig(&inputs.ranges, edge.unitig_index).ok_or_else(|| {
                SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone())
            })?;
        let bucket_id = range_id / ranges_per_bucket;
        path_info_records
            .get_mut(bucket_id)
            .ok_or_else(|| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?
            .push(record);
    }

    let workers = threads.max(1).min(path_info_records.len().max(1));
    let mut manifest = if workers == 1 || path_info_records.len() < 2 {
        let mut writers =
            MaterializedStitchedCoordShardWriters::new(coord_dir, 0, max_unitig_bucket_count);
        for (bucket_id, records) in path_info_records.iter().enumerate() {
            materialize_external_stitched_coord_range_group_from_records::<K>(
                inputs,
                max_unitig_bucket_mask,
                ranges_per_bucket,
                bucket_id,
                records,
                &mut writers,
            )?;
        }
        writers.finish()?
    } else {
        let next_bucket = AtomicUsize::new(0);
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for worker_id in 0..workers {
                let next_bucket = &next_bucket;
                let path_info_records = &path_info_records;
                handles.push(scope.spawn(move || {
                    let mut writers = MaterializedStitchedCoordShardWriters::new(
                        coord_dir,
                        worker_id,
                        max_unitig_bucket_count,
                    );
                    loop {
                        let bucket_id = next_bucket.fetch_add(1, Ordering::Relaxed);
                        let Some(records) = path_info_records.get(bucket_id) else {
                            break;
                        };
                        materialize_external_stitched_coord_range_group_from_records::<K>(
                            inputs,
                            max_unitig_bucket_mask,
                            ranges_per_bucket,
                            bucket_id,
                            records,
                            &mut writers,
                        )?;
                    }
                    writers.finish()
                }));
            }

            let mut manifest = Vec::new();
            for handle in handles {
                manifest.extend(
                    handle
                        .join()
                        .map_err(|_| SerialCollationError::WorkerPanic)??,
                );
            }
            Ok::<_, SerialCollationError>(manifest)
        })?
    };
    manifest.extend(phantom_writers.finish()?);
    manifest.sort_by(|left, right| {
        left.bucket_id
            .cmp(&right.bucket_id)
            .then_with(|| left.coord_path.cmp(&right.coord_path))
    });
    Ok(ExternalCppPathInfoMaterialized {
        manifest,
        retained: Vec::new(),
        direct_local_unitigs: 0,
        direct_local_unitigs_complete: false,
    })
}

#[allow(dead_code)]
fn write_external_cpp_path_info_materialized_coord_buckets_legacy<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    coord_dir: &Path,
    threads: usize,
    expansion: &SerialExpansion<K>,
) -> Result<Vec<MaterializedStitchedCoordBucketEntry>, SerialCollationError> {
    if coord_dir.exists() {
        fs::remove_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
            path: coord_dir.to_path_buf(),
            source,
        })?;
    }
    fs::create_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
        path: coord_dir.to_path_buf(),
        source,
    })?;

    let bucket_count = materialized_stitched_coord_bucket_count(threads);
    let bucket_mask = bucket_count - 1;
    let mut unitigs = Vec::with_capacity(inputs.unitig_count());
    let reader = ExternalDiscontinuityReader::open(inputs)?;
    for unitig in reader.iter()? {
        unitigs.push(unitig?);
    }
    let mut label = Vec::new();
    let mut writers = MaterializedStitchedCoordShardWriters::new(coord_dir, 0, bucket_count);

    for edge in &expansion.edges {
        let path_id = u64::try_from(edge.info.path_id.as_u128())
            .map_err(|_| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;
        let unitig_index = u32::try_from(edge.unitig_index)
            .map_err(|_| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;
        let record = StitchedCoordRecord {
            path_id,
            rank: edge.info.rank,
            unitig_index,
            reverse: edge.info.exit_side == Side::Front,
            is_cycle: edge.info.is_cycle,
        };
        let bucket_id = stitched_coord_bucket(path_id, bucket_mask);
        if let Some(phantom) = edge.phantom_unitig {
            label = phantom_unitig_label(phantom);
        } else {
            let unitig = unitigs.get(edge.unitig_index).ok_or_else(|| {
                SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone())
            })?;
            reader.read_label(unitig, &mut label)?;
        }
        writers.write_materialized_record(bucket_id, &record, &label)?;
    }

    writers.finish()
}

#[allow(dead_code)]
fn write_external_cpp_path_info_materialized_coord_buckets_streaming<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    coord_dir: &Path,
    threads: usize,
    expansion: &SerialExpansion<K>,
    final_buckets: &mut FinalUnitigBucketWriters,
) -> Result<ExternalCppPathInfoMaterialized, SerialCollationError> {
    if coord_dir.exists() {
        fs::remove_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
            path: coord_dir.to_path_buf(),
            source,
        })?;
    }
    fs::create_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
        path: coord_dir.to_path_buf(),
        source,
    })?;

    let bucket_count = materialized_stitched_coord_bucket_count(threads);
    let bucket_mask = bucket_count - 1;
    let mut path_info_by_unitig = vec![None::<StitchedCoordRecord>; inputs.unitig_count()];
    let mut writers = MaterializedStitchedCoordShardWriters::new(coord_dir, 0, bucket_count);
    let mut duplicate_path_info = false;

    for edge in &expansion.edges {
        let path_id = u64::try_from(edge.info.path_id.as_u128())
            .map_err(|_| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;
        let unitig_index = u32::try_from(edge.unitig_index)
            .map_err(|_| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;
        let record = StitchedCoordRecord {
            path_id,
            rank: edge.info.rank,
            unitig_index,
            reverse: edge.info.exit_side == Side::Front,
            is_cycle: edge.info.is_cycle,
        };

        if let Some(phantom) = edge.phantom_unitig {
            let label = phantom_unitig_label(phantom);
            let bucket_id = stitched_coord_bucket(path_id, bucket_mask);
            writers.write_materialized_record(bucket_id, &record, &label)?;
            continue;
        }

        let Some(slot) = path_info_by_unitig.get_mut(edge.unitig_index) else {
            return Err(SerialCollationError::MalformedCoordBucket(
                inputs.unitig_path.clone(),
            ));
        };
        if slot.replace(record).is_some() {
            duplicate_path_info = true;
        }
    }

    if duplicate_path_info {
        drop(writers);
        let manifest = write_external_cpp_path_info_materialized_coord_buckets_legacy::<K>(
            inputs, coord_dir, threads, expansion,
        )?;
        return Ok(ExternalCppPathInfoMaterialized {
            manifest,
            retained: Vec::new(),
            direct_local_unitigs: 0,
            direct_local_unitigs_complete: false,
        });
    }

    let unitig_file =
        File::open(&inputs.unitig_path).map_err(|source| SerialCollationError::Io {
            path: inputs.unitig_path.clone(),
            source,
        })?;
    let mut unitigs: ExternalDiscontinuityIter<K> = ExternalDiscontinuityIter {
        input: BufReader::with_capacity(1024 * 1024, unitig_file),
        path: inputs.unitig_path.clone(),
        remaining: inputs.unitigs,
        compact_unitigs: inputs.compact_unitigs,
    };
    let label_file = File::open(&inputs.label_path).map_err(|source| SerialCollationError::Io {
        path: inputs.label_path.clone(),
        source,
    })?;
    let mut labels = BufReader::with_capacity(1024 * 1024, label_file);
    let mut label_pos = 0u64;
    let mut label = Vec::new();

    let mut direct_local_unitigs = 0u64;
    for (unitig_index, record) in path_info_by_unitig.into_iter().enumerate() {
        let unitig = unitigs.next_unitig()?;
        if label_pos != unitig.label_start {
            labels
                .seek(SeekFrom::Start(unitig.label_start))
                .map_err(|source| SerialCollationError::Io {
                    path: inputs.label_path.clone(),
                    source,
                })?;
            label_pos = unitig.label_start;
        }
        label.resize(unitig.label_len as usize, 0);
        labels
            .read_exact(&mut label)
            .map_err(|source| SerialCollationError::Io {
                path: inputs.label_path.clone(),
                source,
            })?;
        label_pos += u64::from(unitig.label_len);

        let Some(record) = record else {
            if unitig.left_exit().is_none() && unitig.right_exit().is_none() {
                let label = canonical_label(label.clone());
                final_buckets.write_label(&label)?;
                direct_local_unitigs += 1;
            }
            continue;
        };
        if usize::try_from(record.unitig_index).ok() != Some(unitig_index) {
            return Err(SerialCollationError::MalformedCoordBucket(
                inputs.unitig_path.clone(),
            ));
        }

        let bucket_id = stitched_coord_bucket(record.path_id, bucket_mask);
        writers.write_materialized_record(bucket_id, &record, &label)?;
    }

    Ok(ExternalCppPathInfoMaterialized {
        manifest: writers.finish()?,
        retained: Vec::new(),
        direct_local_unitigs,
        direct_local_unitigs_complete: true,
    })
}

fn phantom_unitig_label<const K: usize>(endpoint: DiscontinuityEndpoint<K>) -> Vec<u8> {
    let label = endpoint.vertex.to_ascii_string().into_bytes();
    if endpoint.side == Side::Front {
        label
    } else {
        reverse_complement_label(&label)
    }
}

fn emit_external_ordered_stitched_labels_to_final_buckets<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    coord_dir: &Path,
    label_refs: &[ExternalLabelRef],
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    starts: &[usize],
    final_buckets: &mut FinalUnitigBucketWriters,
) -> Result<u64, SerialCollationError> {
    const NO_RECORD_INDEX: u32 = u32::MAX;
    let map_started = Instant::now();
    let mut record_index_by_unitig = vec![NO_RECORD_INDEX; inputs.unitig_count()];
    let mut path_offsets = Vec::with_capacity(starts.len() + 1);
    let mut records = Vec::<StitchedCoordRecord>::new();
    let mut path_records = Vec::new();

    for (path_id, &start) in starts.iter().enumerate() {
        path_offsets.push(records.len());
        path_records.clear();
        walk_simple_stitched_component_coords(
            half_ends,
            join_neighbor,
            start,
            path_id as u64,
            &mut path_records,
        );
        for record in &path_records {
            let unitig_index = record.unitig_index as usize;
            let slot = record_index_by_unitig
                .get_mut(unitig_index)
                .ok_or_else(|| {
                    SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone())
                })?;
            if *slot != NO_RECORD_INDEX {
                return Err(SerialCollationError::MalformedCoordBucket(
                    inputs.unitig_path.clone(),
                ));
            }
            *slot = u32::try_from(records.len()).map_err(|_| {
                SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone())
            })?;
            records.push(*record);
        }
    }
    path_offsets.push(records.len());
    let map_elapsed = map_started.elapsed();

    let offset_started = Instant::now();
    let mut record_label_offsets = Vec::with_capacity(records.len());
    let mut ordered_label_bytes = 0u64;
    for record in &records {
        let label_ref = label_refs
            .get(record.unitig_index as usize)
            .ok_or_else(|| SerialCollationError::MalformedCoordBucket(inputs.label_path.clone()))?;
        record_label_offsets.push(ordered_label_bytes);
        ordered_label_bytes += u64::from(label_ref.label_len);
    }
    let offset_elapsed = offset_started.elapsed();

    let spill_started = Instant::now();
    fs::create_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
        path: coord_dir.to_path_buf(),
        source,
    })?;
    let ordered_label_path = coord_dir.join("ordered-stitch-labels.bin");
    let ordered_label_file =
        File::create(&ordered_label_path).map_err(|source| SerialCollationError::Io {
            path: ordered_label_path.clone(),
            source,
        })?;
    ordered_label_file
        .set_len(ordered_label_bytes)
        .map_err(|source| SerialCollationError::Io {
            path: ordered_label_path.clone(),
            source,
        })?;

    let source_label_file =
        File::open(&inputs.label_path).map_err(|source| SerialCollationError::Io {
            path: inputs.label_path.clone(),
            source,
        })?;
    let mut source_labels = BufReader::with_capacity(1024 * 1024, source_label_file);
    let mut label = Vec::new();
    for (unitig_index, label_ref) in label_refs.iter().copied().enumerate() {
        label.resize(label_ref.label_len as usize, 0);
        source_labels
            .read_exact(&mut label)
            .map_err(|source| SerialCollationError::Io {
                path: inputs.label_path.clone(),
                source,
            })?;
        let record_index = record_index_by_unitig[unitig_index];
        if record_index == NO_RECORD_INDEX {
            continue;
        }
        ordered_label_file
            .write_all_at(&label, record_label_offsets[record_index as usize])
            .map_err(|source| SerialCollationError::Io {
                path: ordered_label_path.clone(),
                source,
            })?;
    }
    drop(ordered_label_file);
    let spill_elapsed = spill_started.elapsed();

    let assemble_started = Instant::now();
    let ordered_label_file =
        File::open(&ordered_label_path).map_err(|source| SerialCollationError::Io {
            path: ordered_label_path.clone(),
            source,
        })?;
    let mut ordered_labels = BufReader::with_capacity(1024 * 1024, ordered_label_file);
    let mut emitted = 0u64;
    let mut unitig_label = Vec::new();
    let mut stitched_label = Vec::new();
    for bounds in path_offsets.windows(2) {
        let start = bounds[0];
        let end = bounds[1];
        if start == end {
            continue;
        }
        stitched_label.clear();
        let is_cycle = records[start].is_cycle;
        for record in &records[start..end] {
            let label_ref = label_refs
                .get(record.unitig_index as usize)
                .ok_or_else(|| {
                    SerialCollationError::MalformedCoordBucket(inputs.label_path.clone())
                })?;
            unitig_label.resize(label_ref.label_len as usize, 0);
            ordered_labels
                .read_exact(&mut unitig_label)
                .map_err(|source| SerialCollationError::Io {
                    path: ordered_label_path.clone(),
                    source,
                })?;
            let mut reverse = record.reverse;
            if !stitched_label.is_empty()
                && !labels_overlap_oriented_fast::<K>(&stitched_label, &unitig_label, reverse)
            {
                let alternate = oriented_label(&unitig_label, !reverse);
                if labels_overlap::<K>(&stitched_label, &alternate) {
                    reverse = !reverse;
                }
            }
            append_or_init_oriented_fast::<K>(&mut stitched_label, &unitig_label, reverse);
        }

        if stitched_label.len() >= K {
            let label = if is_cycle {
                normalize_stitched_cycle::<K>(&stitched_label)
            } else {
                canonical_label(stitched_label.clone())
            };
            final_buckets.write_label(&label)?;
            emitted += 1;
        }
    }
    let assemble_elapsed = assemble_started.elapsed();

    eprintln!(
        "cuttlefish3-rs: ordered external stitch labels: map {:.3}s, offsets {:.3}s, label spill {:.3}s, assemble {:.3}s, records {}, label bytes {}",
        map_elapsed.as_secs_f64(),
        offset_elapsed.as_secs_f64(),
        spill_elapsed.as_secs_f64(),
        assemble_elapsed.as_secs_f64(),
        records.len(),
        ordered_label_bytes
    );
    Ok(emitted)
}

fn emit_external_stitched_labels_to_final_buckets<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    threads: usize,
    label_refs: &[ExternalLabelRef],
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    starts: &[usize],
    final_buckets: &mut FinalUnitigBucketWriters,
) -> Result<u64, SerialCollationError> {
    let workers = threads.max(1).min(starts.len().max(1));
    if workers == 1 || starts.len() < 1024 {
        let label_file =
            File::open(&inputs.label_path).map_err(|source| SerialCollationError::Io {
                path: inputs.label_path.clone(),
                source,
            })?;
        let mut emitted = 0u64;
        for &start in starts {
            if let Some(label) = walk_external_simple_stitched_component_label::<K>(
                &label_file,
                &inputs.label_path,
                label_refs,
                half_ends,
                join_neighbor,
                start,
            )? {
                final_buckets.write_label(&label)?;
                emitted += 1;
            }
        }
        return Ok(emitted);
    }

    const STITCH_LABEL_BATCH: usize = 256;
    let chunk_size = starts.len().div_ceil(workers);
    let (tx, rx) = mpsc::sync_channel::<Result<Vec<Vec<u8>>, SerialCollationError>>(workers * 2);
    let mut emitted = 0u64;
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for chunk in starts.chunks(chunk_size) {
            let tx = tx.clone();
            handles.push(scope.spawn(move || {
                let label_file =
                    File::open(&inputs.label_path).map_err(|source| SerialCollationError::Io {
                        path: inputs.label_path.clone(),
                        source,
                    })?;
                let mut batch = Vec::with_capacity(STITCH_LABEL_BATCH);
                for &start in chunk {
                    if let Some(label) = walk_external_simple_stitched_component_label::<K>(
                        &label_file,
                        &inputs.label_path,
                        label_refs,
                        half_ends,
                        join_neighbor,
                        start,
                    )? {
                        batch.push(label);
                        if batch.len() == STITCH_LABEL_BATCH {
                            let full = std::mem::take(&mut batch);
                            if tx.send(Ok(full)).is_err() {
                                return Ok(());
                            }
                        }
                    }
                }
                if !batch.is_empty() {
                    let _ = tx.send(Ok(batch));
                }
                Ok::<_, SerialCollationError>(())
            }));
        }
        drop(tx);

        let mut first_error = None;
        for result in rx {
            match result {
                Ok(labels) if first_error.is_none() => {
                    for label in labels {
                        if let Err(err) = final_buckets.write_label(&label) {
                            first_error = Some(err);
                            break;
                        }
                        emitted += 1;
                    }
                }
                Ok(_) => {}
                Err(err) if first_error.is_none() => first_error = Some(err),
                Err(_) => {}
            }
        }

        for handle in handles {
            if let Err(err) = handle
                .join()
                .map_err(|_| SerialCollationError::WorkerPanic)?
            {
                if first_error.is_none() {
                    first_error = Some(err);
                }
            }
        }

        if let Some(err) = first_error {
            Err(err)
        } else {
            Ok(())
        }
    })?;

    Ok(emitted)
}

fn walk_external_simple_stitched_component_label<const K: usize>(
    label_file: &File,
    label_path: &Path,
    label_refs: &[ExternalLabelRef],
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    start: usize,
) -> Result<Option<Vec<u8>>, SerialCollationError> {
    let mut current = start;
    let mut label = Vec::new();
    let mut unitig_label = Vec::new();
    let mut is_cycle = false;

    loop {
        let unitig_index = half_ends[current].unitig_index();
        let label_ref = label_refs
            .get(unitig_index)
            .ok_or_else(|| SerialCollationError::MalformedCoordBucket(label_path.to_path_buf()))?;
        read_external_label_ref(label_file, label_path, *label_ref, &mut unitig_label)?;
        let reverse = reverse_for_stitch_node(current);
        let mut append_reverse = reverse;
        if !label.is_empty() && !labels_overlap_oriented_fast::<K>(&label, &unitig_label, reverse) {
            let alternate = oriented_label(&unitig_label, !reverse);
            if labels_overlap::<K>(&label, &alternate) {
                append_reverse = !reverse;
            }
        }
        append_or_init_oriented_fast::<K>(&mut label, &unitig_label, append_reverse);

        let other = current ^ 1;
        let next = join_neighbor[other];
        if next == STITCH_NO_NODE {
            break;
        }
        if stitch_node_index(next) == start {
            is_cycle = true;
            break;
        }
        current = stitch_node_index(next);
    }

    if label.len() < K {
        Ok(None)
    } else if is_cycle {
        Ok(Some(normalize_stitched_cycle::<K>(&label)))
    } else {
        Ok(Some(canonical_label(label)))
    }
}

fn read_external_label_ref(
    label_file: &File,
    label_path: &Path,
    label_ref: ExternalLabelRef,
    scratch: &mut Vec<u8>,
) -> Result<(), SerialCollationError> {
    scratch.resize(label_ref.label_len as usize, 0);
    label_file
        .read_exact_at(scratch, label_ref.label_start)
        .map_err(|source| SerialCollationError::Io {
            path: label_path.to_path_buf(),
            source,
        })
}

fn stitch_discontinuity_paths_impl<const K: usize>(
    inputs: &DiscontinuityInputs<K>,
    skip_unitigs: &[bool],
    threads: usize,
    coord_dir: Option<&Path>,
) -> Result<Vec<Vec<u8>>, SerialCollationError> {
    let debug_stitch = std::env::var_os("CF3_RS_DEBUG_STITCH").is_some();
    let use_adjacency_stitch = std::env::var_os("CF3_RS_ADJACENCY_STITCH").is_some();
    let endpoint_dir = coord_dir.map(|dir| dir.join("endpoints"));
    let mut endpoint_writers = if !use_adjacency_stitch {
        endpoint_dir
            .as_deref()
            .map(|dir| StitchEndpointBucketWriters::create(dir, threads))
            .transpose()?
    } else {
        None
    };
    let build_started = Instant::now();
    let mut half_ends = Vec::<HalfEnd>::new();
    let mut segment_edges = if use_adjacency_stitch {
        Some(Vec::<(usize, usize, usize)>::new())
    } else {
        None
    };
    let mut endpoint_records = if endpoint_writers.is_none() {
        Some(Vec::<StitchEndpointRecord<K>>::new())
    } else {
        None
    };

    for (unitig_index, unitig) in inputs.unitigs.iter().enumerate() {
        if skip_unitigs.get(unitig_index).copied().unwrap_or(false) {
            continue;
        }
        if unitig.left_exit().is_none() && unitig.right_exit().is_none() {
            continue;
        }

        let (left_endpoint, right_endpoint) = endpoints_by_label_end(unitig);
        let left_node = half_ends.len();
        half_ends.push(HalfEnd {
            unitig_index: u32::try_from(unitig_index).expect("unitig index exceeds u32"),
        });
        if let Some(endpoint) = left_endpoint {
            let record = StitchEndpointRecord {
                vertex: endpoint.vertex,
                side: endpoint.side,
                node: stitch_node(left_node),
            };
            if let Some(writers) = &mut endpoint_writers {
                writers.write_record(&record)?;
            } else {
                endpoint_records
                    .as_mut()
                    .expect("in-memory endpoint records are enabled")
                    .push(record);
            }
        }
        let right_node = half_ends.len();
        half_ends.push(HalfEnd {
            unitig_index: u32::try_from(unitig_index).expect("unitig index exceeds u32"),
        });
        if let Some(endpoint) = right_endpoint {
            let record = StitchEndpointRecord {
                vertex: endpoint.vertex,
                side: endpoint.side,
                node: stitch_node(right_node),
            };
            if let Some(writers) = &mut endpoint_writers {
                writers.write_record(&record)?;
            } else {
                endpoint_records
                    .as_mut()
                    .expect("in-memory endpoint records are enabled")
                    .push(record);
            }
        }
        if let Some(segment_edges) = &mut segment_edges {
            segment_edges.push((left_node, right_node, unitig_index));
        }
    }
    let half_end_elapsed = build_started.elapsed();

    if debug_stitch {
        eprintln!("stitch inputs: {} unitigs", inputs.unitigs.len());
        for (index, unitig) in inputs.unitigs.iter().enumerate() {
            let (left_endpoint, right_endpoint) = endpoints_by_label_end(unitig);
            eprintln!(
                "unitig {index}: label={} left_exit={:?} right_exit={:?} label_left={:?} label_right={:?}",
                String::from_utf8_lossy(unitig.label(inputs)),
                unitig.left_exit().map(debug_endpoint),
                unitig.right_exit().map(debug_endpoint),
                left_endpoint.map(debug_endpoint),
                right_endpoint.map(debug_endpoint),
            );
        }
    }

    if use_adjacency_stitch {
        return stitch_discontinuity_paths_with_adjacency(
            inputs,
            threads,
            coord_dir,
            debug_stitch,
            half_ends,
            segment_edges.expect("adjacency stitching requires segment edges"),
            endpoint_records.expect("adjacency stitching uses in-memory endpoint records"),
            half_end_elapsed,
        );
    }

    let mut join_neighbor = vec![STITCH_NO_NODE; half_ends.len()];

    let endpoint_sort_started = Instant::now();
    let endpoint_manifest = endpoint_writers
        .map(StitchEndpointBucketWriters::finish)
        .transpose()?;
    let mut endpoint_records = endpoint_records.unwrap_or_default();
    if endpoint_manifest.is_none() {
        endpoint_records.sort_unstable_by_key(|record| record.vertex.as_u128());
    }
    let endpoint_sort_elapsed = endpoint_sort_started.elapsed();
    let endpoint_join_started = Instant::now();
    if let Some(manifest) = &endpoint_manifest {
        join_neighbors_from_endpoint_buckets::<K>(manifest, &mut join_neighbor, threads)?;
    } else {
        join_neighbors_from_sorted_endpoints(&endpoint_records, &mut join_neighbor);
    }
    let endpoint_join_elapsed = endpoint_join_started.elapsed();
    drop(endpoint_records);
    trim_process_allocations();

    if debug_stitch {
        for (node, half_end) in half_ends.iter().enumerate() {
            let (left_endpoint, right_endpoint) =
                endpoints_by_label_end(&inputs.unitigs[half_end.unitig_index()]);
            let label_end = label_end_for_node(node);
            let endpoint = match label_end {
                LabelEnd::Left => left_endpoint,
                LabelEnd::Right => right_endpoint,
            };
            eprintln!(
                "node {node}: unitig={} end={:?} endpoint={:?} join={}",
                half_end.unitig_index(),
                label_end,
                endpoint.map(debug_endpoint),
                join_neighbor[node],
            );
        }
    }

    let component_started = Instant::now();
    let component_starts = stitch_component_starts_from_neighbors(&half_ends, &join_neighbor);
    let component_elapsed = component_started.elapsed();
    let walk_started = Instant::now();
    let mut unitigs = if let Some(coord_dir) = coord_dir {
        let manifest = write_stitched_coord_buckets(
            coord_dir,
            threads,
            &half_ends,
            &join_neighbor,
            &component_starts,
        )?;
        reduce_stitched_coord_bucket_files::<K>(inputs, &manifest, threads)?
    } else {
        stitch_simple_components_labels_with_threads(
            inputs,
            &half_ends,
            &join_neighbor,
            &component_starts,
            threads,
        )
    };
    let walk_elapsed = walk_started.elapsed();

    let sort_started = Instant::now();
    unitigs = bucketed_maximal_unitig_reduce(unitigs, threads);
    let sort_elapsed = sort_started.elapsed();
    eprintln!(
        "cuttlefish3-rs: stitch detail: half-ends {:.3}s, endpoint sort {:.3}s, endpoint join {:.3}s, components {:.3}s, {} walk {:.3}s, bucket reduce {:.3}s",
        half_end_elapsed.as_secs_f64(),
        endpoint_sort_elapsed.as_secs_f64(),
        endpoint_join_elapsed.as_secs_f64(),
        component_elapsed.as_secs_f64(),
        "component",
        walk_elapsed.as_secs_f64(),
        sort_elapsed.as_secs_f64()
    );
    if debug_stitch {
        for unitig in &unitigs {
            eprintln!("stitched: {}", String::from_utf8_lossy(unitig));
        }
    }
    Ok(unitigs)
}

fn stitch_discontinuity_paths_with_adjacency<const K: usize>(
    inputs: &DiscontinuityInputs<K>,
    threads: usize,
    coord_dir: Option<&Path>,
    debug_stitch: bool,
    half_ends: Vec<HalfEnd>,
    segment_edges: Vec<(usize, usize, usize)>,
    mut endpoint_records: Vec<StitchEndpointRecord<K>>,
    half_end_elapsed: Duration,
) -> Result<Vec<Vec<u8>>, SerialCollationError> {
    let mut adjacency = vec![StitchAdjacencyList::default(); half_ends.len()];
    for &(left_node, right_node, unitig_index) in &segment_edges {
        push_stitch_edge(&mut adjacency, left_node, right_node, Some(unitig_index));
    }

    let endpoint_sort_started = Instant::now();
    endpoint_records.sort_unstable_by_key(|record| record.vertex.as_u128());
    let endpoint_sort_elapsed = endpoint_sort_started.elapsed();
    let endpoint_join_started = Instant::now();
    let mut start = 0;
    while start < endpoint_records.len() {
        let vertex = endpoint_records[start].vertex;
        let mut end = start + 1;
        while end < endpoint_records.len() && endpoint_records[end].vertex == vertex {
            end += 1;
        }

        let mut ends = StitchVertexEnds::default();
        for record in &endpoint_records[start..end] {
            ends.push(record.side, record.node);
        }
        if ends.front_len == 1 && ends.back_len == 1 {
            push_stitch_edge(
                &mut adjacency,
                stitch_node_index(ends.fronts[0]),
                stitch_node_index(ends.backs[0]),
                None,
            );
        } else if ends.total_len == 2 {
            push_stitch_edge(
                &mut adjacency,
                stitch_node_index(ends.nodes[0]),
                stitch_node_index(ends.nodes[1]),
                None,
            );
        }
        start = end;
    }
    let endpoint_join_elapsed = endpoint_join_started.elapsed();
    drop(endpoint_records);
    trim_process_allocations();

    if debug_stitch {
        for (node, half_end) in half_ends.iter().enumerate() {
            let (left_endpoint, right_endpoint) =
                endpoints_by_label_end(&inputs.unitigs[half_end.unitig_index()]);
            let label_end = label_end_for_node(node);
            let endpoint = match label_end {
                LabelEnd::Left => left_endpoint,
                LabelEnd::Right => right_endpoint,
            };
            let edges = adjacency[node]
                .iter()
                .map(|edge| format!("{}:{:?}", edge.to, edge.unitig_index))
                .collect::<Vec<_>>()
                .join(",");
            eprintln!(
                "node {node}: unitig={} end={:?} endpoint={:?} edges=[{}]",
                half_end.unitig_index(),
                label_end,
                endpoint.map(debug_endpoint),
                edges,
            );
        }
    }

    let mut unitigs = Vec::new();
    let component_started = Instant::now();
    let component_starts = stitch_component_starts_from_adjacency(&adjacency);
    let component_elapsed = component_started.elapsed();
    let walk_started = Instant::now();
    let used_simple_components = component_starts.is_some();
    if let Some(starts) = component_starts {
        unitigs = if let Some(coord_dir) = coord_dir {
            let manifest = write_stitched_coord_buckets_from_adjacency(
                coord_dir, threads, &adjacency, &starts,
            )?;
            reduce_stitched_coord_bucket_files::<K>(inputs, &manifest, threads)?
        } else {
            let records = map_simple_stitched_components_with_threads_from_adjacency(
                &adjacency, &starts, threads,
            );
            reduce_stitched_coord_buckets::<K>(inputs, records, threads)
        };
    } else {
        let mut visited_segments = vec![false; inputs.unitigs.len()];

        for start in (0..half_ends.len()).filter(|&node| adjacency[node].len() == 1) {
            let Some(path) = walk_stitched_path(inputs, &adjacency, &mut visited_segments, start)
            else {
                continue;
            };
            push_stitched_path_label::<K>(&mut unitigs, path);
        }

        for &(left_node, _, unitig_index) in &segment_edges {
            if visited_segments[unitig_index] {
                continue;
            }
            let Some(path) =
                walk_stitched_path(inputs, &adjacency, &mut visited_segments, left_node)
            else {
                continue;
            };
            push_stitched_path_label::<K>(&mut unitigs, path);
        }
    }
    let walk_elapsed = walk_started.elapsed();

    let sort_started = Instant::now();
    unitigs = bucketed_maximal_unitig_reduce(unitigs, threads);
    let sort_elapsed = sort_started.elapsed();
    eprintln!(
        "cuttlefish3-rs: stitch detail: half-ends {:.3}s, endpoint sort {:.3}s, endpoint join {:.3}s, components {:.3}s, {} walk {:.3}s, bucket reduce {:.3}s",
        half_end_elapsed.as_secs_f64(),
        endpoint_sort_elapsed.as_secs_f64(),
        endpoint_join_elapsed.as_secs_f64(),
        component_elapsed.as_secs_f64(),
        if used_simple_components {
            "component"
        } else {
            "fallback"
        },
        walk_elapsed.as_secs_f64(),
        sort_elapsed.as_secs_f64()
    );
    if debug_stitch {
        for unitig in &unitigs {
            eprintln!("stitched: {}", String::from_utf8_lossy(unitig));
        }
    }
    Ok(unitigs)
}

fn join_neighbors_from_endpoint_buckets<const K: usize>(
    manifest: &[StitchEndpointBucketEntry],
    join_neighbor: &mut [u32],
    threads: usize,
) -> Result<(), SerialCollationError> {
    let workers = threads.max(1).min(manifest.len().max(1));
    if workers == 1 || manifest.len() < 2 {
        for entry in manifest {
            let mut records = read_stitch_endpoint_bucket::<K>(entry)?;
            records.sort_unstable_by_key(|record| record.vertex.as_u128());
            join_neighbors_from_sorted_endpoint_assignments(&records, join_neighbor);
        }
        return Ok(());
    }

    let next_entry = AtomicUsize::new(0);
    let join_writer = JoinNeighborWriter::new(join_neighbor);
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for _ in 0..workers {
            let next_entry = &next_entry;
            let join_writer = &join_writer;
            handles.push(scope.spawn(move || {
                loop {
                    let entry_idx = next_entry.fetch_add(1, Ordering::Relaxed);
                    let Some(entry) = manifest.get(entry_idx) else {
                        break;
                    };
                    let mut records = read_stitch_endpoint_bucket::<K>(entry)?;
                    records.sort_unstable_by_key(|record| record.vertex.as_u128());
                    join_writer.write_sorted_endpoint_assignments(&records);
                }
                Ok::<_, SerialCollationError>(())
            }));
        }

        for handle in handles {
            handle
                .join()
                .map_err(|_| SerialCollationError::WorkerPanic)??;
        }
        Ok::<_, SerialCollationError>(())
    })
}

struct JoinNeighborWriter {
    ptr: *mut u32,
    len: usize,
}

// Endpoint buckets are partitioned by vertex, so each stitch node's join slot is
// assigned by at most one worker. The wrapper lets workers fill those disjoint
// slots without materializing a second global assignment vector.
unsafe impl Sync for JoinNeighborWriter {}

impl JoinNeighborWriter {
    fn new(join_neighbor: &mut [u32]) -> Self {
        Self {
            ptr: join_neighbor.as_mut_ptr(),
            len: join_neighbor.len(),
        }
    }

    fn write_sorted_endpoint_assignments<const K: usize>(
        &self,
        endpoint_records: &[StitchEndpointRecord<K>],
    ) {
        let mut start = 0;
        while start < endpoint_records.len() {
            let vertex = endpoint_records[start].vertex;
            let mut end = start + 1;
            while end < endpoint_records.len() && endpoint_records[end].vertex == vertex {
                end += 1;
            }

            let mut ends = StitchVertexEnds::default();
            for record in &endpoint_records[start..end] {
                ends.push(record.side, record.node);
            }
            if ends.front_len == 1 && ends.back_len == 1 {
                self.write(ends.fronts[0], ends.backs[0]);
                self.write(ends.backs[0], ends.fronts[0]);
            } else if ends.total_len == 2 {
                self.write(ends.nodes[0], ends.nodes[1]);
                self.write(ends.nodes[1], ends.nodes[0]);
            }
            start = end;
        }
    }

    fn write(&self, node: u32, neighbor: u32) {
        let index = stitch_node_index(node);
        debug_assert!(index < self.len);
        if index < self.len {
            unsafe {
                *self.ptr.add(index) = neighbor;
            }
        }
    }
}

fn read_stitch_endpoint_bucket<const K: usize>(
    entry: &StitchEndpointBucketEntry,
) -> Result<Vec<StitchEndpointRecord<K>>, SerialCollationError> {
    let file = File::open(&entry.path).map_err(|source| SerialCollationError::Io {
        path: entry.path.clone(),
        source,
    })?;
    let mut input = BufReader::with_capacity(1024 * 1024, file);
    let mut records = Vec::with_capacity(entry.records as usize);
    let mut vertex_bytes = [0u8; 16];
    let mut side_byte = [0u8; 1];
    let mut node_bytes = [0u8; 4];
    for _ in 0..entry.records {
        input
            .read_exact(&mut vertex_bytes)
            .and_then(|_| input.read_exact(&mut side_byte))
            .and_then(|_| input.read_exact(&mut node_bytes))
            .map_err(|source| SerialCollationError::Io {
                path: entry.path.clone(),
                source,
            })?;
        let side = match side_byte[0] {
            0 => Side::Front,
            1 => Side::Back,
            _ => {
                return Err(SerialCollationError::MalformedCoordBucket(
                    entry.path.clone(),
                ));
            }
        };
        records.push(StitchEndpointRecord {
            vertex: Kmer::from_bits(u128::from_le_bytes(vertex_bytes)),
            side,
            node: u32::from_le_bytes(node_bytes),
        });
    }

    let mut trailing = [0u8; 1];
    if input
        .read(&mut trailing)
        .map_err(|source| SerialCollationError::Io {
            path: entry.path.clone(),
            source,
        })?
        != 0
    {
        return Err(SerialCollationError::MalformedCoordBucket(
            entry.path.clone(),
        ));
    }

    Ok(records)
}

fn join_neighbors_from_sorted_endpoints<const K: usize>(
    endpoint_records: &[StitchEndpointRecord<K>],
    join_neighbor: &mut [u32],
) {
    join_neighbors_from_sorted_endpoint_assignments(endpoint_records, join_neighbor);
}

fn join_neighbors_from_sorted_endpoint_assignments<const K: usize>(
    endpoint_records: &[StitchEndpointRecord<K>],
    join_neighbor: &mut [u32],
) {
    let mut start = 0;
    while start < endpoint_records.len() {
        let vertex = endpoint_records[start].vertex;
        let mut end = start + 1;
        while end < endpoint_records.len() && endpoint_records[end].vertex == vertex {
            end += 1;
        }

        let mut ends = StitchVertexEnds::default();
        for record in &endpoint_records[start..end] {
            ends.push(record.side, record.node);
        }
        if ends.front_len == 1 && ends.back_len == 1 {
            join_neighbor[stitch_node_index(ends.fronts[0])] = ends.backs[0];
            join_neighbor[stitch_node_index(ends.backs[0])] = ends.fronts[0];
        } else if ends.total_len == 2 {
            join_neighbor[stitch_node_index(ends.nodes[0])] = ends.nodes[1];
            join_neighbor[stitch_node_index(ends.nodes[1])] = ends.nodes[0];
        }
        start = end;
    }
}

fn stitch_component_starts_from_adjacency(adjacency: &[StitchAdjacencyList]) -> Option<Vec<usize>> {
    let mut visited = vec![false; adjacency.len()];
    let mut starts = Vec::new();

    for node in 0..adjacency.len() {
        if adjacency[node].len() > 2 {
            return None;
        }
    }

    for node in 0..adjacency.len() {
        if visited[node] || adjacency[node].len() != 1 {
            continue;
        }
        starts.push(node);
        mark_linear_stitch_component(adjacency, node, &mut visited);
    }

    for node in 0..adjacency.len() {
        if visited[node] || adjacency[node].len() == 0 {
            continue;
        }
        starts.push(node);
        mark_linear_stitch_component(adjacency, node, &mut visited);
    }

    Some(starts)
}

fn mark_linear_stitch_component(
    adjacency: &[StitchAdjacencyList],
    start: usize,
    visited: &mut [bool],
) {
    let mut previous = usize::MAX;
    let mut current = start;

    loop {
        if visited[current] {
            break;
        }
        visited[current] = true;

        let next = adjacency[current]
            .iter()
            .map(|edge| edge.to)
            .find(|&next| next != previous && !visited[next]);
        let Some(next) = next else {
            break;
        };
        previous = current;
        current = next;
    }
}

fn stitch_component_starts_from_neighbors(
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
) -> Vec<usize> {
    debug_assert_eq!(half_ends.len(), join_neighbor.len());
    debug_assert_eq!(half_ends.len() % 2, 0);
    let mut visited = vec![0u8; half_ends.len().div_ceil(2)];
    let mut starts = Vec::new();

    for node in 0..join_neighbor.len() {
        if join_neighbor[node] != STITCH_NO_NODE {
            continue;
        }
        let unitig_pair = node >> 1;
        if visited[unitig_pair] != 0 {
            continue;
        }
        mark_neighbor_component(join_neighbor, node, &mut visited);
        starts.push(node);
    }

    for node in (0..join_neighbor.len()).step_by(2) {
        let unitig_pair = node >> 1;
        if visited[unitig_pair] != 0 {
            continue;
        }
        mark_neighbor_component(join_neighbor, node, &mut visited);
        starts.push(node);
    }

    starts
}

fn mark_neighbor_component(join_neighbor: &[u32], start: usize, visited: &mut [u8]) {
    let mut current = start;
    loop {
        let unitig_pair = current >> 1;
        if visited[unitig_pair] != 0 {
            break;
        }
        visited[unitig_pair] = 1;

        let other = current ^ 1;
        let next = join_neighbor[other];
        if next == STITCH_NO_NODE || stitch_node_index(next) == start {
            break;
        }
        current = stitch_node_index(next);
    }
}

fn push_stitched_path_label<const K: usize>(unitigs: &mut Vec<Vec<u8>>, path: StitchedPath) {
    let label = path.label;
    if label.len() >= K {
        if path.is_cycle {
            unitigs.push(normalize_stitched_cycle::<K>(&label));
        } else {
            unitigs.push(canonical_label(label));
        }
    }
}

fn stitch_simple_components_labels_with_threads<const K: usize>(
    inputs: &DiscontinuityInputs<K>,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    starts: &[usize],
    threads: usize,
) -> Vec<Vec<u8>> {
    let workers = threads.max(1).min(starts.len().max(1));
    if workers == 1 || starts.len() < 1024 {
        let mut unitigs = Vec::new();
        for &start in starts {
            if let Some(label) =
                walk_simple_stitched_component_label(inputs, half_ends, join_neighbor, start)
            {
                unitigs.push(label);
            }
        }
        return unitigs;
    }

    let chunk_size = starts.len().div_ceil(workers);
    let worker_unitigs = std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for chunk in starts.chunks(chunk_size) {
            handles.push(scope.spawn(move || {
                let mut unitigs = Vec::new();
                for &start in chunk {
                    if let Some(label) = walk_simple_stitched_component_label(
                        inputs,
                        half_ends,
                        join_neighbor,
                        start,
                    ) {
                        unitigs.push(label);
                    }
                }
                unitigs
            }));
        }

        let mut worker_unitigs = Vec::new();
        for handle in handles {
            worker_unitigs.push(handle.join().expect("stitch label worker panicked"));
        }
        worker_unitigs
    });

    let total = worker_unitigs.iter().map(Vec::len).sum();
    let mut unitigs = Vec::with_capacity(total);
    for mut worker in worker_unitigs {
        unitigs.append(&mut worker);
    }
    unitigs
}

fn walk_simple_stitched_component_label<const K: usize>(
    inputs: &DiscontinuityInputs<K>,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    start: usize,
) -> Option<Vec<u8>> {
    let mut current = start;
    let mut label = Vec::new();
    let mut is_cycle = false;

    loop {
        let unitig = &inputs.unitigs[half_ends[current].unitig_index()];
        let unitig_label = unitig.label(inputs);
        let reverse = reverse_for_stitch_node(current);
        let mut append_reverse = reverse;
        if !label.is_empty() && !labels_overlap_oriented_fast::<K>(&label, unitig_label, reverse) {
            let alternate = oriented_label(unitig_label, !reverse);
            if labels_overlap::<K>(&label, &alternate) {
                append_reverse = !reverse;
            }
        }
        append_or_init_oriented_fast::<K>(&mut label, unitig_label, append_reverse);

        let other = current ^ 1;
        let next = join_neighbor[other];
        if next == STITCH_NO_NODE {
            break;
        }
        if stitch_node_index(next) == start {
            is_cycle = true;
            break;
        }
        current = stitch_node_index(next);
    }

    if label.len() < K {
        None
    } else if is_cycle {
        Some(normalize_stitched_cycle::<K>(&label))
    } else {
        Some(canonical_label(label))
    }
}

fn map_simple_stitched_components_with_threads_from_adjacency(
    adjacency: &[StitchAdjacencyList],
    starts: &[usize],
    threads: usize,
) -> Vec<StitchedCoordRecord> {
    let workers = threads.max(1).min(starts.len().max(1));
    if workers == 1 || starts.len() < 1024 {
        let mut records = Vec::new();
        for (path_id, &start) in starts.iter().enumerate() {
            walk_simple_stitched_component_coords_from_adjacency(
                adjacency,
                start,
                path_id as u64,
                &mut records,
            );
        }
        return records;
    }

    let chunk_size = starts.len().div_ceil(workers);
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for (chunk_index, chunk) in starts.chunks(chunk_size).enumerate() {
            let base_path_id = (chunk_index * chunk_size) as u64;
            handles.push(scope.spawn(move || {
                let mut records = Vec::new();
                for (offset, &start) in chunk.iter().enumerate() {
                    walk_simple_stitched_component_coords_from_adjacency(
                        adjacency,
                        start,
                        base_path_id + offset as u64,
                        &mut records,
                    );
                }
                records
            }));
        }

        let mut records = Vec::with_capacity(starts.len());
        for handle in handles {
            records.extend(handle.join().expect("stitch mapper worker panicked"));
        }
        records
    })
}

fn write_stitched_coord_buckets(
    coord_dir: &Path,
    threads: usize,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    starts: &[usize],
) -> Result<Vec<StitchedCoordBucketEntry>, SerialCollationError> {
    if coord_dir.exists() {
        fs::remove_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
            path: coord_dir.to_path_buf(),
            source,
        })?;
    }
    fs::create_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
        path: coord_dir.to_path_buf(),
        source,
    })?;

    let bucket_count = stitched_coord_bucket_count(threads);
    let bucket_mask = bucket_count - 1;
    let workers = threads.max(1).min(starts.len().max(1));

    let mut manifest = if workers == 1 || starts.len() < 1024 {
        write_neighbor_stitched_coord_shards(
            coord_dir,
            0,
            bucket_count,
            bucket_mask,
            half_ends,
            join_neighbor,
            starts,
            0,
        )?
    } else {
        let chunk_size = starts.len().div_ceil(workers);
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for (chunk_index, chunk) in starts.chunks(chunk_size).enumerate() {
                let base_path_id = (chunk_index * chunk_size) as u64;
                handles.push(scope.spawn(move || {
                    write_neighbor_stitched_coord_shards(
                        coord_dir,
                        chunk_index,
                        bucket_count,
                        bucket_mask,
                        half_ends,
                        join_neighbor,
                        chunk,
                        base_path_id,
                    )
                }));
            }

            let mut manifest = Vec::new();
            for handle in handles {
                manifest.extend(
                    handle
                        .join()
                        .map_err(|_| SerialCollationError::WorkerPanic)??,
                );
            }
            Ok::<_, SerialCollationError>(manifest)
        })?
    };
    manifest.sort_by(|left, right| {
        left.bucket_id
            .cmp(&right.bucket_id)
            .then_with(|| left.path.cmp(&right.path))
    });

    Ok(manifest)
}

fn write_stitched_coord_buckets_from_adjacency(
    coord_dir: &Path,
    threads: usize,
    adjacency: &[StitchAdjacencyList],
    starts: &[usize],
) -> Result<Vec<StitchedCoordBucketEntry>, SerialCollationError> {
    if coord_dir.exists() {
        fs::remove_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
            path: coord_dir.to_path_buf(),
            source,
        })?;
    }
    fs::create_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
        path: coord_dir.to_path_buf(),
        source,
    })?;

    let bucket_count = stitched_coord_bucket_count(threads);
    let bucket_mask = bucket_count - 1;
    let workers = threads.max(1).min(starts.len().max(1));

    let mut manifest = if workers == 1 || starts.len() < 1024 {
        write_adjacency_stitched_coord_shards(
            coord_dir,
            0,
            bucket_count,
            bucket_mask,
            adjacency,
            starts,
            0,
        )?
    } else {
        let chunk_size = starts.len().div_ceil(workers);
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for (chunk_index, chunk) in starts.chunks(chunk_size).enumerate() {
                let base_path_id = (chunk_index * chunk_size) as u64;
                handles.push(scope.spawn(move || {
                    write_adjacency_stitched_coord_shards(
                        coord_dir,
                        chunk_index,
                        bucket_count,
                        bucket_mask,
                        adjacency,
                        chunk,
                        base_path_id,
                    )
                }));
            }

            let mut manifest = Vec::new();
            for handle in handles {
                manifest.extend(
                    handle
                        .join()
                        .map_err(|_| SerialCollationError::WorkerPanic)??,
                );
            }
            Ok::<_, SerialCollationError>(manifest)
        })?
    };
    manifest.sort_by(|left, right| {
        left.bucket_id
            .cmp(&right.bucket_id)
            .then_with(|| left.path.cmp(&right.path))
    });
    Ok(manifest)
}

fn write_neighbor_stitched_coord_shards(
    coord_dir: &Path,
    worker_id: usize,
    bucket_count: usize,
    bucket_mask: usize,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    starts: &[usize],
    base_path_id: u64,
) -> Result<Vec<StitchedCoordBucketEntry>, SerialCollationError> {
    let mut writers = StitchedCoordShardWriters::new(coord_dir, worker_id, bucket_count);
    let mut records = Vec::new();
    for (offset, &start) in starts.iter().enumerate() {
        records.clear();
        let path_id = base_path_id + offset as u64;
        walk_simple_stitched_component_coords(
            half_ends,
            join_neighbor,
            start,
            path_id,
            &mut records,
        );
        writers.write_path_records(stitched_coord_bucket(path_id, bucket_mask), &records)?;
    }
    writers.finish()
}

fn write_adjacency_stitched_coord_shards(
    coord_dir: &Path,
    worker_id: usize,
    bucket_count: usize,
    bucket_mask: usize,
    adjacency: &[StitchAdjacencyList],
    starts: &[usize],
    base_path_id: u64,
) -> Result<Vec<StitchedCoordBucketEntry>, SerialCollationError> {
    let mut writers = StitchedCoordShardWriters::new(coord_dir, worker_id, bucket_count);
    let mut records = Vec::new();
    for (offset, &start) in starts.iter().enumerate() {
        records.clear();
        let path_id = base_path_id + offset as u64;
        walk_simple_stitched_component_coords_from_adjacency(
            adjacency,
            start,
            path_id,
            &mut records,
        );
        writers.write_path_records(stitched_coord_bucket(path_id, bucket_mask), &records)?;
    }
    writers.finish()
}

fn write_materialized_stitched_coord_buckets<const K: usize>(
    inputs: &DiscontinuityInputs<K>,
    coord_dir: &Path,
    threads: usize,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    starts: &[usize],
) -> Result<Vec<MaterializedStitchedCoordBucketEntry>, SerialCollationError> {
    if coord_dir.exists() {
        fs::remove_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
            path: coord_dir.to_path_buf(),
            source,
        })?;
    }
    fs::create_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
        path: coord_dir.to_path_buf(),
        source,
    })?;

    let bucket_count = materialized_stitched_coord_bucket_count(threads);
    let bucket_mask = bucket_count - 1;
    let workers = threads.max(1).min(starts.len().max(1));

    let mut manifest = if workers == 1 || starts.len() < 1024 {
        write_materialized_neighbor_stitched_coord_shards(
            inputs,
            coord_dir,
            0,
            bucket_count,
            bucket_mask,
            half_ends,
            join_neighbor,
            starts,
            0,
        )?
    } else {
        let chunk_size = starts.len().div_ceil(workers);
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for (chunk_index, chunk) in starts.chunks(chunk_size).enumerate() {
                let base_path_id = (chunk_index * chunk_size) as u64;
                handles.push(scope.spawn(move || {
                    write_materialized_neighbor_stitched_coord_shards(
                        inputs,
                        coord_dir,
                        chunk_index,
                        bucket_count,
                        bucket_mask,
                        half_ends,
                        join_neighbor,
                        chunk,
                        base_path_id,
                    )
                }));
            }

            let mut manifest = Vec::new();
            for handle in handles {
                manifest.extend(
                    handle
                        .join()
                        .map_err(|_| SerialCollationError::WorkerPanic)??,
                );
            }
            Ok::<_, SerialCollationError>(manifest)
        })?
    };
    manifest.sort_by(|left, right| {
        left.bucket_id
            .cmp(&right.bucket_id)
            .then_with(|| left.coord_path.cmp(&right.coord_path))
    });
    Ok(manifest)
}

fn write_external_materialized_stitched_coord_buckets_from_neighbors<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    coord_dir: &Path,
    threads: usize,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
) -> Result<Vec<MaterializedStitchedCoordBucketEntry>, SerialCollationError> {
    let estimated_path_info_bytes = (half_ends.len() / 2)
        .checked_mul(std::mem::size_of::<StitchedCoordRecord>())
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(coord_dir.to_path_buf()))?;
    let inmem_limit = in_memory_stitch_path_info_limit();
    if estimated_path_info_bytes <= inmem_limit {
        return write_external_materialized_stitched_coord_buckets_from_neighbor_path_info::<K>(
            inputs,
            coord_dir,
            threads,
            half_ends,
            join_neighbor,
            estimated_path_info_bytes,
        );
    }

    eprintln!(
        "cuttlefish3-rs: external stitch path-info estimate {:.1} MiB exceeds in-memory limit {:.1} MiB; using disk path-info buckets",
        estimated_path_info_bytes as f64 / (1024.0 * 1024.0),
        inmem_limit as f64 / (1024.0 * 1024.0)
    );

    if std::env::var_os("CF3_RS_NEIGHBOR_DISK_PATH_INFO").is_none() {
        let component_starts = stitch_component_starts_from_neighbors(half_ends, join_neighbor);
        return write_external_materialized_stitched_coord_buckets(
            inputs,
            coord_dir,
            threads,
            half_ends,
            join_neighbor,
            &component_starts,
        );
    }

    if coord_dir.exists() {
        fs::remove_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
            path: coord_dir.to_path_buf(),
            source,
        })?;
    }
    fs::create_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
        path: coord_dir.to_path_buf(),
        source,
    })?;

    let path_info_started = Instant::now();
    let path_info_manifest = write_external_stitched_path_info_buckets_from_neighbors(
        inputs,
        coord_dir,
        threads,
        half_ends,
        join_neighbor,
    )?;
    let path_info_elapsed = path_info_started.elapsed();

    let materialize_started = Instant::now();
    let manifest = materialize_external_stitched_coord_buckets_from_path_info::<K>(
        inputs,
        coord_dir,
        threads,
        &path_info_manifest,
    )?;
    eprintln!(
        "cuttlefish3-rs: external neighbor disk path-info {:.3}s, label materialize {:.3}s",
        path_info_elapsed.as_secs_f64(),
        materialize_started.elapsed().as_secs_f64()
    );
    Ok(manifest)
}

fn write_external_materialized_stitched_coord_buckets<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    coord_dir: &Path,
    threads: usize,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    starts: &[usize],
) -> Result<Vec<MaterializedStitchedCoordBucketEntry>, SerialCollationError> {
    let estimated_path_info_bytes = (half_ends.len() / 2)
        .checked_mul(std::mem::size_of::<StitchedCoordRecord>())
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(coord_dir.to_path_buf()))?;
    let inmem_limit = in_memory_stitch_path_info_limit();
    if estimated_path_info_bytes <= inmem_limit {
        return write_external_materialized_stitched_coord_buckets_in_memory::<K>(
            inputs,
            coord_dir,
            threads,
            half_ends,
            join_neighbor,
            starts,
            estimated_path_info_bytes,
        );
    }

    eprintln!(
        "cuttlefish3-rs: external stitch path-info estimate {:.1} MiB exceeds in-memory limit {:.1} MiB; using disk path-info buckets",
        estimated_path_info_bytes as f64 / (1024.0 * 1024.0),
        inmem_limit as f64 / (1024.0 * 1024.0)
    );

    if coord_dir.exists() {
        fs::remove_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
            path: coord_dir.to_path_buf(),
            source,
        })?;
    }
    fs::create_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
        path: coord_dir.to_path_buf(),
        source,
    })?;

    let path_info_started = Instant::now();
    let path_info_manifest = write_external_stitched_path_info_buckets::<K>(
        inputs,
        coord_dir,
        threads,
        half_ends,
        join_neighbor,
        starts,
    )?;
    let path_info_elapsed = path_info_started.elapsed();

    let materialize_started = Instant::now();
    let manifest = materialize_external_stitched_coord_buckets_from_path_info::<K>(
        inputs,
        coord_dir,
        threads,
        &path_info_manifest,
    )?;
    eprintln!(
        "cuttlefish3-rs: external materialized stitch path-info {:.3}s, label materialize {:.3}s",
        path_info_elapsed.as_secs_f64(),
        materialize_started.elapsed().as_secs_f64()
    );
    Ok(manifest)
}

#[allow(clippy::too_many_arguments)]
fn write_external_materialized_stitched_coord_buckets_in_memory<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    coord_dir: &Path,
    threads: usize,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    starts: &[usize],
    estimated_path_info_bytes: usize,
) -> Result<Vec<MaterializedStitchedCoordBucketEntry>, SerialCollationError> {
    if coord_dir.exists() {
        fs::remove_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
            path: coord_dir.to_path_buf(),
            source,
        })?;
    }
    fs::create_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
        path: coord_dir.to_path_buf(),
        source,
    })?;

    let ranges_per_bucket = ranges_per_path_info_bucket(inputs.ranges.len(), threads.max(1));
    let path_info_bucket_count = inputs.ranges.len().div_ceil(ranges_per_bucket).max(1);
    let path_info_started = Instant::now();
    let records_by_range_bucket = collect_external_stitched_path_info_by_range_bucket(
        inputs,
        threads,
        half_ends,
        join_neighbor,
        starts,
        ranges_per_bucket,
        path_info_bucket_count,
    )?;
    let path_info_elapsed = path_info_started.elapsed();

    let materialize_started = Instant::now();
    let manifest = materialize_external_stitched_coord_buckets_from_memory::<K>(
        inputs,
        coord_dir,
        threads,
        ranges_per_bucket,
        records_by_range_bucket,
    )?;
    eprintln!(
        "cuttlefish3-rs: external in-memory stitch path-info {:.3}s ({:.1} MiB est), label materialize {:.3}s",
        path_info_elapsed.as_secs_f64(),
        estimated_path_info_bytes as f64 / (1024.0 * 1024.0),
        materialize_started.elapsed().as_secs_f64()
    );
    Ok(manifest)
}

fn write_external_materialized_stitched_coord_buckets_from_neighbor_path_info<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    coord_dir: &Path,
    threads: usize,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    estimated_path_info_bytes: usize,
) -> Result<Vec<MaterializedStitchedCoordBucketEntry>, SerialCollationError> {
    if coord_dir.exists() {
        fs::remove_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
            path: coord_dir.to_path_buf(),
            source,
        })?;
    }
    fs::create_dir_all(coord_dir).map_err(|source| SerialCollationError::Io {
        path: coord_dir.to_path_buf(),
        source,
    })?;

    let ranges_per_bucket = ranges_per_path_info_bucket(inputs.ranges.len(), threads.max(1));
    let path_info_bucket_count = inputs.ranges.len().div_ceil(ranges_per_bucket).max(1);
    let path_info_started = Instant::now();
    let records_by_range_bucket = collect_external_stitched_path_info_from_neighbors(
        inputs,
        half_ends,
        join_neighbor,
        ranges_per_bucket,
        path_info_bucket_count,
    )?;
    let path_info_elapsed = path_info_started.elapsed();

    let materialize_started = Instant::now();
    let manifest = materialize_external_stitched_coord_buckets_from_memory::<K>(
        inputs,
        coord_dir,
        threads,
        ranges_per_bucket,
        records_by_range_bucket,
    )?;
    eprintln!(
        "cuttlefish3-rs: external neighbor path-info {:.3}s ({:.1} MiB est), label materialize {:.3}s",
        path_info_elapsed.as_secs_f64(),
        estimated_path_info_bytes as f64 / (1024.0 * 1024.0),
        materialize_started.elapsed().as_secs_f64()
    );
    Ok(manifest)
}

fn in_memory_stitch_path_info_limit() -> usize {
    if std::env::var_os("CF3_RS_DISABLE_INMEM_STITCH_PATH_INFO").is_some() {
        return 0;
    }
    std::env::var("CF3_RS_INMEM_STITCH_PATH_INFO_BYTES")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(DEFAULT_INMEM_STITCH_PATH_INFO_LIMIT)
}

fn collect_external_stitched_path_info_by_range_bucket<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    threads: usize,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    starts: &[usize],
    ranges_per_bucket: usize,
    path_info_bucket_count: usize,
) -> Result<Vec<Vec<StitchedCoordRecord>>, SerialCollationError> {
    let workers = threads.max(1).min(starts.len().max(1));
    if workers == 1 || starts.len() < 1024 {
        let mut records_by_bucket = empty_stitched_record_buckets(path_info_bucket_count);
        collect_external_stitched_path_info_range_into_memory(
            &inputs.ranges,
            &inputs.unitig_path,
            ranges_per_bucket,
            half_ends,
            join_neighbor,
            starts,
            0,
            &mut records_by_bucket,
        )?;
        return Ok(records_by_bucket);
    }

    let chunk_size = starts.len().div_ceil(workers);
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for (chunk_index, chunk) in starts.chunks(chunk_size).enumerate() {
            let base_path_id = (chunk_index * chunk_size) as u64;
            handles.push(scope.spawn(move || {
                let mut records_by_bucket = empty_stitched_record_buckets(path_info_bucket_count);
                collect_external_stitched_path_info_range_into_memory(
                    &inputs.ranges,
                    &inputs.unitig_path,
                    ranges_per_bucket,
                    half_ends,
                    join_neighbor,
                    chunk,
                    base_path_id,
                    &mut records_by_bucket,
                )?;
                Ok::<_, SerialCollationError>(records_by_bucket)
            }));
        }

        let mut merged = empty_stitched_record_buckets(path_info_bucket_count);
        for handle in handles {
            let mut worker_buckets = handle
                .join()
                .map_err(|_| SerialCollationError::WorkerPanic)??;
            for (bucket_id, bucket) in worker_buckets.iter_mut().enumerate() {
                merged[bucket_id].append(bucket);
            }
        }
        Ok::<_, SerialCollationError>(merged)
    })
}

fn empty_stitched_record_buckets(bucket_count: usize) -> Vec<Vec<StitchedCoordRecord>> {
    let mut buckets = Vec::with_capacity(bucket_count);
    buckets.resize_with(bucket_count, Vec::new);
    buckets
}

fn collect_external_stitched_path_info_from_neighbors<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    ranges_per_bucket: usize,
    path_info_bucket_count: usize,
) -> Result<Vec<Vec<StitchedCoordRecord>>, SerialCollationError> {
    debug_assert_eq!(half_ends.len(), join_neighbor.len());
    debug_assert_eq!(half_ends.len() % 2, 0);
    let mut visited = vec![0u8; half_ends.len().div_ceil(2)];
    let mut records_by_bucket = empty_stitched_record_buckets(path_info_bucket_count);
    let mut path_records = Vec::new();
    let unitig_bucket = unitig_path_info_bucket_map(
        &inputs.ranges,
        ranges_per_bucket,
        path_info_bucket_count,
        inputs.unitig_count(),
        &inputs.unitig_path,
    )?;
    let mut path_id = 0u64;

    for node in 0..join_neighbor.len() {
        if join_neighbor[node] != STITCH_NO_NODE {
            continue;
        }
        let unitig_pair = node >> 1;
        if visited[unitig_pair] != 0 {
            continue;
        }
        path_records.clear();
        walk_simple_stitched_component_coords_marked(
            half_ends,
            join_neighbor,
            node,
            path_id,
            &mut visited,
            &mut path_records,
        );
        if !path_records.is_empty() {
            push_stitched_records_to_range_buckets(
                &unitig_bucket,
                &inputs.unitig_path,
                &path_records,
                &mut records_by_bucket,
            )?;
            path_id += 1;
        }
    }

    for node in (0..join_neighbor.len()).step_by(2) {
        let unitig_pair = node >> 1;
        if visited[unitig_pair] != 0 {
            continue;
        }
        path_records.clear();
        walk_simple_stitched_component_coords_marked(
            half_ends,
            join_neighbor,
            node,
            path_id,
            &mut visited,
            &mut path_records,
        );
        if !path_records.is_empty() {
            push_stitched_records_to_range_buckets(
                &unitig_bucket,
                &inputs.unitig_path,
                &path_records,
                &mut records_by_bucket,
            )?;
            path_id += 1;
        }
    }

    Ok(records_by_bucket)
}

fn unitig_path_info_bucket_map(
    ranges: &[ExternalLocalUnitigRange],
    ranges_per_bucket: usize,
    path_info_bucket_count: usize,
    unitig_count: usize,
    malformed_path: &Path,
) -> Result<Vec<u32>, SerialCollationError> {
    let bucket_limit = u32::try_from(path_info_bucket_count)
        .map_err(|_| SerialCollationError::MalformedCoordBucket(malformed_path.to_path_buf()))?;
    let mut unitig_bucket = vec![bucket_limit; unitig_count];
    for (range_id, range) in ranges.iter().enumerate() {
        let bucket_id = u32::try_from(range_id / ranges_per_bucket).map_err(|_| {
            SerialCollationError::MalformedCoordBucket(malformed_path.to_path_buf())
        })?;
        let end = range
            .start_unitig
            .checked_add(range.unitigs)
            .ok_or_else(|| {
                SerialCollationError::MalformedCoordBucket(malformed_path.to_path_buf())
            })?;
        if end > unitig_bucket.len() {
            return Err(SerialCollationError::MalformedCoordBucket(
                malformed_path.to_path_buf(),
            ));
        }
        unitig_bucket[range.start_unitig..end].fill(bucket_id);
    }
    Ok(unitig_bucket)
}

fn walk_simple_stitched_component_coords_marked(
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    start: usize,
    path_id: u64,
    visited: &mut [u8],
    records: &mut Vec<StitchedCoordRecord>,
) {
    let record_start = records.len();
    let mut current = start;
    let mut rank = 0u64;
    let mut is_cycle = false;

    loop {
        let unitig_pair = current >> 1;
        if visited[unitig_pair] != 0 {
            break;
        }
        visited[unitig_pair] = 1;

        let unitig_index = half_ends[current].unitig_index;
        let reverse = reverse_for_stitch_node(current);
        records.push(StitchedCoordRecord {
            path_id,
            rank,
            unitig_index,
            reverse,
            is_cycle: false,
        });
        rank += 1;

        let other = current ^ 1;
        let next = join_neighbor[other];
        if next == STITCH_NO_NODE {
            break;
        }
        if stitch_node_index(next) == start {
            is_cycle = true;
            break;
        }
        current = stitch_node_index(next);
    }

    if rank == 0 {
        records.truncate(record_start);
        return;
    }

    if is_cycle {
        for record in &mut records[record_start..] {
            record.is_cycle = true;
        }
    }
}

fn push_stitched_records_to_range_buckets(
    unitig_bucket: &[u32],
    malformed_path: &Path,
    records: &[StitchedCoordRecord],
    records_by_bucket: &mut [Vec<StitchedCoordRecord>],
) -> Result<(), SerialCollationError> {
    for &record in records {
        let bucket_id = unitig_bucket
            .get(record.unitig_index as usize)
            .copied()
            .ok_or_else(|| {
                SerialCollationError::MalformedCoordBucket(malformed_path.to_path_buf())
            })?;
        let bucket_id = usize::try_from(bucket_id).map_err(|_| {
            SerialCollationError::MalformedCoordBucket(malformed_path.to_path_buf())
        })?;
        if bucket_id >= records_by_bucket.len() {
            return Err(SerialCollationError::MalformedCoordBucket(
                malformed_path.to_path_buf(),
            ));
        }
        records_by_bucket[bucket_id].push(record);
    }
    Ok(())
}

fn write_stitched_records_to_range_bucket_writers_by_range(
    ranges: &[ExternalLocalUnitigRange],
    ranges_per_bucket: usize,
    malformed_path: &Path,
    records: &[StitchedCoordRecord],
    writers: &mut StitchedCoordShardWriters<'_>,
) -> Result<(), SerialCollationError> {
    for &record in records {
        let range_id = external_range_id_for_unitig(ranges, record.unitig_index as usize)
            .ok_or_else(|| {
                SerialCollationError::MalformedCoordBucket(malformed_path.to_path_buf())
            })?;
        let bucket_id = range_id / ranges_per_bucket;
        writers.write_record(bucket_id, record)?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn collect_external_stitched_path_info_range_into_memory(
    ranges: &[ExternalLocalUnitigRange],
    malformed_path: &Path,
    ranges_per_bucket: usize,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    starts: &[usize],
    base_path_id: u64,
    records_by_bucket: &mut [Vec<StitchedCoordRecord>],
) -> Result<(), SerialCollationError> {
    let mut records = Vec::new();
    for (offset, &start) in starts.iter().enumerate() {
        records.clear();
        let path_id = base_path_id + offset as u64;
        walk_simple_stitched_component_coords(
            half_ends,
            join_neighbor,
            start,
            path_id,
            &mut records,
        );
        for &record in &records {
            let range_id = external_range_id_for_unitig(ranges, record.unitig_index as usize)
                .ok_or_else(|| {
                    SerialCollationError::MalformedCoordBucket(malformed_path.to_path_buf())
                })?;
            let bucket_id = range_id / ranges_per_bucket;
            records_by_bucket[bucket_id].push(record);
        }
    }
    Ok(())
}

fn materialize_external_stitched_coord_buckets_from_memory<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    coord_dir: &Path,
    threads: usize,
    ranges_per_bucket: usize,
    records_by_range_bucket: Vec<Vec<StitchedCoordRecord>>,
) -> Result<Vec<MaterializedStitchedCoordBucketEntry>, SerialCollationError> {
    let bucket_count = materialized_stitched_coord_bucket_count(threads);
    let bucket_mask = bucket_count - 1;
    let non_empty_buckets = records_by_range_bucket
        .iter()
        .enumerate()
        .filter_map(|(bucket_id, records)| (!records.is_empty()).then_some(bucket_id))
        .collect::<Vec<_>>();

    let workers = threads.max(1).min(non_empty_buckets.len().max(1));
    let mut manifest = if workers == 1 || non_empty_buckets.len() < 2 {
        let mut writers = MaterializedStitchedCoordShardWriters::new(coord_dir, 0, bucket_count);
        for bucket_id in non_empty_buckets {
            materialize_external_stitched_coord_range_group_from_records::<K>(
                inputs,
                bucket_mask,
                ranges_per_bucket,
                bucket_id,
                &records_by_range_bucket[bucket_id],
                &mut writers,
            )?;
        }
        writers.finish()?
    } else {
        let next_group = AtomicUsize::new(0);
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for worker_id in 0..workers {
                let next_group = &next_group;
                let non_empty_buckets = &non_empty_buckets;
                let records_by_range_bucket = &records_by_range_bucket;
                handles.push(scope.spawn(move || {
                    let mut writers = MaterializedStitchedCoordShardWriters::new(
                        coord_dir,
                        worker_id,
                        bucket_count,
                    );
                    loop {
                        let group_idx = next_group.fetch_add(1, Ordering::Relaxed);
                        let Some(&bucket_id) = non_empty_buckets.get(group_idx) else {
                            break;
                        };
                        materialize_external_stitched_coord_range_group_from_records::<K>(
                            inputs,
                            bucket_mask,
                            ranges_per_bucket,
                            bucket_id,
                            &records_by_range_bucket[bucket_id],
                            &mut writers,
                        )?;
                    }
                    writers.finish()
                }));
            }

            let mut manifest = Vec::new();
            for handle in handles {
                manifest.extend(
                    handle
                        .join()
                        .map_err(|_| SerialCollationError::WorkerPanic)??,
                );
            }
            Ok::<_, SerialCollationError>(manifest)
        })?
    };

    manifest.sort_by(|left, right| {
        left.bucket_id
            .cmp(&right.bucket_id)
            .then_with(|| left.coord_path.cmp(&right.coord_path))
    });
    Ok(manifest)
}

#[allow(clippy::too_many_arguments)]
fn write_external_stitched_path_info_buckets<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    coord_dir: &Path,
    threads: usize,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    starts: &[usize],
) -> Result<Vec<StitchedCoordBucketEntry>, SerialCollationError> {
    let path_info_dir = coord_dir.join("path-info");
    fs::create_dir_all(&path_info_dir).map_err(|source| SerialCollationError::Io {
        path: path_info_dir.clone(),
        source,
    })?;
    let ranges_per_bucket = ranges_per_path_info_bucket(inputs.ranges.len(), threads.max(1));
    let path_info_bucket_count = inputs.ranges.len().div_ceil(ranges_per_bucket).max(1);
    let workers = threads.max(1).min(starts.len().max(1));

    let mut manifest = if workers == 1 || starts.len() < 1024 {
        write_external_stitched_path_info_range(
            &path_info_dir,
            0,
            &inputs.ranges,
            ranges_per_bucket,
            path_info_bucket_count,
            half_ends,
            join_neighbor,
            starts,
            0,
        )?
    } else {
        let chunk_size = starts.len().div_ceil(workers);
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for (chunk_index, chunk) in starts.chunks(chunk_size).enumerate() {
                let base_path_id = (chunk_index * chunk_size) as u64;
                let worker_dir = path_info_dir.clone();
                handles.push(scope.spawn(move || {
                    write_external_stitched_path_info_range(
                        &worker_dir,
                        chunk_index,
                        &inputs.ranges,
                        ranges_per_bucket,
                        path_info_bucket_count,
                        half_ends,
                        join_neighbor,
                        chunk,
                        base_path_id,
                    )
                }));
            }

            let mut manifest = Vec::new();
            for handle in handles {
                manifest.extend(
                    handle
                        .join()
                        .map_err(|_| SerialCollationError::WorkerPanic)??,
                );
            }
            Ok::<_, SerialCollationError>(manifest)
        })?
    };

    manifest.sort_by(|left, right| {
        left.bucket_id
            .cmp(&right.bucket_id)
            .then_with(|| left.path.cmp(&right.path))
    });
    Ok(manifest)
}

fn write_external_stitched_path_info_buckets_from_neighbors<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    coord_dir: &Path,
    threads: usize,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
) -> Result<Vec<StitchedCoordBucketEntry>, SerialCollationError> {
    debug_assert_eq!(half_ends.len(), join_neighbor.len());
    debug_assert_eq!(half_ends.len() % 2, 0);

    let path_info_dir = coord_dir.join("path-info");
    fs::create_dir_all(&path_info_dir).map_err(|source| SerialCollationError::Io {
        path: path_info_dir.clone(),
        source,
    })?;
    let ranges_per_bucket = ranges_per_path_info_bucket(inputs.ranges.len(), threads.max(1));
    let path_info_bucket_count = inputs.ranges.len().div_ceil(ranges_per_bucket).max(1);
    let mut writers = StitchedCoordShardWriters::new(&path_info_dir, 0, path_info_bucket_count);
    let mut visited = vec![0u8; half_ends.len().div_ceil(2)];
    let mut path_records = Vec::new();
    let mut path_id = 0u64;

    for node in 0..join_neighbor.len() {
        if join_neighbor[node] != STITCH_NO_NODE {
            continue;
        }
        let unitig_pair = node >> 1;
        if visited[unitig_pair] != 0 {
            continue;
        }
        path_records.clear();
        walk_simple_stitched_component_coords_marked(
            half_ends,
            join_neighbor,
            node,
            path_id,
            &mut visited,
            &mut path_records,
        );
        if !path_records.is_empty() {
            write_stitched_records_to_range_bucket_writers_by_range(
                &inputs.ranges,
                ranges_per_bucket,
                &inputs.unitig_path,
                &path_records,
                &mut writers,
            )?;
            path_id += 1;
        }
    }

    for node in (0..join_neighbor.len()).step_by(2) {
        let unitig_pair = node >> 1;
        if visited[unitig_pair] != 0 {
            continue;
        }
        path_records.clear();
        walk_simple_stitched_component_coords_marked(
            half_ends,
            join_neighbor,
            node,
            path_id,
            &mut visited,
            &mut path_records,
        );
        if !path_records.is_empty() {
            write_stitched_records_to_range_bucket_writers_by_range(
                &inputs.ranges,
                ranges_per_bucket,
                &inputs.unitig_path,
                &path_records,
                &mut writers,
            )?;
            path_id += 1;
        }
    }

    let mut manifest = writers.finish()?;
    manifest.sort_by(|left, right| {
        left.bucket_id
            .cmp(&right.bucket_id)
            .then_with(|| left.path.cmp(&right.path))
    });
    Ok(manifest)
}

fn write_external_stitched_path_info_range(
    path_info_dir: &Path,
    worker_id: usize,
    ranges: &[ExternalLocalUnitigRange],
    ranges_per_bucket: usize,
    path_info_bucket_count: usize,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    starts: &[usize],
    base_path_id: u64,
) -> Result<Vec<StitchedCoordBucketEntry>, SerialCollationError> {
    let mut writers =
        StitchedCoordShardWriters::new(path_info_dir, worker_id, path_info_bucket_count);
    let mut records = Vec::new();
    for (offset, &start) in starts.iter().enumerate() {
        records.clear();
        let path_id = base_path_id + offset as u64;
        walk_simple_stitched_component_coords(
            half_ends,
            join_neighbor,
            start,
            path_id,
            &mut records,
        );
        for record in &records {
            let range_id = external_range_id_for_unitig(ranges, record.unitig_index as usize)
                .ok_or_else(|| {
                    SerialCollationError::MalformedCoordBucket(path_info_dir.to_path_buf())
                })?;
            let bucket_id = range_id / ranges_per_bucket;
            writers.write_record(bucket_id, *record)?;
        }
    }
    writers.finish()
}

fn ranges_per_path_info_bucket(range_count: usize, workers: usize) -> usize {
    let target_buckets = range_count.div_ceil(STITCH_PATH_INFO_BUCKET_TARGET).max(1);
    // Range-bucket writers are shared across workers, so the descriptor limit
    // applies to the global bucket set rather than independently per worker.
    // Dividing it by worker count produced multi-gigabyte sort batches at
    // HumGut scale.
    let _ = workers;
    let bucket_count = target_buckets.min(MAX_OPEN_STITCH_PATH_INFO_WRITERS).max(1);
    range_count.div_ceil(bucket_count).max(1)
}

fn external_range_id_for_unitig(
    ranges: &[ExternalLocalUnitigRange],
    unitig_index: usize,
) -> Option<usize> {
    let range_id = ranges.partition_point(|range| range.start_unitig <= unitig_index);
    let range_id = range_id.checked_sub(1)?;
    let range = ranges.get(range_id)?;
    (unitig_index < range.start_unitig + range.unitigs).then_some(range_id)
}

struct ExternalRangeIndex {
    page_starts: Vec<usize>,
}

impl ExternalRangeIndex {
    const PAGE_SHIFT: usize = 16;

    fn new(ranges: &[ExternalLocalUnitigRange]) -> Self {
        let unitigs = ranges
            .last()
            .map_or(0, |range| range.start_unitig + range.unitigs);
        let page_count = unitigs.div_ceil(1 << Self::PAGE_SHIFT);
        let mut page_starts = Vec::with_capacity(page_count);
        let mut range_id = 0;
        for page in 0..page_count {
            let unitig = page << Self::PAGE_SHIFT;
            while range_id + 1 < ranges.len() && ranges[range_id + 1].start_unitig <= unitig {
                range_id += 1;
            }
            page_starts.push(range_id);
        }
        Self { page_starts }
    }

    #[inline(always)]
    fn find(&self, ranges: &[ExternalLocalUnitigRange], unitig_index: usize) -> Option<usize> {
        let mut range_id = *self.page_starts.get(unitig_index >> Self::PAGE_SHIFT)?;
        while range_id + 1 < ranges.len() && ranges[range_id + 1].start_unitig <= unitig_index {
            range_id += 1;
        }
        let range = ranges.get(range_id)?;
        (unitig_index < range.start_unitig + range.unitigs).then_some(range_id)
    }
}

fn materialize_external_stitched_coord_buckets_from_path_info<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    coord_dir: &Path,
    threads: usize,
    path_info_manifest: &[StitchedCoordBucketEntry],
) -> Result<Vec<MaterializedStitchedCoordBucketEntry>, SerialCollationError> {
    let bucket_count = materialized_stitched_coord_bucket_count(threads);
    let bucket_mask = bucket_count - 1;
    let ranges_per_bucket = ranges_per_path_info_bucket(inputs.ranges.len(), threads.max(1));

    let mut group_bounds = Vec::new();
    let mut group_start = 0;
    while group_start < path_info_manifest.len() {
        let bucket_id = path_info_manifest[group_start].bucket_id;
        let mut group_end = group_start + 1;
        while group_end < path_info_manifest.len()
            && path_info_manifest[group_end].bucket_id == bucket_id
        {
            group_end += 1;
        }
        group_bounds.push((bucket_id, group_start, group_end));
        group_start = group_end;
    }

    let workers = threads.max(1).min(group_bounds.len().max(1));
    let mut manifest = if workers == 1 || group_bounds.len() < 2 {
        let mut writers = MaterializedStitchedCoordShardWriters::new(coord_dir, 0, bucket_count);
        for &(bucket_id, start, end) in &group_bounds {
            materialize_external_stitched_coord_range_group::<K>(
                inputs,
                bucket_mask,
                ranges_per_bucket,
                bucket_id,
                &path_info_manifest[start..end],
                &mut writers,
            )?;
        }
        writers.finish()?
    } else {
        let next_group = AtomicUsize::new(0);
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for worker_id in 0..workers {
                let next_group = &next_group;
                let group_bounds = &group_bounds;
                handles.push(scope.spawn(move || {
                    let mut writers = MaterializedStitchedCoordShardWriters::new(
                        coord_dir,
                        worker_id,
                        bucket_count,
                    );
                    loop {
                        let group_idx = next_group.fetch_add(1, Ordering::Relaxed);
                        let Some(&(bucket_id, start, end)) = group_bounds.get(group_idx) else {
                            break;
                        };
                        materialize_external_stitched_coord_range_group::<K>(
                            inputs,
                            bucket_mask,
                            ranges_per_bucket,
                            bucket_id,
                            &path_info_manifest[start..end],
                            &mut writers,
                        )?;
                    }
                    writers.finish()
                }));
            }

            let mut manifest = Vec::new();
            for handle in handles {
                manifest.extend(
                    handle
                        .join()
                        .map_err(|_| SerialCollationError::WorkerPanic)??,
                );
            }
            Ok::<_, SerialCollationError>(manifest)
        })?
    };

    manifest.sort_by(|left, right| {
        left.bucket_id
            .cmp(&right.bucket_id)
            .then_with(|| left.coord_path.cmp(&right.coord_path))
    });
    Ok(manifest)
}

fn materialize_external_stitched_coord_range_group<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    bucket_mask: usize,
    ranges_per_bucket: usize,
    path_info_bucket_id: usize,
    path_info_entries: &[StitchedCoordBucketEntry],
    writers: &mut MaterializedStitchedCoordShardWriters<'_>,
) -> Result<(), SerialCollationError> {
    if path_info_entries.is_empty() {
        return Ok(());
    }
    let range_start = path_info_bucket_id
        .checked_mul(ranges_per_bucket)
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;
    let range_end = (range_start + ranges_per_bucket).min(inputs.ranges.len());
    let Some(first_range) = inputs.ranges.get(range_start) else {
        return Err(SerialCollationError::MalformedCoordBucket(
            inputs.unitig_path.clone(),
        ));
    };
    let last_range = inputs
        .ranges
        .get(range_end.saturating_sub(1))
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;
    let start_unitig = first_range.start_unitig;
    let end_unitig = last_range.start_unitig + last_range.unitigs;
    let unitig_span = end_unitig
        .checked_sub(start_unitig)
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;

    let mut path_info: Vec<Option<StitchedCoordRecord>> = vec![None; unitig_span];
    for entry in path_info_entries {
        for record in read_stitched_coord_bucket_file(entry)? {
            let unitig_index = record.unitig_index as usize;
            if unitig_index < start_unitig || unitig_index >= end_unitig {
                return Err(SerialCollationError::MalformedCoordBucket(
                    entry.path.clone(),
                ));
            }
            let slot = &mut path_info[unitig_index - start_unitig];
            if slot.replace(record).is_some() {
                return Err(SerialCollationError::MalformedCoordBucket(
                    entry.path.clone(),
                ));
            }
        }
    }

    let unitig_file =
        File::open(&inputs.unitig_path).map_err(|source| SerialCollationError::Io {
            path: inputs.unitig_path.clone(),
            source,
        })?;
    let mut unitig_input = BufReader::with_capacity(1024 * 1024, unitig_file);
    unitig_input
        .seek(SeekFrom::Start(
            (start_unitig * external_unitig_record_len::<K>(inputs.compact_unitigs)) as u64,
        ))
        .map_err(|source| SerialCollationError::Io {
            path: inputs.unitig_path.clone(),
            source,
        })?;
    let label_file = File::open(&inputs.label_path).map_err(|source| SerialCollationError::Io {
        path: inputs.label_path.clone(),
        source,
    })?;
    let mut label_input = BufReader::with_capacity(1024 * 1024, label_file);
    label_input
        .seek(SeekFrom::Start(first_range.label_start))
        .map_err(|source| SerialCollationError::Io {
            path: inputs.label_path.clone(),
            source,
        })?;
    let mut label = Vec::new();
    let mut discard = vec![0u8; 64 * 1024];

    for record in &path_info {
        let unitig: DiscontinuityUnitig<K> = read_discontinuity_unitig_from_reader(
            &mut unitig_input,
            &inputs.unitig_path,
            inputs.compact_unitigs,
        )?;
        if let Some(record) = record {
            label.resize(unitig.label_len as usize, 0);
            label_input
                .read_exact(&mut label)
                .map_err(|source| SerialCollationError::Io {
                    path: inputs.label_path.clone(),
                    source,
                })?;
            let bucket_id = stitched_coord_bucket(record.path_id, bucket_mask);
            writers.write_materialized_record(bucket_id, record, &label)?;
        } else {
            read_and_discard_exact(
                &mut label_input,
                &inputs.label_path,
                &mut discard,
                unitig.label_len as u64,
            )?;
        }
    }
    Ok(())
}

fn materialize_external_stitched_coord_range_group_from_records<const K: usize>(
    inputs: &ExternalDiscontinuityInputs<K>,
    bucket_mask: usize,
    ranges_per_bucket: usize,
    path_info_bucket_id: usize,
    path_info_records: &[StitchedCoordRecord],
    writers: &mut MaterializedStitchedCoordShardWriters<'_>,
) -> Result<(), SerialCollationError> {
    if path_info_records.is_empty() {
        return Ok(());
    }
    let range_start = path_info_bucket_id
        .checked_mul(ranges_per_bucket)
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;
    let range_end = (range_start + ranges_per_bucket).min(inputs.ranges.len());
    let Some(first_range) = inputs.ranges.get(range_start) else {
        return Err(SerialCollationError::MalformedCoordBucket(
            inputs.unitig_path.clone(),
        ));
    };
    let last_range = inputs
        .ranges
        .get(range_end.saturating_sub(1))
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;
    let start_unitig = first_range.start_unitig;
    let end_unitig = last_range.start_unitig + last_range.unitigs;
    let unitig_span = end_unitig
        .checked_sub(start_unitig)
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;

    let mut path_info = path_info_records.to_vec();
    let empty_record = StitchedCoordRecord {
        path_id: 0,
        rank: 0,
        unitig_index: u32::MAX,
        reverse: false,
        is_cycle: false,
    };
    let mut path_info_by_unitig = vec![empty_record; unitig_span];
    let expected_record_count = path_info.len();
    for record in path_info.drain(..) {
        let unitig_index = record.unitig_index as usize;
        if unitig_index < start_unitig || unitig_index >= end_unitig {
            return Err(SerialCollationError::MalformedCoordBucket(
                inputs.unitig_path.clone(),
            ));
        }
        let slot = &mut path_info_by_unitig[unitig_index - start_unitig];
        if slot.unitig_index != u32::MAX {
            return Err(SerialCollationError::MalformedCoordBucket(
                inputs.unitig_path.clone(),
            ));
        }
        *slot = record;
    }

    let unitig_file =
        File::open(&inputs.unitig_path).map_err(|source| SerialCollationError::Io {
            path: inputs.unitig_path.clone(),
            source,
        })?;
    let mut unitig_input = BufReader::with_capacity(1024 * 1024, unitig_file);
    unitig_input
        .seek(SeekFrom::Start(
            (start_unitig * external_unitig_record_len::<K>(inputs.compact_unitigs)) as u64,
        ))
        .map_err(|source| SerialCollationError::Io {
            path: inputs.unitig_path.clone(),
            source,
        })?;
    let label_file = File::open(&inputs.label_path).map_err(|source| SerialCollationError::Io {
        path: inputs.label_path.clone(),
        source,
    })?;
    let mut label_input = BufReader::with_capacity(1024 * 1024, label_file);
    label_input
        .seek(SeekFrom::Start(first_range.label_start))
        .map_err(|source| SerialCollationError::Io {
            path: inputs.label_path.clone(),
            source,
        })?;
    let mut label = Vec::new();
    let mut discard = vec![0u8; 64 * 1024];
    let mut record_count = 0usize;
    for unitig_index in start_unitig..end_unitig {
        let unitig: DiscontinuityUnitig<K> = read_discontinuity_unitig_from_reader(
            &mut unitig_input,
            &inputs.unitig_path,
            inputs.compact_unitigs,
        )?;
        let dense_record = &path_info_by_unitig[unitig_index - start_unitig];
        let record = (dense_record.unitig_index != u32::MAX).then_some(dense_record);
        if let Some(record) = record {
            label.resize(unitig.label_len as usize, 0);
            label_input
                .read_exact(&mut label)
                .map_err(|source| SerialCollationError::Io {
                    path: inputs.label_path.clone(),
                    source,
                })?;
            let bucket_id = stitched_coord_bucket(record.path_id, bucket_mask);
            writers.write_materialized_record(bucket_id, record, &label)?;
            record_count += 1;
        } else {
            read_and_discard_exact(
                &mut label_input,
                &inputs.label_path,
                &mut discard,
                unitig.label_len as u64,
            )?;
        }
    }
    if record_count != expected_record_count {
        return Err(SerialCollationError::MalformedCoordBucket(
            inputs.unitig_path.clone(),
        ));
    }
    debug_assert_eq!(unitig_span, end_unitig - start_unitig);
    Ok(())
}

fn map_external_cpp_path_info_range_bucket<const K: usize, F, W>(
    inputs: &ExternalDiscontinuityInputs<K>,
    bucket_mask: usize,
    ranges_per_bucket: usize,
    path_info_bucket_id: usize,
    path_info_records: &[StitchedCoordRecord],
    writers: &mut W,
    emit_direct_local: &mut F,
) -> Result<(), SerialCollationError>
where
    F: FnMut(FinalUnitigRecord) -> Result<(), SerialCollationError>,
    W: MaterializedRecordSink,
{
    map_external_cpp_path_info_range_bucket_owned(
        inputs,
        bucket_mask,
        ranges_per_bucket,
        path_info_bucket_id,
        path_info_records.to_vec(),
        writers,
        emit_direct_local,
    )
}

fn map_external_cpp_path_info_range_bucket_owned<const K: usize, F, W>(
    inputs: &ExternalDiscontinuityInputs<K>,
    bucket_mask: usize,
    ranges_per_bucket: usize,
    path_info_bucket_id: usize,
    mut path_info: Vec<StitchedCoordRecord>,
    writers: &mut W,
    emit_direct_local: &mut F,
) -> Result<(), SerialCollationError>
where
    F: FnMut(FinalUnitigRecord) -> Result<(), SerialCollationError>,
    W: MaterializedRecordSink,
{
    let range_start = path_info_bucket_id
        .checked_mul(ranges_per_bucket)
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;
    let range_end = (range_start + ranges_per_bucket).min(inputs.ranges.len());
    let Some(first_range) = inputs.ranges.get(range_start) else {
        if path_info.is_empty() {
            return Ok(());
        }
        return Err(SerialCollationError::MalformedCoordBucket(
            inputs.unitig_path.clone(),
        ));
    };
    let last_range = inputs
        .ranges
        .get(range_end.saturating_sub(1))
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;
    let start_unitig = first_range.start_unitig;
    let end_unitig = last_range.start_unitig + last_range.unitigs;
    let unitig_span = end_unitig
        .checked_sub(start_unitig)
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(inputs.unitig_path.clone()))?;

    path_info.sort_unstable_by_key(|record| record.unitig_index);
    for pair in path_info.windows(2) {
        if pair[0].unitig_index == pair[1].unitig_index {
            return Err(SerialCollationError::MalformedCoordBucket(
                inputs.unitig_path.clone(),
            ));
        }
    }
    for record in &path_info {
        let unitig_index = record.unitig_index as usize;
        if unitig_index < start_unitig || unitig_index >= end_unitig {
            return Err(SerialCollationError::MalformedCoordBucket(
                inputs.unitig_path.clone(),
            ));
        }
    }

    let unitig_file =
        File::open(&inputs.unitig_path).map_err(|source| SerialCollationError::Io {
            path: inputs.unitig_path.clone(),
            source,
        })?;
    let mut unitig_input = BufReader::with_capacity(1024 * 1024, unitig_file);
    unitig_input
        .seek(SeekFrom::Start(
            (start_unitig * external_unitig_record_len::<K>(inputs.compact_unitigs)) as u64,
        ))
        .map_err(|source| SerialCollationError::Io {
            path: inputs.unitig_path.clone(),
            source,
        })?;
    let label_file = File::open(&inputs.label_path).map_err(|source| SerialCollationError::Io {
        path: inputs.label_path.clone(),
        source,
    })?;
    let mut label_input = BufReader::with_capacity(1024 * 1024, label_file);
    label_input
        .seek(SeekFrom::Start(first_range.label_start))
        .map_err(|source| SerialCollationError::Io {
            path: inputs.label_path.clone(),
            source,
        })?;
    let mut label = Vec::new();
    let mut discard = vec![0u8; 64 * 1024];
    let mut record_idx = 0usize;
    let mut color_reader = inputs
        .color_runs
        .as_ref()
        .map(|sidecar| sidecar.reader_at(first_range.color_start))
        .transpose()?;
    let mut colors = Vec::new();
    let mut current_range = range_start;

    for unitig_index in start_unitig..end_unitig {
        if inputs
            .ranges
            .get(current_range)
            .is_some_and(|range| range.start_unitig == unitig_index)
        {
            let range = &inputs.ranges[current_range];
            label_input
                .seek(SeekFrom::Start(range.label_start))
                .map_err(|source| SerialCollationError::Io {
                    path: inputs.label_path.clone(),
                    source,
                })?;
            color_reader = inputs
                .color_runs
                .as_ref()
                .map(|sidecar| sidecar.reader_at(range.color_start))
                .transpose()?;
            current_range += 1;
        }
        let unitig: DiscontinuityUnitig<K> = read_discontinuity_unitig_from_reader(
            &mut unitig_input,
            &inputs.unitig_path,
            inputs.compact_unitigs,
        )?;
        let has_colors = if let Some(reader) = color_reader.as_mut() {
            reader.read_next_into(&mut colors)?;
            true
        } else {
            false
        };
        let record = path_info
            .get(record_idx)
            .filter(|record| record.unitig_index as usize == unitig_index);
        if let Some(record) = record {
            label.resize(unitig.label_len as usize, 0);
            label_input
                .read_exact(&mut label)
                .map_err(|source| SerialCollationError::Io {
                    path: inputs.label_path.clone(),
                    source,
                })?;
            let bucket_id = stitched_coord_bucket(record.path_id, bucket_mask);
            write_external_materialized_record(
                writers,
                bucket_id,
                record,
                &label,
                has_colors.then_some(colors.as_slice()),
            )?;
            record_idx += 1;
        } else if unitig.left_exit().is_none() && unitig.right_exit().is_none() {
            label.resize(unitig.label_len as usize, 0);
            label_input
                .read_exact(&mut label)
                .map_err(|source| SerialCollationError::Io {
                    path: inputs.label_path.clone(),
                    source,
                })?;
            let reverse = reverse_complement_is_less(&label);
            let label = if reverse {
                reverse_complement_label(&label)
            } else {
                label.clone()
            };
            let mut direct_colors = if has_colors {
                std::mem::take(&mut colors)
            } else {
                Vec::new()
            };
            if reverse && !direct_colors.is_empty() {
                direct_colors =
                    reverse_color_runs(&direct_colors, (unitig.label_len as usize - K + 1) as u32);
            }
            emit_direct_local(FinalUnitigRecord {
                label,
                colors: direct_colors,
            })?;
        } else {
            read_and_discard_exact(
                &mut label_input,
                &inputs.label_path,
                &mut discard,
                unitig.label_len as u64,
            )?;
        }
    }
    if record_idx != path_info.len() {
        return Err(SerialCollationError::MalformedCoordBucket(
            inputs.unitig_path.clone(),
        ));
    }
    debug_assert_eq!(unitig_span, end_unitig - start_unitig);
    Ok(())
}

fn write_external_materialized_record<W: MaterializedRecordSink>(
    writers: &mut W,
    bucket_id: usize,
    record: &StitchedCoordRecord,
    label: &[u8],
    colors: Option<&[UnitigColor]>,
) -> Result<(), SerialCollationError> {
    if let Some(colors) = colors {
        writers.write_materialized_colored_record(bucket_id, record, label, colors)
    } else {
        writers.write_materialized_record(bucket_id, record, label)
    }
}

fn read_and_discard_exact(
    input: &mut BufReader<File>,
    path: &Path,
    scratch: &mut [u8],
    mut len: u64,
) -> Result<(), SerialCollationError> {
    while len != 0 {
        let chunk_len = scratch.len().min(len as usize);
        input
            .read_exact(&mut scratch[..chunk_len])
            .map_err(|source| SerialCollationError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        len -= chunk_len as u64;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_materialized_neighbor_stitched_coord_shards<const K: usize>(
    inputs: &DiscontinuityInputs<K>,
    coord_dir: &Path,
    worker_id: usize,
    bucket_count: usize,
    bucket_mask: usize,
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    starts: &[usize],
    base_path_id: u64,
) -> Result<Vec<MaterializedStitchedCoordBucketEntry>, SerialCollationError> {
    let mut writers =
        MaterializedStitchedCoordShardWriters::new(coord_dir, worker_id, bucket_count);
    let mut records = Vec::new();
    for (offset, &start) in starts.iter().enumerate() {
        records.clear();
        let path_id = base_path_id + offset as u64;
        walk_simple_stitched_component_coords(
            half_ends,
            join_neighbor,
            start,
            path_id,
            &mut records,
        );
        writers.write_path_records(
            inputs,
            stitched_coord_bucket(path_id, bucket_mask),
            &records,
        )?;
    }
    writers.finish()
}

struct MaterializedStitchedCoordShardWriters<'a> {
    coord_dir: &'a Path,
    worker_id: usize,
    writers: Vec<Option<MaterializedStitchedCoordShardWriter>>,
    open_writers: usize,
}

trait MaterializedRecordSink {
    fn write_materialized_record(
        &mut self,
        bucket_id: usize,
        record: &StitchedCoordRecord,
        label: &[u8],
    ) -> Result<(), SerialCollationError>;

    fn write_materialized_colored_record(
        &mut self,
        bucket_id: usize,
        record: &StitchedCoordRecord,
        label: &[u8],
        colors: &[UnitigColor],
    ) -> Result<(), SerialCollationError>;
}

impl MaterializedRecordSink for MaterializedStitchedCoordShardWriters<'_> {
    fn write_materialized_record(
        &mut self,
        bucket_id: usize,
        record: &StitchedCoordRecord,
        label: &[u8],
    ) -> Result<(), SerialCollationError> {
        MaterializedStitchedCoordShardWriters::write_materialized_record(
            self, bucket_id, record, label,
        )
    }

    fn write_materialized_colored_record(
        &mut self,
        bucket_id: usize,
        record: &StitchedCoordRecord,
        label: &[u8],
        colors: &[UnitigColor],
    ) -> Result<(), SerialCollationError> {
        self.ensure_writer(bucket_id)?;
        self.writers[bucket_id]
            .as_mut()
            .expect("bucket writer was just created")
            .write_colored_record(record, label, colors)
    }
}

struct SharedMaterializedWriters<'a> {
    coord_dir: &'a Path,
    writers: Vec<Mutex<Option<MaterializedStitchedCoordShardWriter>>>,
    retained: Vec<Mutex<Vec<PendingMaterializedBucket>>>,
    keep_open: bool,
}

#[derive(Default)]
struct PendingMaterializedBucket {
    records: Vec<LoadedMaterializedStitchedCoordRecord>,
    labels: Vec<u8>,
    colors: Vec<UnitigColor>,
}

impl<'a> SharedMaterializedWriters<'a> {
    fn new(coord_dir: &'a Path, bucket_count: usize) -> Self {
        Self {
            coord_dir,
            writers: (0..bucket_count).map(|_| Mutex::new(None)).collect(),
            retained: (0..bucket_count).map(|_| Mutex::new(Vec::new())).collect(),
            keep_open: raise_open_file_limit(),
        }
    }

    fn retain_batch(
        &self,
        bucket_id: usize,
        batch: &mut PendingMaterializedBucket,
    ) -> Result<(), SerialCollationError> {
        if batch.records.is_empty() {
            return Ok(());
        }
        self.retained[bucket_id]
            .lock()
            .map_err(|_| SerialCollationError::WorkerPanic)?
            .push(std::mem::take(batch));
        Ok(())
    }

    fn write_batch(
        &self,
        bucket_id: usize,
        batch: &mut PendingMaterializedBucket,
    ) -> Result<(), SerialCollationError> {
        if batch.records.is_empty() {
            return Ok(());
        }
        let mut writer = self.writers[bucket_id]
            .lock()
            .map_err(|_| SerialCollationError::WorkerPanic)?;
        if writer.is_none() {
            *writer = Some(MaterializedStitchedCoordShardWriter::create(
                self.coord_dir,
                0,
                bucket_id,
            )?);
        }
        let writer = writer.as_mut().expect("global writer was just created");
        writer.ensure_open()?;
        writer.write_pending_batch(batch)?;
        if !self.keep_open {
            writer.flush_and_close()?;
        }
        Ok(())
    }

    fn write_colored_record(
        &self,
        bucket_id: usize,
        record: &StitchedCoordRecord,
        label: &[u8],
        colors: &[UnitigColor],
    ) -> Result<(), SerialCollationError> {
        let mut writer = self.writers[bucket_id]
            .lock()
            .map_err(|_| SerialCollationError::WorkerPanic)?;
        if writer.is_none() {
            *writer = Some(MaterializedStitchedCoordShardWriter::create(
                self.coord_dir,
                0,
                bucket_id,
            )?);
        }
        let writer = writer.as_mut().expect("global writer was just created");
        writer.ensure_open()?;
        writer.write_colored_record(record, label, colors)?;
        if !self.keep_open {
            writer.flush_and_close()?;
        }
        Ok(())
    }

    fn finish(self) -> Result<MaterializedStitchedCoordBuckets, SerialCollationError> {
        let mut manifest = Vec::new();
        for writer in self.writers {
            if let Some(writer) = writer
                .into_inner()
                .map_err(|_| SerialCollationError::WorkerPanic)?
            {
                manifest.push(writer.finish()?);
            }
        }
        let retained = self
            .retained
            .into_iter()
            .map(|bucket| {
                bucket
                    .into_inner()
                    .map_err(|_| SerialCollationError::WorkerPanic)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(MaterializedStitchedCoordBuckets { manifest, retained })
    }
}

struct MaterializedStitchedCoordBuckets {
    manifest: Vec<MaterializedStitchedCoordBucketEntry>,
    retained: Vec<Vec<PendingMaterializedBucket>>,
}

struct SharedMaterializedBatch<'a, 'b> {
    shared: &'a SharedMaterializedWriters<'b>,
    buckets: Vec<PendingMaterializedBucket>,
}

impl<'a, 'b> SharedMaterializedBatch<'a, 'b> {
    // Match C++ Unitig_Coord_Bucket_Concurrent: keep worker-local collation
    // tails small so the 1024-way map does not retain gigabytes before reduce.
    const FLUSH_BYTES: usize = 8 * 1024;

    fn new(shared: &'a SharedMaterializedWriters<'b>, bucket_count: usize) -> Self {
        Self {
            shared,
            buckets: (0..bucket_count)
                .map(|_| PendingMaterializedBucket::default())
                .collect(),
        }
    }

    fn finish(mut self) -> Result<(), SerialCollationError> {
        for bucket_id in 0..self.buckets.len() {
            self.shared
                .retain_batch(bucket_id, &mut self.buckets[bucket_id])?;
        }
        Ok(())
    }

    fn finish_bucket(&mut self, bucket_id: usize) -> Result<(), SerialCollationError> {
        self.shared
            .write_batch(bucket_id, &mut self.buckets[bucket_id])?;
        Ok(())
    }

    #[inline]
    fn uncolored_bucket_ready(bucket: &PendingMaterializedBucket) -> bool {
        bucket.records.len() >= Self::FLUSH_BYTES / STITCH_COORD_RECORD_LEN as usize
            && bucket.labels.len() >= Self::FLUSH_BYTES
    }

    #[inline]
    fn colored_bucket_ready(bucket: &PendingMaterializedBucket) -> bool {
        Self::uncolored_bucket_ready(bucket)
            && bucket.colors.len() * std::mem::size_of::<UnitigColor>() >= Self::FLUSH_BYTES
    }
}

impl MaterializedRecordSink for SharedMaterializedBatch<'_, '_> {
    fn write_materialized_record(
        &mut self,
        bucket_id: usize,
        record: &StitchedCoordRecord,
        label: &[u8],
    ) -> Result<(), SerialCollationError> {
        let bucket = &mut self.buckets[bucket_id];
        let label_offset = u32::try_from(bucket.labels.len()).map_err(|_| {
            SerialCollationError::MalformedCoordBucket(self.shared.coord_dir.to_path_buf())
        })?;
        bucket.labels.extend_from_slice(label);
        bucket
            .records
            .push(LoadedMaterializedStitchedCoordRecord::new(
                record.path_id,
                record.rank,
                label_offset,
                label.len() as u32,
                record.reverse,
                record.is_cycle,
                u32::MAX,
                0,
            ));
        if Self::uncolored_bucket_ready(bucket) {
            self.shared
                .write_batch(bucket_id, &mut self.buckets[bucket_id])?;
        }
        Ok(())
    }

    fn write_materialized_colored_record(
        &mut self,
        bucket_id: usize,
        record: &StitchedCoordRecord,
        label: &[u8],
        colors: &[UnitigColor],
    ) -> Result<(), SerialCollationError> {
        let color_count = u32::try_from(colors.len()).map_err(|_| {
            SerialCollationError::MalformedCoordBucket(self.shared.coord_dir.to_path_buf())
        })?;
        let bucket = &mut self.buckets[bucket_id];
        let label_offset = u32::try_from(bucket.labels.len()).map_err(|_| {
            SerialCollationError::MalformedCoordBucket(self.shared.coord_dir.to_path_buf())
        })?;
        let color_start = bucket.colors.len() as u32;
        bucket.labels.extend_from_slice(label);
        bucket.colors.extend_from_slice(colors);
        bucket
            .records
            .push(LoadedMaterializedStitchedCoordRecord::new(
                record.path_id,
                record.rank,
                label_offset,
                label.len() as u32,
                record.reverse,
                record.is_cycle,
                color_start,
                color_count,
            ));
        if Self::colored_bucket_ready(bucket) {
            self.finish_bucket(bucket_id)?;
        }
        Ok(())
    }
}

impl<'a> MaterializedStitchedCoordShardWriters<'a> {
    fn new(coord_dir: &'a Path, worker_id: usize, bucket_count: usize) -> Self {
        let mut writers = Vec::with_capacity(bucket_count);
        writers.resize_with(bucket_count, || None);
        Self {
            coord_dir,
            worker_id,
            writers,
            open_writers: 0,
        }
    }

    fn write_path_records<const K: usize>(
        &mut self,
        inputs: &DiscontinuityInputs<K>,
        bucket_id: usize,
        records: &[StitchedCoordRecord],
    ) -> Result<(), SerialCollationError> {
        if records.is_empty() {
            return Ok(());
        }
        self.ensure_writer(bucket_id)?;
        let writer = self.writers[bucket_id]
            .as_mut()
            .expect("bucket writer was just created");
        writer.ensure_open()?;
        for record in records {
            let unitig = &inputs.unitigs[record.unitig_index as usize];
            writer.write_record(record, unitig.label(inputs))?;
        }
        Ok(())
    }

    fn ensure_writer(&mut self, bucket_id: usize) -> Result<(), SerialCollationError> {
        if self.writers[bucket_id].is_none() {
            self.evict_writer_if_needed(bucket_id)?;
            self.writers[bucket_id] = Some(MaterializedStitchedCoordShardWriter::create(
                self.coord_dir,
                self.worker_id,
                bucket_id,
            )?);
            self.open_writers += 1;
        } else if self
            .writers
            .get(bucket_id)
            .and_then(Option::as_ref)
            .is_some_and(|writer| !writer.is_open())
        {
            self.evict_writer_if_needed(bucket_id)?;
            self.writers[bucket_id]
                .as_mut()
                .expect("checked that writer exists")
                .ensure_open()?;
            self.open_writers += 1;
        }
        Ok(())
    }

    fn evict_writer_if_needed(
        &mut self,
        requested_bucket_id: usize,
    ) -> Result<(), SerialCollationError> {
        if self.open_writers < MAX_OPEN_MATERIALIZED_STITCH_WRITERS_PER_SHARD {
            return Ok(());
        }
        let evict_bucket_id = self
            .writers
            .iter()
            .enumerate()
            .find_map(|(bucket_id, writer)| {
                (bucket_id != requested_bucket_id
                    && writer
                        .as_ref()
                        .is_some_and(MaterializedStitchedCoordShardWriter::is_open))
                .then_some(bucket_id)
            })
            .unwrap_or(requested_bucket_id);
        if let Some(writer) = self.writers[evict_bucket_id].as_mut()
            && writer.is_open()
        {
            writer.flush_and_close()?;
            self.open_writers -= 1;
        }
        Ok(())
    }

    fn write_materialized_record(
        &mut self,
        bucket_id: usize,
        record: &StitchedCoordRecord,
        label: &[u8],
    ) -> Result<(), SerialCollationError> {
        self.ensure_writer(bucket_id)?;
        self.writers[bucket_id]
            .as_mut()
            .expect("bucket writer was just created")
            .write_record(record, label)
    }

    fn finish(self) -> Result<Vec<MaterializedStitchedCoordBucketEntry>, SerialCollationError> {
        let mut manifest = Vec::new();
        for writer in self.writers.into_iter().flatten() {
            manifest.push(writer.finish()?);
        }
        Ok(manifest)
    }
}

struct MaterializedStitchedCoordShardWriter {
    bucket_id: usize,
    coord_path: PathBuf,
    label_path: PathBuf,
    color_path: PathBuf,
    coord_out: Option<BufWriter<File>>,
    label_out: Option<BufWriter<File>>,
    color_out: Option<BufWriter<File>>,
    record_buffer: Vec<u8>,
    records: u64,
    label_bytes: u64,
    color_runs: u64,
}

#[inline]
fn unitig_colors_as_bytes(colors: &[UnitigColor]) -> &[u8] {
    // UnitigColor is repr(transparent) over u64; materialized buckets use the
    // same native little-endian private layout as their coordinate records.
    unsafe {
        std::slice::from_raw_parts(colors.as_ptr().cast::<u8>(), std::mem::size_of_val(colors))
    }
}

impl MaterializedStitchedCoordShardWriter {
    fn create(
        coord_dir: &Path,
        worker_id: usize,
        bucket_id: usize,
    ) -> Result<Self, SerialCollationError> {
        let coord_path = coord_dir.join(format!("{bucket_id:05}.{worker_id:03}.mcoord"));
        let label_path = coord_dir.join(format!("{bucket_id:05}.{worker_id:03}.mlabel"));
        let color_path = coord_dir.join(format!("{bucket_id:05}.{worker_id:03}.mcolor"));
        let coord_file = File::create(&coord_path).map_err(|source| SerialCollationError::Io {
            path: coord_path.clone(),
            source,
        })?;
        let label_file = File::create(&label_path).map_err(|source| SerialCollationError::Io {
            path: label_path.clone(),
            source,
        })?;
        let mut coord_out =
            BufWriter::with_capacity(MATERIALIZED_STITCH_COORD_SHARD_WRITE_BUFFER, coord_file);
        coord_out
            .write_all(MATERIALIZED_STITCH_COORD_MAGIC)
            .and_then(|_| coord_out.write_all(&(bucket_id as u64).to_le_bytes()))
            .and_then(|_| coord_out.write_all(&0u64.to_le_bytes()))
            .and_then(|_| coord_out.write_all(&0u64.to_le_bytes()))
            .map_err(|source| SerialCollationError::Io {
                path: coord_path.clone(),
                source,
            })?;
        Ok(Self {
            bucket_id,
            coord_path,
            label_path,
            color_path,
            coord_out: Some(coord_out),
            label_out: Some(BufWriter::with_capacity(1024 * 1024, label_file)),
            color_out: None,
            record_buffer: Vec::with_capacity(STITCH_COORD_RECORD_WRITE_BUFFER),
            records: 0,
            label_bytes: 0,
            color_runs: 0,
        })
    }

    fn is_open(&self) -> bool {
        self.coord_out.is_some() || self.label_out.is_some() || self.color_out.is_some()
    }

    fn ensure_open(&mut self) -> Result<(), SerialCollationError> {
        if self.coord_out.is_none() {
            let coord_file = OpenOptions::new()
                .append(true)
                .open(&self.coord_path)
                .map_err(|source| SerialCollationError::Io {
                    path: self.coord_path.clone(),
                    source,
                })?;
            self.coord_out = Some(BufWriter::with_capacity(
                MATERIALIZED_STITCH_COORD_SHARD_WRITE_BUFFER,
                coord_file,
            ));
        }
        if self.label_out.is_none() {
            let label_file = OpenOptions::new()
                .append(true)
                .open(&self.label_path)
                .map_err(|source| SerialCollationError::Io {
                    path: self.label_path.clone(),
                    source,
                })?;
            self.label_out = Some(BufWriter::with_capacity(1024 * 1024, label_file));
        }
        if self.color_runs != 0 && self.color_out.is_none() {
            let color_file = OpenOptions::new()
                .append(true)
                .open(&self.color_path)
                .map_err(|source| SerialCollationError::Io {
                    path: self.color_path.clone(),
                    source,
                })?;
            self.color_out = Some(BufWriter::with_capacity(1024 * 1024, color_file));
        }
        Ok(())
    }

    fn write_record(
        &mut self,
        record: &StitchedCoordRecord,
        label: &[u8],
    ) -> Result<(), SerialCollationError> {
        self.write_record_with_color_index(record, label, u32::MAX, 0)
    }

    fn write_record_with_color_index(
        &mut self,
        record: &StitchedCoordRecord,
        label: &[u8],
        color_index: u32,
        color_count: u32,
    ) -> Result<(), SerialCollationError> {
        let label_len = u32::try_from(label.len())
            .map_err(|_| SerialCollationError::MalformedCoordBucket(self.coord_path.clone()))?;
        self.record_buffer
            .extend_from_slice(&encoded_materialized_stitched_coord_record(
                MaterializedStitchedCoordRecord {
                    path_id: record.path_id,
                    rank: record.rank,
                    label_offset: self.label_bytes,
                    label_len,
                    reverse: record.reverse,
                    is_cycle: record.is_cycle,
                    color_index,
                    color_count,
                },
            ));
        if self.record_buffer.len() >= STITCH_COORD_RECORD_WRITE_BUFFER {
            self.flush_record_buffer()?;
        }
        self.label_out
            .as_mut()
            .expect("label writer is open")
            .write_all(label)
            .map_err(|source| SerialCollationError::Io {
                path: self.label_path.clone(),
                source,
            })?;
        self.records += 1;
        self.label_bytes += u64::from(label_len);
        Ok(())
    }

    fn write_colored_record(
        &mut self,
        record: &StitchedCoordRecord,
        label: &[u8],
        colors: &[UnitigColor],
    ) -> Result<(), SerialCollationError> {
        let count = u32::try_from(colors.len())
            .map_err(|_| SerialCollationError::MalformedCoordBucket(self.color_path.clone()))?;
        if self.color_out.is_none() {
            let file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(&self.color_path)
                .map_err(|source| SerialCollationError::Io {
                    path: self.color_path.clone(),
                    source,
                })?;
            self.color_out = Some(BufWriter::with_capacity(1024 * 1024, file));
        }
        let color_out = self
            .color_out
            .as_mut()
            .expect("color writer was just created");
        color_out
            .write_all(unitig_colors_as_bytes(colors))
            .map_err(|source| SerialCollationError::Io {
                path: self.color_path.clone(),
                source,
            })?;
        let color_index = u32::try_from(self.color_runs)
            .map_err(|_| SerialCollationError::MalformedCoordBucket(self.color_path.clone()))?;
        if color_index >= 0x3fff_ffff {
            return Err(SerialCollationError::MalformedCoordBucket(
                self.color_path.clone(),
            ));
        }
        self.color_runs += u64::from(count);
        self.write_record_with_color_index(record, label, color_index, count)
    }

    fn write_encoded_colored_record(
        &mut self,
        record: &StitchedCoordRecord,
        label: &[u8],
        count: u32,
        encoded_colors: &[u8],
    ) -> Result<(), SerialCollationError> {
        if encoded_colors.len() != count as usize * std::mem::size_of::<u64>() {
            return Err(SerialCollationError::MalformedCoordBucket(
                self.color_path.clone(),
            ));
        }
        if self.color_out.is_none() {
            let file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(&self.color_path)
                .map_err(|source| SerialCollationError::Io {
                    path: self.color_path.clone(),
                    source,
                })?;
            self.color_out = Some(BufWriter::with_capacity(1024 * 1024, file));
        }
        let color_out = self
            .color_out
            .as_mut()
            .expect("color writer was just created");
        color_out
            .write_all(encoded_colors)
            .map_err(|source| SerialCollationError::Io {
                path: self.color_path.clone(),
                source,
            })?;
        let color_index = u32::try_from(self.color_runs)
            .map_err(|_| SerialCollationError::MalformedCoordBucket(self.color_path.clone()))?;
        if color_index >= 0x3fff_ffff {
            return Err(SerialCollationError::MalformedCoordBucket(
                self.color_path.clone(),
            ));
        }
        self.color_runs += u64::from(count);
        self.write_record_with_color_index(record, label, color_index, count)
    }

    fn write_pending_batch(
        &mut self,
        batch: &mut PendingMaterializedBucket,
    ) -> Result<(), SerialCollationError> {
        let color_base = u32::try_from(self.color_runs)
            .map_err(|_| SerialCollationError::MalformedCoordBucket(self.color_path.clone()))?;
        for pending in &batch.records {
            let color_index = if pending.color_start == u32::MAX {
                u32::MAX
            } else {
                color_base.checked_add(pending.color_start).ok_or_else(|| {
                    SerialCollationError::MalformedCoordBucket(self.color_path.clone())
                })?
            };
            if color_index != u32::MAX && color_index >= 0x3fff_ffff {
                return Err(SerialCollationError::MalformedCoordBucket(
                    self.color_path.clone(),
                ));
            }
            self.record_buffer
                .extend_from_slice(&encoded_materialized_stitched_coord_record(
                    MaterializedStitchedCoordRecord {
                        path_id: pending.path_id,
                        rank: u64::from(pending.rank),
                        label_offset: self.label_bytes + u64::from(pending.label_offset),
                        label_len: u32::from(pending.label_len),
                        reverse: pending.reverse(),
                        is_cycle: pending.is_cycle(),
                        color_index,
                        color_count: pending.color_count(),
                    },
                ));
        }

        self.label_out
            .as_mut()
            .expect("label writer is open")
            .write_all(&batch.labels)
            .map_err(|source| SerialCollationError::Io {
                path: self.label_path.clone(),
                source,
            })?;
        if !batch.colors.is_empty() {
            if self.color_out.is_none() {
                let file = OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(&self.color_path)
                    .map_err(|source| SerialCollationError::Io {
                        path: self.color_path.clone(),
                        source,
                    })?;
                self.color_out = Some(BufWriter::with_capacity(1024 * 1024, file));
            }
            let color_out = self
                .color_out
                .as_mut()
                .expect("color writer was just created");
            color_out
                .write_all(unitig_colors_as_bytes(&batch.colors))
                .map_err(|source| SerialCollationError::Io {
                    path: self.color_path.clone(),
                    source,
                })?;
        }
        self.records += batch.records.len() as u64;
        self.label_bytes += batch.labels.len() as u64;
        self.color_runs += batch.colors.len() as u64;
        batch.records.clear();
        batch.labels.clear();
        batch.colors.clear();
        if self.record_buffer.len() >= STITCH_COORD_RECORD_WRITE_BUFFER {
            self.flush_record_buffer()?;
        }
        Ok(())
    }

    fn flush_record_buffer(&mut self) -> Result<(), SerialCollationError> {
        if self.record_buffer.is_empty() {
            return Ok(());
        }
        self.coord_out
            .as_mut()
            .expect("coord writer is open")
            .write_all(&self.record_buffer)
            .map_err(|source| SerialCollationError::Io {
                path: self.coord_path.clone(),
                source,
            })?;
        self.record_buffer.clear();
        Ok(())
    }

    fn flush_and_close(&mut self) -> Result<(), SerialCollationError> {
        self.flush_record_buffer()?;
        if let Some(mut coord_out) = self.coord_out.take() {
            coord_out
                .flush()
                .map_err(|source| SerialCollationError::Io {
                    path: self.coord_path.clone(),
                    source,
                })?;
        }
        if let Some(mut label_out) = self.label_out.take() {
            label_out
                .flush()
                .map_err(|source| SerialCollationError::Io {
                    path: self.label_path.clone(),
                    source,
                })?;
        }
        if let Some(mut color_out) = self.color_out.take() {
            color_out
                .flush()
                .map_err(|source| SerialCollationError::Io {
                    path: self.color_path.clone(),
                    source,
                })?;
        }
        Ok(())
    }

    fn finish(mut self) -> Result<MaterializedStitchedCoordBucketEntry, SerialCollationError> {
        self.flush_and_close()?;
        let mut coord_out = BufWriter::with_capacity(
            MATERIALIZED_STITCH_COORD_SHARD_WRITE_BUFFER,
            OpenOptions::new()
                .write(true)
                .open(&self.coord_path)
                .map_err(|source| SerialCollationError::Io {
                    path: self.coord_path.clone(),
                    source,
                })?,
        );
        coord_out
            .seek(SeekFrom::Start(16))
            .map_err(|source| SerialCollationError::Io {
                path: self.coord_path.clone(),
                source,
            })?;
        coord_out
            .write_all(&self.records.to_le_bytes())
            .and_then(|_| coord_out.write_all(&self.label_bytes.to_le_bytes()))
            .and_then(|_| coord_out.flush())
            .map_err(|source| SerialCollationError::Io {
                path: self.coord_path.clone(),
                source,
            })?;
        Ok(MaterializedStitchedCoordBucketEntry {
            bucket_id: self.bucket_id,
            records: self.records,
            label_bytes: self.label_bytes,
            coord_path: self.coord_path,
            label_path: self.label_path,
            color_path: (self.color_runs != 0).then_some(self.color_path),
            color_runs: self.color_runs,
        })
    }
}

struct StitchedCoordShardWriters<'a> {
    coord_dir: &'a Path,
    worker_id: usize,
    writers: Vec<Option<StitchedCoordShardWriter>>,
}

struct ConcurrentStitchedCoordWriters<'a> {
    coord_dir: &'a Path,
    writers: Vec<Mutex<Option<StitchedCoordShardWriter>>>,
    worker_buffers: Vec<UnsafeCell<Vec<Option<Vec<StitchedCoordRecord>>>>>,
}

// A Rayon worker exclusively accesses the buffer at its worker index. The final
// slot is used by the serial expansion work, and finish runs after the pool joins.
unsafe impl Sync for ConcurrentStitchedCoordWriters<'_> {}

impl<'a> ConcurrentStitchedCoordWriters<'a> {
    fn new(coord_dir: &'a Path, bucket_count: usize, workers: usize) -> Self {
        Self {
            coord_dir,
            writers: (0..bucket_count).map(|_| Mutex::new(None)).collect(),
            worker_buffers: (0..workers.max(1) + 1)
                .map(|_| {
                    let mut buckets = Vec::with_capacity(bucket_count);
                    buckets.resize_with(bucket_count, || None);
                    UnsafeCell::new(buckets)
                })
                .collect(),
        }
    }

    fn flush_records(
        &self,
        bucket_id: usize,
        records: &[StitchedCoordRecord],
    ) -> Result<(), SerialCollationError> {
        let mut writer = self.writers[bucket_id]
            .lock()
            .map_err(|_| SerialCollationError::WorkerPanic)?;
        if writer.is_none() {
            *writer = Some(StitchedCoordShardWriter::create(
                self.coord_dir,
                0,
                bucket_id,
            )?);
        }
        let writer = writer.as_mut().expect("bucket writer was just created");
        for &record in records {
            writer.write_record(record)?;
        }
        Ok(())
    }

    fn write_path_records(
        &self,
        bucket_id: usize,
        records: &[StitchedCoordRecord],
    ) -> Result<(), SerialCollationError> {
        if records.is_empty() {
            return Ok(());
        }
        let serial_slot = self.worker_buffers.len() - 1;
        let worker_id = rayon::current_thread_index()
            .filter(|&worker_id| worker_id < serial_slot)
            .unwrap_or(serial_slot);
        // SAFETY: each Rayon worker has a unique index in the expansion pool,
        // and non-pool calls are serial and use the dedicated final slot.
        let buckets = unsafe { &mut *self.worker_buffers[worker_id].get() };
        let capacity = EDGE_PATH_INFO_WORKER_BUFFER
            .div_ceil(std::mem::size_of::<StitchedCoordRecord>())
            .max(1);
        let buffer = buckets[bucket_id].get_or_insert_with(|| Vec::with_capacity(capacity));
        let mut remaining = records;
        while !remaining.is_empty() {
            let take = (capacity - buffer.len()).min(remaining.len());
            buffer.extend_from_slice(&remaining[..take]);
            remaining = &remaining[take..];
            if buffer.len() == capacity {
                self.flush_records(bucket_id, buffer)?;
                buffer.clear();
            }
        }
        Ok(())
    }

    fn write_record(
        &self,
        bucket_id: usize,
        record: StitchedCoordRecord,
    ) -> Result<(), SerialCollationError> {
        self.write_path_records(bucket_id, std::slice::from_ref(&record))
    }

    fn finish(
        self,
        pool: &ThreadPool,
    ) -> Result<Vec<StitchedCoordBucketEntry>, SerialCollationError> {
        let bucket_count = self.writers.len();
        let mut pending_by_bucket = (0..bucket_count)
            .map(|_| Vec::<Vec<StitchedCoordRecord>>::new())
            .collect::<Vec<_>>();
        for worker in self.worker_buffers {
            for (bucket_id, records) in worker.into_inner().into_iter().enumerate() {
                if let Some(records) = records
                    && !records.is_empty()
                {
                    pending_by_bucket[bucket_id].push(records);
                }
            }
        }
        let writers = self
            .writers
            .into_iter()
            .map(|writer| {
                writer
                    .into_inner()
                    .map_err(|_| SerialCollationError::WorkerPanic)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let coord_dir = self.coord_dir;
        let entries = pool.install(|| {
            writers
                .into_par_iter()
                .zip(pending_by_bucket.into_par_iter())
                .enumerate()
                .map(|(bucket_id, (writer, pending))| {
                    if writer.is_none() && pending.is_empty() {
                        return Ok(None);
                    }
                    let mut writer = match writer {
                        Some(writer) => writer,
                        None => StitchedCoordShardWriter::create(coord_dir, 0, bucket_id)?,
                    };
                    for records in pending {
                        for record in records {
                            writer.write_record(record)?;
                        }
                    }
                    writer.finish().map(Some)
                })
                .collect::<Result<Vec<_>, SerialCollationError>>()
        })?;
        Ok(entries.into_iter().flatten().collect())
    }
}

impl<'a> StitchedCoordShardWriters<'a> {
    fn new(coord_dir: &'a Path, worker_id: usize, bucket_count: usize) -> Self {
        let mut writers = Vec::with_capacity(bucket_count);
        writers.resize_with(bucket_count, || None);
        Self {
            coord_dir,
            worker_id,
            writers,
        }
    }

    fn write_path_records(
        &mut self,
        bucket_id: usize,
        records: &[StitchedCoordRecord],
    ) -> Result<(), SerialCollationError> {
        if records.is_empty() {
            return Ok(());
        }
        if self.writers[bucket_id].is_none() {
            self.writers[bucket_id] = Some(StitchedCoordShardWriter::create(
                self.coord_dir,
                self.worker_id,
                bucket_id,
            )?);
        }
        let writer = self.writers[bucket_id]
            .as_mut()
            .expect("bucket writer was just created");
        for &record in records {
            writer.write_record(record)?;
        }
        Ok(())
    }

    fn write_record(
        &mut self,
        bucket_id: usize,
        record: StitchedCoordRecord,
    ) -> Result<(), SerialCollationError> {
        if self.writers[bucket_id].is_none() {
            self.writers[bucket_id] = Some(StitchedCoordShardWriter::create(
                self.coord_dir,
                self.worker_id,
                bucket_id,
            )?);
        }
        self.writers[bucket_id]
            .as_mut()
            .expect("bucket writer was just created")
            .write_record(record)
    }

    fn finish(self) -> Result<Vec<StitchedCoordBucketEntry>, SerialCollationError> {
        let mut manifest = Vec::new();
        for writer in self.writers.into_iter().flatten() {
            manifest.push(writer.finish()?);
        }
        Ok(manifest)
    }
}

struct StitchedCoordShardWriter {
    bucket_id: usize,
    path: PathBuf,
    out: BufWriter<File>,
    record_buffer: Vec<u8>,
    records: u64,
}

impl StitchedCoordShardWriter {
    fn create(
        coord_dir: &Path,
        worker_id: usize,
        bucket_id: usize,
    ) -> Result<Self, SerialCollationError> {
        let path = coord_dir.join(format!("{bucket_id:05}.{worker_id:03}.scb"));
        let file = File::create(&path).map_err(|source| SerialCollationError::Io {
            path: path.clone(),
            source,
        })?;
        let mut out = BufWriter::with_capacity(STITCH_COORD_SHARD_WRITE_BUFFER, file);
        out.write_all(STITCH_COORD_MAGIC)
            .map_err(|source| SerialCollationError::Io {
                path: path.clone(),
                source,
            })?;
        out.write_all(&(bucket_id as u64).to_le_bytes())
            .map_err(|source| SerialCollationError::Io {
                path: path.clone(),
                source,
            })?;
        out.write_all(&0u64.to_le_bytes())
            .map_err(|source| SerialCollationError::Io {
                path: path.clone(),
                source,
            })?;
        out.write_all(&[0u8; 8])
            .map_err(|source| SerialCollationError::Io {
                path: path.clone(),
                source,
            })?;
        Ok(Self {
            bucket_id,
            path,
            out,
            record_buffer: Vec::with_capacity(STITCH_COORD_RECORD_WRITE_BUFFER),
            records: 0,
        })
    }

    fn write_record(&mut self, record: StitchedCoordRecord) -> Result<(), SerialCollationError> {
        self.record_buffer
            .extend_from_slice(&encoded_stitched_coord_record(record));
        self.records += 1;
        if self.record_buffer.len() >= STITCH_COORD_RECORD_WRITE_BUFFER {
            self.flush_record_buffer()?;
        }
        Ok(())
    }

    fn flush_record_buffer(&mut self) -> Result<(), SerialCollationError> {
        if self.record_buffer.is_empty() {
            return Ok(());
        }
        self.out
            .write_all(&self.record_buffer)
            .map_err(|source| SerialCollationError::Io {
                path: self.path.clone(),
                source,
            })?;
        self.record_buffer.clear();
        Ok(())
    }

    fn finish(mut self) -> Result<StitchedCoordBucketEntry, SerialCollationError> {
        self.flush_record_buffer()?;
        self.out
            .flush()
            .map_err(|source| SerialCollationError::Io {
                path: self.path.clone(),
                source,
            })?;
        self.out
            .seek(SeekFrom::Start(16))
            .map_err(|source| SerialCollationError::Io {
                path: self.path.clone(),
                source,
            })?;
        self.out
            .write_all(&self.records.to_le_bytes())
            .map_err(|source| SerialCollationError::Io {
                path: self.path.clone(),
                source,
            })?;
        self.out
            .flush()
            .map_err(|source| SerialCollationError::Io {
                path: self.path.clone(),
                source,
            })?;
        Ok(StitchedCoordBucketEntry {
            bucket_id: self.bucket_id,
            records: self.records,
            path: self.path,
        })
    }
}

fn stitched_coord_bucket_count(threads: usize) -> usize {
    (threads.max(1) * 1024)
        .next_power_of_two()
        .clamp(1024, 16_384)
}

fn stitch_endpoint_bucket_count(threads: usize) -> usize {
    let _ = threads;
    MAX_OPEN_STITCH_ENDPOINT_WRITERS
}

fn materialized_stitched_coord_bucket_count(threads: usize) -> usize {
    let target = (threads.max(1) * 16).next_power_of_two().clamp(32, 64);
    let max_open_buckets = MAX_OPEN_MATERIALIZED_STITCH_WRITERS
        .checked_div(threads.max(1))
        .unwrap_or(1)
        .max(1);
    target.min(floor_power_of_two(max_open_buckets)).max(1)
}

fn floor_power_of_two(value: usize) -> usize {
    if value == 0 {
        0
    } else {
        1usize << (usize::BITS - 1 - value.leading_zeros())
    }
}

fn final_unitig_bucket_count(threads: usize) -> usize {
    (threads.max(1) * 8)
        .next_power_of_two()
        .clamp(128, MAX_OPEN_STITCH_ENDPOINT_WRITERS)
}

fn stitched_coord_bucket(path_id: u64, bucket_mask: usize) -> usize {
    (xxh3_64(&path_id.to_ne_bytes()) as usize) & bucket_mask
}

fn encoded_stitched_coord_record(
    record: StitchedCoordRecord,
) -> [u8; STITCH_PATH_INFO_RECORD_LEN as usize] {
    let mut bytes = [0u8; STITCH_PATH_INFO_RECORD_LEN as usize];
    bytes[..8].copy_from_slice(&record.path_id.to_le_bytes());
    bytes[8..16].copy_from_slice(&record.rank.to_le_bytes());
    bytes[16..20].copy_from_slice(&record.unitig_index.to_le_bytes());
    let mut flags = 0u8;
    if record.reverse {
        flags |= STITCH_COORD_REVERSE_FLAG;
    }
    if record.is_cycle {
        flags |= STITCH_COORD_CYCLE_FLAG;
    }
    bytes[20] = flags;
    bytes
}

fn encoded_materialized_stitched_coord_record(
    record: MaterializedStitchedCoordRecord,
) -> [u8; STITCH_COORD_RECORD_LEN as usize] {
    let mut bytes = [0u8; STITCH_COORD_RECORD_LEN as usize];
    bytes[..8].copy_from_slice(&record.path_id.to_le_bytes());
    debug_assert!(u32::try_from(record.label_offset).is_ok());
    bytes[8..12].copy_from_slice(&(record.label_offset as u32).to_le_bytes());
    bytes[12..16].copy_from_slice(&record.color_index.to_le_bytes());
    bytes[16..18].copy_from_slice(
        &u16::try_from(record.rank)
            .expect("materialized rank fits C++ weight_t")
            .to_le_bytes(),
    );
    bytes[18..20].copy_from_slice(
        &u16::try_from(record.label_len)
            .expect("materialized label length fits C++ uni_len_t")
            .to_le_bytes(),
    );
    bytes[20..22].copy_from_slice(
        &u16::try_from(record.color_count)
            .expect("materialized color count fits C++ uni_len_t")
            .to_le_bytes(),
    );
    let mut flags = 0u16;
    if record.reverse {
        flags |= LoadedMaterializedStitchedCoordRecord::REVERSE_FLAG;
    }
    if record.is_cycle {
        flags |= LoadedMaterializedStitchedCoordRecord::CYCLE_FLAG;
    }
    bytes[22..24].copy_from_slice(&flags.to_le_bytes());
    bytes
}

fn reduce_materialized_stitched_coord_bucket_files_to_final<const K: usize>(
    manifest: &[MaterializedStitchedCoordBucketEntry],
    retained: &[Vec<PendingMaterializedBucket>],
    threads: usize,
    final_buckets: &mut FinalUnitigBucketWriters,
) -> Result<u64, SerialCollationError> {
    if manifest.is_empty() && retained.iter().all(Vec::is_empty) {
        return Ok(0);
    }

    let bucket_count = manifest
        .iter()
        .map(|entry| entry.bucket_id + 1)
        .max()
        .unwrap_or(0)
        .max(retained.len());
    let mut files_by_bucket = vec![Vec::new(); bucket_count];
    for entry in manifest {
        files_by_bucket[entry.bucket_id].push(entry.clone());
    }
    let groups = files_by_bucket
        .into_iter()
        .enumerate()
        .filter(|(bucket_id, files)| {
            !files.is_empty()
                || retained
                    .get(*bucket_id)
                    .is_some_and(|tails| !tails.is_empty())
        })
        .collect::<Vec<_>>();

    let workers = threads.max(1).min(groups.len());
    if final_buckets.direct_output.is_some() && workers > 1 {
        return reduce_materialized_groups_to_direct_fasta::<K>(
            &groups,
            retained,
            workers,
            final_buckets,
        );
    }
    let mut emitted = 0u64;
    if workers == 1 {
        for (bucket_id, files) in &groups {
            let unitigs = reduce_materialized_stitched_coord_bucket_file_group_with_tails::<K>(
                files,
                retained.get(*bucket_id).map_or(&[], Vec::as_slice),
            )?;
            for unitig in unitigs {
                if unitig.colors.is_empty() {
                    final_buckets.write_label(&unitig.label)?;
                } else {
                    final_buckets.write_colored_label(&unitig.label, &unitig.colors)?;
                }
                emitted += 1;
            }
        }
        return Ok(emitted);
    }

    let next_group = AtomicUsize::new(0);
    let (tx, rx) =
        mpsc::sync_channel::<Result<Vec<FinalUnitigRecord>, SerialCollationError>>(workers * 2);
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for _ in 0..workers {
            let tx = tx.clone();
            let next_group = &next_group;
            let groups = &groups;
            handles.push(
                scope.spawn(move || {
                    loop {
                        let group_idx = next_group.fetch_add(1, Ordering::Relaxed);
                        let Some(group) = groups.get(group_idx) else {
                            break;
                        };
                        let result =
                            reduce_materialized_stitched_coord_bucket_file_group_with_tails::<K>(
                                &group.1,
                                retained.get(group.0).map_or(&[], Vec::as_slice),
                            );
                        let should_stop = result.is_err();
                        if tx.send(result).is_err() || should_stop {
                            break;
                        }
                    }
                }),
            );
        }
        drop(tx);

        let mut first_error = None;
        for result in rx {
            match result {
                Ok(unitigs) if first_error.is_none() => {
                    for unitig in unitigs {
                        let result = if unitig.colors.is_empty() {
                            final_buckets.write_label(&unitig.label)
                        } else {
                            final_buckets.write_colored_label(&unitig.label, &unitig.colors)
                        };
                        if let Err(err) = result {
                            first_error = Some(err);
                            break;
                        }
                        emitted += 1;
                    }
                }
                Ok(_) => {}
                Err(err) if first_error.is_none() => first_error = Some(err),
                Err(_) => {}
            }
        }

        for handle in handles {
            handle
                .join()
                .map_err(|_| SerialCollationError::WorkerPanic)?;
        }

        if let Some(err) = first_error {
            Err(err)
        } else {
            Ok(())
        }
    })?;

    Ok(emitted)
}

struct EncodedFinalBatch {
    bytes: Vec<u8>,
    records: u64,
    bases: u64,
}

struct DirectFinalBatchBuilder {
    bytes: Vec<u8>,
    records: u64,
    bases: u64,
}

impl DirectFinalBatchBuilder {
    fn new(_first_record: u64) -> Self {
        Self {
            bytes: Vec::with_capacity(512 * 1024),
            records: 0,
            bases: 0,
        }
    }

    fn push(&mut self, label: &[u8], colors: &[UnitigColor]) {
        self.bytes.extend_from_slice(b">0");
        for color in colors {
            self.bytes.push(b' ');
            append_decimal_u64(&mut self.bytes, color.raw());
        }
        self.bytes.push(b'\n');
        self.bytes.extend_from_slice(label);
        self.bytes.push(b'\n');
        self.records += 1;
        self.bases += label.len() as u64;
    }

    fn finish(self) -> EncodedFinalBatch {
        EncodedFinalBatch {
            bytes: self.bytes,
            records: self.records,
            bases: self.bases,
        }
    }
}

fn reduce_materialized_groups_to_direct_fasta<const K: usize>(
    groups: &[(usize, Vec<MaterializedStitchedCoordBucketEntry>)],
    retained: &[Vec<PendingMaterializedBucket>],
    workers: usize,
    final_buckets: &mut FinalUnitigBucketWriters,
) -> Result<u64, SerialCollationError> {
    const DIRECT_FINAL_BATCH_RECORDS: usize = 16 * 1024;
    let next_group = AtomicUsize::new(0);
    let next_record = AtomicU64::new(final_buckets.direct_record_id_highwater + 1);
    let (output, output_path, initial_offset) = final_buckets.prepare_parallel_direct_output()?;
    let next_offset = AtomicU64::new(initial_offset);
    let emitted = AtomicU64::new(0);
    let emitted_bases = AtomicU64::new(0);
    let load_ns = AtomicU64::new(0);
    let sort_ns = AtomicU64::new(0);
    let assemble_ns = AtomicU64::new(0);
    let send_ns = AtomicU64::new(0);
    let write_started = Instant::now();
    let worker_timings = std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(workers);
        for _ in 0..workers {
            let next_group = &next_group;
            let next_record = &next_record;
            let next_offset = &next_offset;
            let emitted = &emitted;
            let emitted_bases = &emitted_bases;
            let output = &output;
            let output_path = &output_path;
            let load_ns = &load_ns;
            let sort_ns = &sort_ns;
            let assemble_ns = &assemble_ns;
            let send_ns = &send_ns;
            handles.push(scope.spawn(move || {
                let worker_started = Instant::now();
                let mut groups_processed = 0u64;
                let mut max_group_ns = 0u64;
                let mut max_load_ns = 0u64;
                let mut max_sort_ns = 0u64;
                let mut max_assemble_ns = 0u64;
                let mut max_send_ns = 0u64;
                let mut shard = MaterializedStitchedCoordBucket::default();
                loop {
                    let group_idx = next_group.fetch_add(1, Ordering::Relaxed);
                    let Some(group) = groups.get(group_idx) else {
                        break;
                    };
                    let group_started = Instant::now();
                    let started = Instant::now();
                    load_materialized_stitched_coord_bucket_file_group_with_tails_into(
                        &group.1,
                        retained.get(group.0).map_or(&[], Vec::as_slice),
                        &mut shard,
                    )?;
                    let elapsed_ns = started.elapsed().as_nanos() as u64;
                    load_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                    max_load_ns = max_load_ns.max(elapsed_ns);
                    let started = Instant::now();
                    shard
                        .records
                        .sort_unstable_by_key(|record| (record.path_id, record.rank));
                    let elapsed_ns = started.elapsed().as_nanos() as u64;
                    sort_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                    max_sort_ns = max_sort_ns.max(elapsed_ns);
                    let mut batch = None::<DirectFinalBatchBuilder>;
                    let mut write_error = None;
                    let started = Instant::now();
                    reduce_sorted_materialized_stitched_coord_bucket_with::<K, _>(
                        &mut shard.records,
                        &shard.labels,
                        &shard.colors,
                        |label, colors| {
                            if write_error.is_some() {
                                return;
                            }
                            let builder = batch.get_or_insert_with(|| {
                                let first_record = next_record.fetch_add(
                                    DIRECT_FINAL_BATCH_RECORDS as u64,
                                    Ordering::Relaxed,
                                );
                                DirectFinalBatchBuilder::new(first_record)
                            });
                            builder.push(label, colors);
                            if builder.records as usize >= DIRECT_FINAL_BATCH_RECORDS {
                                let write_started = Instant::now();
                                let encoded = batch.take().expect("batch exists").finish();
                                let result = write_encoded_final_batch_at(
                                    output,
                                    output_path,
                                    next_offset,
                                    emitted,
                                    emitted_bases,
                                    encoded,
                                );
                                let elapsed_ns = write_started.elapsed().as_nanos() as u64;
                                send_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                                max_send_ns = max_send_ns.max(elapsed_ns);
                                if let Err(error) = result {
                                    write_error = Some(error);
                                }
                            }
                        },
                    );
                    let elapsed_ns = started.elapsed().as_nanos() as u64;
                    assemble_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                    max_assemble_ns = max_assemble_ns.max(elapsed_ns);
                    if let Some(error) = write_error {
                        return Err(error);
                    }
                    if let Some(batch) = batch {
                        let started = Instant::now();
                        write_encoded_final_batch_at(
                            output,
                            output_path,
                            next_offset,
                            emitted,
                            emitted_bases,
                            batch.finish(),
                        )?;
                        let elapsed_ns = started.elapsed().as_nanos() as u64;
                        send_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                        max_send_ns = max_send_ns.max(elapsed_ns);
                    }
                    groups_processed += 1;
                    max_group_ns = max_group_ns.max(group_started.elapsed().as_nanos() as u64);
                }
                Ok::<_, SerialCollationError>(ReducerWorkerTiming {
                    elapsed_ns: worker_started.elapsed().as_nanos() as u64,
                    groups: groups_processed,
                    max_group_ns,
                    max_load_ns,
                    max_sort_ns,
                    max_assemble_ns,
                    max_send_ns,
                })
            }));
        }
        let mut first_error = None;
        let mut timings = Vec::with_capacity(handles.len());
        for handle in handles {
            let result = handle
                .join()
                .map_err(|_| SerialCollationError::WorkerPanic)?;
            match result {
                Ok(timing) => timings.push(timing),
                Err(err) if first_error.is_none() => first_error = Some(err),
                Err(_) => {}
            }
        }
        if let Some(err) = first_error {
            Err(err)
        } else {
            Ok(timings)
        }
    })?;
    final_buckets.direct_record_id_highwater =
        next_record.load(Ordering::Relaxed).saturating_sub(1);
    let emitted = emitted.load(Ordering::Relaxed);
    final_buckets.direct_records += emitted;
    final_buckets.direct_bases += emitted_bases.load(Ordering::Relaxed);
    eprintln!(
        "cuttlefish3-rs: path-info reduce worker detail: load {:.3}s, sort {:.3}s, assemble/encode {:.3}s, blocked send {:.3}s; reducer wall/write-drain {:.3}s",
        load_ns.load(Ordering::Relaxed) as f64 / 1e9,
        sort_ns.load(Ordering::Relaxed) as f64 / 1e9,
        assemble_ns.load(Ordering::Relaxed) as f64 / 1e9,
        send_ns.load(Ordering::Relaxed) as f64 / 1e9,
        write_started.elapsed().as_secs_f64(),
    );
    if let Some(slowest) = worker_timings.iter().max_by_key(|timing| timing.elapsed_ns) {
        eprintln!(
            "cuttlefish3-rs: path-info reduce critical worker: elapsed {:.3}s, groups {}; largest group {:.3}s (load {:.3}s, sort {:.3}s, assemble/write {:.3}s, write {:.3}s)",
            slowest.elapsed_ns as f64 / 1e9,
            slowest.groups,
            slowest.max_group_ns as f64 / 1e9,
            slowest.max_load_ns as f64 / 1e9,
            slowest.max_sort_ns as f64 / 1e9,
            slowest.max_assemble_ns as f64 / 1e9,
            slowest.max_send_ns as f64 / 1e9,
        );
    }
    Ok(emitted)
}

struct ReducerWorkerTiming {
    elapsed_ns: u64,
    groups: u64,
    max_group_ns: u64,
    max_load_ns: u64,
    max_sort_ns: u64,
    max_assemble_ns: u64,
    max_send_ns: u64,
}

fn write_encoded_final_batch_at(
    output: &File,
    output_path: &Path,
    next_offset: &AtomicU64,
    emitted: &AtomicU64,
    emitted_bases: &AtomicU64,
    batch: EncodedFinalBatch,
) -> Result<(), SerialCollationError> {
    let offset = next_offset.fetch_add(batch.bytes.len() as u64, Ordering::Relaxed);
    output
        .write_all_at(&batch.bytes, offset)
        .map_err(|source| SerialCollationError::Io {
            path: output_path.to_path_buf(),
            source,
        })?;
    emitted.fetch_add(batch.records, Ordering::Relaxed);
    emitted_bases.fetch_add(batch.bases, Ordering::Relaxed);
    Ok(())
}

fn encode_final_unitig_batch(
    unitigs: Vec<FinalUnitigRecord>,
    _first_record: u64,
) -> EncodedFinalBatch {
    let records = unitigs.len() as u64;
    let bases = unitigs.iter().map(|unitig| unitig.label.len() as u64).sum();
    let mut bytes = Vec::with_capacity(bases as usize + unitigs.len() * 32);
    for unitig in unitigs {
        bytes.extend_from_slice(b">0");
        for color in unitig.colors {
            bytes.push(b' ');
            append_decimal_u64(&mut bytes, color.raw());
        }
        bytes.push(b'\n');
        bytes.extend_from_slice(&unitig.label);
        bytes.push(b'\n');
    }
    EncodedFinalBatch {
        bytes,
        records,
        bases,
    }
}

fn reduce_materialized_stitched_coord_bucket_file_group<const K: usize>(
    group: &[MaterializedStitchedCoordBucketEntry],
) -> Result<Vec<FinalUnitigRecord>, SerialCollationError> {
    let mut shard = load_materialized_stitched_coord_bucket_file_group(group)?;
    Ok(reduce_materialized_stitched_coord_bucket::<K>(
        &mut shard.records,
        &shard.labels,
        &shard.colors,
    ))
}

fn reduce_materialized_stitched_coord_bucket_file_group_with_tails<const K: usize>(
    group: &[MaterializedStitchedCoordBucketEntry],
    tails: &[PendingMaterializedBucket],
) -> Result<Vec<FinalUnitigRecord>, SerialCollationError> {
    let mut shard = load_materialized_stitched_coord_bucket_file_group_with_tails(group, tails)?;
    Ok(reduce_materialized_stitched_coord_bucket::<K>(
        &mut shard.records,
        &shard.labels,
        &shard.colors,
    ))
}

fn load_materialized_stitched_coord_bucket_file_group_with_tails(
    group: &[MaterializedStitchedCoordBucketEntry],
    tails: &[PendingMaterializedBucket],
) -> Result<MaterializedStitchedCoordBucket, SerialCollationError> {
    let mut shard = if group.is_empty() {
        MaterializedStitchedCoordBucket {
            records: Vec::new(),
            labels: Vec::new(),
            colors: Vec::new(),
        }
    } else {
        load_materialized_stitched_coord_bucket_file_group(group)?
    };
    for tail in tails {
        append_pending_materialized_bucket(&mut shard, tail)?;
    }
    Ok(shard)
}

fn load_materialized_stitched_coord_bucket_file_group_with_tails_into(
    group: &[MaterializedStitchedCoordBucketEntry],
    tails: &[PendingMaterializedBucket],
    shard: &mut MaterializedStitchedCoordBucket,
) -> Result<(), SerialCollationError> {
    shard.records.clear();
    shard.labels.clear();
    shard.colors.clear();
    let total_records = group.iter().map(|entry| entry.records).sum::<u64>()
        + tails
            .iter()
            .map(|tail| tail.records.len() as u64)
            .sum::<u64>();
    let total_label_bytes = group.iter().map(|entry| entry.label_bytes).sum::<u64>()
        + tails
            .iter()
            .map(|tail| tail.labels.len() as u64)
            .sum::<u64>();
    let total_colors = group.iter().map(|entry| entry.color_runs).sum::<u64>()
        + tails
            .iter()
            .map(|tail| tail.colors.len() as u64)
            .sum::<u64>();
    shard.records.reserve(total_records as usize);
    shard.labels.reserve(total_label_bytes as usize);
    shard.colors.reserve(total_colors as usize);
    for entry in group {
        append_materialized_stitched_coord_bucket_file(entry, shard)?;
        remove_serial_file(&entry.coord_path)?;
        remove_serial_file(&entry.label_path)?;
        if let Some(path) = entry.color_path.as_ref() {
            remove_serial_file(path)?;
        }
    }
    for tail in tails {
        append_pending_materialized_bucket(shard, tail)?;
    }
    Ok(())
}

fn append_pending_materialized_bucket(
    shard: &mut MaterializedStitchedCoordBucket,
    tail: &PendingMaterializedBucket,
) -> Result<(), SerialCollationError> {
    let label_base = u32::try_from(shard.labels.len()).map_err(|_| {
        SerialCollationError::MalformedCoordBucket(PathBuf::from("retained-materialized-bucket"))
    })?;
    let color_base = u32::try_from(shard.colors.len()).map_err(|_| {
        SerialCollationError::MalformedCoordBucket(PathBuf::from("retained-materialized-bucket"))
    })?;
    shard.records.reserve(tail.records.len());
    shard.labels.extend_from_slice(&tail.labels);
    for pending in &tail.records {
        let color_start = if pending.color_start == u32::MAX {
            u32::MAX
        } else {
            color_base.checked_add(pending.color_start).ok_or_else(|| {
                SerialCollationError::MalformedCoordBucket(PathBuf::from(
                    "retained-materialized-bucket",
                ))
            })?
        };
        let label_offset = label_base
            .checked_add(pending.label_offset)
            .ok_or_else(|| {
                SerialCollationError::MalformedCoordBucket(PathBuf::from(
                    "retained-materialized-bucket",
                ))
            })?;
        shard
            .records
            .push(LoadedMaterializedStitchedCoordRecord::new(
                pending.path_id,
                u64::from(pending.rank),
                label_offset,
                u32::from(pending.label_len),
                pending.reverse(),
                pending.is_cycle(),
                color_start,
                pending.color_count(),
            ));
    }
    shard.colors.extend_from_slice(&tail.colors);
    Ok(())
}

fn load_materialized_stitched_coord_bucket_file_group(
    group: &[MaterializedStitchedCoordBucketEntry],
) -> Result<MaterializedStitchedCoordBucket, SerialCollationError> {
    if let [entry] = group {
        let shard = read_materialized_stitched_coord_bucket_file(entry)?;
        remove_serial_file(&entry.coord_path)?;
        remove_serial_file(&entry.label_path)?;
        if let Some(path) = entry.color_path.as_ref() {
            remove_serial_file(path)?;
        }
        return Ok(shard);
    }
    let total_records = group.iter().map(|entry| entry.records).sum::<u64>();
    let total_label_bytes = group.iter().map(|entry| entry.label_bytes).sum::<u64>();
    let mut records = Vec::with_capacity(total_records as usize);
    let mut labels = Vec::with_capacity(total_label_bytes as usize);
    let mut colors = Vec::new();
    for entry in group {
        let label_offset = u32::try_from(labels.len())
            .map_err(|_| SerialCollationError::MalformedCoordBucket(entry.label_path.clone()))?;
        let mut shard = read_materialized_stitched_coord_bucket_file(entry)?;
        remove_serial_file(&entry.coord_path)?;
        remove_serial_file(&entry.label_path)?;
        if let Some(path) = entry.color_path.as_ref() {
            remove_serial_file(path)?;
        }
        let color_offset = u32::try_from(colors.len())
            .map_err(|_| SerialCollationError::MalformedCoordBucket(entry.coord_path.clone()))?;
        for record in &mut shard.records {
            record.label_offset =
                record
                    .label_offset
                    .checked_add(label_offset)
                    .ok_or_else(|| {
                        SerialCollationError::MalformedCoordBucket(entry.label_path.clone())
                    })?;
            if record.color_start != u32::MAX {
                record.color_start =
                    record
                        .color_start
                        .checked_add(color_offset)
                        .ok_or_else(|| {
                            SerialCollationError::MalformedCoordBucket(entry.coord_path.clone())
                        })?;
            }
        }
        records.extend(shard.records);
        labels.extend(shard.labels);
        colors.extend(shard.colors);
    }
    Ok(MaterializedStitchedCoordBucket {
        records,
        labels,
        colors,
    })
}

#[derive(Default)]
struct MaterializedStitchedCoordBucket {
    records: Vec<LoadedMaterializedStitchedCoordRecord>,
    labels: Vec<u8>,
    colors: Vec<UnitigColor>,
}

struct FinalUnitigRecord {
    label: Vec<u8>,
    colors: Vec<UnitigColor>,
}

fn read_materialized_stitched_coord_bucket_file(
    entry: &MaterializedStitchedCoordBucketEntry,
) -> Result<MaterializedStitchedCoordBucket, SerialCollationError> {
    let mut bucket = MaterializedStitchedCoordBucket::default();
    append_materialized_stitched_coord_bucket_file(entry, &mut bucket)?;
    Ok(bucket)
}

fn append_materialized_stitched_coord_bucket_file(
    entry: &MaterializedStitchedCoordBucketEntry,
    bucket: &mut MaterializedStitchedCoordBucket,
) -> Result<(), SerialCollationError> {
    let mut file = File::open(&entry.coord_path).map_err(|source| SerialCollationError::Io {
        path: entry.coord_path.clone(),
        source,
    })?;
    let actual_len = file
        .metadata()
        .map_err(|source| SerialCollationError::Io {
            path: entry.coord_path.clone(),
            source,
        })?
        .len();
    let mut header = [0u8; STITCH_COORD_HEADER_LEN as usize];
    file.read_exact(&mut header)
        .map_err(|source| SerialCollationError::Io {
            path: entry.coord_path.clone(),
            source,
        })?;
    if &header[..8] != MATERIALIZED_STITCH_COORD_MAGIC {
        return Err(SerialCollationError::MalformedCoordBucket(
            entry.coord_path.clone(),
        ));
    }
    let bucket_id = u64::from_le_bytes(header[8..16].try_into().expect("bucket ID")) as usize;
    let records = u64::from_le_bytes(header[16..24].try_into().expect("record count"));
    let label_bytes = u64::from_le_bytes(header[24..32].try_into().expect("label bytes"));
    if bucket_id != entry.bucket_id || records != entry.records || label_bytes != entry.label_bytes
    {
        return Err(SerialCollationError::MalformedCoordBucket(
            entry.coord_path.clone(),
        ));
    }
    let expected_len = STITCH_COORD_HEADER_LEN
        .checked_add(
            records
                .checked_mul(STITCH_COORD_RECORD_LEN)
                .ok_or_else(|| {
                    SerialCollationError::MalformedCoordBucket(entry.coord_path.clone())
                })?,
        )
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(entry.coord_path.clone()))?;
    if actual_len != expected_len {
        return Err(SerialCollationError::MalformedCoordBucket(
            entry.coord_path.clone(),
        ));
    }
    let record_base = bucket.records.len();
    let label_base = u32::try_from(bucket.labels.len())
        .map_err(|_| SerialCollationError::MalformedCoordBucket(entry.label_path.clone()))?;
    bucket.records.reserve(records as usize);
    // The private bucket format is the native little-endian in-memory layout.
    let coord_bytes = unsafe {
        std::slice::from_raw_parts_mut(
            bucket
                .records
                .spare_capacity_mut()
                .as_mut_ptr()
                .cast::<u8>(),
            records as usize * std::mem::size_of::<LoadedMaterializedStitchedCoordRecord>(),
        )
    };
    file.read_exact(coord_bytes)
        .map_err(|source| SerialCollationError::Io {
            path: entry.coord_path.clone(),
            source,
        })?;
    // SAFETY: read_exact initialized every byte of each native POD record.
    unsafe { bucket.records.set_len(record_base + records as usize) };

    let mut label_file =
        File::open(&entry.label_path).map_err(|source| SerialCollationError::Io {
            path: entry.label_path.clone(),
            source,
        })?;
    let actual_label_len = label_file
        .metadata()
        .map_err(|source| SerialCollationError::Io {
            path: entry.label_path.clone(),
            source,
        })?
        .len();
    if actual_label_len != entry.label_bytes {
        return Err(SerialCollationError::MalformedCoordBucket(
            entry.label_path.clone(),
        ));
    }
    let label_record_base = bucket.labels.len();
    bucket.labels.reserve(entry.label_bytes as usize);
    let uninitialized_labels = unsafe {
        std::slice::from_raw_parts_mut(
            bucket.labels.spare_capacity_mut().as_mut_ptr().cast::<u8>(),
            entry.label_bytes as usize,
        )
    };
    label_file
        .read_exact(uninitialized_labels)
        .map_err(|source| SerialCollationError::Io {
            path: entry.label_path.clone(),
            source,
        })?;
    // SAFETY: read_exact initialized the entire requested label buffer.
    unsafe {
        bucket
            .labels
            .set_len(label_record_base + entry.label_bytes as usize)
    };

    let color_base = u32::try_from(bucket.colors.len())
        .map_err(|_| SerialCollationError::MalformedCoordBucket(entry.coord_path.clone()))?;
    bucket.colors.reserve(entry.color_runs as usize);
    if let Some(color_path) = &entry.color_path {
        let mut color_file = File::open(color_path).map_err(|source| SerialCollationError::Io {
            path: color_path.clone(),
            source,
        })?;
        let expected_color_bytes = entry.color_runs * std::mem::size_of::<UnitigColor>() as u64;
        if color_file
            .metadata()
            .map_err(|source| SerialCollationError::Io {
                path: color_path.clone(),
                source,
            })?
            .len()
            != expected_color_bytes
        {
            return Err(SerialCollationError::MalformedCoordBucket(
                color_path.clone(),
            ));
        }
        let color_record_base = bucket.colors.len();
        let uninitialized_colors = unsafe {
            std::slice::from_raw_parts_mut(
                bucket.colors.spare_capacity_mut().as_mut_ptr().cast::<u8>(),
                expected_color_bytes as usize,
            )
        };
        color_file
            .read_exact(uninitialized_colors)
            .map_err(|source| SerialCollationError::Io {
                path: color_path.clone(),
                source,
            })?;
        // SAFETY: UnitigColor is transparent over u64 and read_exact
        // initialized every byte of each native private-format record.
        unsafe {
            bucket
                .colors
                .set_len(color_record_base + entry.color_runs as usize)
        };
    }

    for record in &mut bucket.records[record_base..] {
        record.label_offset = record
            .label_offset
            .checked_add(label_base)
            .ok_or_else(|| SerialCollationError::MalformedCoordBucket(entry.label_path.clone()))?;
        if record.color_start != u32::MAX {
            let color_end = u64::from(record.color_start) + u64::from(record.color_count());
            if color_end > entry.color_runs {
                return Err(SerialCollationError::MalformedCoordBucket(
                    entry.coord_path.clone(),
                ));
            }
            record.color_start = color_base.checked_add(record.color_start).ok_or_else(|| {
                SerialCollationError::MalformedCoordBucket(entry.coord_path.clone())
            })?;
        }
    }

    Ok(())
}

fn decoded_materialized_stitched_coord_record(
    bytes: &[u8],
    path: &Path,
) -> Result<MaterializedStitchedCoordRecord, SerialCollationError> {
    if bytes.len() != STITCH_COORD_RECORD_LEN as usize {
        return Err(SerialCollationError::MalformedCoordBucket(
            path.to_path_buf(),
        ));
    }
    let path_id = u64::from_le_bytes(bytes[..8].try_into().expect("u64 path_id field"));
    let label_offset = u32::from_le_bytes(bytes[8..12].try_into().expect("u32 label offset field"));
    let color_index = u32::from_le_bytes(bytes[12..16].try_into().expect("u32 color index field"));
    let rank = u16::from_le_bytes(bytes[16..18].try_into().expect("u16 rank field"));
    let label_len = u16::from_le_bytes(bytes[18..20].try_into().expect("u16 label length field"));
    let color_count = u16::from_le_bytes(bytes[20..22].try_into().expect("u16 color count field"));
    let flags = u16::from_le_bytes(bytes[22..24].try_into().expect("u16 flags field"));
    Ok(MaterializedStitchedCoordRecord {
        path_id,
        rank: u64::from(rank),
        label_offset: u64::from(label_offset),
        label_len: u32::from(label_len),
        reverse: flags & LoadedMaterializedStitchedCoordRecord::REVERSE_FLAG != 0,
        is_cycle: flags & LoadedMaterializedStitchedCoordRecord::CYCLE_FLAG != 0,
        color_index,
        color_count: u32::from(color_count),
    })
}

fn reduce_materialized_stitched_coord_bucket<const K: usize>(
    records: &mut [LoadedMaterializedStitchedCoordRecord],
    labels: &[u8],
    color_runs: &[UnitigColor],
) -> Vec<FinalUnitigRecord> {
    let mut unitigs = Vec::new();
    reduce_materialized_stitched_coord_bucket_with::<K, _>(
        records,
        labels,
        color_runs,
        |label, colors| {
            unitigs.push(FinalUnitigRecord {
                label: label.to_vec(),
                colors: colors.to_vec(),
            });
        },
    );
    unitigs
}

fn reduce_materialized_stitched_coord_bucket_with<const K: usize, F>(
    records: &mut [LoadedMaterializedStitchedCoordRecord],
    labels: &[u8],
    color_runs: &[UnitigColor],
    emit: F,
) where
    F: FnMut(&[u8], &[UnitigColor]),
{
    if records.is_empty() {
        return;
    }

    if !materialized_stitched_coord_records_are_ordered(records) {
        records.sort_by_key(|record| (record.path_id, record.rank));
    }
    reduce_sorted_materialized_stitched_coord_bucket_with::<K, _>(
        records, labels, color_runs, emit,
    );
}

fn reduce_sorted_materialized_stitched_coord_bucket_with<const K: usize, F>(
    records: &mut [LoadedMaterializedStitchedCoordRecord],
    labels: &[u8],
    color_runs: &[UnitigColor],
    mut emit: F,
) where
    F: FnMut(&[u8], &[UnitigColor]),
{
    if records.is_empty() {
        return;
    }
    let mut label = Vec::new();
    let mut colors = Vec::new();
    let mut start = 0;
    while start < records.len() {
        let path_id = records[start].path_id;
        let is_cycle = records[start].is_cycle();
        let mut end = start + 1;
        while end < records.len() && records[end].path_id == path_id {
            end += 1;
        }

        label.clear();
        colors.clear();
        if end - start == 2 && !is_cycle && records[start].rank == 0 && records[start + 1].rank == 0
        {
            for (idx, record) in records[start..end].iter().enumerate() {
                let label_start = record.label_offset as usize;
                let label_end = label_start + record.label_len as usize;
                let unitig_label = &labels[label_start..label_end];
                let reverse = if idx == 0 {
                    record.reverse()
                } else {
                    !record.reverse()
                };
                append_materialized_colors::<K>(
                    &mut colors,
                    label.len(),
                    record,
                    color_runs,
                    reverse,
                );
                append_or_init_oriented_fast::<K>(&mut label, unitig_label, reverse);
            }
        } else {
            for record in &records[start..end] {
                let label_start = record.label_offset as usize;
                let label_end = label_start + record.label_len as usize;
                let unitig_label = &labels[label_start..label_end];
                append_materialized_colors::<K>(
                    &mut colors,
                    label.len(),
                    record,
                    color_runs,
                    record.reverse(),
                );
                append_or_init_oriented_fast::<K>(&mut label, unitig_label, record.reverse());
            }
        }

        if label.len() >= K {
            let colored = !colors.is_empty();
            if is_cycle && colored {
                label.pop();
                let vertex_count = label.len().saturating_sub(K - 1) as u32;
                while colors
                    .last()
                    .is_some_and(|run| run.offset() >= vertex_count)
                {
                    colors.pop();
                }
            } else if is_cycle {
                label = normalize_stitched_cycle::<K>(&label);
            } else {
                let reverse = reverse_complement_is_less(&label);
                if reverse && colored {
                    reverse_color_runs_in_place(&mut colors, (label.len() - K + 1) as u32);
                }
                if reverse {
                    reverse_complement_label_in_place(&mut label);
                }
            }
            emit(&label, &colors);
        }
        start = end;
    }
}

fn append_materialized_colors<const K: usize>(
    output: &mut Vec<UnitigColor>,
    output_label_len: usize,
    record: &LoadedMaterializedStitchedCoordRecord,
    color_runs: &[UnitigColor],
    reverse: bool,
) {
    if record.color_start == u32::MAX {
        return;
    }
    let start = record.color_start as usize;
    let Some(runs) = color_runs.get(start..start + record.color_count() as usize) else {
        return;
    };
    let unitig_vertex_count = record.label_len as usize - K + 1;
    if output.is_empty() {
        output.reserve(runs.len());
        if reverse {
            for index in (0..runs.len()).rev() {
                let end = runs
                    .get(index + 1)
                    .map_or(unitig_vertex_count as u32, |next| next.offset());
                output.push(UnitigColor::new(
                    unitig_vertex_count as u32 - end,
                    crate::state::ColorCoordinate::from_u40(runs[index].coordinate()),
                ));
            }
        } else {
            output.extend_from_slice(runs);
        }
    } else {
        let output_vertex_count = output_label_len - K + 1;
        append_color_runs(
            output,
            output_vertex_count as u32,
            runs,
            unitig_vertex_count as u32,
            reverse,
        );
    }
}

fn materialized_stitched_coord_records_are_ordered(
    records: &[LoadedMaterializedStitchedCoordRecord],
) -> bool {
    records
        .windows(2)
        .all(|pair| (pair[0].path_id, pair[0].rank) <= (pair[1].path_id, pair[1].rank))
}

fn reduce_stitched_coord_buckets<const K: usize>(
    inputs: &DiscontinuityInputs<K>,
    records: Vec<StitchedCoordRecord>,
    threads: usize,
) -> Vec<Vec<u8>> {
    if records.is_empty() {
        return Vec::new();
    }

    let bucket_count = (threads.max(1) * 1024)
        .next_power_of_two()
        .clamp(1024, 16_384);
    let bucket_mask = bucket_count - 1;
    let mut buckets = (0..bucket_count).map(|_| Vec::new()).collect::<Vec<_>>();
    for record in records {
        let bucket = (hash_u64(record.path_id, 0) as usize) & bucket_mask;
        buckets[bucket].push(record);
    }

    let workers = threads.max(1).min(bucket_count);
    let mut reduced_buckets = Vec::new();
    if workers == 1 {
        for mut bucket in buckets {
            reduced_buckets.push(reduce_stitched_coord_bucket::<K>(
                inputs,
                &mut bucket,
                false,
            ));
        }
    } else {
        let chunk_size = bucket_count.div_ceil(workers);
        reduced_buckets = std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for chunk in buckets.chunks_mut(chunk_size) {
                handles.push(scope.spawn(move || {
                    let mut chunk_unitigs = Vec::new();
                    for bucket in chunk {
                        chunk_unitigs
                            .extend(reduce_stitched_coord_bucket::<K>(inputs, bucket, false));
                    }
                    chunk_unitigs
                }));
            }

            let mut reduced = Vec::new();
            for handle in handles {
                reduced.push(
                    handle
                        .join()
                        .expect("stitched coordinate reducer worker panicked"),
                );
            }
            reduced
        });
    }

    let total = reduced_buckets.iter().map(Vec::len).sum();
    let mut unitigs = Vec::with_capacity(total);
    for mut bucket in reduced_buckets {
        unitigs.append(&mut bucket);
    }
    unitigs
}

fn reduce_stitched_coord_bucket_files<const K: usize>(
    inputs: &DiscontinuityInputs<K>,
    manifest: &[StitchedCoordBucketEntry],
    threads: usize,
) -> Result<Vec<Vec<u8>>, SerialCollationError> {
    if manifest.is_empty() {
        return Ok(Vec::new());
    }

    let mut groups: Vec<Vec<StitchedCoordBucketEntry>> = Vec::new();
    for entry in manifest {
        if groups
            .last()
            .and_then(|group| group.first())
            .is_some_and(|first| first.bucket_id == entry.bucket_id)
        {
            groups
                .last_mut()
                .expect("checked that a final group exists")
                .push(entry.clone());
        } else {
            groups.push(vec![entry.clone()]);
        }
    }

    let workers = threads.max(1).min(groups.len());
    let reduced_buckets = if workers == 1 {
        let mut reduced = Vec::with_capacity(groups.len());
        for group in &groups {
            reduced.push(reduce_stitched_coord_bucket_file_group::<K>(inputs, group)?);
        }
        reduced
    } else {
        let next_group = AtomicUsize::new(0);
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for _ in 0..workers {
                let next_group = &next_group;
                let groups = &groups;
                handles.push(scope.spawn(move || {
                    let mut local = Vec::new();
                    loop {
                        let group_idx = next_group.fetch_add(1, Ordering::Relaxed);
                        let Some(group) = groups.get(group_idx) else {
                            break;
                        };
                        local.extend(reduce_stitched_coord_bucket_file_group::<K>(inputs, group)?);
                    }
                    Ok::<_, SerialCollationError>(local)
                }));
            }

            let mut reduced = Vec::new();
            for handle in handles {
                reduced.push(
                    handle
                        .join()
                        .map_err(|_| SerialCollationError::WorkerPanic)??,
                );
            }
            Ok::<_, SerialCollationError>(reduced)
        })?
    };

    let total = reduced_buckets.iter().map(Vec::len).sum();
    let mut unitigs = Vec::with_capacity(total);
    for mut bucket in reduced_buckets {
        unitigs.append(&mut bucket);
    }
    Ok(unitigs)
}

fn reduce_stitched_coord_bucket_file_group<const K: usize>(
    inputs: &DiscontinuityInputs<K>,
    group: &[StitchedCoordBucketEntry],
) -> Result<Vec<Vec<u8>>, SerialCollationError> {
    let total_records = group.iter().map(|entry| entry.records).sum::<u64>();
    let mut records = Vec::with_capacity(total_records as usize);
    for entry in group {
        records.extend(read_stitched_coord_bucket_file(entry)?);
    }
    Ok(reduce_stitched_coord_bucket::<K>(
        inputs,
        &mut records,
        true,
    ))
}

fn read_stitched_coord_bucket_file(
    entry: &StitchedCoordBucketEntry,
) -> Result<Vec<StitchedCoordRecord>, SerialCollationError> {
    let file = File::open(&entry.path).map_err(|source| SerialCollationError::Io {
        path: entry.path.clone(),
        source,
    })?;
    let mut input = BufReader::with_capacity(1024 * 1024, file);
    let mut magic = [0u8; 8];
    read_exact_coord(&mut input, &entry.path, &mut magic)?;
    if &magic != STITCH_COORD_MAGIC {
        return Err(SerialCollationError::MalformedCoordBucket(
            entry.path.clone(),
        ));
    }

    let bucket_id = read_u64_coord(&mut input, &entry.path)? as usize;
    let records = read_u64_coord(&mut input, &entry.path)?;
    let _reserved = read_u64_coord(&mut input, &entry.path)?;
    if bucket_id != entry.bucket_id || records != entry.records {
        return Err(SerialCollationError::MalformedCoordBucket(
            entry.path.clone(),
        ));
    }

    let expected_len = STITCH_COORD_HEADER_LEN
        .checked_add(
            records
                .checked_mul(STITCH_PATH_INFO_RECORD_LEN)
                .ok_or_else(|| SerialCollationError::MalformedCoordBucket(entry.path.clone()))?,
        )
        .ok_or_else(|| SerialCollationError::MalformedCoordBucket(entry.path.clone()))?;
    let actual_len = input
        .get_ref()
        .metadata()
        .map_err(|source| SerialCollationError::Io {
            path: entry.path.clone(),
            source,
        })?
        .len();
    if actual_len != expected_len {
        return Err(SerialCollationError::MalformedCoordBucket(
            entry.path.clone(),
        ));
    }

    let mut coord_bytes = vec![0u8; (records * STITCH_PATH_INFO_RECORD_LEN) as usize];
    input
        .read_exact(&mut coord_bytes)
        .map_err(|source| SerialCollationError::Io {
            path: entry.path.clone(),
            source,
        })?;
    let mut out = Vec::with_capacity(records as usize);
    for chunk in coord_bytes.chunks_exact(STITCH_PATH_INFO_RECORD_LEN as usize) {
        out.push(decoded_stitched_coord_record(chunk, &entry.path)?);
    }
    Ok(out)
}

fn read_stitched_coord_bucket_group(
    group: &[StitchedCoordBucketEntry],
) -> Result<Vec<StitchedCoordRecord>, SerialCollationError> {
    if let [entry] = group {
        return read_stitched_coord_bucket_file(entry);
    }
    let total_records = group.iter().map(|entry| entry.records).sum::<u64>();
    let mut records = Vec::with_capacity(total_records as usize);
    for entry in group {
        records.extend(read_stitched_coord_bucket_file(entry)?);
    }
    Ok(records)
}

fn read_stitched_coord_bucket_group_dense(
    group: &[StitchedCoordBucketEntry],
    unitigs: usize,
    unitig_path: &Path,
    dense: &mut Vec<DenseLocalPathInfo>,
) -> Result<(), SerialCollationError> {
    dense.clear();
    dense.resize(unitigs, DenseLocalPathInfo::EMPTY);
    for entry in group {
        let file = File::open(&entry.path).map_err(|source| SerialCollationError::Io {
            path: entry.path.clone(),
            source,
        })?;
        let mut input = BufReader::with_capacity(1024 * 1024, file);
        let mut magic = [0u8; 8];
        read_exact_coord(&mut input, &entry.path, &mut magic)?;
        let bucket_id = read_u64_coord(&mut input, &entry.path)? as usize;
        let records = read_u64_coord(&mut input, &entry.path)?;
        let _reserved = read_u64_coord(&mut input, &entry.path)?;
        let expected_len = STITCH_COORD_HEADER_LEN
            .checked_add(
                records
                    .checked_mul(STITCH_PATH_INFO_RECORD_LEN)
                    .ok_or_else(|| {
                        SerialCollationError::MalformedCoordBucket(entry.path.clone())
                    })?,
            )
            .ok_or_else(|| SerialCollationError::MalformedCoordBucket(entry.path.clone()))?;
        let actual_len = input
            .get_ref()
            .metadata()
            .map_err(|source| SerialCollationError::Io {
                path: entry.path.clone(),
                source,
            })?
            .len();
        if &magic != STITCH_COORD_MAGIC
            || bucket_id != entry.bucket_id
            || records != entry.records
            || actual_len != expected_len
        {
            return Err(SerialCollationError::MalformedCoordBucket(
                entry.path.clone(),
            ));
        }

        let mut bytes = vec![0u8; (records * STITCH_PATH_INFO_RECORD_LEN) as usize];
        read_exact_coord(&mut input, &entry.path, &mut bytes)?;
        for record in bytes.chunks_exact(STITCH_PATH_INFO_RECORD_LEN as usize) {
            let path_id = u64::from_le_bytes(record[..8].try_into().expect("path ID"));
            let rank = u64::from_le_bytes(record[8..16].try_into().expect("path rank"));
            let unitig_index =
                u32::from_le_bytes(record[16..20].try_into().expect("local unitig index")) as usize;
            let flags = record[20];
            let Some(slot) = dense.get_mut(unitig_index) else {
                return Err(SerialCollationError::MalformedCoordBucket(
                    unitig_path.to_path_buf(),
                ));
            };
            if slot.rank_and_flags != u64::MAX || record[21..].iter().any(|&byte| byte != 0) {
                return Err(SerialCollationError::MalformedCoordBucket(
                    entry.path.clone(),
                ));
            }
            *slot = DenseLocalPathInfo {
                path_id,
                rank_and_flags: (rank << 2)
                    | u64::from((flags & STITCH_COORD_REVERSE_FLAG) != 0)
                    | (u64::from((flags & STITCH_COORD_CYCLE_FLAG) != 0) << 1),
            };
        }
    }
    Ok(())
}

fn decoded_stitched_coord_record(
    bytes: &[u8],
    path: &Path,
) -> Result<StitchedCoordRecord, SerialCollationError> {
    if bytes.len() != STITCH_PATH_INFO_RECORD_LEN as usize {
        return Err(SerialCollationError::MalformedCoordBucket(
            path.to_path_buf(),
        ));
    }
    let path_id = u64::from_le_bytes(bytes[..8].try_into().expect("u64 path_id field"));
    let rank = u64::from_le_bytes(bytes[8..16].try_into().expect("u64 rank field"));
    let unitig_index = u32::from_le_bytes(bytes[16..20].try_into().expect("u32 unitig field"));
    let flags = bytes[20];
    if bytes[21..].iter().any(|&byte| byte != 0) {
        return Err(SerialCollationError::MalformedCoordBucket(
            path.to_path_buf(),
        ));
    }
    Ok(StitchedCoordRecord {
        path_id,
        rank,
        unitig_index,
        reverse: (flags & STITCH_COORD_REVERSE_FLAG) != 0,
        is_cycle: (flags & STITCH_COORD_CYCLE_FLAG) != 0,
    })
}

fn read_u64_coord(input: &mut BufReader<File>, path: &Path) -> Result<u64, SerialCollationError> {
    let mut bytes = [0u8; 8];
    read_exact_coord(input, path, &mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

fn read_exact_coord(
    input: &mut BufReader<File>,
    path: &Path,
    bytes: &mut [u8],
) -> Result<(), SerialCollationError> {
    input
        .read_exact(bytes)
        .map_err(|source| SerialCollationError::Io {
            path: path.to_path_buf(),
            source,
        })
}

fn reduce_stitched_coord_bucket<const K: usize>(
    inputs: &DiscontinuityInputs<K>,
    records: &mut Vec<StitchedCoordRecord>,
    sort_records: bool,
) -> Vec<Vec<u8>> {
    if records.is_empty() {
        return Vec::new();
    }

    if sort_records {
        records.sort_unstable_by_key(|record| (record.path_id, record.rank));
    }
    debug_assert!(
        records
            .windows(2)
            .all(|pair| (pair[0].path_id, pair[0].rank) <= (pair[1].path_id, pair[1].rank))
    );

    let mut unitigs = Vec::new();
    let mut start = 0;
    while start < records.len() {
        let path_id = records[start].path_id;
        let is_cycle = records[start].is_cycle;
        let mut end = start + 1;
        while end < records.len() && records[end].path_id == path_id {
            end += 1;
        }

        let mut label = Vec::new();
        for record in &records[start..end] {
            let unitig = &inputs.unitigs[record.unitig_index as usize];
            let unitig_label = unitig.label(inputs);
            let mut reverse = record.reverse;
            if !label.is_empty()
                && !labels_overlap_oriented_fast::<K>(&label, unitig_label, reverse)
            {
                let alternate = oriented_label(unitig_label, !reverse);
                if labels_overlap::<K>(&label, &alternate) {
                    reverse = !reverse;
                }
            }
            append_or_init_oriented_fast::<K>(&mut label, unitig_label, reverse);
        }

        if label.len() >= K {
            let label = if is_cycle {
                normalize_stitched_cycle::<K>(&label)
            } else {
                canonical_label(label)
            };
            unitigs.push(label);
        }

        start = end;
    }

    unitigs
}

fn debug_endpoint<const K: usize>(endpoint: DiscontinuityEndpoint<K>) -> (String, Side) {
    (endpoint.vertex.to_ascii_string(), endpoint.side)
}

fn endpoints_by_label_end<const K: usize>(
    unitig: &DiscontinuityUnitig<K>,
) -> (
    Option<DiscontinuityEndpoint<K>>,
    Option<DiscontinuityEndpoint<K>>,
) {
    (unitig.left_exit(), unitig.right_exit())
}

fn push_stitch_edge(
    adjacency: &mut [StitchAdjacencyList],
    first: usize,
    second: usize,
    unitig_index: Option<usize>,
) {
    adjacency[first].push(StitchAdjacency {
        to: second,
        unitig_index,
    });
    adjacency[second].push(StitchAdjacency {
        to: first,
        unitig_index,
    });
}

fn walk_stitched_path<const K: usize>(
    inputs: &DiscontinuityInputs<K>,
    adjacency: &[StitchAdjacencyList],
    visited_segments: &mut [bool],
    start: usize,
) -> Option<StitchedPath> {
    let mut label = Vec::new();
    let mut prev = None;
    let mut current = start;
    let mut emitted_segment = false;
    let mut is_cycle = false;

    loop {
        let next = adjacency[current]
            .iter()
            .copied()
            .find(|edge| {
                Some(edge.to) != prev
                    && edge
                        .unitig_index
                        .is_none_or(|unitig_index| !visited_segments[unitig_index])
            })
            .or_else(|| {
                adjacency[current].iter().copied().find(|edge| {
                    edge.unitig_index
                        .is_some_and(|unitig_index| !visited_segments[unitig_index])
                })
            });

        let Some(edge) = next else {
            break;
        };

        if let Some(unitig_index) = edge.unitig_index {
            if visited_segments[unitig_index] {
                break;
            }
            visited_segments[unitig_index] = true;
            let unitig = &inputs.unitigs[unitig_index];
            let unitig_label = unitig.label(inputs);
            let reverse = reverse_for_stitch_node(current);
            let mut append_reverse = reverse;
            if !label.is_empty()
                && !labels_overlap_oriented_fast::<K>(&label, unitig_label, reverse)
            {
                let alternate = oriented_label(unitig_label, !reverse);
                if labels_overlap::<K>(&label, &alternate) {
                    append_reverse = !reverse;
                }
            }
            append_or_init_oriented_fast::<K>(&mut label, unitig_label, append_reverse);
            emitted_segment = true;
        }

        prev = Some(current);
        current = edge.to;
        if emitted_segment && current == start {
            is_cycle = true;
            break;
        }
    }

    emitted_segment.then_some(StitchedPath { label, is_cycle })
}

fn walk_simple_stitched_component_coords(
    half_ends: &[HalfEnd],
    join_neighbor: &[u32],
    start: usize,
    path_id: u64,
    records: &mut Vec<StitchedCoordRecord>,
) {
    let record_start = records.len();
    let mut current = start;
    let mut rank = 0u64;
    let mut is_cycle = false;

    loop {
        let unitig_index = half_ends[current].unitig_index;
        let reverse = reverse_for_stitch_node(current);
        records.push(StitchedCoordRecord {
            path_id,
            rank,
            unitig_index,
            reverse,
            is_cycle: false,
        });
        rank += 1;

        let other = current ^ 1;
        let next = join_neighbor[other];
        if next == STITCH_NO_NODE {
            break;
        }
        if stitch_node_index(next) == start {
            is_cycle = true;
            break;
        }
        current = stitch_node_index(next);
    }

    if rank == 0 {
        records.truncate(record_start);
        return;
    }

    if is_cycle {
        for record in &mut records[record_start..] {
            record.is_cycle = true;
        }
    }
}

fn walk_simple_stitched_component_coords_from_adjacency(
    adjacency: &[StitchAdjacencyList],
    start: usize,
    path_id: u64,
    records: &mut Vec<StitchedCoordRecord>,
) {
    let record_start = records.len();
    let mut prev = None;
    let mut current = start;
    let mut rank = 0u64;
    let mut emitted_segment = false;
    let mut is_cycle = false;

    loop {
        let next = adjacency[current]
            .iter()
            .copied()
            .find(|edge| Some(edge.to) != prev);

        let Some(edge) = next else {
            break;
        };

        if let Some(unitig_index) = edge.unitig_index {
            let reverse = reverse_for_stitch_node(current);
            records.push(StitchedCoordRecord {
                path_id,
                rank,
                unitig_index: u32::try_from(unitig_index).expect("unitig index exceeds u32"),
                reverse,
                is_cycle: false,
            });
            rank += 1;
            emitted_segment = true;
        }

        prev = Some(current);
        current = edge.to;
        if emitted_segment && current == start {
            is_cycle = true;
            break;
        }
    }

    if !emitted_segment {
        records.truncate(record_start);
        return;
    }

    if is_cycle {
        for record in &mut records[record_start..] {
            record.is_cycle = true;
        }
    }
}

fn labels_overlap<const K: usize>(left: &[u8], right: &[u8]) -> bool {
    left.len() >= K && right.len() >= K && left[left.len() - K..] == right[..K]
}

fn labels_overlap_oriented_fast<const K: usize>(left: &[u8], right: &[u8], reverse: bool) -> bool {
    if left.len() < K || right.len() < K {
        return false;
    }

    let suffix = &left[left.len() - K..];
    if reverse {
        suffix
            .iter()
            .enumerate()
            .all(|(i, &base)| base == complement_ascii(right[right.len() - 1 - i]))
    } else {
        suffix == &right[..K]
    }
}

fn append_or_init_oriented_fast<const K: usize>(label: &mut Vec<u8>, next: &[u8], reverse: bool) {
    let start = if label.is_empty() { 0 } else { K };
    if !reverse {
        label.extend_from_slice(&next[start..]);
        return;
    }

    label.reserve(next.len().saturating_sub(start));
    for &base in next[..next.len() - start].iter().rev() {
        label.push(complement_ascii(base));
    }
}

#[inline]
fn reverse_complement_label_in_place(label: &mut [u8]) {
    let mut left = 0;
    let mut right = label.len();
    while left < right {
        right -= 1;
        if left == right {
            label[left] = complement_ascii(label[left]);
            break;
        }
        let left_base = label[left];
        label[left] = complement_ascii(label[right]);
        label[right] = complement_ascii(left_base);
        left += 1;
    }
}

fn normalize_stitched_cycle<const K: usize>(label: &[u8]) -> Vec<u8> {
    let mut graph = StitchCycleGraph::<K>::new();
    if graph.add_label(label).is_err() {
        return canonical_cycle_label(label.to_vec());
    }

    let mut labels = graph.contract();
    labels.sort_by_key(|candidate| std::cmp::Reverse(candidate.len()));
    labels
        .into_iter()
        .next()
        .map(canonical_label)
        .unwrap_or_else(|| canonical_cycle_label(label.to_vec()))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct StitchCycleEdge<const K: usize> {
    from: Kmer<K>,
    to: Kmer<K>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct StitchCycleKmer<const K: usize> {
    observed: Kmer<K>,
}

impl<const K: usize> StitchCycleKmer<K> {
    fn new(observed: Kmer<K>) -> Self {
        Self { observed }
    }

    fn canonical(self) -> Kmer<K> {
        self.observed.canonical()
    }

    fn in_canonical_form(self) -> bool {
        self.observed.is_canonical()
    }

    fn entrance_side(self) -> Side {
        if self.in_canonical_form() {
            Side::Front
        } else {
            Side::Back
        }
    }

    fn roll_forward(self, base: Base) -> Self {
        Self::new(self.observed.roll_forward(base))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StitchCycleWalk<const K: usize> {
    label: Vec<u8>,
    anchor: Kmer<K>,
}

impl<const K: usize> StitchCycleWalk<K> {
    fn init(v: StitchCycleKmer<K>) -> Self {
        Self {
            label: v.observed.to_ascii_string().into_bytes(),
            anchor: v.canonical(),
        }
    }

    fn extend(&mut self, v: StitchCycleKmer<K>, base: Base) -> bool {
        if v.canonical() == self.anchor {
            return false;
        }

        self.label.push(base.to_ascii());
        true
    }
}

struct StitchCycleGraph<const K: usize> {
    vertices: FastHashMap<Kmer<K>, VertexState>,
    unique_edges: BTreeSet<StitchCycleEdge<K>>,
}

impl<const K: usize> StitchCycleGraph<K> {
    fn new() -> Self {
        Self {
            vertices: FastHashMap::default(),
            unique_edges: BTreeSet::new(),
        }
    }

    fn add_label(&mut self, label: &[u8]) -> Result<(), ()> {
        if label.len() < K {
            return Ok(());
        }

        let last_vertex_offset = label.len() - K;
        let mut prev = None;
        for offset in 0..=last_vertex_offset {
            let directed = StitchCycleKmer::new(
                Kmer::<K>::from_ascii(&label[offset..offset + K]).map_err(|_| ())?,
            );
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
                self.unique_edges.insert(StitchCycleEdge {
                    from,
                    to: canonical,
                });
            }
            prev = Some(canonical);
        }

        Ok(())
    }

    fn contract(&mut self) -> Vec<Vec<u8>> {
        let mut unitigs = Vec::new();
        let mut vertices = self.vertices.keys().copied().collect::<Vec<_>>();
        vertices.sort_unstable();

        for v_hat in vertices {
            let Some(state) = self.vertices.get(&v_hat).copied() else {
                continue;
            };
            if state.is_visited() || state.is_isolated(1) {
                continue;
            }

            unitigs.push(self.extract_unitig(v_hat));
        }

        unitigs
    }

    fn extract_unitig(&mut self, v_hat: Kmer<K>) -> Vec<u8> {
        let (back_walk, back_is_cycle) = self.walk_unitig(v_hat, Side::Back);
        if back_is_cycle {
            return back_walk.label;
        }

        let (front_walk, _) = self.walk_unitig(v_hat, Side::Front);
        let mut label = reverse_complement_label(&front_walk.label);
        label.extend_from_slice(&back_walk.label[K..]);
        label
    }

    fn walk_unitig(&mut self, v_hat: Kmer<K>, start_side: Side) -> (StitchCycleWalk<K>, bool) {
        let icc_return_side = start_side.inverse();
        let mut v = if start_side == Side::Back {
            StitchCycleKmer::new(v_hat)
        } else {
            StitchCycleKmer::new(v_hat.reverse_complement())
        };
        let mut side = start_side;
        let mut walk = StitchCycleWalk::init(v);

        loop {
            let canonical = v.canonical();
            let Some(state) = self.vertices.get(&canonical).copied() else {
                return (walk, false);
            };
            let Some(vertex_state) = self.vertices.get_mut(&canonical) else {
                return (walk, false);
            };
            vertex_state.mark_visited();

            let mut edge = state.edge_at(side, 1);
            if edge == Base::N || edge == Base::E {
                return (walk, false);
            }

            if side == Side::Front {
                edge = edge.complement();
            }
            v = v.roll_forward(edge);

            let Some(next_state) = self.vertices.get(&v.canonical()).copied() else {
                return (walk, false);
            };
            side = v.entrance_side();
            if next_state.is_branching_side(side, 1) {
                return (walk, false);
            }
            if next_state.is_visited() {
                return (walk, v.canonical() == v_hat && side == icc_return_side);
            }

            if !walk.extend(v, edge) {
                return (walk, false);
            }
            side = side.inverse();
        }
    }
}

fn canonical_cycle_label(label: Vec<u8>) -> Vec<u8> {
    if label.is_empty() {
        return label;
    }
    let forward = minimal_rotation(&label);
    let reverse = minimal_rotation(&reverse_complement_label(&label));
    if reverse < forward { reverse } else { forward }
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct DiagonalOtherEnd<const K: usize> {
    vertex: Kmer<K>,
    side_at_vertex: Side,
    side_at_current: Side,
    weight: u64,
    unitig_index: usize,
    unitig_exit_side: Side,
    is_phi: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PartitionOtherEnd<const K: usize> {
    endpoint: MatrixEndpoint<K>,
    side_at_current: Side,
    weight: u64,
    in_same_part: bool,
    processed: bool,
}

struct PathInfoSlot<const K: usize> {
    tag: AtomicU64,
    key: UnsafeCell<MaybeUninit<Kmer<K>>>,
    value: UnsafeCell<MaybeUninit<PathInfo<K>>>,
}

unsafe impl<const K: usize> Sync for PathInfoSlot<K> {}

struct ConcurrentPathInfoTable<const K: usize> {
    slots: Vec<PathInfoSlot<K>>,
    mask: usize,
}

struct CompactPathInfoSlot {
    key: AtomicU64,
    path_id: AtomicU64,
    rank_and_flags: AtomicU64,
    epoch: AtomicU8,
}

const _: () = assert!(std::mem::size_of::<CompactPathInfoSlot>() == 32);

#[derive(Clone, Copy)]
struct CompactExpansionPathInfo {
    path_id: u64,
    rank_and_flags: u64,
}

impl CompactExpansionPathInfo {
    #[inline(always)]
    fn new(path_id: u64, rank: u64, exit_side: Side, is_cycle: bool) -> Self {
        let flags = u64::from(exit_side == Side::Back) | (u64::from(is_cycle) << 1);
        Self {
            path_id,
            rank_and_flags: (rank << 2) | flags,
        }
    }

    #[inline(always)]
    fn rank(self) -> u64 {
        self.rank_and_flags >> 2
    }

    #[inline(always)]
    fn exit_side(self) -> Side {
        if self.rank_and_flags & 1 == 0 {
            Side::Front
        } else {
            Side::Back
        }
    }

    #[inline(always)]
    fn is_cycle(self) -> bool {
        self.rank_and_flags & 2 != 0
    }
}

struct CompactPathInfoTable {
    slots: Vec<CompactPathInfoSlot>,
    mask: usize,
    key_mask: u64,
    generation: AtomicU8,
}

impl CompactPathInfoTable {
    fn with_max_entries(max_entries: usize, k: usize) -> Self {
        debug_assert!(k <= 31);
        let capacity = max_entries
            // Match C++'s expansion map sizing. Expansion performs several
            // lookups per edge, so its 50% maximum load is materially faster
            // than the denser contraction-table sizing at scale.
            .saturating_mul(2)
            .next_power_of_two()
            .max(8);
        let key_mask = (1u64 << (2 * k)) - 1;
        Self {
            slots: (0..capacity)
                .map(|_| CompactPathInfoSlot {
                    key: AtomicU64::new(0),
                    path_id: AtomicU64::new(0),
                    rank_and_flags: AtomicU64::new(0),
                    epoch: AtomicU8::new(0),
                })
                .collect(),
            mask: capacity - 1,
            key_mask,
            generation: AtomicU8::new(0),
        }
    }

    fn clear(&self) {
        let current = self.generation.load(Ordering::Relaxed);
        let mut next = current.wrapping_add(1);
        if next == 0 || next == u8::MAX {
            self.slots
                .par_iter()
                .for_each(|slot| slot.epoch.store(0, Ordering::Relaxed));
            next = 1;
        }
        self.generation.store(next, Ordering::Relaxed);
    }

    fn insert<const K: usize>(&self, vertex: Kmer<K>, value: PathInfo<K>) -> bool {
        let flags = u8::from(value.exit_side == Side::Back) | (u8::from(value.is_cycle) << 1);
        self.insert_raw(
            vertex.as_u128() as u64,
            value.path_id.as_u128() as u64,
            value.rank,
            flags,
        )
    }

    #[inline(always)]
    fn insert_compact<const K: usize>(
        &self,
        vertex: Kmer<K>,
        value: CompactExpansionPathInfo,
    ) -> bool {
        self.insert_raw(
            vertex.as_u128() as u64,
            value.path_id,
            value.rank(),
            value.rank_and_flags as u8 & 3,
        )
    }

    #[inline(always)]
    fn insert_raw(&self, key: u64, path_id: u64, rank: u64, flags: u8) -> bool {
        debug_assert!(rank <= (u64::MAX >> 2));
        debug_assert!(flags < 4);
        self.insert_packed(key, path_id, (rank << 2) | u64::from(flags))
    }

    #[inline(always)]
    fn insert_packed(&self, key: u64, path_id: u64, rank_and_flags: u64) -> bool {
        const BUSY: u8 = u8::MAX;
        debug_assert_eq!(key & !self.key_mask, 0);
        let generation = self.generation.load(Ordering::Relaxed);
        let mut idx = wyhash_u64(key, 0) as usize & self.mask;
        loop {
            let slot = &self.slots[idx];
            let observed_epoch = slot.epoch.load(Ordering::Acquire);
            if observed_epoch == generation {
                if slot.key.load(Ordering::Relaxed) == key {
                    return false;
                }
                idx = (idx + 1) & self.mask;
                continue;
            }
            if observed_epoch == BUSY {
                std::hint::spin_loop();
                continue;
            }
            if slot
                .epoch
                .compare_exchange(observed_epoch, BUSY, Ordering::AcqRel, Ordering::Acquire)
                .is_ok()
            {
                // Parallel loading only inserts values; lookups begin after the
                // worker pool joins. Later diagonal inserts are single-threaded.
                slot.key.store(key, Ordering::Relaxed);
                slot.path_id.store(path_id, Ordering::Relaxed);
                slot.rank_and_flags.store(rank_and_flags, Ordering::Relaxed);
                slot.epoch.store(generation, Ordering::Release);
                return true;
            }
        }
    }

    #[inline(always)]
    fn get<const K: usize>(&self, vertex: Kmer<K>) -> Option<PathInfo<K>> {
        let value = self.get_compact(vertex)?;
        Some(PathInfo {
            path_id: Kmer::from_bits(value.path_id as u128),
            rank: value.rank(),
            exit_side: value.exit_side(),
            is_cycle: value.is_cycle(),
        })
    }

    #[inline(always)]
    fn get_compact<const K: usize>(&self, vertex: Kmer<K>) -> Option<CompactExpansionPathInfo> {
        self.get_compact_raw(vertex.as_u128() as u64)
    }

    #[inline(always)]
    fn get_compact_raw(&self, key: u64) -> Option<CompactExpansionPathInfo> {
        let generation = self.generation.load(Ordering::Relaxed);
        let mut idx = wyhash_u64(key, 0) as usize & self.mask;
        loop {
            let slot = &self.slots[idx];
            if slot.epoch.load(Ordering::Acquire) != generation {
                return None;
            }
            if slot.key.load(Ordering::Relaxed) == key {
                let path_id = slot.path_id.load(Ordering::Relaxed);
                let rank_and_flags = slot.rank_and_flags.load(Ordering::Relaxed);
                return Some(CompactExpansionPathInfo {
                    path_id,
                    rank_and_flags,
                });
            }
            idx = (idx + 1) & self.mask;
        }
    }
}

enum ExpansionPathInfoTable<const K: usize> {
    Compact(CompactPathInfoTable),
    Wide(ConcurrentPathInfoTable<K>),
}

impl<const K: usize> ExpansionPathInfoTable<K> {
    fn with_max_entries(max_entries: usize) -> Self {
        if K <= 31 {
            Self::Compact(CompactPathInfoTable::with_max_entries(max_entries, K))
        } else {
            Self::Wide(ConcurrentPathInfoTable::with_max_entries(max_entries))
        }
    }

    fn clear(&self) {
        match self {
            Self::Compact(table) => table.clear(),
            Self::Wide(table) => table.clear(),
        }
    }

    #[inline(always)]
    fn insert(&self, vertex: Kmer<K>, value: PathInfo<K>) -> bool {
        match self {
            Self::Compact(table) => table.insert(vertex, value),
            Self::Wide(table) => table.insert(vertex, value),
        }
    }

    #[inline(always)]
    fn get(&self, vertex: Kmer<K>) -> Option<PathInfo<K>> {
        match self {
            Self::Compact(table) => table.get(vertex),
            Self::Wide(table) => table.get(vertex),
        }
    }

    #[inline(always)]
    fn get_compact(&self, vertex: Kmer<K>) -> Option<CompactExpansionPathInfo> {
        match self {
            Self::Compact(table) => table.get_compact(vertex),
            Self::Wide(_) => unreachable!("compact path info is only used for k <= 31"),
        }
    }

    #[inline(always)]
    fn insert_compact(&self, vertex: Kmer<K>, value: CompactExpansionPathInfo) -> bool {
        match self {
            Self::Compact(table) => table.insert_compact(vertex, value),
            Self::Wide(_) => unreachable!("compact path info is only used for k <= 31"),
        }
    }

    #[inline(always)]
    fn insert_compact_record(&self, record: CompactVertexPathInfoRecord) -> bool {
        match self {
            Self::Compact(table) => {
                table.insert_packed(record.vertex, record.path_id, record.rank_and_flags)
            }
            Self::Wide(_) => unreachable!("compact path info is only used for k <= 31"),
        }
    }

    fn contains_key(&self, vertex: &Kmer<K>) -> bool {
        self.get(*vertex).is_some()
    }

    fn capacity(&self) -> usize {
        match self {
            Self::Compact(table) => table.slots.len(),
            Self::Wide(table) => table.slots.len(),
        }
    }

    fn slot_size(&self) -> usize {
        match self {
            Self::Compact(_) => std::mem::size_of::<CompactPathInfoSlot>(),
            Self::Wide(_) => std::mem::size_of::<PathInfoSlot<K>>(),
        }
    }

    #[inline(always)]
    fn insert_encoded(&self, bytes: &[u8]) -> bool {
        match self {
            Self::Compact(table) => {
                if discontinuity_edge_kmer_bytes::<K>() != 8 {
                    let record = decoded_vertex_path_info_record::<K>(bytes);
                    return table.insert(record.vertex, record.info);
                }
                let mut key = [0u8; 8];
                key.copy_from_slice(&bytes[..8]);
                let mut path_id = [0u8; 8];
                path_id.copy_from_slice(&bytes[8..16]);
                let mut rank = [0u8; 8];
                rank.copy_from_slice(&bytes[16..24]);
                table.insert_packed(
                    u64::from_le_bytes(key),
                    u64::from_le_bytes(path_id),
                    u64::from_le_bytes(rank),
                )
            }
            Self::Wide(table) => {
                let record = decoded_vertex_path_info_record::<K>(bytes);
                table.insert(record.vertex, record.info)
            }
        }
    }
}

impl<const K: usize> ConcurrentPathInfoTable<K> {
    const EMPTY: u64 = 0;
    const BUSY: u64 = 1;

    fn with_max_entries(max_entries: usize) -> Self {
        let capacity = max_entries.saturating_mul(2).next_power_of_two().max(8);
        Self {
            slots: (0..capacity)
                .map(|_| PathInfoSlot {
                    tag: AtomicU64::new(Self::EMPTY),
                    key: UnsafeCell::new(MaybeUninit::uninit()),
                    value: UnsafeCell::new(MaybeUninit::uninit()),
                })
                .collect(),
            mask: capacity - 1,
        }
    }

    #[inline]
    fn tag(vertex: Kmer<K>) -> u64 {
        vertex.hash64(0).wrapping_add(2).max(2)
    }

    fn clear(&self) {
        self.slots
            .par_iter()
            .for_each(|slot| slot.tag.store(Self::EMPTY, Ordering::Relaxed));
    }

    fn insert(&self, vertex: Kmer<K>, value: PathInfo<K>) -> bool {
        let tag = Self::tag(vertex);
        let mut idx = vertex.hash64(0) as usize & self.mask;
        loop {
            let slot = &self.slots[idx];
            let observed = slot.tag.load(Ordering::Acquire);
            if observed == Self::EMPTY
                && slot
                    .tag
                    .compare_exchange(
                        Self::EMPTY,
                        Self::BUSY,
                        Ordering::Acquire,
                        Ordering::Relaxed,
                    )
                    .is_ok()
            {
                unsafe {
                    (*slot.key.get()).write(vertex);
                    (*slot.value.get()).write(value);
                }
                slot.tag.store(tag, Ordering::Release);
                return true;
            }
            if observed == Self::BUSY {
                std::hint::spin_loop();
                continue;
            }
            if observed == tag && unsafe { (*slot.key.get()).assume_init() == vertex } {
                return false;
            }
            idx = (idx + 1) & self.mask;
        }
    }

    fn get(&self, vertex: Kmer<K>) -> Option<PathInfo<K>> {
        let tag = Self::tag(vertex);
        let mut idx = vertex.hash64(0) as usize & self.mask;
        loop {
            let slot = &self.slots[idx];
            let observed = slot.tag.load(Ordering::Acquire);
            if observed == Self::EMPTY {
                return None;
            }
            if observed == Self::BUSY {
                std::hint::spin_loop();
                continue;
            }
            if observed == tag && unsafe { (*slot.key.get()).assume_init() == vertex } {
                return Some(unsafe { (*slot.value.get()).assume_init() });
            }
            idx = (idx + 1) & self.mask;
        }
    }

    #[inline]
    fn contains_key(&self, vertex: &Kmer<K>) -> bool {
        self.get(*vertex).is_some()
    }
}

struct OwnedPartitionTables<const K: usize> {
    maps: Vec<FastHashMap<Kmer<K>, PartitionOtherEnd<K>>>,
    mask: usize,
}

impl<const K: usize> OwnedPartitionTables<K> {
    fn new(threads: usize) -> Self {
        let count = threads.max(1).next_power_of_two();
        Self {
            maps: (0..count).map(|_| FastHashMap::default()).collect(),
            mask: count - 1,
        }
    }

    #[inline]
    fn owner(&self, vertex: Kmer<K>) -> usize {
        partition_column_vertex_shard(vertex, self.mask)
    }

    fn clear(&mut self) {
        for map in &mut self.maps {
            map.clear();
        }
    }

    fn contains_key(&self, vertex: &Kmer<K>) -> bool {
        self.maps[self.owner(*vertex)].contains_key(vertex)
    }

    fn get(&self, vertex: &Kmer<K>) -> Option<&PartitionOtherEnd<K>> {
        self.maps[self.owner(*vertex)].get(vertex)
    }

    fn get_mut(&mut self, vertex: &Kmer<K>) -> Option<&mut PartitionOtherEnd<K>> {
        let owner = self.owner(*vertex);
        self.maps[owner].get_mut(vertex)
    }

    fn insert(&mut self, vertex: Kmer<K>, end: PartitionOtherEnd<K>) {
        let owner = self.owner(vertex);
        self.maps[owner].insert(vertex, end);
    }
}

#[derive(Clone, Copy)]
struct CompactPartitionOtherEnd {
    endpoint: u64,
    weight: u16,
    flags: u8,
}

impl CompactPartitionOtherEnd {
    const PHI: u8 = 1 << 4;

    fn pack<const K: usize>(end: PartitionOtherEnd<K>) -> Self {
        let (endpoint, endpoint_side, phi) = match end.endpoint {
            MatrixEndpoint::Phi => (0, Side::Front, true),
            MatrixEndpoint::Vertex(endpoint) => {
                (endpoint.vertex.as_u128() as u64, endpoint.side, false)
            }
        };
        let mut flags = u8::from(endpoint_side == Side::Back)
            | (u8::from(end.side_at_current == Side::Back) << 1)
            | (u8::from(end.in_same_part) << 2)
            | (u8::from(end.processed) << 3);
        if phi {
            flags |= Self::PHI;
        }
        Self {
            endpoint,
            weight: u16::try_from(end.weight).expect("partition edge weight fits u16"),
            flags,
        }
    }

    fn unpack<const K: usize>(self) -> PartitionOtherEnd<K> {
        let side = |bit| {
            if self.flags & bit == 0 {
                Side::Front
            } else {
                Side::Back
            }
        };
        PartitionOtherEnd {
            endpoint: if self.flags & Self::PHI != 0 {
                MatrixEndpoint::Phi
            } else {
                MatrixEndpoint::Vertex(DiscontinuityEndpoint {
                    vertex: Kmer::from_bits(self.endpoint as u128),
                    side: side(1),
                })
            },
            side_at_current: side(1 << 1),
            weight: u64::from(self.weight),
            in_same_part: self.flags & (1 << 2) != 0,
            processed: self.flags & (1 << 3) != 0,
        }
    }
}

struct AtomicPartitionSlot {
    key: AtomicU64,
    value: UnsafeCell<MaybeUninit<CompactPartitionOtherEnd>>,
}

struct FlatPartitionTable<const K: usize> {
    keys: Vec<u64>,
    values: Vec<MaybeUninit<PartitionOtherEnd<K>>>,
    occupied: Vec<usize>,
    mask: usize,
}

impl<const K: usize> FlatPartitionTable<K> {
    const EMPTY: u64 = u64::MAX;

    fn with_max_entries(max_entries: usize) -> Self {
        let capacity = max_entries
            .saturating_mul(4)
            .div_ceil(3)
            .next_power_of_two()
            .max(8);
        Self {
            keys: vec![Self::EMPTY; capacity],
            values: (0..capacity).map(|_| MaybeUninit::uninit()).collect(),
            occupied: Vec::with_capacity(max_entries),
            mask: capacity - 1,
        }
    }

    fn clear(&mut self) {
        for idx in self.occupied.drain(..) {
            self.keys[idx] = Self::EMPTY;
        }
    }

    #[inline]
    fn find_index(&self, vertex: Kmer<K>) -> Result<usize, usize> {
        let key = vertex.as_u128() as u64;
        let mut idx = vertex.hash64(0) as usize & self.mask;
        loop {
            let observed = self.keys[idx];
            if observed == key {
                return Ok(idx);
            }
            if observed == Self::EMPTY {
                return Err(idx);
            }
            idx = (idx + 1) & self.mask;
        }
    }

    fn absorb(
        &mut self,
        vertex: Kmer<K>,
        incoming: PartitionOtherEnd<K>,
        partition: usize,
        edges: &mut Vec<DiscontinuityEdge<K>>,
        meta_vertices: &mut Vec<SerialMetaVertex<K>>,
    ) {
        match self.find_index(vertex) {
            Ok(idx) => {
                let existing = unsafe { self.values[idx].assume_init() };
                absorb_partition_collision_outputs(
                    vertex,
                    incoming,
                    existing,
                    partition,
                    edges,
                    meta_vertices,
                );
                let mut updated = existing;
                if updated.in_same_part {
                    updated = incoming;
                }
                updated.processed = true;
                self.values[idx].write(updated);
            }
            Err(idx) => {
                self.keys[idx] = vertex.as_u128() as u64;
                self.values[idx].write(incoming);
                self.occupied.push(idx);
            }
        }
    }

    fn get(&self, vertex: Kmer<K>) -> Option<PartitionOtherEnd<K>> {
        self.find_index(vertex)
            .ok()
            .map(|idx| unsafe { self.values[idx].assume_init() })
    }

    fn insert(&mut self, vertex: Kmer<K>, value: PartitionOtherEnd<K>) {
        match self.find_index(vertex) {
            Ok(idx) => {
                self.values[idx].write(value);
            }
            Err(idx) => {
                self.keys[idx] = vertex.as_u128() as u64;
                self.values[idx].write(value);
                self.occupied.push(idx);
            }
        }
    }

    fn mark_processed(&mut self, vertex: Kmer<K>) {
        let idx = self.find_index(vertex).expect("partition endpoint exists");
        let mut value = unsafe { self.values[idx].assume_init() };
        value.processed = true;
        self.values[idx].write(value);
    }
}

unsafe impl Sync for AtomicPartitionSlot {}

struct AtomicPartitionTable<const K: usize> {
    slots: Vec<AtomicPartitionSlot>,
    mask: usize,
}

impl<const K: usize> AtomicPartitionTable<K> {
    const EMPTY: u64 = u64::MAX;
    const LOCKED: u64 = 1 << 63;
    const HASH_SEED: u64 = 0xAAAA_AAAA_5555_5555;

    #[inline(always)]
    fn index(&self, vertex: Kmer<K>) -> usize {
        vertex.hash64(Self::HASH_SEED) as usize & self.mask
    }

    #[allow(dead_code)]
    fn with_max_entries(max_entries: usize) -> Self {
        let capacity = max_entries
            .saturating_mul(4)
            .div_ceil(3)
            .next_power_of_two()
            .max(8);
        Self {
            slots: (0..capacity)
                .map(|_| AtomicPartitionSlot {
                    key: AtomicU64::new(Self::EMPTY),
                    value: UnsafeCell::new(MaybeUninit::uninit()),
                })
                .collect(),
            mask: capacity - 1,
        }
    }

    fn clear(&self, threads: usize, pool: &ThreadPool) {
        let workers = threads.max(1).min(self.slots.len());
        let chunk = self.slots.len().div_ceil(workers);
        pool.install(|| {
            self.slots.par_chunks(chunk).for_each(|slots| {
                for slot in slots {
                    slot.key.store(Self::EMPTY, Ordering::Relaxed);
                }
            })
        });
    }

    fn absorb(
        &self,
        vertex: Kmer<K>,
        incoming: PartitionOtherEnd<K>,
        partition: usize,
        edges: &mut Vec<DiscontinuityEdge<K>>,
        meta_vertices: &mut Vec<SerialMetaVertex<K>>,
    ) {
        let key = vertex.as_u128() as u64;
        debug_assert!(K <= 31 && key != Self::EMPTY);
        let mut idx = self.index(vertex);
        loop {
            let slot = &self.slots[idx];
            let observed = slot.key.load(Ordering::Acquire);
            if observed == key {
                if slot
                    .key
                    .compare_exchange_weak(
                        key,
                        key | Self::LOCKED,
                        Ordering::Acquire,
                        Ordering::Relaxed,
                    )
                    .is_err()
                {
                    std::hint::spin_loop();
                    continue;
                }
                let existing = unsafe { (*slot.value.get()).assume_init() }.unpack::<K>();
                absorb_partition_collision_outputs(
                    vertex,
                    incoming,
                    existing,
                    partition,
                    edges,
                    meta_vertices,
                );
                let mut updated = existing;
                if updated.in_same_part {
                    updated = incoming;
                }
                updated.processed = true;
                unsafe { (*slot.value.get()).write(CompactPartitionOtherEnd::pack(updated)) };
                slot.key.store(key, Ordering::Release);
                return;
            }
            if observed == Self::EMPTY {
                if slot
                    .key
                    .compare_exchange_weak(
                        Self::EMPTY,
                        Self::LOCKED,
                        Ordering::Acquire,
                        Ordering::Relaxed,
                    )
                    .is_err()
                {
                    std::hint::spin_loop();
                    continue;
                }
                unsafe { (*slot.value.get()).write(CompactPartitionOtherEnd::pack(incoming)) };
                slot.key.store(key, Ordering::Release);
                return;
            }
            if observed & Self::LOCKED != 0 {
                std::hint::spin_loop();
                continue;
            }
            idx = (idx + 1) & self.mask;
        }
    }

    fn absorb_prepared(
        &self,
        vertex: Kmer<K>,
        incoming: PartitionOtherEnd<K>,
        partition: usize,
        vertex_partitions: usize,
        edges: &mut Vec<PreparedBlockedEdge>,
        meta_vertices: &mut Vec<SerialMetaVertex<K>>,
    ) {
        let key = vertex.as_u128() as u64;
        debug_assert!(K <= 31 && key != Self::EMPTY);
        let mut idx = self.index(vertex);
        loop {
            let slot = &self.slots[idx];
            let observed = slot.key.load(Ordering::Acquire);
            if observed == key {
                if slot
                    .key
                    .compare_exchange_weak(
                        key,
                        key | Self::LOCKED,
                        Ordering::Acquire,
                        Ordering::Relaxed,
                    )
                    .is_err()
                {
                    std::hint::spin_loop();
                    continue;
                }
                let existing = unsafe { (*slot.value.get()).assume_init() }.unpack::<K>();
                if existing.endpoint.is_phi() && incoming.endpoint.is_phi() {
                    meta_vertices.push(two_weight_meta_vertex(
                        vertex,
                        partition,
                        incoming.side_at_current,
                        incoming.weight,
                        existing.weight,
                        false,
                    ));
                } else if !existing.in_same_part {
                    edges.push(prepare_existing_blocked_edge(
                        join_other_ends(
                            incoming.endpoint,
                            existing.endpoint,
                            incoming.weight + existing.weight,
                        ),
                        vertex_partitions,
                    ));
                }
                let mut updated = existing;
                if updated.in_same_part {
                    updated = incoming;
                }
                updated.processed = true;
                unsafe { (*slot.value.get()).write(CompactPartitionOtherEnd::pack(updated)) };
                slot.key.store(key, Ordering::Release);
                return;
            }
            if observed == Self::EMPTY {
                if slot
                    .key
                    .compare_exchange_weak(
                        Self::EMPTY,
                        Self::LOCKED,
                        Ordering::Acquire,
                        Ordering::Relaxed,
                    )
                    .is_err()
                {
                    std::hint::spin_loop();
                    continue;
                }
                unsafe { (*slot.value.get()).write(CompactPartitionOtherEnd::pack(incoming)) };
                slot.key.store(key, Ordering::Release);
                return;
            }
            if observed & Self::LOCKED != 0 {
                std::hint::spin_loop();
                continue;
            }
            idx = (idx + 1) & self.mask;
        }
    }

    fn to_fast_map(&self) -> FastHashMap<Kmer<K>, PartitionOtherEnd<K>> {
        let occupied = self
            .slots
            .iter()
            .filter(|slot| slot.key.load(Ordering::Relaxed) != Self::EMPTY)
            .count();
        let mut map = FastHashMap::with_capacity_and_hasher(occupied, FastBuildHasher::default());
        for slot in &self.slots {
            let key = slot.key.load(Ordering::Relaxed);
            if key != Self::EMPTY {
                let value = unsafe { (*slot.value.get()).assume_init() }.unpack::<K>();
                map.insert(Kmer::from_bits(key as u128), value);
            }
        }
        map
    }

    fn get(&self, vertex: Kmer<K>) -> Option<PartitionOtherEnd<K>> {
        let key = vertex.as_u128() as u64;
        let mut idx = self.index(vertex);
        loop {
            let slot = &self.slots[idx];
            let observed = slot.key.load(Ordering::Acquire);
            if observed == key {
                return Some(unsafe { (*slot.value.get()).assume_init() }.unpack::<K>());
            }
            if observed == Self::EMPTY {
                return None;
            }
            idx = (idx + 1) & self.mask;
        }
    }

    fn insert_serial(&self, vertex: Kmer<K>, value: PartitionOtherEnd<K>) {
        let key = vertex.as_u128() as u64;
        let mut idx = self.index(vertex);
        loop {
            let slot = &self.slots[idx];
            let observed = slot.key.load(Ordering::Relaxed);
            if observed == key {
                unsafe { (*slot.value.get()).write(CompactPartitionOtherEnd::pack(value)) };
                return;
            }
            if observed == Self::EMPTY {
                unsafe { (*slot.value.get()).write(CompactPartitionOtherEnd::pack(value)) };
                slot.key.store(key, Ordering::Release);
                return;
            }
            idx = (idx + 1) & self.mask;
        }
    }

    fn mark_processed(&self, vertex: Kmer<K>) {
        let key = vertex.as_u128() as u64;
        let mut idx = self.index(vertex);
        loop {
            let slot = &self.slots[idx];
            let observed = slot.key.load(Ordering::Relaxed);
            if observed == key {
                let mut value = unsafe { (*slot.value.get()).assume_init() }.unpack::<K>();
                value.processed = true;
                unsafe { (*slot.value.get()).write(CompactPartitionOtherEnd::pack(value)) };
                return;
            }
            assert_ne!(observed, Self::EMPTY, "partition endpoint disappeared");
            idx = (idx + 1) & self.mask;
        }
    }
}

fn absorb_partition_other_end<const K: usize>(
    vertex: Kmer<K>,
    incoming: PartitionOtherEnd<K>,
    partition: usize,
    table: &mut FastHashMap<Kmer<K>, PartitionOtherEnd<K>>,
    edges: &mut Vec<DiscontinuityEdge<K>>,
    meta_vertices: &mut Vec<SerialMetaVertex<K>>,
) {
    match table.get_mut(&vertex) {
        None => {
            table.insert(vertex, incoming);
        }
        Some(existing) => {
            if existing.endpoint.is_phi() && incoming.endpoint.is_phi() {
                meta_vertices.push(two_weight_meta_vertex(
                    vertex,
                    partition,
                    incoming.side_at_current,
                    incoming.weight,
                    existing.weight,
                    false,
                ));
            } else if existing.in_same_part {
                *existing = incoming;
            } else {
                edges.push(join_other_ends(
                    incoming.endpoint,
                    existing.endpoint,
                    incoming.weight + existing.weight,
                ));
            }
            existing.processed = true;
        }
    }
}

fn absorb_partition_collision_outputs<const K: usize>(
    vertex: Kmer<K>,
    incoming: PartitionOtherEnd<K>,
    existing: PartitionOtherEnd<K>,
    partition: usize,
    edges: &mut Vec<DiscontinuityEdge<K>>,
    meta_vertices: &mut Vec<SerialMetaVertex<K>>,
) {
    if existing.endpoint.is_phi() && incoming.endpoint.is_phi() {
        meta_vertices.push(two_weight_meta_vertex(
            vertex,
            partition,
            incoming.side_at_current,
            incoming.weight,
            existing.weight,
            false,
        ));
    } else if !existing.in_same_part {
        edges.push(join_other_ends(
            incoming.endpoint,
            existing.endpoint,
            incoming.weight + existing.weight,
        ));
    }
}

fn endpoint_in_partition<const K: usize>(
    matrix: &SerialEdgeMatrix<K>,
    edge: &DiscontinuityEdge<K>,
    partition: usize,
) -> Option<(MatrixEndpoint<K>, DiscontinuityEndpoint<K>)> {
    match (edge.first, edge.second) {
        (lower, MatrixEndpoint::Vertex(current))
            if matrix.partition(MatrixEndpoint::Vertex(current)) == partition =>
        {
            Some((lower, current))
        }
        (MatrixEndpoint::Vertex(current), lower)
            if matrix.partition(MatrixEndpoint::Vertex(current)) == partition =>
        {
            Some((lower, current))
        }
        _ => None,
    }
}

fn join_other_ends<const K: usize>(
    first: MatrixEndpoint<K>,
    second: MatrixEndpoint<K>,
    weight: u64,
) -> DiscontinuityEdge<K> {
    join_other_ends_with_phantom(first, second, weight, None)
}

fn join_other_ends_with_phantom<const K: usize>(
    first: MatrixEndpoint<K>,
    second: MatrixEndpoint<K>,
    weight: u64,
    phantom_unitig: Option<DiscontinuityEndpoint<K>>,
) -> DiscontinuityEdge<K> {
    DiscontinuityEdge {
        first,
        second,
        weight,
        unitig_bucket: 0,
        unitig_index: 0,
        unitig_exit_side: Side::Back,
        phantom_unitig,
        swapped: false,
    }
}

fn two_weight_meta_vertex<const K: usize>(
    vertex: Kmer<K>,
    partition: usize,
    side: Side,
    front_weight: u64,
    back_weight: u64,
    is_cycle: bool,
) -> SerialMetaVertex<K> {
    SerialMetaVertex {
        vertex,
        partition,
        entry_side: Side::Back,
        weight: if side == Side::Front {
            front_weight
        } else {
            back_weight
        },
        is_cycle,
    }
}

fn max_endpoint_partition<const K: usize>(
    matrix: &SerialEdgeMatrix<K>,
    first: MatrixEndpoint<K>,
    second: MatrixEndpoint<K>,
) -> usize {
    matrix.partition(first).max(matrix.partition(second))
}

fn edge_matrix_partition<const K: usize>(
    vertex_partitions: usize,
    endpoint: MatrixEndpoint<K>,
) -> usize {
    match endpoint {
        MatrixEndpoint::Phi => 0,
        MatrixEndpoint::Vertex(endpoint) => {
            ((endpoint.vertex.hash64(0) as usize) & (vertex_partitions - 1)) + 1
        }
    }
}

fn edge_matrix_row_col<const K: usize>(
    vertex_partitions: usize,
    first: MatrixEndpoint<K>,
    second: MatrixEndpoint<K>,
) -> (usize, usize) {
    let first_partition = edge_matrix_partition(vertex_partitions, first);
    let second_partition = edge_matrix_partition(vertex_partitions, second);
    if first_partition <= second_partition {
        (first_partition, second_partition)
    } else {
        (second_partition, first_partition)
    }
}

fn endpoint_sort_key<const K: usize>(endpoint: MatrixEndpoint<K>) -> (u8, u128, u8) {
    match endpoint {
        MatrixEndpoint::Phi => (0, 0, 0),
        MatrixEndpoint::Vertex(endpoint) => (1, endpoint.vertex.as_u128(), endpoint.side as u8),
    }
}

fn endpoint_id<const K: usize>(
    ids: &mut FastHashMap<EndpointKey<K>, usize>,
    adjacency: &mut Vec<Vec<(usize, u64, bool)>>,
    endpoint: MatrixEndpoint<K>,
) -> usize {
    let next = ids.len();
    *ids.entry(EndpointKey(endpoint)).or_insert_with(|| {
        adjacency.push(Vec::new());
        next
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SerialEdgeMatrixError {
    InvalidPartitionCount(usize),
}

impl std::fmt::Display for SerialEdgeMatrixError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidPartitionCount(count) => write!(
                f,
                "vertex partition count must be a non-zero power of two, got {count}"
            ),
        }
    }
}

impl std::error::Error for SerialEdgeMatrixError {}

pub fn emit_uncolored_discontinuity_inputs<const K: usize>(
    bucket_dir: impl AsRef<Path>,
    cutoff: u32,
) -> Result<DiscontinuityInputs<K>, DiscontinuityInputError> {
    emit_uncolored_discontinuity_inputs_with_threads::<K>(bucket_dir, cutoff, 1)
}

pub fn emit_uncolored_discontinuity_inputs_with_threads<const K: usize>(
    bucket_dir: impl AsRef<Path>,
    cutoff: u32,
    threads: usize,
) -> Result<DiscontinuityInputs<K>, DiscontinuityInputError> {
    emit_uncolored_discontinuity_inputs_with_threads_impl::<K>(bucket_dir, cutoff, threads, None)
}

pub fn emit_uncolored_discontinuity_inputs_with_threads_in_dir<const K: usize>(
    bucket_dir: impl AsRef<Path>,
    cutoff: u32,
    threads: usize,
    label_path: impl AsRef<Path>,
) -> Result<DiscontinuityInputs<K>, DiscontinuityInputError> {
    emit_uncolored_discontinuity_inputs_with_threads_impl::<K>(
        bucket_dir,
        cutoff,
        threads,
        Some(label_path.as_ref()),
    )
}

pub fn emit_uncolored_external_discontinuity_inputs_with_threads_in_dir<const K: usize>(
    bucket_dir: impl AsRef<Path>,
    cutoff: u32,
    threads: usize,
    label_path: impl AsRef<Path>,
) -> Result<ExternalDiscontinuityInputs<K>, DiscontinuityInputError> {
    if cutoff == 0 {
        return Err(DiscontinuityInputError::InvalidCutoff);
    }
    if threads == 0 {
        return Err(DiscontinuityInputError::InvalidThreadCount);
    }
    let entries = read_manifest(bucket_dir.as_ref())?;
    contract_local_subgraphs_into_external_inputs::<K>(
        &entries,
        cutoff,
        threads,
        label_path.as_ref(),
        None,
    )
}

pub fn emit_colored_external_discontinuity_inputs_with_threads_in_dir<const K: usize>(
    bucket_dir: impl AsRef<Path>,
    cutoff: u32,
    threads: usize,
    label_path: impl AsRef<Path>,
    color_path: impl AsRef<Path>,
) -> Result<ExternalDiscontinuityInputs<K>, DiscontinuityInputError> {
    if cutoff == 0 {
        return Err(DiscontinuityInputError::InvalidCutoff);
    }
    if threads == 0 {
        return Err(DiscontinuityInputError::InvalidThreadCount);
    }
    let entries = read_manifest(bucket_dir.as_ref())?;
    contract_local_subgraphs_into_external_inputs::<K>(
        &entries,
        cutoff,
        threads,
        label_path.as_ref(),
        Some(color_path.as_ref()),
    )
}

fn emit_uncolored_discontinuity_inputs_with_threads_impl<const K: usize>(
    bucket_dir: impl AsRef<Path>,
    cutoff: u32,
    threads: usize,
    label_path: Option<&Path>,
) -> Result<DiscontinuityInputs<K>, DiscontinuityInputError> {
    if cutoff == 0 {
        return Err(DiscontinuityInputError::InvalidCutoff);
    }
    if threads == 0 {
        return Err(DiscontinuityInputError::InvalidThreadCount);
    }

    let entries = read_manifest(bucket_dir.as_ref())?;
    if let Some(label_path) = label_path {
        let external = contract_local_subgraphs_into_external_inputs::<K>(
            &entries, cutoff, threads, label_path, None,
        )?;
        return external_inputs_to_memory_inputs(external);
    }

    let outputs = contract_local_subgraphs::<K>(&entries, cutoff, threads)?;
    let mut inputs = DiscontinuityInputs::empty(DiscontinuityInputStats {
        input_buckets: entries.len(),
        ..DiscontinuityInputStats::default()
    });

    for output in outputs {
        inputs.stats.weak_superkmers += output.weak_superkmers;
        let label_offset = inputs.labels.len() as u64;
        inputs.labels.extend_from_slice(&output.labels);
        for mut unitig in output.unitigs {
            inputs.stats.local_unitigs += 1;
            inputs.stats.discontinuity_exits +=
                u64::from(unitig.left_exit().is_some()) + u64::from(unitig.right_exit().is_some());
            inputs.stats.unitig_bases += unitig.label_len as u64;
            unitig.label_start += label_offset;
            inputs.unitigs.push(unitig);
        }
    }

    Ok(inputs)
}

fn external_inputs_to_memory_inputs<const K: usize>(
    external: ExternalDiscontinuityInputs<K>,
) -> Result<DiscontinuityInputs<K>, DiscontinuityInputError> {
    let unitig_file =
        File::open(&external.unitig_path).map_err(|source| DiscontinuityInputError::Io {
            path: external.unitig_path.clone(),
            source,
        })?;
    let mut unitig_input = BufReader::with_capacity(1024 * 1024, unitig_file);
    let mut unitigs = Vec::with_capacity(external.unitigs);
    for _ in 0..external.unitigs {
        unitigs.push(read_discontinuity_unitig_from_reader_for_input(
            &mut unitig_input,
            &external.unitig_path,
        )?);
    }

    let mut label_file =
        File::open(&external.label_path).map_err(|source| DiscontinuityInputError::Io {
            path: external.label_path.clone(),
            source,
        })?;
    let label_len = label_file
        .metadata()
        .map_err(|source| DiscontinuityInputError::Io {
            path: external.label_path.clone(),
            source,
        })?
        .len() as usize;
    let mut labels = vec![0u8; label_len];
    label_file
        .read_exact(&mut labels)
        .map_err(|source| DiscontinuityInputError::Io {
            path: external.label_path.clone(),
            source,
        })?;

    Ok(DiscontinuityInputs {
        unitigs,
        labels,
        stats: external.stats,
    })
}

fn contract_local_subgraphs_into_external_inputs<const K: usize>(
    entries: &[BucketManifestEntry],
    cutoff: u32,
    threads: usize,
    label_path: &Path,
    color_path: Option<&Path>,
) -> Result<ExternalDiscontinuityInputs<K>, DiscontinuityInputError> {
    // C++ hands both colored and uncolored local contractions to expansion in
    // max-unitig coordinate buckets. The representation does not depend on
    // colors; keeping uncolored on the legacy global unitig/range stream adds
    // a full random-access materialization at scale.
    let compact_unitigs = std::env::var_os("CF3_RS_ENDPOINT_STITCH").is_none();
    let mut inputs = DiscontinuityInputs::empty(DiscontinuityInputStats {
        input_buckets: entries.len(),
        ..DiscontinuityInputStats::default()
    });

    if let Some(parent) = label_path.parent() {
        fs::create_dir_all(parent).map_err(|source| DiscontinuityInputError::Io {
            path: parent.to_path_buf(),
            source,
        })?;
    }
    let file = File::create(label_path).map_err(|source| DiscontinuityInputError::Io {
        path: label_path.to_path_buf(),
        source,
    })?;
    let mut labels = BufWriter::with_capacity(8 * 1024 * 1024, file);
    let unitig_path = unitig_table_path_for_labels(label_path);
    let unitig_file = File::create(&unitig_path).map_err(|source| DiscontinuityInputError::Io {
        path: unitig_path.clone(),
        source,
    })?;
    let mut unitigs = BufWriter::with_capacity(8 * 1024 * 1024, unitig_file);
    let trivial_path = label_path.with_extension("trivial.fa");
    let trivial_file =
        File::create(&trivial_path).map_err(|source| DiscontinuityInputError::Io {
            path: trivial_path.clone(),
            source,
        })?;
    let mut trivial_output = BufWriter::with_capacity(8 * 1024 * 1024, trivial_file);
    let mut trivial_unitigs = 0u64;
    let mut trivial_bases = 0u64;
    let matrix_dir = unitig_path.with_file_name(format!(
        "{}.edge-matrix",
        unitig_path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("unitigs")
    ));
    let mut edge_matrix = BlockedEdgeMatrix::create(&matrix_dir, DEFAULT_VERTEX_PARTITIONS)
        .map_err(serial_collation_to_input_error)?;
    let mut label_offset = 0u64;
    let mut ranges = Vec::new();
    let expected_colors = entries
        .iter()
        .map(|entry| entry.records)
        .sum::<u64>()
        .div_ceil(16)
        .min(48 * 1024 * 1024) as usize;
    let color_repository = color_path
        .map(|path| {
            ConcurrentColorRepository::create(
                path.with_extension("color-repository"),
                threads.min(256),
                expected_colors.max(entries.len() * 8),
            )
        })
        .transpose()?;
    let mut groups = local_bucket_groups(entries);
    groups.sort_by_key(|group| std::cmp::Reverse(group.records));
    let workers = threads.min(groups.len().max(1));
    let local_unitig_bucket_count: usize = if compact_unitigs && workers > 1 && groups.len() >= 1024
    {
        1024
    } else {
        0
    };
    let local_unitig_bucket_dir = unitig_path.with_extension("local-unitig-buckets");
    let mut local_unitig_writers = if local_unitig_bucket_count == 0 {
        None
    } else {
        let _ = raise_open_file_limit();
        if local_unitig_bucket_dir.exists() {
            fs::remove_dir_all(&local_unitig_bucket_dir).map_err(|source| {
                DiscontinuityInputError::Io {
                    path: local_unitig_bucket_dir.clone(),
                    source,
                }
            })?;
        }
        fs::create_dir_all(&local_unitig_bucket_dir).map_err(|source| {
            DiscontinuityInputError::Io {
                path: local_unitig_bucket_dir.clone(),
                source,
            }
        })?;
        Some(
            (1..=local_unitig_bucket_count)
                .map(|bucket_id| {
                    LocalUnitigBucketWriter::create(
                        &local_unitig_bucket_dir,
                        bucket_id as u16,
                        color_path.is_some(),
                    )
                    .map(Mutex::new)
                })
                .collect::<Result<Vec<_>, _>>()?,
        )
    };
    let mut color_run_writer = if workers == 1 {
        color_path.map(ColorRunSidecarWriter::create).transpose()?
    } else {
        None
    };
    let concurrent_color_runs = if workers > 1 && local_unitig_bucket_count == 0 {
        color_path
            .map(ConcurrentColorRunSidecarWriter::create)
            .transpose()?
    } else {
        None
    };
    let started = Instant::now();
    if groups.len() >= 1024 {
        eprintln!(
            "cuttlefish3-rs: contracting {} local subgraph(s) from {} bucket file(s) with {} worker(s)",
            groups.len(),
            entries.len(),
            workers
        );
    }

    let mut build_elapsed = Duration::default();
    let mut contract_elapsed = Duration::default();
    let mut sink_io_elapsed = Duration::default();
    let mut sink_color_elapsed = Duration::default();
    let mut sink_edge_elapsed = Duration::default();

    if workers == 1 {
        let mut reusable_vertices = None;
        for (offset, group) in groups.iter().enumerate() {
            let output = contract_local_subgraph::<K>(
                group,
                cutoff,
                color_repository.as_ref(),
                &mut reusable_vertices,
            )?;
            append_local_contraction_output_to_external_inputs(
                &mut trivial_output,
                &trivial_path,
                &mut trivial_unitigs,
                &mut trivial_bases,
                output,
                &mut inputs,
                &mut labels,
                label_path,
                &mut unitigs,
                &unitig_path,
                &mut label_offset,
                &mut ranges,
                &mut edge_matrix,
                &mut build_elapsed,
                &mut contract_elapsed,
                color_repository.as_ref(),
                color_run_writer.as_mut(),
                threads,
                compact_unitigs,
            )?;
            report_local_contraction_progress(offset + 1, groups.len(), started);
        }
    } else {
        let total_groups = groups.len();
        let next_group = AtomicUsize::new(0);
        let completed = Arc::new(AtomicUsize::new(0));
        let sink = ConcurrentLocalOutputSink::<K> {
            labels: labels
                .get_ref()
                .try_clone()
                .map_err(|source| DiscontinuityInputError::Io {
                    path: label_path.to_path_buf(),
                    source,
                })?,
            label_path,
            unitigs: unitigs.get_ref().try_clone().map_err(|source| {
                DiscontinuityInputError::Io {
                    path: unitig_path.clone(),
                    source,
                }
            })?,
            unitig_path: &unitig_path,
            next_label: AtomicU64::new(0),
            next_unitig: AtomicUsize::new(0),
            weak_superkmers: AtomicU64::new(0),
            discontinuity_exits: AtomicU64::new(0),
            unitig_bases: AtomicU64::new(0),
            build_nanos: AtomicU64::new(0),
            contract_nanos: AtomicU64::new(0),
            sink_io_nanos: AtomicU64::new(0),
            sink_color_nanos: AtomicU64::new(0),
            sink_edge_nanos: AtomicU64::new(0),
            ranges: Mutex::new(Vec::with_capacity(groups.len())),
            edge_writers: ConcurrentBlockedEdgeWriters::new(&edge_matrix),
            color_repository: color_repository.as_ref(),
            color_runs: concurrent_color_runs.as_ref(),
            local_unitig_writers: local_unitig_writers.as_deref(),
            compact_unitigs,
            trivial_output: trivial_output.get_ref().try_clone().map_err(|source| {
                DiscontinuityInputError::Io {
                    path: trivial_path.clone(),
                    source,
                }
            })?,
            trivial_path: &trivial_path,
            next_trivial_byte: AtomicU64::new(0),
            trivial_unitigs: AtomicU64::new(0),
            trivial_bases: AtomicU64::new(0),
            marker: PhantomData,
        };
        let worker_result = std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for worker_id in 0..workers {
                let completed = Arc::clone(&completed);
                let next_group = &next_group;
                let groups = &groups;
                let sink = &sink;
                handles.push(scope.spawn(move || {
                    let mut reusable_vertices = None;
                    let writers_per_worker = local_unitig_bucket_count.div_ceil(workers);
                    let bucket_start = worker_id * writers_per_worker;
                    let bucket_end =
                        (bucket_start + writers_per_worker).min(local_unitig_bucket_count);
                    let mut next_bucket = bucket_start;
                    loop {
                        let group_idx = next_group.fetch_add(1, Ordering::Relaxed);
                        let Some(group) = groups.get(group_idx) else {
                            break;
                        };
                        let output = contract_local_subgraph::<K>(
                            group,
                            cutoff,
                            sink.color_repository,
                            &mut reusable_vertices,
                        )
                        .and_then(|output| {
                            let bucket_id = if bucket_start < bucket_end {
                                let bucket = next_bucket;
                                next_bucket += 1;
                                if next_bucket == bucket_end {
                                    next_bucket = bucket_start;
                                }
                                (bucket + 1) as u16
                            } else {
                                0
                            };
                            sink.write(output, bucket_id)
                        });
                        let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
                        report_local_contraction_progress(done, total_groups, started);
                        let should_stop = output.is_err();
                        if should_stop {
                            return output;
                        }
                    }
                    Ok(())
                }));
            }
            let mut first_error = None;
            for handle in handles {
                let result = handle
                    .join()
                    .map_err(|_| DiscontinuityInputError::WorkerPanic)?;
                if let Err(err) = result
                    && first_error.is_none()
                {
                    first_error = Some(err);
                }
            }

            if let Some(err) = first_error {
                Err(err)
            } else {
                Ok(())
            }
        });
        worker_result?;
        sink.edge_writers
            .finish_into(&mut edge_matrix)
            .map_err(serial_collation_to_input_error)?;
        inputs.stats.weak_superkmers = sink.weak_superkmers.load(Ordering::Relaxed);
        inputs.stats.local_unitigs = sink.next_unitig.load(Ordering::Relaxed) as u64;
        trivial_unitigs = sink.trivial_unitigs.load(Ordering::Relaxed);
        trivial_bases = sink.trivial_bases.load(Ordering::Relaxed);
        inputs.stats.local_unitigs += trivial_unitigs;
        inputs.stats.discontinuity_exits = sink.discontinuity_exits.load(Ordering::Relaxed);
        inputs.stats.unitig_bases = sink.unitig_bases.load(Ordering::Relaxed);
        inputs.stats.unitig_bases += trivial_bases;
        build_elapsed = Duration::from_nanos(sink.build_nanos.load(Ordering::Relaxed));
        contract_elapsed = Duration::from_nanos(sink.contract_nanos.load(Ordering::Relaxed));
        sink_io_elapsed = Duration::from_nanos(sink.sink_io_nanos.load(Ordering::Relaxed));
        sink_color_elapsed = Duration::from_nanos(sink.sink_color_nanos.load(Ordering::Relaxed));
        sink_edge_elapsed = Duration::from_nanos(sink.sink_edge_nanos.load(Ordering::Relaxed));
        ranges = sink
            .ranges
            .into_inner()
            .map_err(|_| DiscontinuityInputError::WorkerPanic)?;
        ranges.sort_by_key(|range| range.start_unitig);
    }

    eprintln!(
        "cuttlefish3-rs: local worker time: bucket read/build {:.3}s, unitig walk {:.3}s",
        build_elapsed.as_secs_f64(),
        contract_elapsed.as_secs_f64()
    );
    if workers > 1 {
        eprintln!(
            "cuttlefish3-rs: local sink worker time: label/unitig I/O {:.3}s, color resolve/write {:.3}s, edge/range emission {:.3}s",
            sink_io_elapsed.as_secs_f64(),
            sink_color_elapsed.as_secs_f64(),
            sink_edge_elapsed.as_secs_f64(),
        );
    }

    labels
        .flush()
        .map_err(|source| DiscontinuityInputError::Io {
            path: label_path.to_path_buf(),
            source,
        })?;
    unitigs
        .flush()
        .map_err(|source| DiscontinuityInputError::Io {
            path: unitig_path.clone(),
            source,
        })?;
    trivial_output
        .flush()
        .map_err(|source| DiscontinuityInputError::Io {
            path: trivial_path.clone(),
            source,
        })?;
    drop(labels);
    drop(unitigs);
    edge_matrix
        .flush_all_with_threads(threads)
        .map_err(serial_collation_to_input_error)?;
    let color_runs = match (color_run_writer, concurrent_color_runs) {
        (Some(writer), None) => Some(writer.finish()?),
        (None, Some(writer)) => Some(writer.finish()?),
        (None, None) => None,
        (Some(_), Some(_)) => unreachable!("color writers are mutually exclusive"),
    };
    let local_unitig_buckets = local_unitig_writers
        .take()
        .map(|writers| finish_local_unitig_writers(writers, workers))
        .transpose()?;
    let color_repository = color_repository
        .map(|repository| repository.finish())
        .transpose()?;
    let trivial_fasta = if trivial_unitigs == 0 {
        let _ = fs::remove_file(&trivial_path);
        None
    } else {
        Some((trivial_path, trivial_unitigs, trivial_bases))
    };
    Ok(ExternalDiscontinuityInputs {
        unitig_path,
        label_path: label_path.to_path_buf(),
        // Trivial uncolored unitigs are already in the direct FASTA artifact.
        unitigs: inputs.stats.local_unitigs.saturating_sub(trivial_unitigs) as usize,
        compact_unitigs,
        ranges,
        edge_matrix: Some(edge_matrix),
        color_runs,
        local_unitig_buckets,
        trivial_fasta,
        color_repository,
        stats: inputs.stats,
    })
}

fn unitig_table_path_for_labels(label_path: &Path) -> PathBuf {
    let Some(file_name) = label_path.file_name().and_then(|name| name.to_str()) else {
        return label_path.with_extension("unitigs");
    };
    let unitig_name = if let Some(prefix) = file_name.strip_suffix("labels") {
        format!("{prefix}unitigs")
    } else {
        format!("{file_name}.unitigs")
    };
    label_path.with_file_name(unitig_name)
}

fn append_local_contraction_output_to_external_inputs<const K: usize>(
    trivial_output: &mut BufWriter<File>,
    trivial_path: &Path,
    trivial_unitigs: &mut u64,
    trivial_bases: &mut u64,
    output: LocalContractionOutput<K>,
    inputs: &mut DiscontinuityInputs<K>,
    labels: &mut BufWriter<File>,
    label_path: &Path,
    unitigs: &mut BufWriter<File>,
    unitig_path: &Path,
    label_offset: &mut u64,
    ranges: &mut Vec<ExternalLocalUnitigRange>,
    edge_matrix: &mut BlockedEdgeMatrix<K>,
    build_elapsed: &mut Duration,
    contract_elapsed: &mut Duration,
    _color_repository: Option<&ConcurrentColorRepository>,
    mut color_run_writer: Option<&mut ColorRunSidecarWriter>,
    _threads: usize,
    compact_unitigs: bool,
) -> Result<(), DiscontinuityInputError> {
    trivial_output
        .write_all(&output.trivial_fasta)
        .map_err(|source| DiscontinuityInputError::Io {
            path: trivial_path.to_path_buf(),
            source,
        })?;
    *trivial_unitigs += output.trivial_unitigs;
    *trivial_bases += output.trivial_bases;
    inputs.stats.local_unitigs += output.trivial_unitigs;
    inputs.stats.unitig_bases += output.trivial_bases;
    inputs.stats.weak_superkmers += output.weak_superkmers;
    let start_unitig = usize::try_from(inputs.stats.local_unitigs).unwrap_or(usize::MAX);
    let output_unitigs = output.unitigs.len();
    let output_label_len = output.labels.len() as u64;
    let color_start = color_run_writer
        .as_deref()
        .map(ColorRunSidecarWriter::position)
        .unwrap_or(0);
    if output_unitigs != 0 {
        ranges.push(ExternalLocalUnitigRange {
            start_unitig,
            unitigs: output_unitigs,
            label_start: *label_offset,
            label_len: output_label_len,
            color_start,
        });
    }
    labels
        .write_all(&output.labels)
        .map_err(|source| DiscontinuityInputError::Io {
            path: label_path.to_path_buf(),
            source,
        })?;
    edge_matrix
        .add_prepared_edges(&output.matrix_edges, start_unitig)
        .map_err(serial_collation_to_input_error)?;
    let mut output_colors = output.color_runs.map(Vec::into_iter);
    for mut unitig in output.unitigs {
        inputs.stats.local_unitigs += 1;
        inputs.stats.discontinuity_exits +=
            u64::from(unitig.left_exit().is_some()) + u64::from(unitig.right_exit().is_some());
        inputs.stats.unitig_bases += unitig.label_len as u64;
        unitig.label_start += *label_offset;
        write_discontinuity_unitig_record(unitigs, unitig_path, &unitig, compact_unitigs)?;
        if let Some(writer) = color_run_writer.as_deref_mut() {
            let runs = output_colors
                .as_mut()
                .and_then(Iterator::next)
                .ok_or(DiscontinuityInputError::MissingColorRuns)?;
            writer.write_unitig(&runs)?;
        }
    }
    if output_colors
        .as_mut()
        .is_some_and(|colors| colors.next().is_some())
    {
        return Err(DiscontinuityInputError::MissingColorRuns);
    }
    *label_offset += output.labels.len() as u64;
    *build_elapsed += output.build_elapsed;
    *contract_elapsed += output.contract_elapsed;
    Ok(())
}

fn write_discontinuity_unitig_record<const K: usize>(
    out: &mut BufWriter<File>,
    path: &Path,
    unitig: &DiscontinuityUnitig<K>,
    compact: bool,
) -> Result<(), DiscontinuityInputError> {
    if compact {
        let mut bytes = [0u8; 8];
        bytes[..4].copy_from_slice(&unitig.label_len.to_le_bytes());
        bytes[4] = unitig.flags;
        return out
            .write_all(&bytes)
            .map_err(|source| DiscontinuityInputError::Io {
                path: path.to_path_buf(),
                source,
            });
    }
    let mut stored = MaybeUninit::<DiscontinuityUnitig<K>>::zeroed();
    unsafe {
        let ptr = stored.as_mut_ptr();
        (*ptr).label_start = unitig.label_start;
        (*ptr).left_vertex = unitig.left_vertex;
        (*ptr).right_vertex = unitig.right_vertex;
        (*ptr).label_len = unitig.label_len;
        (*ptr).flags = unitig.flags;
        let stored = stored.assume_init();
        let bytes = std::slice::from_raw_parts(
            (&stored as *const DiscontinuityUnitig<K>).cast::<u8>(),
            std::mem::size_of::<DiscontinuityUnitig<K>>(),
        );
        out.write_all(bytes)
            .map_err(|source| DiscontinuityInputError::Io {
                path: path.to_path_buf(),
                source,
            })?;
    }
    Ok(())
}

struct LocalContractionOutput<const K: usize> {
    index: usize,
    weak_superkmers: u64,
    build_elapsed: Duration,
    contract_elapsed: Duration,
    unitigs: Vec<DiscontinuityUnitig<K>>,
    labels: Vec<u8>,
    matrix_edges: Vec<PreparedBlockedEdge>,
    color_runs: Option<Vec<Vec<UnitigColor>>>,
    trivial_fasta: Vec<u8>,
    trivial_unitigs: u64,
    trivial_bases: u64,
}

struct ConcurrentLocalOutputSink<'a, const K: usize> {
    labels: File,
    label_path: &'a Path,
    unitigs: File,
    unitig_path: &'a Path,
    next_label: AtomicU64,
    next_unitig: AtomicUsize,
    weak_superkmers: AtomicU64,
    discontinuity_exits: AtomicU64,
    unitig_bases: AtomicU64,
    build_nanos: AtomicU64,
    contract_nanos: AtomicU64,
    sink_io_nanos: AtomicU64,
    sink_color_nanos: AtomicU64,
    sink_edge_nanos: AtomicU64,
    ranges: Mutex<Vec<ExternalLocalUnitigRange>>,
    edge_writers: ConcurrentBlockedEdgeWriters,
    color_repository: Option<&'a ConcurrentColorRepository>,
    color_runs: Option<&'a ConcurrentColorRunSidecarWriter>,
    local_unitig_writers: Option<&'a [Mutex<LocalUnitigBucketWriter>]>,
    compact_unitigs: bool,
    trivial_output: File,
    trivial_path: &'a Path,
    next_trivial_byte: AtomicU64,
    trivial_unitigs: AtomicU64,
    trivial_bases: AtomicU64,
    marker: PhantomData<[(); K]>,
}

impl<'a, const K: usize> ConcurrentLocalOutputSink<'a, K> {
    fn write(
        &self,
        mut output: LocalContractionOutput<K>,
        unitig_bucket: u16,
    ) -> Result<(), DiscontinuityInputError> {
        let io_started = Instant::now();
        if !output.trivial_fasta.is_empty() {
            let offset = self
                .next_trivial_byte
                .fetch_add(output.trivial_fasta.len() as u64, Ordering::Relaxed);
            self.trivial_output
                .write_all_at(&output.trivial_fasta, offset)
                .map_err(|source| DiscontinuityInputError::Io {
                    path: self.trivial_path.to_path_buf(),
                    source,
                })?;
            self.trivial_unitigs
                .fetch_add(output.trivial_unitigs, Ordering::Relaxed);
            self.trivial_bases
                .fetch_add(output.trivial_bases, Ordering::Relaxed);
        }
        let output_unitigs = output.unitigs.len();
        let unitig_base = self
            .next_unitig
            .fetch_add(output_unitigs, Ordering::Relaxed);
        let label_base = self
            .next_label
            .fetch_add(output.labels.len() as u64, Ordering::Relaxed);

        let bucket_local = unitig_bucket != 0;
        if !bucket_local {
            self.labels
                .write_all_at(&output.labels, label_base)
                .map_err(|source| DiscontinuityInputError::Io {
                    path: self.label_path.to_path_buf(),
                    source,
                })?;
        }

        let mut encoded_unitigs = Vec::with_capacity(
            output_unitigs * external_unitig_record_len::<K>(self.compact_unitigs),
        );
        let mut exits = 0u64;
        let mut bases = 0u64;
        for unitig in &mut output.unitigs {
            exits +=
                u64::from(unitig.left_exit().is_some()) + u64::from(unitig.right_exit().is_some());
            bases += u64::from(unitig.label_len);
            if !bucket_local {
                unitig.label_start += label_base;
                append_encoded_discontinuity_unitig_record(
                    &mut encoded_unitigs,
                    unitig,
                    self.compact_unitigs,
                );
            }
        }
        let unitig_offset = unitig_base
            .checked_mul(external_unitig_record_len::<K>(self.compact_unitigs))
            .ok_or_else(|| DiscontinuityInputError::Io {
                path: self.unitig_path.to_path_buf(),
                source: std::io::Error::from(std::io::ErrorKind::FileTooLarge),
            })? as u64;
        if !bucket_local {
            self.unitigs
                .write_all_at(&encoded_unitigs, unitig_offset)
                .map_err(|source| DiscontinuityInputError::Io {
                    path: self.unitig_path.to_path_buf(),
                    source,
                })?;
        }
        self.sink_io_nanos.fetch_add(
            io_started.elapsed().as_nanos().min(u128::from(u64::MAX)) as u64,
            Ordering::Relaxed,
        );

        let color_started = Instant::now();
        let resolved_colors = if self.color_repository.is_some() {
            let runs = output
                .color_runs
                .take()
                .ok_or(DiscontinuityInputError::MissingColorRuns)?;
            if runs.len() != output_unitigs {
                return Err(DiscontinuityInputError::MissingColorRuns);
            }
            Some(runs)
        } else {
            None
        };
        let color_start = if let Some(writer) = self.color_runs {
            writer.write_unitigs(
                resolved_colors
                    .as_deref()
                    .ok_or(DiscontinuityInputError::MissingColorRuns)?,
            )?
        } else {
            0
        };
        let local_unitig_base = if unitig_bucket != 0 {
            let writers = self
                .local_unitig_writers
                .ok_or(DiscontinuityInputError::WorkerPanic)?;
            let writer = writers
                .get(unitig_bucket as usize - 1)
                .ok_or(DiscontinuityInputError::WorkerPanic)?;
            writer
                .lock()
                .map_err(|_| DiscontinuityInputError::WorkerPanic)?
                .write::<K>(&output.labels, &output.unitigs, resolved_colors.as_deref())?
        } else {
            unitig_base
        };
        self.sink_color_nanos.fetch_add(
            color_started.elapsed().as_nanos().min(u128::from(u64::MAX)) as u64,
            Ordering::Relaxed,
        );

        let edge_started = Instant::now();
        self.edge_writers
            .add_prepared_edges::<K>(&mut output.matrix_edges, local_unitig_base, unitig_bucket)
            .map_err(serial_collation_to_input_error)?;
        if output_unitigs != 0 && !bucket_local {
            self.ranges
                .lock()
                .map_err(|_| DiscontinuityInputError::WorkerPanic)?
                .push(ExternalLocalUnitigRange {
                    start_unitig: unitig_base,
                    unitigs: output_unitigs,
                    label_start: label_base,
                    label_len: output.labels.len() as u64,
                    color_start,
                });
        }
        self.sink_edge_nanos.fetch_add(
            edge_started.elapsed().as_nanos().min(u128::from(u64::MAX)) as u64,
            Ordering::Relaxed,
        );
        self.weak_superkmers
            .fetch_add(output.weak_superkmers, Ordering::Relaxed);
        self.discontinuity_exits.fetch_add(exits, Ordering::Relaxed);
        self.unitig_bases.fetch_add(bases, Ordering::Relaxed);
        self.build_nanos.fetch_add(
            output.build_elapsed.as_nanos().min(u128::from(u64::MAX)) as u64,
            Ordering::Relaxed,
        );
        self.contract_nanos.fetch_add(
            output.contract_elapsed.as_nanos().min(u128::from(u64::MAX)) as u64,
            Ordering::Relaxed,
        );
        Ok(())
    }
}

fn append_encoded_discontinuity_unitig_record<const K: usize>(
    output: &mut Vec<u8>,
    unitig: &DiscontinuityUnitig<K>,
    compact: bool,
) {
    if compact {
        output.extend_from_slice(&unitig.label_len.to_le_bytes());
        output.push(unitig.flags);
        output.extend_from_slice(&[0; 3]);
        return;
    }
    let mut stored = MaybeUninit::<DiscontinuityUnitig<K>>::zeroed();
    unsafe {
        let ptr = stored.as_mut_ptr();
        (*ptr).label_start = unitig.label_start;
        (*ptr).left_vertex = unitig.left_vertex;
        (*ptr).right_vertex = unitig.right_vertex;
        (*ptr).label_len = unitig.label_len;
        (*ptr).flags = unitig.flags;
        let stored = stored.assume_init();
        output.extend_from_slice(std::slice::from_raw_parts(
            (&stored as *const DiscontinuityUnitig<K>).cast::<u8>(),
            std::mem::size_of::<DiscontinuityUnitig<K>>(),
        ));
    }
}

#[derive(Debug, Clone)]
struct LocalBucketGroup {
    index: usize,
    graph_id: usize,
    records: u64,
    entries: Vec<BucketManifestEntry>,
}

fn contract_local_subgraphs<const K: usize>(
    entries: &[BucketManifestEntry],
    cutoff: u32,
    threads: usize,
) -> Result<Vec<LocalContractionOutput<K>>, DiscontinuityInputError> {
    let mut groups = local_bucket_groups(entries);
    groups.sort_by_key(|group| std::cmp::Reverse(group.records));
    let workers = threads.min(groups.len().max(1));
    let started = Instant::now();
    if groups.len() >= 1024 {
        eprintln!(
            "cuttlefish3-rs: contracting {} local subgraph(s) from {} bucket file(s) with {} worker(s)",
            groups.len(),
            entries.len(),
            workers
        );
    }
    if workers == 1 {
        let mut outputs = Vec::with_capacity(groups.len());
        let mut reusable_vertices = None;
        for (offset, group) in groups.iter().enumerate() {
            outputs.push(contract_local_subgraph::<K>(
                group,
                cutoff,
                None,
                &mut reusable_vertices,
            )?);
            report_local_contraction_progress(offset + 1, groups.len(), started);
        }
        return Ok(outputs);
    }

    let total_groups = groups.len();
    let next_group = AtomicUsize::new(0);
    let completed = Arc::new(AtomicUsize::new(0));
    let mut outputs = std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for _ in 0..workers {
            let completed = Arc::clone(&completed);
            let total_groups = total_groups;
            let next_group = &next_group;
            let groups = &groups;
            handles.push(scope.spawn(move || {
                let mut chunk_outputs = Vec::new();
                let mut reusable_vertices = None;
                loop {
                    let group_idx = next_group.fetch_add(1, Ordering::Relaxed);
                    let Some(group) = groups.get(group_idx) else {
                        break;
                    };
                    chunk_outputs.push(contract_local_subgraph::<K>(
                        group,
                        cutoff,
                        None,
                        &mut reusable_vertices,
                    )?);
                    let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
                    report_local_contraction_progress(done, total_groups, started);
                }
                Ok::<_, DiscontinuityInputError>(chunk_outputs)
            }));
        }

        let mut outputs = Vec::with_capacity(entries.len());
        for handle in handles {
            outputs.extend(
                handle
                    .join()
                    .map_err(|_| DiscontinuityInputError::WorkerPanic)??,
            );
        }
        Ok::<_, DiscontinuityInputError>(outputs)
    })?;

    outputs.sort_by_key(|output| output.index);
    let build_elapsed = outputs
        .iter()
        .map(|output| output.build_elapsed)
        .sum::<Duration>();
    let contract_elapsed = outputs
        .iter()
        .map(|output| output.contract_elapsed)
        .sum::<Duration>();
    eprintln!(
        "cuttlefish3-rs: local worker time: bucket read/build {:.3}s, unitig walk {:.3}s",
        build_elapsed.as_secs_f64(),
        contract_elapsed.as_secs_f64()
    );
    Ok(outputs)
}

fn local_bucket_groups(entries: &[BucketManifestEntry]) -> Vec<LocalBucketGroup> {
    let mut by_graph = BTreeMap::<usize, Vec<BucketManifestEntry>>::new();
    for entry in entries {
        by_graph
            .entry(entry.graph_id)
            .or_default()
            .push(entry.clone());
    }

    by_graph
        .into_iter()
        .enumerate()
        .map(|(index, (graph_id, entries))| {
            let records = entries.iter().map(|entry| entry.records).sum();
            LocalBucketGroup {
                index,
                graph_id,
                records,
                entries,
            }
        })
        .collect()
}

fn report_local_contraction_progress(done: usize, total: usize, started: Instant) {
    if total < 1024 {
        return;
    }
    if done == total || done % 1024 == 0 {
        eprintln!(
            "cuttlefish3-rs: contracted {done}/{total} local subgraph bucket(s) in {:.1}s",
            started.elapsed().as_secs_f64()
        );
        report_process_memory(&format!("local contraction progress {done}/{total}"));
    }
}

fn report_discontinuity_contraction_progress(done: usize, total: usize, started: Instant) {
    if total < 16 {
        return;
    }
    if done == total || done % 16 == 0 {
        eprintln!(
            "cuttlefish3-rs: contracted {done}/{total} discontinuity partition(s) in {:.1}s",
            started.elapsed().as_secs_f64()
        );
    }
}

fn contract_local_subgraph<const K: usize>(
    group: &LocalBucketGroup,
    cutoff: u32,
    color_repository: Option<&ConcurrentColorRepository>,
    reusable_vertices: &mut Option<LocalVertexMap<K>>,
) -> Result<LocalContractionOutput<K>, DiscontinuityInputError> {
    let build_start = Instant::now();
    let mut subgraph = LocalSubgraph::<K>::from_manifest_entries_reusing(
        &group.entries,
        cutoff,
        reusable_vertices.take(),
    )?;
    let build_elapsed = build_start.elapsed();
    let weak_superkmers = subgraph.stats.weak_superkmers;
    let graph_id = subgraph.graph_id;
    debug_assert_eq!(graph_id, group.graph_id);
    let mut inputs = DiscontinuityInputs::empty(DiscontinuityInputStats::default());

    let contract_start = Instant::now();
    let mut color_runs = None;
    let local_unitigs = if subgraph.colored {
        let repository = color_repository.ok_or(DiscontinuityInputError::MissingColorRuns)?;
        let (unitigs, runs) = subgraph.contract_colored_resolved(
            &group.entries,
            repository,
            group.index % repository.worker_count(),
        )?;
        color_runs = Some(runs);
        unitigs
    } else {
        subgraph.contract_compact()?
    };
    let mut trivial_fasta = Vec::new();
    let mut trivial_unitigs = 0u64;
    let mut trivial_bases = 0u64;
    for unitig in local_unitigs {
        let left_exit = unitig
            .left_exit
            .map(|(vertex, side)| DiscontinuityEndpoint { vertex, side });
        let right_exit = unitig
            .right_exit
            .map(|(vertex, side)| DiscontinuityEndpoint { vertex, side });

        if !subgraph.colored && left_exit.is_none() && right_exit.is_none() {
            let label = canonical_label(unitig.label);
            trivial_fasta.extend_from_slice(b">0\n");
            trivial_fasta.extend_from_slice(&label);
            trivial_fasta.push(b'\n');
            trivial_unitigs += 1;
            trivial_bases += label.len() as u64;
            continue;
        }

        inputs.push_unitig(OwnedDiscontinuityUnitig {
            graph_id,
            label: unitig.label,
            left_exit,
            right_exit,
            is_cycle: unitig.is_cycle,
        });
    }
    let contract_elapsed = contract_start.elapsed();
    let matrix_edges = inputs
        .unitigs
        .iter()
        .enumerate()
        .filter_map(|(unitig_index, unitig)| {
            prepare_unitig_blocked_edge(unitig, unitig_index, DEFAULT_VERTEX_PARTITIONS)
        })
        .collect();

    let output = LocalContractionOutput {
        index: group.index,
        weak_superkmers,
        build_elapsed,
        contract_elapsed,
        unitigs: inputs.unitigs,
        labels: inputs.labels,
        matrix_edges,
        color_runs,
        trivial_fasta,
        trivial_unitigs,
        trivial_bases,
    };
    *reusable_vertices = Some(subgraph.into_vertex_map());
    if !keep_intermediates() {
        for entry in &group.entries {
            match fs::remove_file(&entry.path) {
                Ok(()) => {}
                Err(source) if source.kind() == std::io::ErrorKind::NotFound => {}
                Err(source) => {
                    return Err(DiscontinuityInputError::Io {
                        path: entry.path.clone(),
                        source,
                    });
                }
            }
        }
    }
    Ok(output)
}

#[derive(Debug)]
pub enum DiscontinuityInputError {
    Bucket(BucketError),
    LocalSubgraph(LocalSubgraphError),
    Color(ColorError),
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    InvalidCutoff,
    InvalidThreadCount,
    WorkerPanic,
    MissingColorRuns,
}

impl From<BucketError> for DiscontinuityInputError {
    fn from(value: BucketError) -> Self {
        Self::Bucket(value)
    }
}

impl From<LocalSubgraphError> for DiscontinuityInputError {
    fn from(value: LocalSubgraphError) -> Self {
        Self::LocalSubgraph(value)
    }
}

impl From<ColorError> for DiscontinuityInputError {
    fn from(value: ColorError) -> Self {
        Self::Color(value)
    }
}

impl std::fmt::Display for DiscontinuityInputError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Bucket(err) => write!(f, "{err}"),
            Self::LocalSubgraph(err) => write!(f, "{err}"),
            Self::Color(err) => write!(f, "{err}"),
            Self::Io { path, source } => write!(f, "{}: {source}", path.display()),
            Self::InvalidCutoff => write!(f, "discontinuity input cutoff must be at least 1"),
            Self::InvalidThreadCount => {
                write!(f, "discontinuity input thread count must be at least 1")
            }
            Self::WorkerPanic => write!(f, "local subgraph worker thread panicked"),
            Self::MissingColorRuns => write!(f, "colored local unitig is missing color runs"),
        }
    }
}

impl std::error::Error for DiscontinuityInputError {}

#[cfg(test)]
mod materialized_record_tests {
    use super::*;

    #[test]
    fn compact_discontinuity_edge_uses_cpp_widths() {
        let edge = DiscontinuityEdge::<31> {
            first: MatrixEndpoint::Vertex(DiscontinuityEndpoint {
                vertex: Kmer::from_bits(0x1234),
                side: Side::Back,
            }),
            second: MatrixEndpoint::Vertex(DiscontinuityEndpoint {
                vertex: Kmer::from_bits(0x5678),
                side: Side::Front,
            }),
            weight: 65_534,
            unitig_bucket: 1_023,
            unitig_index: u32::MAX as usize - 1,
            unitig_exit_side: Side::Back,
            phantom_unitig: None,
            swapped: true,
        };

        assert_eq!(discontinuity_edge_record_len::<31>(), 25);
        assert_eq!(blocked_edge_unitig_offset::<31>(), 18);
        let mut encoded = encode_discontinuity_edge(&edge);
        assert_eq!(decode_discontinuity_edge::<31>(&encoded[..25]), edge);

        set_encoded_edge_unitig_bucket::<31>(&mut encoded, 511);
        assert_eq!(
            decode_discontinuity_edge::<31>(&encoded[..25]).unitig_bucket,
            511
        );
    }

    #[test]
    fn shared_materialized_batch_matches_cpp_buffer_thresholds() {
        let mut bucket = PendingMaterializedBucket::default();
        bucket.records.resize_with(
            SharedMaterializedBatch::FLUSH_BYTES / STITCH_COORD_RECORD_LEN as usize,
            || LoadedMaterializedStitchedCoordRecord::new(0, 0, 0, 31, false, false, u32::MAX, 0),
        );
        assert!(!SharedMaterializedBatch::uncolored_bucket_ready(&bucket));

        bucket
            .labels
            .resize(SharedMaterializedBatch::FLUSH_BYTES, 0);
        assert!(SharedMaterializedBatch::uncolored_bucket_ready(&bucket));
        assert!(!SharedMaterializedBatch::colored_bucket_ready(&bucket));

        bucket.colors.resize(
            SharedMaterializedBatch::FLUSH_BYTES / std::mem::size_of::<UnitigColor>(),
            UnitigColor::new(0, crate::state::ColorCoordinate::from_u40(0)),
        );
        assert!(SharedMaterializedBatch::colored_bucket_ready(&bucket));
    }

    #[test]
    fn stitched_path_info_record_uses_compact_external_layout() {
        let record = StitchedCoordRecord {
            path_id: 17,
            rank: 23,
            unitig_index: u32::MAX - 1,
            reverse: true,
            is_cycle: true,
        };
        let encoded = encoded_stitched_coord_record(record);
        assert_eq!(encoded.len(), 24);
        assert_eq!(
            decoded_stitched_coord_record(&encoded, Path::new("test.scb")).unwrap(),
            record
        );
    }

    #[test]
    fn in_place_reverse_complement_matches_allocating_version() {
        for label in [b"".as_slice(), b"A", b"AC", b"ACG", b"AACGTT"] {
            let mut in_place = label.to_vec();
            reverse_complement_label_in_place(&mut in_place);
            assert_eq!(in_place, reverse_complement_label(label));
        }
    }

    #[test]
    fn compact_path_info_table_survives_collisions_and_generation_wraps() {
        let table = CompactPathInfoTable::with_max_entries(6, 31);
        assert_eq!(std::mem::size_of::<CompactPathInfoSlot>(), 32);
        assert_eq!(table.slots.len(), 16);

        let mut first_by_bucket = vec![None; table.slots.len()];
        let (first, second) = (0u64..)
            .find_map(|key| {
                let bucket = wyhash_u64(key, 0) as usize & table.mask;
                if let Some(first) = first_by_bucket[bucket] {
                    Some((first, key))
                } else {
                    first_by_bucket[bucket] = Some(key);
                    None
                }
            })
            .unwrap();

        for generation in 0..300u64 {
            table.clear();
            assert!(table.get::<31>(Kmer::from_bits(first as u128)).is_none());
            let first_info = PathInfo {
                path_id: Kmer::from_bits((100 + generation) as u128),
                rank: 10 + generation,
                exit_side: Side::Front,
                is_cycle: false,
            };
            let second_info = PathInfo {
                path_id: Kmer::from_bits((200 + generation) as u128),
                rank: 20 + generation,
                exit_side: Side::Back,
                is_cycle: true,
            };
            assert!(table.insert::<31>(Kmer::from_bits(first as u128), first_info));
            assert!(table.insert::<31>(Kmer::from_bits(second as u128), second_info));
            assert!(!table.insert::<31>(Kmer::from_bits(first as u128), first_info));
            assert_eq!(
                table.get::<31>(Kmer::from_bits(first as u128)),
                Some(first_info)
            );
            assert_eq!(
                table.get::<31>(Kmer::from_bits(second as u128)),
                Some(second_info)
            );
        }
    }

    #[test]
    fn materialized_color_index_uses_thirty_bits() {
        for color_index in [0, 0xff_ffff, 0x100_0000, 0x3fff_fffe, u32::MAX] {
            let record = MaterializedStitchedCoordRecord {
                path_id: 17,
                rank: 23,
                label_offset: 41,
                label_len: 59,
                reverse: true,
                is_cycle: true,
                color_index,
                color_count: 61,
            };
            let encoded = encoded_materialized_stitched_coord_record(record);
            assert_eq!(encoded.len(), 24);
            let decoded = decoded_materialized_stitched_coord_record(
                &encoded,
                Path::new("materialized-record-test"),
            )
            .unwrap();
            assert_eq!(decoded, record);
        }
    }

    #[test]
    fn appending_materialized_shards_rebases_labels_and_colors() {
        let directory = std::env::temp_dir().join(format!(
            "cf3-materialized-append-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&directory).unwrap();

        let make_entry = |worker_id, path_id, label: &[u8], color_offset| {
            let mut writer =
                MaterializedStitchedCoordShardWriter::create(&directory, worker_id, 7).unwrap();
            writer
                .write_colored_record(
                    &StitchedCoordRecord {
                        path_id,
                        rank: worker_id as u64,
                        unitig_index: 0,
                        reverse: false,
                        is_cycle: false,
                    },
                    label,
                    &[UnitigColor::new(
                        color_offset,
                        crate::state::ColorCoordinate::from_u40(path_id),
                    )],
                )
                .unwrap();
            writer.finish().unwrap()
        };

        let first = make_entry(0, 11, b"ACGT", 3);
        let second = make_entry(1, 29, b"TTAA", 5);
        assert_eq!(
            std::fs::metadata(first.color_path.as_ref().unwrap())
                .unwrap()
                .len(),
            std::mem::size_of::<UnitigColor>() as u64
        );
        assert_eq!(
            std::fs::metadata(second.color_path.as_ref().unwrap())
                .unwrap()
                .len(),
            std::mem::size_of::<UnitigColor>() as u64
        );
        let mut bucket = MaterializedStitchedCoordBucket::default();
        append_materialized_stitched_coord_bucket_file(&first, &mut bucket).unwrap();
        append_materialized_stitched_coord_bucket_file(&second, &mut bucket).unwrap();

        assert_eq!(bucket.records.len(), 2);
        assert_eq!(bucket.labels, b"ACGTTTAA");
        assert_eq!(bucket.records[0].label_offset, 0);
        assert_eq!(bucket.records[1].label_offset, 4);
        assert_eq!(bucket.records[0].color_start, 0);
        assert_eq!(bucket.records[1].color_start, 1);
        assert_eq!(bucket.records[0].color_count(), 1);
        assert_eq!(bucket.records[1].color_count(), 1);
        assert_eq!(
            bucket.colors[0].raw(),
            UnitigColor::new(3, crate::state::ColorCoordinate::from_u40(11),).raw()
        );
        assert_eq!(
            bucket.colors[1].raw(),
            UnitigColor::new(5, crate::state::ColorCoordinate::from_u40(29),).raw()
        );

        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn paged_external_range_index_matches_binary_search() {
        let mut ranges = Vec::new();
        let mut start = 0;
        for (index, unitigs) in [1, 17, 65_535, 2, 131_073, 4096, 70_000]
            .into_iter()
            .enumerate()
        {
            ranges.push(ExternalLocalUnitigRange {
                start_unitig: start,
                unitigs,
                label_start: index as u64,
                label_len: 0,
                color_start: 0,
            });
            start += unitigs;
        }
        let index = ExternalRangeIndex::new(&ranges);
        for unitig in 0..start {
            assert_eq!(
                index.find(&ranges, unitig),
                external_range_id_for_unitig(&ranges, unitig)
            );
        }
        assert_eq!(index.find(&ranges, start), None);
    }
}

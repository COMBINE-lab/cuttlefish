use crate::buckets::{BucketEmitStats, SharedBucketEmitter, SharedBucketSink};
use crate::dna::{ascii_base_bits, valid_ascii_base_bits};
use crate::hash::hash_u64;
use crate::input::{
    BorrowedSequenceFragment, InputError, expand_input_paths, parse_fragments,
    parse_fragments_borrowed,
};
use crate::params::{BuildParams, ParamError};
use std::path::{Path, PathBuf};
use std::sync::{Arc, mpsc};
use std::time::{Duration, Instant};

const PARTITION_FRAGMENT_BATCH: usize = 16 * 1024;
const PARTITION_BATCH_BASES: usize = 8 * 1024 * 1024;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WeakSuperKmer {
    pub graph_id: usize,
    pub offset: usize,
    pub len: usize,
    pub source_id: Option<u32>,
    pub left_discontinuous: bool,
    pub right_discontinuous: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PartitionStats {
    pub input_files: usize,
    pub records: u64,
    pub fragments: u64,
    pub fragment_bases: u64,
    pub weak_superkmers: u64,
    pub weak_superkmer_bases: u64,
    pub graph_histogram: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PartitionEmissionStats {
    pub partition: PartitionStats,
    pub buckets: BucketEmitStats,
    pub parse_elapsed: Duration,
    pub worker_elapsed: Duration,
    pub bucket_flush_elapsed: Duration,
    pub bucket_flushes: u64,
    pub bucket_finish_elapsed: Duration,
}

struct PartitionWorkerOutput {
    stats: PartitionStats,
    elapsed: Duration,
}

struct PartitionFragmentBatch {
    fragments: Vec<PartitionBatchFragment>,
    bases: Vec<u8>,
}

struct PartitionBatchFragment {
    source_id: u32,
    offset: usize,
    len: usize,
}

impl PartitionStats {
    pub fn new(graph_count: usize) -> Self {
        Self {
            input_files: 0,
            records: 0,
            fragments: 0,
            fragment_bases: 0,
            weak_superkmers: 0,
            weak_superkmer_bases: 0,
            graph_histogram: vec![0; graph_count],
        }
    }

    #[inline]
    pub fn non_empty_graphs(&self) -> usize {
        self.graph_histogram
            .iter()
            .filter(|&&count| count > 0)
            .count()
    }

    #[inline]
    pub fn max_graph_superkmers(&self) -> u64 {
        self.graph_histogram.iter().copied().max().unwrap_or(0)
    }

    fn merge_from(&mut self, other: &Self) {
        self.input_files += other.input_files;
        self.records += other.records;
        self.fragments += other.fragments;
        self.fragment_bases += other.fragment_bases;
        self.weak_superkmers += other.weak_superkmers;
        self.weak_superkmer_bases += other.weak_superkmer_bases;
        for (dst, src) in self
            .graph_histogram
            .iter_mut()
            .zip(other.graph_histogram.iter())
        {
            *dst += *src;
        }
    }
}

pub fn emit_weak_superkmer_buckets<const K: usize>(
    params: &BuildParams,
    graph_count: usize,
) -> Result<PartitionEmissionStats, PartitionRunError> {
    params.validate()?;

    let paths = expand_input_paths(params)?;
    if params.color {
        return emit_colored_weak_superkmer_buckets::<K>(params, graph_count, &paths);
    }
    let mut stats = PartitionStats::new(graph_count);
    let workers = params.threads.max(1);
    let mut batch_txs = Vec::with_capacity(workers);
    let mut batch_rxs = Vec::with_capacity(workers);
    for _ in 0..workers {
        let (tx, rx) = mpsc::sync_channel::<PartitionFragmentBatch>(4);
        batch_txs.push(tx);
        batch_rxs.push(rx);
    }
    let sink = SharedBucketSink::create(params, graph_count)?;
    let mut worker_elapsed = Duration::ZERO;
    let parse_started = Instant::now();
    let mut parse_elapsed = Duration::ZERO;

    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for rx in batch_rxs {
            let params = params;
            let sink = Arc::clone(&sink);
            handles.push(scope.spawn(move || {
                let mut worker_stats = PartitionStats::new(graph_count);
                let mut worker_elapsed = Duration::ZERO;
                let mut buckets = sink.emitter();
                loop {
                    let batch = rx.recv();
                    let Ok(batch) = batch else {
                        break;
                    };
                    for fragment in &batch.fragments {
                        let seq = &batch.bases[fragment.offset..fragment.offset + fragment.len];
                        let fragment_started = Instant::now();
                        emit_fragment_seq_weak_superkmer_buckets::<K>(
                            params,
                            graph_count,
                            fragment.source_id,
                            seq,
                            &mut buckets,
                            &mut worker_stats,
                        )?;
                        worker_elapsed += fragment_started.elapsed();
                    }
                }
                buckets.finish()?;
                Ok::<_, PartitionRunError>(PartitionWorkerOutput {
                    stats: worker_stats,
                    elapsed: worker_elapsed,
                })
            }));
        }

        let mut producer_handles = Vec::new();
        for (path_idx, path) in paths.iter().enumerate() {
            let source_id =
                u32::try_from(path_idx + 1).map_err(|_| PartitionRunError::TooManySources)?;
            let txs = batch_txs.clone();
            producer_handles.push(scope.spawn(move || {
                let mut batch = PartitionFragmentBatch::new();
                let mut next_worker = path_idx % txs.len();
                let records = parse_fragments_borrowed(path, source_id, K + 1, |fragment| {
                    batch.push(fragment);
                    if batch.fragments.len() >= PARTITION_FRAGMENT_BATCH
                        || batch.bases.len() >= PARTITION_BATCH_BASES
                    {
                        txs[next_worker]
                            .send(std::mem::take(&mut batch))
                            .map_err(|_| {
                                InputError::Partition(PartitionError::WorkerDisconnected)
                            })?;
                        next_worker = (next_worker + 1) % txs.len();
                    }
                    Ok(())
                })?;
                if !batch.fragments.is_empty() {
                    txs[next_worker]
                        .send(batch)
                        .map_err(|_| InputError::Partition(PartitionError::WorkerDisconnected))?;
                }
                Ok::<_, InputError>(records)
            }));
        }
        drop(batch_txs);
        for handle in producer_handles {
            let records = handle
                .join()
                .map_err(|_| PartitionRunError::WorkerPanic)??;
            stats.input_files += 1;
            stats.records += records;
        }
        parse_elapsed = parse_started.elapsed();

        for handle in handles {
            let worker_stats = handle
                .join()
                .map_err(|_| PartitionRunError::WorkerPanic)??;
            worker_elapsed += worker_stats.elapsed;
            stats.merge_from(&worker_stats.stats);
        }
        Ok::<(), PartitionRunError>(())
    })?;

    let (bucket_flushes, bucket_flush_elapsed) = sink.flush_stats();
    let bucket_finish_started = Instant::now();
    let bucket_stats = sink.finish()?;
    let bucket_finish_elapsed = bucket_finish_started.elapsed();

    Ok(PartitionEmissionStats {
        partition: stats,
        buckets: bucket_stats,
        parse_elapsed,
        worker_elapsed,
        bucket_flush_elapsed,
        bucket_flushes,
        bucket_finish_elapsed,
    })
}

fn emit_colored_weak_superkmer_buckets<const K: usize>(
    params: &BuildParams,
    graph_count: usize,
    paths: &[PathBuf],
) -> Result<PartitionEmissionStats, PartitionRunError> {
    const SOURCE_WINDOW_MAX: usize = 32;

    let sink = SharedBucketSink::create(params, graph_count)?;
    let mut stats = PartitionStats::new(graph_count);
    let mut worker_elapsed = Duration::ZERO;
    let mut parse_elapsed = Duration::ZERO;

    for (window_index, window) in paths
        .chunks(params.threads.max(1).min(SOURCE_WINDOW_MAX))
        .enumerate()
    {
        let source_start = window_index * params.threads.max(1).min(SOURCE_WINDOW_MAX);
        let workers = params.threads.max(1).min(window.len().max(1));
        let mut batch_txs = Vec::with_capacity(workers);
        let mut batch_rxs = Vec::with_capacity(workers);
        for _ in 0..workers {
            let (tx, rx) = mpsc::sync_channel::<PartitionFragmentBatch>(4);
            batch_txs.push(tx);
            batch_rxs.push(rx);
        }

        let parse_started = Instant::now();
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for rx in batch_rxs {
                let sink = Arc::clone(&sink);
                handles.push(scope.spawn(move || {
                    let mut worker_stats = PartitionStats::new(graph_count);
                    let mut elapsed = Duration::ZERO;
                    let mut buckets = sink.emitter();
                    while let Ok(batch) = rx.recv() {
                        for fragment in &batch.fragments {
                            let seq = &batch.bases[fragment.offset..fragment.offset + fragment.len];
                            let started = Instant::now();
                            emit_fragment_seq_weak_superkmer_buckets::<K>(
                                params,
                                graph_count,
                                fragment.source_id,
                                seq,
                                &mut buckets,
                                &mut worker_stats,
                            )?;
                            elapsed += started.elapsed();
                        }
                    }
                    buckets.finish()?;
                    Ok::<_, PartitionRunError>(PartitionWorkerOutput {
                        stats: worker_stats,
                        elapsed,
                    })
                }));
            }

            let mut producers = Vec::new();
            for (offset, path) in window.iter().enumerate() {
                let source_id = u32::try_from(source_start + offset + 1)
                    .map_err(|_| PartitionRunError::TooManySources)?;
                let txs = batch_txs.clone();
                producers.push(scope.spawn(move || {
                    let mut batch = PartitionFragmentBatch::new();
                    let mut next_worker = offset % txs.len();
                    let records =
                        parse_fragments_borrowed(path, source_id, K + 1, |fragment| {
                            batch.push(fragment);
                            if batch.fragments.len() >= PARTITION_FRAGMENT_BATCH
                                || batch.bases.len() >= PARTITION_BATCH_BASES
                            {
                                txs[next_worker].send(std::mem::take(&mut batch)).map_err(
                                    |_| InputError::Partition(PartitionError::WorkerDisconnected),
                                )?;
                                next_worker = (next_worker + 1) % txs.len();
                            }
                            Ok(())
                        })?;
                    if !batch.fragments.is_empty() {
                        txs[next_worker].send(batch).map_err(|_| {
                            InputError::Partition(PartitionError::WorkerDisconnected)
                        })?;
                    }
                    Ok::<_, InputError>(records)
                }));
            }
            drop(batch_txs);

            for producer in producers {
                stats.records += producer
                    .join()
                    .map_err(|_| PartitionRunError::WorkerPanic)??;
                stats.input_files += 1;
            }
            for handle in handles {
                let output = handle
                    .join()
                    .map_err(|_| PartitionRunError::WorkerPanic)??;
                worker_elapsed += output.elapsed;
                stats.merge_from(&output.stats);
            }
            Ok::<_, PartitionRunError>(())
        })?;
        parse_elapsed += parse_started.elapsed();

        let source_min =
            u32::try_from(source_start + 1).map_err(|_| PartitionRunError::TooManySources)?;
        let source_max = u32::try_from(source_start + window.len())
            .map_err(|_| PartitionRunError::TooManySources)?;
        sink.flush_colored_window(source_min, source_max)?;
    }

    let (bucket_flushes, bucket_flush_elapsed) = sink.flush_stats();
    let bucket_finish_started = Instant::now();
    let bucket_stats = sink.finish()?;
    let bucket_finish_elapsed = bucket_finish_started.elapsed();
    Ok(PartitionEmissionStats {
        partition: stats,
        buckets: bucket_stats,
        parse_elapsed,
        worker_elapsed,
        bucket_flush_elapsed,
        bucket_flushes,
        bucket_finish_elapsed,
    })
}

impl PartitionFragmentBatch {
    fn new() -> Self {
        Self {
            fragments: Vec::with_capacity(PARTITION_FRAGMENT_BATCH),
            bases: Vec::with_capacity(PARTITION_BATCH_BASES.min(PARTITION_FRAGMENT_BATCH * 256)),
        }
    }

    fn push(&mut self, fragment: BorrowedSequenceFragment<'_>) {
        let offset = self.bases.len();
        self.bases.extend_from_slice(fragment.seq);
        self.fragments.push(PartitionBatchFragment {
            source_id: fragment.source_id,
            offset,
            len: fragment.seq.len(),
        });
    }
}

impl Default for PartitionFragmentBatch {
    fn default() -> Self {
        Self::new()
    }
}

fn emit_fragment_seq_weak_superkmer_buckets<const K: usize>(
    params: &BuildParams,
    graph_count: usize,
    source_id: u32,
    seq: &[u8],
    buckets: &mut SharedBucketEmitter,
    stats: &mut PartitionStats,
) -> Result<(), PartitionRunError> {
    stats.fragments += 1;
    stats.fragment_bases += seq.len() as u64;
    for_each_valid_weak_superkmer::<K, PartitionRunError, _>(
        seq,
        params.minimizer_len as usize,
        graph_count,
        params.color.then_some(source_id),
        |sk| {
            buckets.add_valid(&sk, sk.sequence(seq))?;
            stats.weak_superkmers += 1;
            stats.weak_superkmer_bases += sk.len as u64;
            stats.graph_histogram[sk.graph_id] += 1;
            Ok(())
        },
    )?;
    Ok(())
}

pub fn partition_inputs<const K: usize>(
    params: &BuildParams,
    graph_count: usize,
) -> Result<PartitionStats, PartitionRunError> {
    params.validate()?;

    let paths = expand_input_paths(params)?;
    let mut stats = PartitionStats::new(graph_count);
    let workers = params.threads.min(paths.len().max(1));
    let chunk_len = paths.len().div_ceil(workers);

    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for (chunk_idx, chunk) in paths.chunks(chunk_len).enumerate() {
            let params = params;
            let start_idx = chunk_idx * chunk_len;
            handles.push(scope.spawn(move || {
                let mut worker_stats = PartitionStats::new(graph_count);
                for (offset, path) in chunk.iter().enumerate() {
                    let source_id = u32::try_from(start_idx + offset + 1)
                        .map_err(|_| PartitionRunError::TooManySources)?;
                    let path_stats = partition_path::<K>(params, graph_count, source_id, path)?;
                    worker_stats.merge_from(&path_stats);
                }
                Ok::<PartitionStats, PartitionRunError>(worker_stats)
            }));
        }

        for handle in handles {
            let worker_stats = handle
                .join()
                .map_err(|_| PartitionRunError::WorkerPanic)??;
            stats.merge_from(&worker_stats);
        }
        Ok::<(), PartitionRunError>(())
    })?;

    Ok(stats)
}

fn partition_path<const K: usize>(
    params: &BuildParams,
    graph_count: usize,
    source_id: u32,
    path: &Path,
) -> Result<PartitionStats, PartitionRunError> {
    let mut stats = PartitionStats::new(graph_count);
    stats.input_files = 1;

    let records = parse_fragments(path, source_id, K + 1, |fragment| {
        stats.fragments += 1;
        stats.fragment_bases += fragment.seq.len() as u64;
        for_each_valid_weak_superkmer::<K, InputError, _>(
            &fragment.seq,
            params.minimizer_len as usize,
            graph_count,
            params.color.then_some(fragment.source_id),
            |sk| {
                stats.weak_superkmers += 1;
                stats.weak_superkmer_bases += sk.len as u64;
                stats.graph_histogram[sk.graph_id] += 1;
                Ok(())
            },
        )?;
        Ok(())
    })?;
    stats.records = records;

    Ok(stats)
}

impl WeakSuperKmer {
    #[inline]
    pub fn sequence<'a>(&self, fragment: &'a [u8]) -> &'a [u8] {
        &fragment[self.offset..self.offset + self.len]
    }
}

pub fn partition_fragment<const K: usize>(
    fragment: &[u8],
    minimizer_len: usize,
    graph_count: usize,
    source_id: Option<u32>,
) -> Result<Vec<WeakSuperKmer>, PartitionError> {
    let mut out = Vec::new();
    for_each_weak_superkmer::<K, PartitionError, _>(
        fragment,
        minimizer_len,
        graph_count,
        source_id,
        |sk| {
            out.push(sk);
            Ok(())
        },
    )?;
    Ok(out)
}

fn for_each_weak_superkmer<const K: usize, E, F>(
    fragment: &[u8],
    minimizer_len: usize,
    graph_count: usize,
    source_id: Option<u32>,
    mut visit: F,
) -> Result<(), E>
where
    E: From<PartitionError>,
    F: FnMut(WeakSuperKmer) -> Result<(), E>,
{
    for_each_weak_superkmer_impl::<K, E, _, true>(
        fragment,
        minimizer_len,
        graph_count,
        source_id,
        &mut visit,
    )
}

fn for_each_valid_weak_superkmer<const K: usize, E, F>(
    fragment: &[u8],
    minimizer_len: usize,
    graph_count: usize,
    source_id: Option<u32>,
    mut visit: F,
) -> Result<(), E>
where
    E: From<PartitionError>,
    F: FnMut(WeakSuperKmer) -> Result<(), E>,
{
    for_each_weak_superkmer_impl::<K, E, _, false>(
        fragment,
        minimizer_len,
        graph_count,
        source_id,
        &mut visit,
    )
}

fn for_each_weak_superkmer_impl<const K: usize, E, F, const CHECK_BASES: bool>(
    fragment: &[u8],
    minimizer_len: usize,
    graph_count: usize,
    source_id: Option<u32>,
    visit: &mut F,
) -> Result<(), E>
where
    E: From<PartitionError>,
    F: FnMut(WeakSuperKmer) -> Result<(), E>,
{
    if K < 3 || K % 2 == 0 || K > 63 {
        return Err(PartitionError::InvalidK(K).into());
    }
    if minimizer_len == 0 || minimizer_len >= K || minimizer_len > 32 {
        return Err(PartitionError::InvalidMinimizerLen(minimizer_len).into());
    }
    if !graph_count.is_power_of_two() {
        return Err(PartitionError::GraphCountNotPowerOfTwo(graph_count).into());
    }
    if fragment.len() < K + 1 {
        return Ok(());
    }

    let window_k = K - 1;
    let max_super_km1_len = 2 * (K - 1) - minimizer_len;
    let mut cur_off = 0usize;
    let mut km1_idx = 0usize;
    let mut cur_g = 0usize;
    let mut prev_g = graph_count;
    let mut initialized = false;
    let mut hash_idx = 0usize;

    for_each_km1_window_hash::<K, E, _, CHECK_BASES>(fragment, minimizer_len, |next_hash| {
        if !initialized {
            cur_g = graph_id(next_hash, graph_count);
            initialized = true;
            hash_idx += 1;
            return Ok(());
        }

        let len = km1_idx + window_k;
        let next_g = graph_id(next_hash, graph_count);
        km1_idx += 1;
        if next_g != cur_g || len == max_super_km1_len {
            let next_off = cur_off + km1_idx;
            let left_joined = cur_off > 0;
            let right_joined = true;
            visit(WeakSuperKmer {
                graph_id: cur_g,
                offset: cur_off - usize::from(left_joined),
                len: usize::from(left_joined) + len + usize::from(right_joined),
                source_id,
                left_discontinuous: left_joined && prev_g != cur_g,
                right_discontinuous: next_g != cur_g,
            })?;

            debug_assert_eq!(next_off, hash_idx);
            cur_off = next_off;
            prev_g = cur_g;
            cur_g = next_g;
            km1_idx = 0;
        }
        hash_idx += 1;
        Ok(())
    })?;

    if !initialized {
        return Ok(());
    }

    let len = fragment.len() - cur_off;
    let left_joined = cur_off > 0;
    visit(WeakSuperKmer {
        graph_id: cur_g,
        offset: cur_off - usize::from(left_joined),
        len: usize::from(left_joined) + len,
        source_id,
        left_discontinuous: left_joined && prev_g != cur_g,
        right_discontinuous: false,
    })?;

    Ok(())
}

fn for_each_km1_window_hash<const K: usize, E, F, const CHECK_BASES: bool>(
    fragment: &[u8],
    minimizer_len: usize,
    mut visit: F,
) -> Result<(), E>
where
    E: From<PartitionError>,
    F: FnMut(u64) -> Result<(), E>,
{
    let window_k = K - 1;
    if fragment.len() < window_k {
        return Ok(());
    }

    let lmers_per_window = window_k - minimizer_len + 1;
    let ring_size = lmers_per_window - 1;
    let mut hashes = [0u64; K];
    let mut prefix_min = [0u64; K];

    let mut fwd = 0u64;
    let mut rev = 0u64;
    for idx in 0..minimizer_len {
        let base_bits = partition_base_bits::<CHECK_BASES, E>(fragment[idx])?;
        fwd = (fwd << 2) | base_bits as u64;
        rev |= ((base_bits ^ 0b11) as u64) << (2 * idx);
    }

    let mask = if minimizer_len == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * minimizer_len)) - 1
    };
    let rev_high_shift = 2 * (minimizer_len - 1);

    let mut pivot = ring_size;
    hashes[pivot] = canonical_lmer_hash(fwd, rev);
    for idx in minimizer_len..window_k {
        let next_bits = partition_base_bits::<CHECK_BASES, E>(fragment[idx])?;
        fwd = ((fwd << 2) | next_bits as u64) & mask;
        rev = (rev >> 2) | (((next_bits ^ 0b11) as u64) << rev_high_shift);
        pivot -= 1;
        hashes[pivot] = canonical_lmer_hash(fwd, rev);
    }

    reset_min_window(&hashes, &mut prefix_min, ring_size, &mut pivot);
    let mut suffix_min = u64::MAX;
    visit(prefix_min[pivot].min(suffix_min))?;

    for &base in &fragment[window_k..] {
        let next_bits = partition_base_bits::<CHECK_BASES, E>(base)?;
        fwd = ((fwd << 2) | next_bits as u64) & mask;
        rev = (rev >> 2) | (((next_bits ^ 0b11) as u64) << rev_high_shift);

        hashes[pivot] = canonical_lmer_hash(fwd, rev);
        suffix_min = suffix_min.min(hashes[pivot]);
        if pivot > 0 {
            pivot -= 1;
        } else {
            reset_min_window(&hashes, &mut prefix_min, ring_size, &mut pivot);
            suffix_min = u64::MAX;
        }
        visit(prefix_min[pivot].min(suffix_min))?;
    }

    Ok(())
}

#[inline(always)]
fn reset_min_window(hashes: &[u64], prefix_min: &mut [u64], ring_size: usize, pivot: &mut usize) {
    prefix_min[0] = hashes[0];
    for idx in 1..=ring_size {
        prefix_min[idx] = prefix_min[idx - 1].min(hashes[idx]);
    }
    *pivot = ring_size;
}

#[inline(always)]
fn canonical_lmer_hash(fwd: u64, rev: u64) -> u64 {
    hash_u64(fwd, 0).min(hash_u64(rev, 0))
}

#[inline(always)]
fn partition_base_bits<const CHECK_BASES: bool, E>(base: u8) -> Result<u8, E>
where
    E: From<PartitionError>,
{
    if CHECK_BASES {
        ascii_base_bits(base)
            .ok_or(PartitionError::InvalidBase(base))
            .map_err(E::from)
    } else {
        debug_assert!(ascii_base_bits(base).is_some());
        Ok(valid_ascii_base_bits(base))
    }
}

#[inline]
fn graph_id(hash: u64, graph_count: usize) -> usize {
    hash as usize & (graph_count - 1)
}

#[inline]
pub fn source_hash(source_id: u32) -> u64 {
    hash_u64(source_id as u64, 0)
}

#[derive(Debug)]
pub enum PartitionError {
    InvalidK(usize),
    InvalidMinimizerLen(usize),
    GraphCountNotPowerOfTwo(usize),
    InvalidBase(u8),
    Kmer(crate::kmer::KmerError),
    WorkerDisconnected,
}

impl From<crate::kmer::KmerError> for PartitionError {
    fn from(value: crate::kmer::KmerError) -> Self {
        Self::Kmer(value)
    }
}

impl std::fmt::Display for PartitionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidK(k) => write!(f, "invalid k for partitioning: {k}"),
            Self::InvalidMinimizerLen(l) => write!(f, "invalid minimizer length: {l}"),
            Self::GraphCountNotPowerOfTwo(c) => {
                write!(f, "graph count must be a power of two: {c}")
            }
            Self::InvalidBase(b) => {
                write!(f, "invalid base in partition fragment: '{}'", *b as char)
            }
            Self::Kmer(err) => write!(f, "{err}"),
            Self::WorkerDisconnected => write!(f, "partition worker disconnected"),
        }
    }
}

impl std::error::Error for PartitionError {}

impl From<PartitionError> for InputError {
    fn from(value: PartitionError) -> Self {
        Self::Partition(value)
    }
}

#[derive(Debug)]
pub enum PartitionRunError {
    Params(ParamError),
    Input(InputError),
    Partition(PartitionError),
    Bucket(crate::buckets::BucketError),
    TooManySources,
    WorkerPanic,
}

impl From<ParamError> for PartitionRunError {
    fn from(value: ParamError) -> Self {
        Self::Params(value)
    }
}

impl From<InputError> for PartitionRunError {
    fn from(value: InputError) -> Self {
        Self::Input(value)
    }
}

impl From<PartitionError> for PartitionRunError {
    fn from(value: PartitionError) -> Self {
        Self::Partition(value)
    }
}

impl From<crate::buckets::BucketError> for PartitionRunError {
    fn from(value: crate::buckets::BucketError) -> Self {
        Self::Bucket(value)
    }
}

impl std::fmt::Display for PartitionRunError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Params(err) => write!(f, "{err}"),
            Self::Input(err) => write!(f, "{err}"),
            Self::Partition(err) => write!(f, "{err}"),
            Self::Bucket(err) => write!(f, "{err}"),
            Self::TooManySources => write!(f, "too many input sources for 32-bit source IDs"),
            Self::WorkerPanic => write!(f, "partition worker thread panicked"),
        }
    }
}

impl std::error::Error for PartitionRunError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn emits_covering_weak_superkmers() {
        let seq = b"ACGTACGTACGTACGTACGTACGT";
        let parts = partition_fragment::<7>(seq, 3, 128, Some(1)).unwrap();
        assert!(!parts.is_empty());
        assert!(parts.iter().all(|p| p.len >= 7));
        assert_eq!(parts[0].offset, 0);
        assert_eq!(
            parts.last().unwrap().offset + parts.last().unwrap().len,
            seq.len()
        );
        assert!(parts.iter().all(|p| p.source_id == Some(1)));
    }

    #[test]
    fn rejects_non_power_of_two_graph_count() {
        let err = partition_fragment::<7>(b"ACGTACGT", 3, 127, None).unwrap_err();
        assert!(matches!(err, PartitionError::GraphCountNotPowerOfTwo(127)));
    }
}

//! Weak-super-k-mer partitioning and bucket emission.
//!
//! The production partitioner dynamically schedules size-sorted inputs, scans
//! packed rolling windows, and drains worker-local 128-subgraph atlas chunks to
//! [`crate::buckets`]. Graph IDs are stable across colored and uncolored paths.

use crate::buckets::{BucketEmitStats, SharedBucketEmitter, SharedBucketSink};
use crate::dna::{ascii_base_bits, valid_ascii_base_bits};
use crate::hash::{hash_u64, wyhash_u64};
use crate::input::{
    BorrowedSequenceFragment, InputError, expand_input_paths, parse_fragments,
    parse_fragments_borrowed, parse_fragments_borrowed_with,
};
use crate::params::{BuildParams, ParamError};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, mpsc};
use std::time::{Duration, Instant};

const PARTITION_FRAGMENT_BATCH: usize = 16 * 1024;
const PARTITION_BATCH_BASES: usize = 8 * 1024 * 1024;

/// Fragment-batch geometry for the streamed partitioner.
///
/// Batches are deliberately much smaller than [`PARTITION_BATCH_BASES`] so that
/// a bounded pool keeps every worker fed without holding gigabytes in flight.
const STREAM_BATCH_BASES: usize = 1024 * 1024;
const STREAM_BATCH_FRAGMENTS: usize = 8 * 1024;
/// Upper bound on the bytes held by the streamed partitioner's batch pool.
const STREAM_POOL_BYTES: usize = 256 * 1024 * 1024;

#[derive(Debug, Clone, PartialEq, Eq)]
/// A weak super-k-mer assigned to one local subgraph.
///
/// `offset` and `len` address a fragment supplied to [`partition_fragment`].
pub struct WeakSuperKmer {
    pub graph_id: usize,
    pub offset: usize,
    pub len: usize,
    pub source_id: Option<u32>,
    pub left_discontinuous: bool,
    pub right_discontinuous: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// Aggregate parser and weak-super-k-mer partitioning counts.
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
/// Partition counts, external bucket manifest, and phase timings.
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
        // Colored vertices retain their previous source in a packed 21-bit
        // field. Reject an oversized source set here, with the input count in
        // hand, rather than letting a contraction worker trip an assertion.
        if paths.len() > crate::state::VertexState::MAX_SOURCE_ID as usize {
            return Err(PartitionRunError::TooManySources);
        }
        return emit_colored_weak_superkmer_buckets::<K>(params, graph_count, &paths);
    }
    if std::env::var_os("CF3_RS_LEGACY_UNCOLORED_PARTITION").is_none() {
        return emit_uncolored_windowed_weak_superkmer_buckets::<K>(params, graph_count, &paths);
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
                        emit_fragment_seq_weak_superkmer_buckets::<K, PartitionRunError>(
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
    populate_emitted_graph_histogram(&mut stats, &bucket_stats)?;

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

/// Returns the reader count for the streamed partitioner, or `None` when
/// file-level parallelism already fills the requested worker count.
///
/// Reference workloads present thousands of sources, so one whole file per
/// worker saturates the machine with no batch copy at all. Read workloads
/// present a handful of very large files, where that mapping leaves nearly
/// every core idle; those inputs are streamed instead.
fn streamed_reader_count(params: &BuildParams, files: usize) -> Option<usize> {
    if std::env::var_os("CF3_RS_DIRECT_PARTITION").is_some() {
        return None;
    }
    // `usize::MAX` asks for the memory bound only, dropping the input-file clamp
    // that the direct path needs because it maps one whole file per worker.
    let workers = params.partition_workers(usize::MAX);
    if files == 0 || files >= workers {
        return None;
    }
    Some(files)
}

fn emit_uncolored_windowed_weak_superkmer_buckets<const K: usize>(
    params: &BuildParams,
    graph_count: usize,
    paths: &[PathBuf],
) -> Result<PartitionEmissionStats, PartitionRunError> {
    if let Some(readers) = streamed_reader_count(params, paths.len()) {
        return emit_uncolored_streamed_weak_superkmer_buckets::<K>(
            params,
            graph_count,
            paths,
            readers,
        );
    }
    emit_uncolored_direct_weak_superkmer_buckets::<K>(params, graph_count, paths)
}

/// Streams fragment batches from per-file readers to a full pool of workers.
///
/// Decompression and record assembly stay on the reader threads while every
/// requested worker runs the minimizer scan and bucket packing, mirroring the
/// C++ partitioner's chunk queue. Uncolored emission has no source ordering
/// requirement, so batches may be consumed in any order.
fn emit_uncolored_streamed_weak_superkmer_buckets<const K: usize>(
    params: &BuildParams,
    graph_count: usize,
    paths: &[PathBuf],
    readers: usize,
) -> Result<PartitionEmissionStats, PartitionRunError> {
    let sink = SharedBucketSink::create(params, graph_count)?;
    let mut stats = PartitionStats::new(graph_count);
    let mut worker_elapsed = Duration::ZERO;
    let workers = params.partition_workers(usize::MAX);
    let buffers = (2 * (workers + readers)).min(STREAM_POOL_BYTES / STREAM_BATCH_BASES);
    eprintln!(
        "cuttlefish3-rs: uncolored partition streaming {readers} reader(s) into {workers} worker(s), {buffers} batch buffer(s)"
    );

    // Readers decompress; only block-structured input can use more than one
    // thread for it, so the budget is shared evenly and costs nothing otherwise.
    let inflate_workers = (workers / readers.max(1)).max(1);
    let pool = BatchPool::new(buffers, readers);
    let next_source = AtomicUsize::new(0);
    let mut work_order = (0..paths.len()).collect::<Vec<_>>();
    work_order.sort_unstable_by_key(|&offset| {
        std::cmp::Reverse(paths[offset].metadata().map_or(0, |meta| meta.len()))
    });

    let started = Instant::now();
    let mut emitters = Vec::with_capacity(workers);
    std::thread::scope(|scope| {
        let mut worker_handles = Vec::with_capacity(workers);
        for _ in 0..workers {
            let sink = Arc::clone(&sink);
            let pool = &pool;
            worker_handles.push(scope.spawn(move || {
                let mut worker_stats = PartitionStats::new(graph_count);
                let mut elapsed = Duration::ZERO;
                let mut buckets = sink.deferred_uncolored_emitter();
                let result = (|| -> Result<(), PartitionRunError> {
                    while let Some(batch) = pool.take_ready() {
                        let fragment_started = Instant::now();
                        let outcome = (|| -> Result<(), PartitionRunError> {
                            for fragment in &batch.fragments {
                                let seq =
                                    &batch.bases[fragment.offset..fragment.offset + fragment.len];
                                emit_fragment_seq_weak_superkmer_buckets::<K, PartitionRunError>(
                                    params,
                                    graph_count,
                                    fragment.source_id,
                                    seq,
                                    &mut buckets,
                                    &mut worker_stats,
                                )?;
                            }
                            Ok(())
                        })();
                        elapsed += fragment_started.elapsed();
                        pool.release(batch);
                        outcome?;
                    }
                    Ok(())
                })();
                if result.is_err() {
                    pool.fail();
                }
                result.map(|()| {
                    (
                        PartitionWorkerOutput {
                            stats: worker_stats,
                            elapsed,
                        },
                        buckets,
                    )
                })
            }));
        }

        let mut reader_handles = Vec::with_capacity(readers);
        for _ in 0..readers {
            let pool = &pool;
            let next_source = &next_source;
            let work_order = &work_order;
            reader_handles.push(scope.spawn(move || {
                let result = (|| -> Result<(u64, usize), PartitionRunError> {
                    let mut records = 0u64;
                    let mut input_files = 0usize;
                    loop {
                        let order_idx = next_source.fetch_add(1, Ordering::Relaxed);
                        let Some(&offset) = work_order.get(order_idx) else {
                            break;
                        };
                        let source_id = u32::try_from(offset + 1)
                            .map_err(|_| PartitionRunError::TooManySources)?;
                        let mut batch = pool
                            .take_free()
                            .ok_or(PartitionRunError::Partition(
                                PartitionError::WorkerDisconnected,
                            ))?;
                        let parsed = parse_fragments_borrowed_with(
                            &paths[offset],
                            source_id,
                            K + 1,
                            inflate_workers,
                            |fragment| {
                                batch.push(fragment);
                                if batch.is_full() {
                                    let replacement = pool.take_free().ok_or(
                                        InputError::Partition(PartitionError::WorkerDisconnected),
                                    )?;
                                    pool.publish(std::mem::replace(&mut batch, replacement));
                                }
                                Ok(())
                            },
                        );
                        // Publish whatever this source produced before deciding
                        // how to treat a failure, so no buffer is leaked.
                        if batch.fragments.is_empty() {
                            pool.release(batch);
                        } else {
                            pool.publish(batch);
                        }
                        match parsed {
                            Ok(parsed) => {
                                records += parsed;
                                input_files += 1;
                            }
                            Err(error) => {
                                handle_source_failure(params, &paths[offset], error)?;
                            }
                        }
                    }
                    Ok((records, input_files))
                })();
                pool.reader_done();
                if result.is_err() {
                    pool.fail();
                }
                result
            }));
        }

        for handle in reader_handles {
            let (records, input_files) = handle
                .join()
                .map_err(|_| PartitionRunError::WorkerPanic)??;
            stats.records += records;
            stats.input_files += input_files;
        }
        for handle in worker_handles {
            let (output, buckets) = handle
                .join()
                .map_err(|_| PartitionRunError::WorkerPanic)??;
            worker_elapsed += output.elapsed;
            stats.merge_from(&output.stats);
            emitters.push(buckets);
        }
        Ok::<(), PartitionRunError>(())
    })?;
    let parse_elapsed = started.elapsed();
    sink.flush_uncolored_emitters(emitters)?;

    let (bucket_flushes, bucket_flush_elapsed) = sink.flush_stats();
    let bucket_finish_started = Instant::now();
    let bucket_stats = sink.finish()?;
    let bucket_finish_elapsed = bucket_finish_started.elapsed();
    populate_emitted_graph_histogram(&mut stats, &bucket_stats)?;
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

fn emit_uncolored_direct_weak_superkmer_buckets<const K: usize>(
    params: &BuildParams,
    graph_count: usize,
    paths: &[PathBuf],
) -> Result<PartitionEmissionStats, PartitionRunError> {
    let sink = SharedBucketSink::create(params, graph_count)?;
    let mut stats = PartitionStats::new(graph_count);
    let mut worker_elapsed = Duration::ZERO;
    let mut parse_elapsed = Duration::ZERO;

    for (source_start, window) in [(0, paths)] {
        let started = Instant::now();
        let mut work_order = (0..window.len()).collect::<Vec<_>>();
        work_order.sort_unstable_by_key(|&offset| {
            std::cmp::Reverse(window[offset].metadata().map_or(0, |meta| meta.len()))
        });
        let next_source = AtomicUsize::new(0);
        let workers = params.partition_workers(window.len());
        eprintln!("cuttlefish3-rs: uncolored partition using {workers} worker(s)");
        let mut emitters = Vec::with_capacity(workers);

        std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(workers);
            for _ in 0..workers {
                let sink = Arc::clone(&sink);
                let next_source = &next_source;
                let work_order = &work_order;
                handles.push(scope.spawn(move || {
                    let mut worker_stats = PartitionStats::new(graph_count);
                    let mut elapsed = Duration::ZERO;
                    let mut records = 0u64;
                    let mut input_files = 0usize;
                    let mut buckets = sink.deferred_uncolored_emitter();
                    loop {
                        let order_idx = next_source.fetch_add(1, Ordering::Relaxed);
                        let Some(&offset) = work_order.get(order_idx) else {
                            break;
                        };
                        let source_id = u32::try_from(source_start + offset + 1)
                            .map_err(|_| PartitionRunError::TooManySources)?;
                        match parse_fragments_borrowed(
                            &window[offset],
                            source_id,
                            K + 1,
                            |fragment| {
                                let fragment_started = Instant::now();
                                emit_fragment_seq_weak_superkmer_buckets::<K, InputError>(
                                    params,
                                    graph_count,
                                    fragment.source_id,
                                    fragment.seq,
                                    &mut buckets,
                                    &mut worker_stats,
                                )?;
                                elapsed += fragment_started.elapsed();
                                Ok(())
                            },
                        ) {
                            Ok(parsed) => {
                                records += parsed;
                                input_files += 1;
                            }
                            Err(error) => {
                                handle_source_failure(params, &window[offset], error)?;
                            }
                        }
                    }
                    Ok::<_, PartitionRunError>((
                        records,
                        input_files,
                        PartitionWorkerOutput {
                            stats: worker_stats,
                            elapsed,
                        },
                        buckets,
                    ))
                }));
            }
            for handle in handles {
                let (records, input_files, output, buckets) = handle
                    .join()
                    .map_err(|_| PartitionRunError::WorkerPanic)??;
                stats.records += records;
                stats.input_files += input_files;
                worker_elapsed += output.elapsed;
                stats.merge_from(&output.stats);
                emitters.push(buckets);
            }
            Ok::<_, PartitionRunError>(())
        })?;
        parse_elapsed += started.elapsed();
        sink.flush_uncolored_emitters(emitters)?;
    }

    let (bucket_flushes, bucket_flush_elapsed) = sink.flush_stats();
    let bucket_finish_started = Instant::now();
    let bucket_stats = sink.finish()?;
    let bucket_finish_elapsed = bucket_finish_started.elapsed();
    populate_emitted_graph_histogram(&mut stats, &bucket_stats)?;
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
    let sink = SharedBucketSink::create(params, graph_count)?;
    let mut stats = PartitionStats::new(graph_count);
    let mut worker_elapsed = Duration::ZERO;
    let mut parse_elapsed = Duration::ZERO;

    // The reference partitioner exposes the complete size-sorted source set to
    // its dynamic scheduler and drains each worker after every source.
    for (source_start, window) in [(0, paths)] {
        let parse_started = Instant::now();
        let mut work_order = (0..window.len()).collect::<Vec<_>>();
        work_order.sort_unstable_by_key(|&offset| {
            std::cmp::Reverse(window[offset].metadata().map_or(0, |meta| meta.len()))
        });
        let next_source = AtomicUsize::new(0);
        let workers = params.partition_workers(window.len());
        eprintln!("cuttlefish3-rs: colored partition using {workers} worker(s)");
        let mut emitters = Vec::with_capacity(workers);
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for _ in 0..workers {
                let sink = Arc::clone(&sink);
                let next_source = &next_source;
                let work_order = &work_order;
                handles.push(scope.spawn(move || {
                    let mut worker_stats = PartitionStats::new(graph_count);
                    let mut elapsed = Duration::ZERO;
                    let mut input_files = 0usize;
                    let mut records = 0u64;
                    let mut buckets = sink.emitter();
                    loop {
                        let order_idx = next_source.fetch_add(1, Ordering::Relaxed);
                        let Some(&offset) = work_order.get(order_idx) else {
                            break;
                        };
                        let source_id = u32::try_from(source_start + offset + 1)
                            .map_err(|_| PartitionRunError::TooManySources)?;
                        match parse_fragments_borrowed(
                            &window[offset],
                            source_id,
                            K + 1,
                            |fragment| {
                                let started = Instant::now();
                                emit_fragment_seq_weak_superkmer_buckets::<K, InputError>(
                                    params,
                                    graph_count,
                                    fragment.source_id,
                                    fragment.seq,
                                    &mut buckets,
                                    &mut worker_stats,
                                )?;
                                elapsed += started.elapsed();
                                Ok(())
                            },
                        ) {
                            Ok(parsed) => {
                                records += parsed;
                                input_files += 1;
                            }
                            Err(error) => {
                                handle_source_failure(params, &window[offset], error)?;
                            }
                        }
                        buckets.flush_colored_worker_if_required()?;
                    }
                    Ok::<_, PartitionRunError>((
                        records,
                        input_files,
                        PartitionWorkerOutput {
                            stats: worker_stats,
                            elapsed,
                        },
                        buckets,
                    ))
                }));
            }
            for handle in handles {
                let (records, input_files, output, buckets) = handle
                    .join()
                    .map_err(|_| PartitionRunError::WorkerPanic)??;
                stats.records += records;
                stats.input_files += input_files;
                worker_elapsed += output.elapsed;
                stats.merge_from(&output.stats);
                emitters.push(buckets);
            }
            Ok::<_, PartitionRunError>(())
        })?;
        parse_elapsed += parse_started.elapsed();

        sink.flush_colored_emitters(emitters)?;
    }

    let (bucket_flushes, bucket_flush_elapsed) = sink.flush_stats();
    let bucket_finish_started = Instant::now();
    let bucket_stats = sink.finish()?;
    let bucket_finish_elapsed = bucket_finish_started.elapsed();
    populate_emitted_graph_histogram(&mut stats, &bucket_stats)?;
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

    fn with_capacity(fragments: usize, bases: usize) -> Self {
        Self {
            fragments: Vec::with_capacity(fragments),
            bases: Vec::with_capacity(bases),
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

    #[inline]
    fn is_full(&self) -> bool {
        self.fragments.len() >= STREAM_BATCH_FRAGMENTS || self.bases.len() >= STREAM_BATCH_BASES
    }

    fn clear(&mut self) {
        self.fragments.clear();
        self.bases.clear();
    }
}

/// A bounded pool of recyclable fragment batches shared by readers and workers.
///
/// The pool is the backpressure mechanism for the streamed partitioner: readers
/// block when every buffer is in flight, and workers block until a reader
/// publishes one. Buffers are recycled rather than reallocated, so the steady
/// state performs no allocation.
struct BatchPool {
    inner: std::sync::Mutex<BatchPoolInner>,
    ready_available: std::sync::Condvar,
    free_available: std::sync::Condvar,
}

struct BatchPoolInner {
    ready: std::collections::VecDeque<PartitionFragmentBatch>,
    free: Vec<PartitionFragmentBatch>,
    readers: usize,
    failed: bool,
}

impl BatchPool {
    fn new(buffers: usize, readers: usize) -> Self {
        let buffers = buffers.max(1);
        let mut free = Vec::with_capacity(buffers);
        for _ in 0..buffers {
            free.push(PartitionFragmentBatch::with_capacity(
                STREAM_BATCH_FRAGMENTS,
                STREAM_BATCH_BASES,
            ));
        }
        Self {
            inner: std::sync::Mutex::new(BatchPoolInner {
                ready: std::collections::VecDeque::with_capacity(buffers),
                free,
                readers,
                failed: false,
            }),
            ready_available: std::sync::Condvar::new(),
            free_available: std::sync::Condvar::new(),
        }
    }

    /// Claims an empty buffer, blocking until one is recycled.
    ///
    /// Returns `None` once a worker has failed, so readers stop early instead of
    /// blocking forever on a pool nobody is draining.
    fn take_free(&self) -> Option<PartitionFragmentBatch> {
        let mut inner = self.inner.lock().ok()?;
        loop {
            if inner.failed {
                return None;
            }
            if let Some(batch) = inner.free.pop() {
                return Some(batch);
            }
            inner = self.free_available.wait(inner).ok()?;
        }
    }

    fn publish(&self, batch: PartitionFragmentBatch) {
        if let Ok(mut inner) = self.inner.lock() {
            inner.ready.push_back(batch);
            self.ready_available.notify_one();
        }
    }

    /// Claims a filled buffer, blocking until one is published.
    ///
    /// Returns `None` when every reader has finished and the queue has drained.
    fn take_ready(&self) -> Option<PartitionFragmentBatch> {
        let mut inner = self.inner.lock().ok()?;
        loop {
            if let Some(batch) = inner.ready.pop_front() {
                return Some(batch);
            }
            if inner.readers == 0 || inner.failed {
                return None;
            }
            inner = self.ready_available.wait(inner).ok()?;
        }
    }

    fn release(&self, mut batch: PartitionFragmentBatch) {
        batch.clear();
        if let Ok(mut inner) = self.inner.lock() {
            inner.free.push(batch);
            self.free_available.notify_one();
        }
    }

    fn reader_done(&self) {
        if let Ok(mut inner) = self.inner.lock() {
            inner.readers -= 1;
            if inner.readers == 0 {
                self.ready_available.notify_all();
            }
        }
    }

    /// Aborts the pool so both sides unblock after a worker or reader error.
    fn fail(&self) {
        if let Ok(mut inner) = self.inner.lock() {
            inner.failed = true;
        }
        self.ready_available.notify_all();
        self.free_available.notify_all();
    }
}

impl Default for PartitionFragmentBatch {
    fn default() -> Self {
        Self::new()
    }
}

/// Reports a source that failed to parse, or propagates the failure.
///
/// With `--skip-unreadable` a corrupt input in a large corpus costs a warning
/// rather than the whole build. Records already emitted from the failed source
/// are retained, and its position in the input list is preserved so colored
/// source assignments do not shift.
fn handle_source_failure(
    params: &BuildParams,
    path: &Path,
    error: InputError,
) -> Result<(), PartitionRunError> {
    if !params.skip_unreadable {
        return Err(error.into());
    }
    eprintln!(
        "cuttlefish3-rs: skipping unreadable input {}: {error}",
        path.display()
    );
    Ok(())
}

/// Whether to skip the minimizer scan entirely (diagnostic only).
fn parse_only_diagnostic() -> bool {
    static PARSE_ONLY: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *PARSE_ONLY.get_or_init(|| std::env::var_os("CF3_RS_PARSE_ONLY").is_some())
}

/// Whether to skip record packing and bucket writes (diagnostic only).
fn scan_only_diagnostic() -> bool {
    static SCAN_ONLY: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *SCAN_ONLY.get_or_init(|| std::env::var_os("CF3_RS_SCAN_ONLY").is_some())
}

fn emit_fragment_seq_weak_superkmer_buckets<const K: usize, E>(
    params: &BuildParams,
    graph_count: usize,
    source_id: u32,
    seq: &[u8],
    buckets: &mut SharedBucketEmitter,
    stats: &mut PartitionStats,
) -> Result<(), E>
where
    E: From<PartitionError> + From<crate::buckets::BucketError>,
{
    stats.fragments += 1;
    stats.fragment_bases += seq.len() as u64;
    // Diagnostic: `CF3_RS_SCAN_ONLY` performs the minimizer scan without packing
    // or writing records, isolating scan cost from emission cost.
    let scan_only = scan_only_diagnostic();
    // Diagnostic: `CF3_RS_PARSE_ONLY` also skips the minimizer scan, leaving only
    // line assembly and fragment splitting, so the two can be costed apart.
    if parse_only_diagnostic() {
        return Ok(());
    }
    for_each_valid_weak_superkmer::<K, E, _>(
        seq,
        params.minimizer_len as usize,
        graph_count,
        params.color.then_some(source_id),
        |sk| {
            if !scan_only {
                buckets.add_valid(&sk, sk.sequence(seq))?;
            }
            stats.weak_superkmers += 1;
            stats.weak_superkmer_bases += sk.len as u64;
            Ok(())
        },
    )?;
    Ok(())
}

fn populate_emitted_graph_histogram(
    stats: &mut PartitionStats,
    buckets: &BucketEmitStats,
) -> Result<(), PartitionRunError> {
    stats.graph_histogram.fill(0);
    for entry in crate::buckets::read_manifest(&buckets.bucket_dir)? {
        stats.graph_histogram[entry.graph_id] = entry.records;
    }
    Ok(())
}

pub fn partition_inputs<const K: usize>(
    params: &BuildParams,
    graph_count: usize,
) -> Result<PartitionStats, PartitionRunError> {
    params.validate()?;

    let paths = expand_input_paths(params)?;
    let mut stats = PartitionStats::new(graph_count);
    let workers = params.partition_workers(paths.len());
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

/// Scans a fragment and reports each weak super-k-mer.
///
/// The rolling minimizer window and the super-k-mer boundary state share one
/// loop, as they do in C++. Splitting them across a per-base visitor closure
/// kept the boundary state in a closure environment, so every base reloaded and
/// stored it through memory and paid a call. `visit` now runs once per weak
/// super-k-mer rather than once per base.
#[inline]
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

    // Boundary state in locals, so the loop needs no closure environment.
    let mut cur_off = 0usize;
    let mut km1_idx = 0usize;
    let mut cur_g = graph_id(prefix_min[pivot].min(suffix_min), graph_count);
    let mut prev_g = graph_count;

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

        let len = km1_idx + window_k;
        let next_g = graph_id(prefix_min[pivot].min(suffix_min), graph_count);
        km1_idx += 1;
        if next_g != cur_g || len == max_super_km1_len {
            let next_off = cur_off + km1_idx;
            let left_joined = cur_off > 0;
            visit(WeakSuperKmer {
                graph_id: cur_g,
                offset: cur_off - usize::from(left_joined),
                len: usize::from(left_joined) + len + 1,
                source_id,
                left_discontinuous: left_joined && prev_g != cur_g,
                right_discontinuous: next_g != cur_g,
            })?;

            cur_off = next_off;
            prev_g = cur_g;
            cur_g = next_g;
            km1_idx = 0;
        }
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

fn reset_min_window(hashes: &[u64], prefix_min: &mut [u64], ring_size: usize, pivot: &mut usize) {
    prefix_min[0] = hashes[0];
    for idx in 1..=ring_size {
        prefix_min[idx] = prefix_min[idx - 1].min(hashes[idx]);
    }
    *pivot = ring_size;
}

#[inline(always)]
fn canonical_lmer_hash(fwd: u64, rev: u64) -> u64 {
    wyhash_u64(fwd, 0).min(wyhash_u64(rev, 0))
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

    #[test]
    fn unreadable_input_aborts_unless_skipping_is_requested() {
        let base = std::env::temp_dir().join(format!("cf3rs-skip-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&base);
        std::fs::create_dir_all(&base).unwrap();
        let good = base.join("good.fa");
        std::fs::write(&good, b">r1\nACGTACGTACGTACGTACGTACGT\n").unwrap();
        // A zero-byte `.gz` is the exact corruption seen in the wild: the file
        // exists and is listed, but no gzip member can be read from it.
        let bad = base.join("bad.fa.gz");
        std::fs::write(&bad, b"").unwrap();

        let build = |skip: bool, threads: usize| {
            let work = base.join(format!("work-{skip}-{threads}"));
            let _ = std::fs::remove_dir_all(&work);
            std::fs::create_dir_all(&work).unwrap();
            let mut params = BuildParams::new(crate::GraphInput::References, "unused".to_string());
            params.seqs.push(good.to_string_lossy().into_owned());
            params.seqs.push(bad.to_string_lossy().into_owned());
            params.k = 7;
            params.minimizer_len = 3;
            params.threads = threads;
            params.skip_unreadable = skip;
            params.work_dir = work.to_string_lossy().into_owned();
            emit_weak_superkmer_buckets::<7>(&params, 128)
        };

        // `threads = 1` takes the direct path, `threads = 8` the streamed one,
        // because the streamed path engages only when files < workers.
        for threads in [1usize, 8] {
            assert!(
                build(false, threads).is_err(),
                "unreadable input must abort by default at {threads} thread(s)"
            );
            let stats = build(true, threads)
                .unwrap_or_else(|e| panic!("skipping should succeed at {threads} thread(s): {e}"));
            assert_eq!(stats.partition.input_files, 1, "only the good source counts");
            assert!(stats.partition.weak_superkmers > 0);
        }

        let _ = std::fs::remove_dir_all(&base);
    }

    #[test]
    fn streamed_partition_matches_direct_partition() {
        // One input file with four workers selects the streamed path; the direct
        // path is the same partitioner mapping that file onto a single worker.
        let manifest = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let input = manifest.join("../../data/refs2.fa");
        assert!(input.exists(), "missing fixture {}", input.display());

        let run = |dir: &std::path::Path| {
            let _ = std::fs::remove_dir_all(dir);
            std::fs::create_dir_all(dir).unwrap();
            let mut params = BuildParams::new(crate::GraphInput::References, "unused".to_string());
            params.seqs.push(input.to_string_lossy().into_owned());
            params.k = 7;
            params.minimizer_len = 3;
            params.threads = 4;
            params.work_dir = dir.to_string_lossy().into_owned();
            let paths = expand_input_paths(&params).unwrap();
            (params, paths)
        };

        let base = std::env::temp_dir().join(format!("cf3rs-stream-{}", std::process::id()));
        let streamed_dir = base.join("streamed");
        let direct_dir = base.join("direct");

        let (params, paths) = run(&streamed_dir);
        assert_eq!(streamed_reader_count(&params, paths.len()), Some(1));
        let streamed =
            emit_uncolored_streamed_weak_superkmer_buckets::<7>(&params, 128, &paths, 1).unwrap();

        let (params, paths) = run(&direct_dir);
        let direct =
            emit_uncolored_direct_weak_superkmer_buckets::<7>(&params, 128, &paths).unwrap();

        assert_eq!(streamed.partition.records, direct.partition.records);
        assert_eq!(streamed.partition.fragments, direct.partition.fragments);
        assert_eq!(
            streamed.partition.fragment_bases,
            direct.partition.fragment_bases
        );
        assert_eq!(
            streamed.partition.weak_superkmers,
            direct.partition.weak_superkmers
        );
        assert_eq!(
            streamed.partition.weak_superkmer_bases,
            direct.partition.weak_superkmer_bases
        );
        // Per-bucket record counts prove every weak super-k-mer reached the same
        // subgraph regardless of which worker packed it.
        assert_eq!(
            streamed.partition.graph_histogram,
            direct.partition.graph_histogram
        );
        assert_eq!(streamed.buckets.bucket_files, direct.buckets.bucket_files);
        assert_eq!(streamed.buckets.bytes_written, direct.buckets.bytes_written);

        let _ = std::fs::remove_dir_all(&base);
    }
}

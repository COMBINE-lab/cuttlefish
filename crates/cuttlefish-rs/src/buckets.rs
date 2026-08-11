//! External weak-super-k-mer bucket storage.
//!
//! Writers use the same 128-subgraph atlas hierarchy as Cuttlefish 3. Worker
//! buffers are drained in source/worker order, and open files are bounded so a
//! build does not require one descriptor per subgraph. The on-disk format is a
//! private, versioned Rust format and may be compressed with LZ4 blocks.

use crate::discontinuity::{current_open_file_count, open_file_limit};
use crate::dna::{Base, ascii_base_bits, valid_ascii_base_bits};
use crate::params::BuildParams;
use crate::partition::WeakSuperKmer;
use std::collections::BTreeMap;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::{
    Arc, Mutex,
    atomic::{AtomicU64, AtomicUsize, Ordering},
};
use std::time::{Duration, Instant};

const MAGIC: &[u8; 8] = b"CF3WSK1\0";
const RECORD_COUNT_OFFSET: u64 = 34;
const HEADER_LEN: u64 = 42;
const COMPRESSED_BLOCK_HEADER_LEN: usize = 12;
const MAX_SOURCE_ID: u32 = 0x1f_ffff;
const MAX_OPEN_BUCKET_WRITERS: usize = 512;
/// Largest record the staging buffers emit: attribute plus four label words.
const MAX_RECORD_BYTES: usize = 4 + 4 * 8;
const MAX_PENDING_BUCKET_BYTES: usize = 1024 * 1024;
// Keep a colored source window coalesced in worker-local atlas buffers. C++
// retains roughly this amount per active worker set; the larger cap avoids
// repeatedly scanning every graph bucket merely to move tiny fragments into
// the shared atlas.
const MAX_TOTAL_PENDING_BYTES: usize = 128 * 1024 * 1024;
const ATLAS_GRAPH_COUNT: usize = 128;
const SUBGRAPH_CHUNK_BYTES: usize = 64 * 1024;

/// Bytes reserved per segment when buckets share container files.
///
/// The segment decouples the write unit from the read and reclaim units. A
/// flush stays at `SUBGRAPH_CHUNK_BYTES` because that is a memory budget --
/// 16,384 buckets times 64 KiB of staging -- while reads become segment-sized
/// and reclaim becomes block-aligned, which is what lets a consumed bucket be
/// punched out in full.
///
/// Sized from measurement on 149,998 Salmonella assemblies: 237.8 GB across
/// 16,385 buckets, mean 14.51 MB, p1 8.80 MB and p99 22.77 MB, so the
/// distribution is tight and the only real cost is the partial final segment
/// each bucket leaves. At 256 KiB that is 2.15 GB, 0.90% of the directory,
/// against 7.3 MB of chain metadata. Larger segments trade waste for fewer
/// reads, and reads were measured not to matter on local storage.
const DEFAULT_BUCKET_SEGMENT_BYTES: u64 = 256 * 1024;

/// Descriptors left for everything downstream of the bucket containers: the
/// edge-matrix containers, local-unitig buckets, stitch writers and the
/// coordinate-bucket fanout, all of which plan against the same budget.
const RESERVED_NON_BUCKET_DESCRIPTORS: usize = 384;

fn bucket_segment_bytes() -> u64 {
    static BYTES: std::sync::OnceLock<u64> = std::sync::OnceLock::new();
    *BYTES.get_or_init(|| {
        std::env::var("CF3_RS_BUCKET_SEGMENT_BYTES")
            .ok()
            .and_then(|value| value.parse::<u64>().ok())
            .filter(|bytes| *bytes >= MAX_RECORD_BYTES as u64 && bytes % 4096 == 0)
            .unwrap_or(DEFAULT_BUCKET_SEGMENT_BYTES)
    })
}

/// The physical files backing the weak-super-k-mer buckets.
///
/// One container per atlas rather than one file per bucket, taking a 16,384
/// bucket build from 16,385 files to 129. Two things make this cheap rather
/// than intricate. Every bucket already writes under its atlas's mutex
/// (`SharedBucketSink::append_bucket`), and a container holds exactly one
/// atlas's buckets, so a container is only ever written by the thread holding
/// that lock and needs no lock of its own. And a flush no longer opens
/// anything: `BucketFile::open_existing` cost an `openat`, seven unbuffered
/// reads to re-read a 42-byte header, a revalidation, an `lseek` and a `close`
/// on *every* 64 KiB flush, which is about eleven syscalls times the 14.4
/// million flushes a full-corpus build performs. A container flush is one
/// `pwrite`.
///
/// The measured prize is smaller than that count suggests -- partitioning's
/// whole system time is 436.7 s of CPU across 64 threads, so the ceiling is
/// a couple of seconds of wall -- and larger somewhere unexpected. XFS
/// speculatively preallocates on extending writes, and with 16,385 repeatedly
/// reopened files it held 332.7 GB for 237.8 GB of data. Writing 128 files
/// instead returns that 94.9 GB.
#[derive(Debug)]
pub struct BucketContainers {
    files: Vec<BucketContainerFile>,
    segment_bytes: u64,
    /// Latched so an unsupported filesystem is reported once, not per bucket.
    punch_unsupported: std::sync::atomic::AtomicBool,
}

#[derive(Debug)]
struct BucketContainerFile {
    path: PathBuf,
    file: File,
    /// Next unreserved byte offset. Atomic for shape rather than contention:
    /// one atlas owns one container.
    cursor: AtomicU64,
}

impl BucketContainers {
    /// Containers a build may hold open, given the descriptor budget.
    ///
    /// One per atlas is the natural choice and what a normal limit allows, but
    /// it must not be a floor: the per-file layout this replaced held no
    /// descriptors between flushes, so a tight `ulimit -n` merely narrowed the
    /// fanout planners rather than failing the build. Sharing a container
    /// between atlases costs nothing -- the segment cursor is atomic, so
    /// concurrent reservations from different atlas locks are already safe --
    /// and keeps that property.
    fn plan_container_count(atlas_count: usize) -> usize {
        let budget = open_file_limit()
            .saturating_sub(current_open_file_count())
            // Local contraction opens the edge-matrix containers and the
            // local-unitig buckets on top of these, and the fanout planners
            // want room of their own.
            .saturating_sub(RESERVED_NON_BUCKET_DESCRIPTORS)
            / 2;
        atlas_count.min(budget.max(1))
    }

    fn create(bucket_dir: &Path, container_count: usize) -> Result<Self, BucketError> {
        let mut files = Vec::with_capacity(container_count);
        for index in 0..container_count {
            let path = bucket_dir.join(format!("{index:05}.wskc"));
            let file = OpenOptions::new()
                .create(true)
                .truncate(true)
                .read(true)
                .write(true)
                .open(&path)
                .map_err(|source| BucketError::Io {
                    path: path.clone(),
                    source,
                })?;
            files.push(BucketContainerFile {
                path,
                file,
                cursor: AtomicU64::new(0),
            });
        }
        Ok(Self {
            files,
            segment_bytes: bucket_segment_bytes(),
            punch_unsupported: std::sync::atomic::AtomicBool::new(false),
        })
    }

    /// Opens the containers a finished manifest names, for reading.
    fn open(
        bucket_dir: &Path,
        container_count: usize,
        segment_bytes: u64,
    ) -> Result<Self, BucketError> {
        let mut files = Vec::with_capacity(container_count);
        for index in 0..container_count {
            let path = bucket_dir.join(format!("{index:05}.wskc"));
            let file = OpenOptions::new()
                .read(true)
                .write(true)
                .open(&path)
                .map_err(|source| BucketError::Io {
                    path: path.clone(),
                    source,
                })?;
            files.push(BucketContainerFile {
                path,
                file,
                cursor: AtomicU64::new(0),
            });
        }
        Ok(Self {
            files,
            segment_bytes,
            punch_unsupported: std::sync::atomic::AtomicBool::new(false),
        })
    }

    #[inline]
    pub fn segment_bytes(&self) -> u64 {
        self.segment_bytes
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.files.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.files.is_empty()
    }

    #[inline]
    fn reserve_segment(&self, container: usize) -> u64 {
        self.files[container]
            .cursor
            .fetch_add(self.segment_bytes, Ordering::Relaxed)
    }

    fn write_at(&self, container: usize, offset: u64, bytes: &[u8]) -> Result<(), BucketError> {
        let file = &self.files[container];
        file.file
            .write_all_at(bytes, offset)
            .map_err(|source| BucketError::Io {
                path: file.path.clone(),
                source,
            })
    }

    fn read_at(&self, container: usize, offset: u64, buf: &mut [u8]) -> Result<(), BucketError> {
        let file = &self.files[container];
        file.file
            .read_exact_at(buf, offset)
            .map_err(|source| BucketError::Io {
                path: file.path.clone(),
                source,
            })
    }

    /// Releases a consumed bucket's segments without disturbing its neighbours.
    ///
    /// This is not optional, which measurement rather than reasoning settled.
    /// The per-file layout unlinked each bucket as local contraction consumed
    /// it, and containers cannot: a container is only droppable once all 128
    /// of its buckets are done. Containers do start about 94 GB below the
    /// per-file layout, because they do not accumulate XFS speculative
    /// preallocation, and the expectation was that this covered it. It does
    /// not -- the work-directory peak moves out of the end of partitioning and
    /// into local contraction, where the containers still hold everything
    /// while local-unitig buckets, labels and the edge matrix accumulate on
    /// top, and peak disk rose 24.5 GB.
    ///
    /// Punching restores the incremental release. Segments are reserved at a
    /// 4 KiB multiple, so a punch frees whole filesystem blocks rather than
    /// leaving partial ones behind -- the one thing raw extents, averaging
    /// 16.1 KiB and unaligned, could not have done. Adjacent segments are
    /// punched in one call.
    pub fn release_segments(&self, container: usize, segments: &[u32]) {
        if segments.is_empty() {
            return;
        }
        let mut ordered = segments.to_vec();
        ordered.sort_unstable();
        let file = &self.files[container];
        let mut start = u64::from(ordered[0]) * self.segment_bytes;
        let mut end = start + self.segment_bytes;
        for &index in &ordered[1..] {
            let offset = u64::from(index) * self.segment_bytes;
            if offset == end {
                end += self.segment_bytes;
                continue;
            }
            self.report_punch(punch_hole(&file.file, start, end - start));
            start = offset;
            end = offset + self.segment_bytes;
        }
        self.report_punch(punch_hole(&file.file, start, end - start));
    }

    /// Says once if the filesystem will not punch holes.
    ///
    /// Reclaim failing is not an error -- the container is unlinked wholesale
    /// at the end regardless -- but it silently costs peak disk, which is most
    /// of what the container layout is for. HFS+ and some network filesystems
    /// do not implement it, and macOS rejects a range that is not aligned to
    /// the filesystem block size. Better to say so than to leave someone
    /// wondering why the work directory is larger than documented.
    fn report_punch(&self, punched: bool) {
        if punched || self.punch_unsupported.swap(true, Ordering::Relaxed) {
            return;
        }
        eprintln!(
            "cuttlefish: this filesystem will not punch holes, so consumed \
             bucket space is held until the build ends; peak disk will be higher"
        );
    }

    pub fn paths(&self) -> impl Iterator<Item = &Path> {
        self.files.iter().map(|file| file.path.as_path())
    }
}

/// Punches `len` bytes at `offset` out of `file`.
///
/// There is no portable interface for this and no crate that abstracts one:
/// Linux spells it `fallocate(FALLOC_FL_PUNCH_HOLE)` and macOS spells it
/// `fcntl(F_PUNCHHOLE)`, and both `nix` and `rustix` gate their `fallocate`
/// behind `target_os = "linux"` with no Apple equivalent offered. So the two
/// are written out here.
///
/// Both interfaces require the range to be filesystem-block aligned -- macOS
/// returns `EINVAL` for a punch that is not a multiple of the block size --
/// which every call here satisfies by construction, because offsets and
/// lengths are whole segments and `bucket_segment_bytes` admits only multiples
/// of 4096. That is a second, independent reason buckets got segments rather
/// than the raw 16.1 KiB extents a flush would otherwise produce.
///
/// Failure is ignored on purpose. A filesystem without hole punching -- an old
/// HFS+ volume, or a network mount -- costs disk rather than correctness,
/// because the container is unlinked wholesale at the end regardless.
#[cfg(target_os = "linux")]
fn punch_hole(file: &File, offset: u64, len: u64) -> bool {
    use std::os::fd::AsRawFd;
    // SAFETY: the descriptor is owned by `file` and outlives the call, and the
    // kernel validates the range itself.
    unsafe {
        libc::fallocate(
            file.as_raw_fd(),
            libc::FALLOC_FL_KEEP_SIZE | libc::FALLOC_FL_PUNCH_HOLE,
            offset as libc::off_t,
            len as libc::off_t,
        ) == 0
    }
}

#[cfg(target_vendor = "apple")]
fn punch_hole(file: &File, offset: u64, len: u64) -> bool {
    use std::os::fd::AsRawFd;
    let punch = libc::fpunchhole_t {
        fp_flags: 0,
        reserved: 0,
        fp_offset: offset as libc::off_t,
        fp_length: len as libc::off_t,
    };
    // SAFETY: as above; `punch` outlives the call and F_PUNCHHOLE reads it as
    // a `*const fpunchhole_t`.
    unsafe { libc::fcntl(file.as_raw_fd(), libc::F_PUNCHHOLE, &punch) == 0 }
}

#[cfg(not(any(target_os = "linux", target_vendor = "apple")))]
fn punch_hole(_file: &File, _offset: u64, _len: u64) -> bool {
    false
}

/// Sequential reader over one bucket's segment chain.
///
/// `BucketReader` uses its source purely as a `Read`, so presenting the chain
/// this way leaves the whole decode path -- record iteration, the compressed
/// block framing, the borrowed-record fast path -- exactly as it was for whole
/// files. Reads stop at each segment boundary, and the caller's `BufReader`
/// hides that.
#[derive(Debug)]
struct SegmentChainReader {
    containers: Arc<BucketContainers>,
    container: usize,
    segments: Vec<u64>,
    len: u64,
    pos: u64,
}

impl Read for SegmentChainReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        if self.pos >= self.len || buf.is_empty() {
            return Ok(0);
        }
        let segment_bytes = self.containers.segment_bytes;
        let index = (self.pos / segment_bytes) as usize;
        let within = self.pos % segment_bytes;
        let room = (segment_bytes - within).min(self.len - self.pos);
        let take = (buf.len() as u64).min(room) as usize;
        let offset = self.segments[index] + within;
        self.containers
            .read_at(self.container, offset, &mut buf[..take])
            .map_err(std::io::Error::other)?;
        self.pos += take as u64;
        Ok(take)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BucketEmitStats {
    pub bucket_dir: PathBuf,
    pub bucket_files: usize,
    pub bytes_written: u64,
}

pub struct BucketEmitter {
    bucket_dir: PathBuf,
    k: u16,
    minimizer_len: u16,
    graph_count: usize,
    colored: bool,
    compress_buckets: bool,
    label_words: usize,
    files: BTreeMap<usize, BucketFileMeta>,
    writers: BTreeMap<usize, BucketFile>,
    pending: Vec<PendingBucket>,
    pending_bytes: usize,
    scratch: CompressionScratch,
}

pub struct SharedBucketSink {
    bucket_dir: PathBuf,
    containers: BucketContainers,
    k: u16,
    minimizer_len: u16,
    graph_count: usize,
    colored: bool,
    compress_buckets: bool,
    label_words: usize,
    workers: usize,
    atlases: Vec<Mutex<SharedBucketAtlas>>,
    flush_calls: AtomicU64,
    flush_nanos: AtomicU64,
}

pub struct SharedBucketEmitter {
    sink: Arc<SharedBucketSink>,
    pending: Vec<PendingBucket>,
    uncolored_pending: Vec<PendingColoredAtlas>,
    colored_pending: Vec<PendingColoredAtlas>,
    pending_bytes: usize,
    deferred_uncolored: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BucketFileMeta {
    path: PathBuf,
    records: u64,
    bytes_written: u64,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
struct PendingBucket {
    records: u64,
    bytes: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
struct PendingColoredAtlas {
    graph_ids: Vec<u16>,
    bytes: Vec<u8>,
}

#[derive(Default)]
struct SharedBucketFileMeta {
    buffer: Vec<u8>,
    buffer_records: u64,
    total_records: u64,
    written_records: u64,
    bytes_written: u64,
    /// Segment indices this bucket owns, in write order.
    segments: Vec<u32>,
    /// Bytes used in the last segment; a flush that overruns it straddles into
    /// a freshly reserved one, which the reader stitches back because it walks
    /// the chain in order.
    segment_used: u64,
}

struct SharedBucketAtlas {
    first_graph_id: usize,
    files: Vec<SharedBucketFileMeta>,
    buffered_bytes: usize,
    /// Shared by every flush through this atlas, which the atlas lock already
    /// serializes, so it costs no contention and saves an allocation per block.
    scratch: CompressionScratch,
}

#[derive(Default)]
struct SharedBucketFlushStats {
    calls: u64,
}

/// Where one bucket's bytes live.
///
/// Both forms are real: the shared production sink writes containers, and
/// `BucketEmitter` -- the serial emitter the tests and the legacy path use --
/// still writes one file per bucket with its header inline. Keeping the enum
/// lets the reader serve both without duplicating the decode path.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BucketLocation {
    /// One whole file, carrying its own 42-byte header.
    File(PathBuf),
    /// A chain of segments inside a shared container. The header that a whole
    /// file would carry lives in the manifest instead, so nothing has to be
    /// re-read and revalidated per flush.
    Container {
        container: usize,
        /// Segment indices in write order; byte offset is index * segment_bytes.
        segments: Vec<u32>,
        /// Logical length, which is the sum of the payload written into those
        /// segments and is generally less than their reserved capacity.
        bytes: u64,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// Manifest entry for one weak-super-k-mer bucket.
pub struct BucketManifestEntry {
    pub graph_id: usize,
    pub records: u64,
    pub location: BucketLocation,
}

impl BucketManifestEntry {
    /// Bytes this bucket occupies, for the longest-bucket-first scheduler.
    ///
    /// Containers make this free. The per-file layout had to `stat` all 16,384
    /// buckets before local contraction could sort them.
    pub fn stored_bytes(&self) -> Result<u64, BucketError> {
        match &self.location {
            BucketLocation::File(path) => {
                fs::metadata(path)
                    .map(|meta| meta.len())
                    .map_err(|source| BucketError::Io {
                        path: path.clone(),
                        source,
                    })
            }
            BucketLocation::Container { bytes, .. } => Ok(*bytes),
        }
    }

    /// The whole-file path, when the bucket is one.
    pub fn file_path(&self) -> Option<&Path> {
        match &self.location {
            BucketLocation::File(path) => Some(path.as_path()),
            BucketLocation::Container { .. } => None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// Versioned parameters decoded from a bucket file header.
pub struct BucketHeader {
    pub k: u16,
    pub minimizer_len: u16,
    pub graph_count: usize,
    pub graph_id: usize,
    pub colored: bool,
    pub compressed: bool,
    pub interleaved_compression: bool,
    pub label_words: usize,
    pub records: u64,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
/// Decoded weak-super-k-mer record with an ASCII label.
pub struct BucketRecord {
    pub graph_id: usize,
    pub len: usize,
    pub source_id: Option<u32>,
    pub left_discontinuous: bool,
    pub right_discontinuous: bool,
    pub label: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
/// Decoded weak-super-k-mer record retaining its packed two-bit label.
pub struct BucketPackedRecord {
    pub graph_id: usize,
    pub len: usize,
    pub source_id: Option<u32>,
    pub left_discontinuous: bool,
    pub right_discontinuous: bool,
    pub words: Vec<u64>,
}

pub(crate) struct BorrowedBucketPackedRecord<'a> {
    pub graph_id: usize,
    pub len: usize,
    pub source_id: Option<u32>,
    pub left_discontinuous: bool,
    pub right_discontinuous: bool,
    pub words: &'a [u64],
}

/// Byte source behind a `BucketReader`.
///
/// The decode path treats its source as a plain sequential `Read`, so a
/// segment chain slots in beside a whole file without any of the record
/// iteration, block framing, or borrowed-record handling needing to know.
enum BucketSource {
    File(File),
    Chain(SegmentChainReader),
}

impl Read for BucketSource {
    #[inline]
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            Self::File(file) => file.read(buf),
            Self::Chain(chain) => chain.read(buf),
        }
    }
}

/// Opens buckets, whichever layout the directory uses.
///
/// A container directory keeps the shared header in its manifest, so this owns
/// both and hands out readers; a whole-file directory carries the header in
/// each bucket and this just opens paths. Threading one of these through the
/// consumers replaces threading a path per bucket.
pub struct BucketStore {
    containers: Option<Arc<BucketContainers>>,
    header: Option<ContainerManifestHeader>,
}

impl BucketStore {
    /// Opens a bucket directory and reads its manifest.
    pub fn open_dir(
        bucket_dir: impl AsRef<Path>,
    ) -> Result<(Self, Vec<BucketManifestEntry>), BucketError> {
        let bucket_dir = bucket_dir.as_ref();
        if let Some((header, entries)) = read_container_manifest(bucket_dir)? {
            let containers =
                BucketContainers::open(bucket_dir, header.container_count, header.segment_bytes)?;
            return Ok((
                Self {
                    containers: Some(Arc::new(containers)),
                    header: Some(header),
                },
                entries,
            ));
        }
        let entries = read_manifest(bucket_dir)?;
        Ok((
            Self {
                containers: None,
                header: None,
            },
            entries,
        ))
    }

    /// A store for a directory that is known to hold whole bucket files.
    pub fn files_only() -> Self {
        Self {
            containers: None,
            header: None,
        }
    }

    pub fn containers(&self) -> Option<&Arc<BucketContainers>> {
        self.containers.as_ref()
    }

    /// Opens one bucket for reading.
    pub fn reader(&self, entry: &BucketManifestEntry) -> Result<BucketReader, BucketError> {
        match &entry.location {
            BucketLocation::File(path) => BucketReader::open(path),
            BucketLocation::Container {
                container,
                segments,
                bytes,
            } => {
                let (Some(containers), Some(header)) = (&self.containers, &self.header) else {
                    return Err(BucketError::MalformedRecord);
                };
                BucketReader::open_chain(
                    Arc::clone(containers),
                    *container,
                    segments.clone(),
                    *bytes,
                    header.bucket_header(entry.graph_id, entry.records),
                )
            }
        }
    }
}

/// Streaming reader for one versioned weak-super-k-mer bucket.
pub struct BucketReader {
    path: PathBuf,
    file: BufReader<BucketSource>,
    header: BucketHeader,
    records_read: u64,
    block_attrs: Vec<u8>,
    block_labels: Vec<u8>,
    block_interleaved: Vec<u8>,
    compressed_block: Vec<u8>,
    block_records: usize,
    block_record: usize,
}

impl BucketReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, BucketError> {
        let path = path.as_ref().to_path_buf();
        let file = File::open(&path).map_err(|source| BucketError::Io {
            path: path.clone(),
            source,
        })?;
        let mut file = BufReader::with_capacity(1024 * 1024, BucketSource::File(file));
        let header = read_header(&mut file, &path)?;

        Ok(Self {
            path,
            file,
            header,
            records_read: 0,
            block_attrs: Vec::new(),
            block_labels: Vec::new(),
            block_interleaved: Vec::new(),
            compressed_block: Vec::new(),
            block_records: 0,
            block_record: 0,
        })
    }

    /// Opens a bucket stored as a segment chain.
    ///
    /// There is no inline header to skip: a container's buckets keep theirs in
    /// the manifest, so the chain starts at the first record byte and the
    /// caller supplies the header it would otherwise have read.
    fn open_chain(
        containers: Arc<BucketContainers>,
        container: usize,
        segments: Vec<u32>,
        bytes: u64,
        header: BucketHeader,
    ) -> Result<Self, BucketError> {
        let segment_bytes = containers.segment_bytes();
        let path = containers.files[container].path.clone();
        let chain = SegmentChainReader {
            containers,
            container,
            segments: segments
                .iter()
                .map(|s| u64::from(*s) * segment_bytes)
                .collect(),
            len: bytes,
            pos: 0,
        };
        Ok(Self {
            path,
            file: BufReader::with_capacity(1024 * 1024, BucketSource::Chain(chain)),
            header,
            records_read: 0,
            block_attrs: Vec::new(),
            block_labels: Vec::new(),
            block_interleaved: Vec::new(),
            compressed_block: Vec::new(),
            block_records: 0,
            block_record: 0,
        })
    }

    #[inline]
    pub fn header(&self) -> &BucketHeader {
        &self.header
    }

    pub fn next_record(&mut self) -> Result<Option<BucketRecord>, BucketError> {
        let mut record = BucketRecord::default();
        if self.next_record_into(&mut record)? {
            Ok(Some(record))
        } else {
            Ok(None)
        }
    }

    pub fn next_record_into(&mut self, record: &mut BucketRecord) -> Result<bool, BucketError> {
        let mut packed = BucketPackedRecord::default();
        if !self.next_packed_record_into(&mut packed)? {
            return Ok(false);
        }

        record.graph_id = packed.graph_id;
        record.len = packed.len;
        record.source_id = packed.source_id;
        record.left_discontinuous = packed.left_discontinuous;
        record.right_discontinuous = packed.right_discontinuous;
        decode_label_into(&packed.words, packed.len, &mut record.label)?;
        Ok(true)
    }

    pub fn next_packed_record_into(
        &mut self,
        record: &mut BucketPackedRecord,
    ) -> Result<bool, BucketError> {
        if self.records_read == self.header.records {
            return Ok(false);
        }

        let mut fixed = [0u8; 4];
        let fixed_len = if self.header.colored { 4 } else { 2 };
        if self.header.compressed {
            if self.block_record == self.block_records {
                self.read_compressed_block()?;
            }
            let start = self.block_record * fixed_len;
            if self.header.interleaved_compression {
                let record_len = record_size(self.header.colored, self.header.label_words);
                let start = self.block_record * record_len;
                fixed[..fixed_len]
                    .copy_from_slice(&self.block_interleaved[start..start + fixed_len]);
            } else {
                fixed[..fixed_len].copy_from_slice(&self.block_attrs[start..start + fixed_len]);
            }
        } else {
            self.file
                .read_exact(&mut fixed[..fixed_len])
                .map_err(|source| BucketError::Io {
                    path: self.path.clone(),
                    source,
                })?;
        }
        let packed_attr = if self.header.colored {
            u32::from_le_bytes(fixed[0..4].try_into().unwrap())
        } else {
            u16::from_le_bytes(fixed[0..2].try_into().unwrap()) as u32
        };

        record.words.clear();
        record.words.resize(self.header.label_words, 0);
        for (word_idx, word) in record.words.iter_mut().enumerate() {
            if self.header.compressed {
                if self.header.interleaved_compression {
                    let record_len = record_size(self.header.colored, self.header.label_words);
                    let start = self.block_record * record_len + fixed_len + word_idx * 8;
                    *word = u64::from_le_bytes(
                        self.block_interleaved[start..start + 8].try_into().unwrap(),
                    );
                } else {
                    let start = (self.block_record * self.header.label_words + word_idx) * 8;
                    *word =
                        u64::from_le_bytes(self.block_labels[start..start + 8].try_into().unwrap());
                }
            } else {
                *word = read_u64(&mut self.file, &self.path)?;
            }
        }

        let len = (packed_attr & 0xff) as usize;
        let source_id = self
            .header
            .colored
            .then_some((packed_attr >> 10) & MAX_SOURCE_ID);
        record.graph_id = self.header.graph_id;
        record.len = len;
        record.source_id = source_id;
        record.left_discontinuous = (packed_attr & (1 << 8)) != 0;
        record.right_discontinuous = (packed_attr & (1 << 9)) != 0;

        self.records_read += 1;
        if self.header.compressed {
            self.block_record += 1;
        }
        Ok(true)
    }

    pub fn try_for_each_packed_record<E, F>(
        &mut self,
        record: &mut BucketPackedRecord,
        mut f: F,
    ) -> Result<(), E>
    where
        E: From<BucketError>,
        F: FnMut(&BucketPackedRecord) -> Result<(), E>,
    {
        let remaining = self.header.records.saturating_sub(self.records_read);
        if remaining == 0 {
            return Ok(());
        }
        if self.header.compressed {
            while self.next_packed_record_into(record).map_err(E::from)? {
                f(record)?;
            }
            return Ok(());
        }
        let record_size = record_size(self.header.colored, self.header.label_words);
        let payload_bytes = remaining
            .checked_mul(record_size as u64)
            .ok_or(BucketError::TooManyRecords)?;
        let mut payload = vec![0u8; payload_bytes as usize];
        self.file
            .read_exact(&mut payload)
            .map_err(|source| BucketError::Io {
                path: self.path.clone(),
                source,
            })?;

        for chunk in payload.chunks_exact(record_size) {
            let (packed_attr, words_start) = if self.header.colored {
                (u32::from_le_bytes(chunk[0..4].try_into().unwrap()), 4)
            } else {
                (
                    u16::from_le_bytes(chunk[0..2].try_into().unwrap()) as u32,
                    2,
                )
            };

            record.words.clear();
            record.words.reserve(self.header.label_words);
            for word_idx in 0..self.header.label_words {
                let start = words_start + word_idx * 8;
                record.words.push(u64::from_le_bytes(
                    chunk[start..start + 8].try_into().unwrap(),
                ));
            }

            let len = (packed_attr & 0xff) as usize;
            let source_id = self
                .header
                .colored
                .then_some((packed_attr >> 10) & MAX_SOURCE_ID);
            record.graph_id = self.header.graph_id;
            record.len = len;
            record.source_id = source_id;
            record.left_discontinuous = (packed_attr & (1 << 8)) != 0;
            record.right_discontinuous = (packed_attr & (1 << 9)) != 0;

            self.records_read += 1;
            f(record)?;
        }
        Ok(())
    }

    pub(crate) fn try_for_each_borrowed_packed_record<E, F>(&mut self, mut f: F) -> Result<(), E>
    where
        E: From<BucketError>,
        F: FnMut(BorrowedBucketPackedRecord<'_>) -> Result<(), E>,
    {
        if !self.header.compressed {
            let mut record = BucketPackedRecord::default();
            return self.try_for_each_packed_record(&mut record, |record| {
                f(BorrowedBucketPackedRecord {
                    graph_id: record.graph_id,
                    len: record.len,
                    source_id: record.source_id,
                    left_discontinuous: record.left_discontinuous,
                    right_discontinuous: record.right_discontinuous,
                    words: &record.words,
                })
            });
        }
        if self.header.label_words > 4 {
            return Err(E::from(BucketError::MalformedRecord));
        }

        let fixed_len = if self.header.colored { 4 } else { 2 };
        while self.records_read < self.header.records {
            if self.block_record == self.block_records {
                self.read_compressed_block().map_err(E::from)?;
            }
            let record_index = self.block_record;
            let (packed_attr, words_bytes) = if self.header.interleaved_compression {
                let record_len = record_size(self.header.colored, self.header.label_words);
                let start = record_index * record_len;
                let attr = if self.header.colored {
                    u32::from_le_bytes(self.block_interleaved[start..start + 4].try_into().unwrap())
                } else {
                    u16::from_le_bytes(self.block_interleaved[start..start + 2].try_into().unwrap())
                        as u32
                };
                (
                    attr,
                    &self.block_interleaved[start + fixed_len..start + record_len],
                )
            } else {
                let attr_start = record_index * fixed_len;
                let attr = if self.header.colored {
                    u32::from_le_bytes(
                        self.block_attrs[attr_start..attr_start + 4]
                            .try_into()
                            .unwrap(),
                    )
                } else {
                    u16::from_le_bytes(
                        self.block_attrs[attr_start..attr_start + 2]
                            .try_into()
                            .unwrap(),
                    ) as u32
                };
                let words_start = record_index * self.header.label_words * 8;
                (
                    attr,
                    &self.block_labels[words_start..words_start + self.header.label_words * 8],
                )
            };
            let mut words = [0u64; 4];
            for (word, bytes) in words[..self.header.label_words]
                .iter_mut()
                .zip(words_bytes.chunks_exact(8))
            {
                *word = u64::from_le_bytes(bytes.try_into().unwrap());
            }
            self.records_read += 1;
            self.block_record += 1;
            f(BorrowedBucketPackedRecord {
                graph_id: self.header.graph_id,
                len: (packed_attr & 0xff) as usize,
                source_id: self
                    .header
                    .colored
                    .then_some((packed_attr >> 10) & MAX_SOURCE_ID),
                left_discontinuous: packed_attr & (1 << 8) != 0,
                right_discontinuous: packed_attr & (1 << 9) != 0,
                words: &words[..self.header.label_words],
            })?;
        }
        Ok(())
    }

    pub fn records(self) -> BucketRecords {
        BucketRecords { reader: self }
    }

    fn read_compressed_block(&mut self) -> Result<(), BucketError> {
        let mut header = [0u8; COMPRESSED_BLOCK_HEADER_LEN];
        self.file
            .read_exact(&mut header)
            .map_err(|source| BucketError::Io {
                path: self.path.clone(),
                source,
            })?;
        let records = u32::from_le_bytes(header[0..4].try_into().unwrap()) as usize;
        let attr_bytes = u32::from_le_bytes(header[4..8].try_into().unwrap()) as usize;
        let label_bytes = u32::from_le_bytes(header[8..12].try_into().unwrap()) as usize;
        if records == 0
            || attr_bytes == 0
            || (!self.header.interleaved_compression && label_bytes == 0)
            || (self.header.interleaved_compression && label_bytes != 0)
        {
            return Err(BucketError::MalformedRecord);
        }
        let remaining = usize::try_from(self.header.records - self.records_read)
            .map_err(|_| BucketError::TooManyRecords)?;
        if records > remaining {
            return Err(BucketError::MalformedRecord);
        }
        self.compressed_block.resize(attr_bytes + label_bytes, 0);
        self.file
            .read_exact(&mut self.compressed_block)
            .map_err(|source| BucketError::Io {
                path: self.path.clone(),
                source,
            })?;
        let fixed_len = if self.header.colored { 4 } else { 2 };
        if self.header.interleaved_compression {
            self.block_interleaved.resize(
                records * record_size(self.header.colored, self.header.label_words),
                0,
            );
            let decoded = lz4_flex::block::decompress_into(
                &self.compressed_block[..attr_bytes],
                &mut self.block_interleaved,
            )
            .map_err(|_| BucketError::MalformedRecord)?;
            if decoded != self.block_interleaved.len() {
                return Err(BucketError::MalformedRecord);
            }
            self.block_records = records;
            self.block_record = 0;
            return Ok(());
        }
        self.block_attrs.resize(records * fixed_len, 0);
        self.block_labels
            .resize(records * self.header.label_words * 8, 0);
        let decoded_attrs = lz4_flex::block::decompress_into(
            &self.compressed_block[..attr_bytes],
            &mut self.block_attrs,
        )
        .map_err(|_| BucketError::MalformedRecord)?;
        let decoded_labels = lz4_flex::block::decompress_into(
            &self.compressed_block[attr_bytes..],
            &mut self.block_labels,
        )
        .map_err(|_| BucketError::MalformedRecord)?;
        if decoded_attrs != self.block_attrs.len() || decoded_labels != self.block_labels.len() {
            return Err(BucketError::MalformedRecord);
        }
        self.block_records = records;
        self.block_record = 0;
        Ok(())
    }
}

/// Iterator over decoded records from a [`BucketReader`].
pub struct BucketRecords {
    reader: BucketReader,
}

impl Iterator for BucketRecords {
    type Item = Result<BucketRecord, BucketError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader.next_record().transpose()
    }
}

const CONTAINER_MANIFEST_MAGIC: &[u8; 8] = b"CF3WSKC1";
const CONTAINER_MANIFEST_NAME: &str = "manifest.bin";

/// The parameters every bucket in a container directory shares.
///
/// This is the 42-byte per-bucket header, hoisted. Under the per-file layout it
/// was written once per bucket and then re-read and revalidated on *every*
/// 64 KiB flush, because a flush reopened the file to append. Stored once for
/// the whole directory, that disappears; the only genuinely per-bucket field
/// was the record count, which the manifest already carried.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContainerManifestHeader {
    pub k: u16,
    pub minimizer_len: u16,
    pub graph_count: usize,
    pub colored: bool,
    pub label_words: usize,
    pub compressed: bool,
    pub interleaved_compression: bool,
    pub segment_bytes: u64,
    pub container_count: usize,
}

impl ContainerManifestHeader {
    fn compression_code(&self) -> u8 {
        if self.interleaved_compression {
            2
        } else {
            u8::from(self.compressed)
        }
    }

    /// The per-bucket header a whole-file reader would have found inline.
    fn bucket_header(&self, graph_id: usize, records: u64) -> BucketHeader {
        BucketHeader {
            k: self.k,
            minimizer_len: self.minimizer_len,
            graph_count: self.graph_count,
            graph_id,
            colored: self.colored,
            label_words: self.label_words,
            compressed: self.compressed,
            interleaved_compression: self.interleaved_compression,
            records,
        }
    }
}

fn write_u32_to(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn write_u64_to(out: &mut Vec<u8>, value: u64) {
    out.extend_from_slice(&value.to_le_bytes());
}

/// Writes the container manifest: shared header, then one record per bucket.
pub fn write_container_manifest(
    bucket_dir: &Path,
    header: &ContainerManifestHeader,
    entries: &[BucketManifestEntry],
) -> Result<(), BucketError> {
    let path = bucket_dir.join(CONTAINER_MANIFEST_NAME);
    let mut out = Vec::with_capacity(64 + entries.len() * 32);
    out.extend_from_slice(CONTAINER_MANIFEST_MAGIC);
    out.extend_from_slice(&header.k.to_le_bytes());
    out.extend_from_slice(&header.minimizer_len.to_le_bytes());
    write_u64_to(&mut out, header.graph_count as u64);
    out.push(u8::from(header.colored));
    out.push(header.label_words as u8);
    out.push(header.compression_code());
    out.push(0);
    write_u64_to(&mut out, header.segment_bytes);
    write_u64_to(&mut out, header.container_count as u64);
    write_u64_to(&mut out, entries.len() as u64);
    for entry in entries {
        let BucketLocation::Container {
            container,
            segments,
            bytes,
        } = &entry.location
        else {
            return Err(BucketError::MalformedManifest(path.clone()));
        };
        write_u64_to(&mut out, entry.graph_id as u64);
        write_u64_to(&mut out, entry.records);
        write_u64_to(&mut out, *bytes);
        write_u32_to(&mut out, *container as u32);
        write_u32_to(&mut out, segments.len() as u32);
        for segment in segments {
            write_u32_to(&mut out, *segment);
        }
    }
    fs::write(&path, &out).map_err(|source| BucketError::Io {
        path: path.clone(),
        source,
    })
}

struct ManifestCursor<'a> {
    bytes: &'a [u8],
    pos: usize,
    path: &'a Path,
}

impl<'a> ManifestCursor<'a> {
    fn take(&mut self, len: usize) -> Result<&'a [u8], BucketError> {
        let end = self
            .pos
            .checked_add(len)
            .filter(|end| *end <= self.bytes.len())
            .ok_or_else(|| BucketError::MalformedManifest(self.path.to_path_buf()))?;
        let slice = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    fn u16(&mut self) -> Result<u16, BucketError> {
        Ok(u16::from_le_bytes(self.take(2)?.try_into().unwrap()))
    }

    fn u32(&mut self) -> Result<u32, BucketError> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }

    fn u64(&mut self) -> Result<u64, BucketError> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
}

/// Reads the container manifest, if this directory has one.
pub fn read_container_manifest(
    bucket_dir: impl AsRef<Path>,
) -> Result<Option<(ContainerManifestHeader, Vec<BucketManifestEntry>)>, BucketError> {
    let path = bucket_dir.as_ref().join(CONTAINER_MANIFEST_NAME);
    let bytes = match fs::read(&path) {
        Ok(bytes) => bytes,
        Err(source) if source.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(source) => {
            return Err(BucketError::Io {
                path: path.clone(),
                source,
            });
        }
    };
    let mut cursor = ManifestCursor {
        bytes: &bytes,
        pos: 0,
        path: &path,
    };
    if cursor.take(8)? != CONTAINER_MANIFEST_MAGIC {
        return Err(BucketError::MalformedManifest(path.clone()));
    }
    let k = cursor.u16()?;
    let minimizer_len = cursor.u16()?;
    let graph_count = cursor.u64()? as usize;
    let flags = cursor.take(4)?;
    let (colored, label_words, compression) = (flags[0], flags[1], flags[2]);
    if colored > 1 || compression > 2 || flags[3] != 0 {
        return Err(BucketError::MalformedManifest(path.clone()));
    }
    let header = ContainerManifestHeader {
        k,
        minimizer_len,
        graph_count,
        colored: colored == 1,
        label_words: label_words as usize,
        compressed: compression != 0,
        interleaved_compression: compression == 2,
        segment_bytes: cursor.u64()?,
        container_count: cursor.u64()? as usize,
    };
    if header.segment_bytes == 0
        || header.label_words as usize != label_word_count(k, minimizer_len)
    {
        return Err(BucketError::MalformedManifest(path.clone()));
    }
    let bucket_count = cursor.u64()? as usize;
    let mut entries = Vec::with_capacity(bucket_count);
    for _ in 0..bucket_count {
        let graph_id = cursor.u64()? as usize;
        let records = cursor.u64()?;
        let bytes_len = cursor.u64()?;
        let container = cursor.u32()? as usize;
        let segment_count = cursor.u32()? as usize;
        if graph_id >= graph_count || container >= header.container_count {
            return Err(BucketError::MalformedManifest(path.clone()));
        }
        let mut segments = Vec::with_capacity(segment_count);
        for _ in 0..segment_count {
            segments.push(cursor.u32()?);
        }
        // The chain must be able to hold the logical length, and must not be
        // more than one segment longer than it needs to be.
        let capacity = segment_count as u64 * header.segment_bytes;
        if bytes_len > capacity || capacity.saturating_sub(bytes_len) >= header.segment_bytes {
            return Err(BucketError::MalformedManifest(path.clone()));
        }
        entries.push(BucketManifestEntry {
            graph_id,
            records,
            location: BucketLocation::Container {
                container,
                segments,
                bytes: bytes_len,
            },
        });
    }
    Ok(Some((header, entries)))
}

/// Reads and validates `manifest.tsv` from a bucket directory.
pub fn read_manifest(
    bucket_dir: impl AsRef<Path>,
) -> Result<Vec<BucketManifestEntry>, BucketError> {
    let bucket_dir = bucket_dir.as_ref();
    let path = bucket_dir.join("manifest.tsv");
    let file = File::open(&path).map_err(|source| BucketError::Io {
        path: path.clone(),
        source,
    })?;
    let mut out = Vec::new();

    for (line_no, line) in BufReader::new(file).lines().enumerate() {
        let line = line.map_err(|source| BucketError::Io {
            path: path.clone(),
            source,
        })?;
        if line_no == 0 {
            if line != "graph_id\trecords\tpath" {
                return Err(BucketError::MalformedManifest(path.clone()));
            }
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }

        let mut fields = line.splitn(3, '\t');
        let graph_id = fields
            .next()
            .ok_or_else(|| BucketError::MalformedManifest(path.clone()))?
            .parse()
            .map_err(|_| BucketError::MalformedManifest(path.clone()))?;
        let records = fields
            .next()
            .ok_or_else(|| BucketError::MalformedManifest(path.clone()))?
            .parse()
            .map_err(|_| BucketError::MalformedManifest(path.clone()))?;
        let bucket_path = fields
            .next()
            .ok_or_else(|| BucketError::MalformedManifest(path.clone()))?;
        out.push(BucketManifestEntry {
            graph_id,
            records,
            location: BucketLocation::File(PathBuf::from(bucket_path)),
        });
    }

    Ok(out)
}

impl BucketEmitter {
    pub fn create(params: &BuildParams, graph_count: usize) -> Result<Self, BucketError> {
        Self::create_in_dir(params, graph_count, bucket_dir(params))
    }

    pub fn create_in_dir(
        params: &BuildParams,
        graph_count: usize,
        bucket_dir: PathBuf,
    ) -> Result<Self, BucketError> {
        if graph_count > u64::MAX as usize {
            return Err(BucketError::GraphCountTooLarge(graph_count));
        }

        if bucket_dir.exists() {
            fs::remove_dir_all(&bucket_dir).map_err(|source| BucketError::Io {
                path: bucket_dir.clone(),
                source,
            })?;
        }
        fs::create_dir_all(&bucket_dir).map_err(|source| BucketError::Io {
            path: bucket_dir.clone(),
            source,
        })?;

        Ok(Self {
            bucket_dir,
            k: params.k,
            minimizer_len: params.minimizer_len,
            graph_count,
            colored: params.color,
            compress_buckets: params.color || params.compress_buckets,
            label_words: label_word_count(params.k, params.minimizer_len),
            files: BTreeMap::new(),
            scratch: CompressionScratch::default(),
            writers: BTreeMap::new(),
            pending: vec![PendingBucket::default(); graph_count],
            pending_bytes: 0,
        })
    }

    pub fn add(&mut self, superkmer: &WeakSuperKmer, seq: &[u8]) -> Result<(), BucketError> {
        if superkmer.graph_id >= self.graph_count {
            return Err(BucketError::InvalidGraphId(superkmer.graph_id));
        }
        if seq.len() > u8::MAX as usize {
            return Err(BucketError::LabelTooLong(seq.len()));
        }

        let attr = if self.colored {
            let source_id = superkmer.source_id.ok_or(BucketError::MissingSourceId)?;
            if source_id > MAX_SOURCE_ID {
                return Err(BucketError::SourceIdTooLarge(source_id));
            }
            pack_colored_attr(
                seq.len(),
                source_id,
                superkmer.left_discontinuous,
                superkmer.right_discontinuous,
            )
        } else {
            pack_uncolored_attr(
                seq.len(),
                superkmer.left_discontinuous,
                superkmer.right_discontinuous,
            )
        };
        let graph_id = superkmer.graph_id;
        let record_len = record_size(self.colored, self.label_words);
        let pending = &mut self.pending[graph_id];
        pending.records = pending
            .records
            .checked_add(1)
            .ok_or(BucketError::TooManyRecords)?;
        append_record(
            &mut pending.bytes,
            attr,
            graph_id,
            seq,
            self.label_words,
            self.colored,
        )?;
        self.pending_bytes += record_len;

        if self.pending[graph_id].bytes.len() >= MAX_PENDING_BUCKET_BYTES {
            self.flush_pending_bucket(graph_id)?;
        } else if self.pending_bytes >= MAX_TOTAL_PENDING_BYTES {
            self.flush_largest_pending_bucket()?;
        }
        Ok(())
    }

    fn flush_largest_pending_bucket(&mut self) -> Result<(), BucketError> {
        let Some((graph_id, _)) = self
            .pending
            .iter()
            .enumerate()
            .max_by_key(|(_, pending)| pending.bytes.len())
        else {
            return Ok(());
        };
        self.flush_pending_bucket(graph_id)
    }

    fn flush_pending_bucket(&mut self, graph_id: usize) -> Result<(), BucketError> {
        if self.pending[graph_id].bytes.is_empty() {
            return Ok(());
        }
        let pending = std::mem::take(&mut self.pending[graph_id]);
        self.pending_bytes -= pending.bytes.len();

        let (records, bytes_written) = {
            self.ensure_writer(graph_id)?;
            let Self {
                writers, scratch, ..
            } = self;
            let writer = writers.get_mut(&graph_id).expect("writer just ensured");
            writer.write_records(&pending.bytes, pending.records, scratch)?
        };
        let meta = self.files.get_mut(&graph_id).unwrap();
        meta.records = records;
        meta.bytes_written = bytes_written;
        Ok(())
    }

    fn ensure_writer(&mut self, graph_id: usize) -> Result<&mut BucketFile, BucketError> {
        if !self.writers.contains_key(&graph_id) {
            self.evict_writer_if_needed(graph_id)?;
            let writer = match self.files.get(&graph_id) {
                Some(meta) => {
                    BucketFile::open_existing(&meta.path, meta.records, meta.bytes_written)?
                }
                None => {
                    let writer = BucketFile::create(
                        &self.bucket_dir,
                        self.k,
                        self.minimizer_len,
                        self.graph_count,
                        graph_id,
                        self.colored,
                        self.label_words,
                        self.compress_buckets,
                    )?;
                    self.files.insert(
                        graph_id,
                        BucketFileMeta {
                            path: writer.path.clone(),
                            records: writer.records,
                            bytes_written: writer.bytes_written,
                        },
                    );
                    writer
                }
            };
            self.writers.insert(graph_id, writer);
        }

        Ok(self.writers.get_mut(&graph_id).unwrap())
    }

    fn evict_writer_if_needed(&mut self, requested_graph_id: usize) -> Result<(), BucketError> {
        if self.writers.len() < MAX_OPEN_BUCKET_WRITERS {
            return Ok(());
        }

        let evict_graph_id = self
            .writers
            .keys()
            .copied()
            .find(|&graph_id| graph_id != requested_graph_id)
            .unwrap_or(requested_graph_id);
        if let Some(mut writer) = self.writers.remove(&evict_graph_id) {
            writer.flush()?;
        }
        Ok(())
    }

    pub fn finish(mut self) -> Result<BucketEmitStats, BucketError> {
        let mut manifest = Vec::new();
        for graph_id in 0..self.pending.len() {
            if !self.pending[graph_id].bytes.is_empty() {
                self.flush_pending_bucket(graph_id)?;
            }
        }

        for (graph_id, meta) in &self.files {
            if let Some(writer) = self.writers.get_mut(graph_id) {
                writer.finish()?;
            } else {
                BucketFile::finish_closed(&meta.path, meta.records)?;
            }
            manifest.push((*graph_id, meta.records, meta.path.clone()));
        }

        write_manifest(&self.bucket_dir, &manifest)?;

        Ok(BucketEmitStats {
            bucket_dir: self.bucket_dir,
            bucket_files: manifest.len(),
            bytes_written: self.files.values().map(|meta| meta.bytes_written).sum(),
        })
    }
}

impl SharedBucketSink {
    pub fn create(params: &BuildParams, graph_count: usize) -> Result<Arc<Self>, BucketError> {
        let bucket_dir = bucket_dir(params);
        if bucket_dir.exists() {
            fs::remove_dir_all(&bucket_dir).map_err(|source| BucketError::Io {
                path: bucket_dir.clone(),
                source,
            })?;
        }
        fs::create_dir_all(&bucket_dir).map_err(|source| BucketError::Io {
            path: bucket_dir.clone(),
            source,
        })?;

        // One container per atlas where descriptors allow. An atlas
        // serializes its own writes with a mutex, so at one-to-one a container
        // has a single writer; when several atlases share one, the atomic
        // segment cursor keeps that safe.
        let atlas_count = graph_count.div_ceil(ATLAS_GRAPH_COUNT);
        let container_count = BucketContainers::plan_container_count(atlas_count);
        if container_count < atlas_count {
            eprintln!(
                "cuttlefish: descriptor budget allows {container_count} bucket container(s) rather than {atlas_count}"
            );
        }
        let containers = BucketContainers::create(&bucket_dir, container_count)?;

        Ok(Arc::new(Self {
            bucket_dir,
            containers,
            k: params.k,
            minimizer_len: params.minimizer_len,
            graph_count,
            colored: params.color,
            compress_buckets: params.color || params.compress_buckets,
            label_words: label_word_count(params.k, params.minimizer_len),
            workers: params.threads.max(1),
            atlases: (0..graph_count.div_ceil(ATLAS_GRAPH_COUNT))
                .map(|atlas_id| {
                    let first_graph_id = atlas_id * ATLAS_GRAPH_COUNT;
                    let file_count = (graph_count - first_graph_id).min(ATLAS_GRAPH_COUNT);
                    Mutex::new(SharedBucketAtlas {
                        first_graph_id,
                        files: (0..file_count)
                            .map(|_| SharedBucketFileMeta::default())
                            .collect(),
                        buffered_bytes: 0,
                        scratch: CompressionScratch::default(),
                    })
                })
                .collect(),
            flush_calls: AtomicU64::new(0),
            flush_nanos: AtomicU64::new(0),
        }))
    }

    pub fn emitter(self: &Arc<Self>) -> SharedBucketEmitter {
        self.make_emitter(false)
    }

    pub fn deferred_uncolored_emitter(self: &Arc<Self>) -> SharedBucketEmitter {
        self.make_emitter(true)
    }

    fn make_emitter(self: &Arc<Self>, deferred_uncolored: bool) -> SharedBucketEmitter {
        SharedBucketEmitter {
            sink: Arc::clone(self),
            pending: if self.colored || deferred_uncolored {
                Vec::new()
            } else {
                vec![PendingBucket::default(); self.graph_count]
            },
            uncolored_pending: if !self.colored && deferred_uncolored {
                (0..self.atlases.len())
                    .map(|_| PendingColoredAtlas::default())
                    .collect()
            } else {
                Vec::new()
            },
            colored_pending: if self.colored {
                (0..self.atlases.len())
                    .map(|_| PendingColoredAtlas::default())
                    .collect()
            } else {
                Vec::new()
            },
            pending_bytes: 0,
            deferred_uncolored,
        }
    }

    pub fn flush_uncolored_emitters(
        &self,
        emitters: Vec<SharedBucketEmitter>,
    ) -> Result<(), BucketError> {
        let started = Instant::now();
        let mut by_atlas = (0..self.atlases.len())
            .map(|_| Mutex::new(Vec::<PendingColoredAtlas>::new()))
            .collect::<Vec<_>>();
        for mut emitter in emitters {
            for (atlas_id, pending) in emitter.uncolored_pending.drain(..).enumerate() {
                if !pending.bytes.is_empty() {
                    by_atlas[atlas_id]
                        .get_mut()
                        .map_err(|_| BucketError::WorkerPanic)?
                        .push(pending);
                }
            }
        }
        let next = AtomicUsize::new(0);
        let workers = self.workers.min(self.atlases.len().max(1));
        let calls = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(workers);
            for _ in 0..workers {
                handles.push(scope.spawn(|| {
                    let mut calls = 0;
                    loop {
                        let atlas_id = next.fetch_add(1, Ordering::Relaxed);
                        let Some(pending) = by_atlas.get(atlas_id) else {
                            break;
                        };
                        let pending = std::mem::take(
                            &mut *pending.lock().map_err(|_| BucketError::WorkerPanic)?,
                        );
                        if pending.is_empty() {
                            continue;
                        }
                        let mut stats = SharedBucketFlushStats::default();
                        for chunk in pending {
                            self.append_uncolored_atlas_inner(atlas_id, chunk, &mut stats)?;
                        }
                        let mut atlas = self.atlases[atlas_id]
                            .lock()
                            .map_err(|_| BucketError::WorkerPanic)?;
                        atlas.flush_all(
                            &self.containers,
                            false,
                            self.label_words,
                            self.compress_buckets,
                            &mut stats,
                        )?;
                        calls += stats.calls;
                    }
                    Ok::<_, BucketError>(calls)
                }));
            }
            let mut calls = 0;
            for handle in handles {
                calls += handle.join().map_err(|_| BucketError::WorkerPanic)??;
            }
            Ok::<_, BucketError>(calls)
        })?;
        self.record_flush_stats(SharedBucketFlushStats { calls }, started.elapsed());
        Ok(())
    }

    fn append_uncolored_atlas(
        &self,
        atlas_id: usize,
        pending: PendingColoredAtlas,
    ) -> Result<(), BucketError> {
        let started = Instant::now();
        let mut stats = SharedBucketFlushStats::default();
        self.append_uncolored_atlas_inner(atlas_id, pending, &mut stats)?;
        self.record_flush_stats(stats, started.elapsed());
        Ok(())
    }

    fn append_uncolored_atlas_inner(
        &self,
        atlas_id: usize,
        pending: PendingColoredAtlas,
        stats: &mut SharedBucketFlushStats,
    ) -> Result<(), BucketError> {
        if pending.bytes.is_empty() {
            return Ok(());
        }
        let record_len = record_size(false, self.label_words);
        if pending.graph_ids.len() * record_len != pending.bytes.len() {
            return Err(BucketError::MalformedRecord);
        }
        let mut atlas = self.atlases[atlas_id]
            .lock()
            .map_err(|_| BucketError::WorkerPanic)?;
        for (&graph_id, record) in pending
            .graph_ids
            .iter()
            .zip(pending.bytes.chunks_exact(record_len))
        {
            let graph_id = usize::from(graph_id);
            let Some(local_graph_id) = graph_id.checked_sub(atlas.first_graph_id) else {
                return Err(BucketError::InvalidGraphId(graph_id));
            };
            let Some(file) = atlas.files.get_mut(local_graph_id) else {
                return Err(BucketError::InvalidGraphId(graph_id));
            };
            file.total_records = file
                .total_records
                .checked_add(1)
                .ok_or(BucketError::TooManyRecords)?;
            file.buffer_records = file
                .buffer_records
                .checked_add(1)
                .ok_or(BucketError::TooManyRecords)?;
            file.buffer.extend_from_slice(record);
        }
        atlas.buffered_bytes += pending.bytes.len();
        for local_graph_id in 0..atlas.files.len() {
            if atlas.files[local_graph_id].buffer.len() >= SUBGRAPH_CHUNK_BYTES {
                atlas.flush_subgraph(
                    local_graph_id,
                    &self.containers,
                    false,
                    self.label_words,
                    self.compress_buckets,
                    stats,
                )?;
            }
        }
        Ok(())
    }

    fn append_bucket(&self, graph_id: usize, pending: PendingBucket) -> Result<(), BucketError> {
        if pending.bytes.is_empty() {
            return Ok(());
        }
        let started = Instant::now();
        let atlas_id = graph_id / ATLAS_GRAPH_COUNT;
        let local_graph_id = graph_id % ATLAS_GRAPH_COUNT;
        let mut atlas = self.atlases[atlas_id]
            .lock()
            .map_err(|_| BucketError::WorkerPanic)?;
        atlas.append_bucket(local_graph_id, pending)?;
        if self.colored {
            return Ok(());
        }
        if atlas.files[local_graph_id].buffer.len() < SUBGRAPH_CHUNK_BYTES {
            return Ok(());
        }

        let mut flush_stats = SharedBucketFlushStats::default();
        atlas.flush_subgraph(
            local_graph_id,
            &self.containers,
            self.colored,
            self.label_words,
            self.compress_buckets,
            &mut flush_stats,
        )?;
        self.record_flush_stats(flush_stats, started.elapsed());
        Ok(())
    }

    fn append_colored_atlas(
        &self,
        atlas_id: usize,
        pending: PendingColoredAtlas,
    ) -> Result<(), BucketError> {
        if pending.bytes.is_empty() {
            return Ok(());
        }
        let record_len = record_size(true, self.label_words);
        if pending.graph_ids.len() * record_len != pending.bytes.len() {
            return Err(BucketError::MalformedRecord);
        }
        let started = Instant::now();
        let mut atlas = self.atlases[atlas_id]
            .lock()
            .map_err(|_| BucketError::WorkerPanic)?;
        for (&graph_id, record) in pending
            .graph_ids
            .iter()
            .zip(pending.bytes.chunks_exact(record_len))
        {
            let graph_id = usize::from(graph_id);
            let Some(local_graph_id) = graph_id.checked_sub(atlas.first_graph_id) else {
                return Err(BucketError::InvalidGraphId(graph_id));
            };
            let Some(file) = atlas.files.get_mut(local_graph_id) else {
                return Err(BucketError::InvalidGraphId(graph_id));
            };
            file.total_records = file
                .total_records
                .checked_add(1)
                .ok_or(BucketError::TooManyRecords)?;
            file.buffer_records = file
                .buffer_records
                .checked_add(1)
                .ok_or(BucketError::TooManyRecords)?;
            file.buffer.extend_from_slice(record);
        }
        atlas.buffered_bytes += pending.bytes.len();
        let mut flush_stats = SharedBucketFlushStats::default();
        for local_graph_id in 0..atlas.files.len() {
            if atlas.files[local_graph_id].buffer.len() >= SUBGRAPH_CHUNK_BYTES {
                atlas.flush_subgraph(
                    local_graph_id,
                    &self.containers,
                    true,
                    self.label_words,
                    self.compress_buckets,
                    &mut flush_stats,
                )?;
            }
        }
        self.record_flush_stats(flush_stats, started.elapsed());
        Ok(())
    }

    pub fn flush_stats(&self) -> (u64, Duration) {
        (
            self.flush_calls.load(Ordering::Relaxed),
            Duration::from_nanos(self.flush_nanos.load(Ordering::Relaxed)),
        )
    }

    pub fn flush_colored_window(
        &self,
        source_min: u32,
        source_max: u32,
    ) -> Result<(), BucketError> {
        if !self.colored || source_min > source_max {
            return Err(BucketError::MalformedRecord);
        }
        let started = Instant::now();
        let record_len = record_size(true, self.label_words);
        let next_atlas = AtomicUsize::new(0);
        let workers = self.workers.min(self.atlases.len().max(1));
        let flush_calls = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(workers);
            for _ in 0..workers {
                handles.push(scope.spawn(|| {
                    let mut stats = SharedBucketFlushStats::default();
                    loop {
                        let atlas_id = next_atlas.fetch_add(1, Ordering::Relaxed);
                        let Some(atlas) = self.atlases.get(atlas_id) else {
                            break;
                        };
                        let mut atlas = atlas.lock().map_err(|_| BucketError::WorkerPanic)?;
                        for local_graph_id in 0..atlas.files.len() {
                            sort_colored_payload_by_source(
                                &mut atlas.files[local_graph_id].buffer,
                                record_len,
                                source_min,
                                source_max,
                            )?;
                            atlas.flush_subgraph(
                                local_graph_id,
                                &self.containers,
                                true,
                                self.label_words,
                                true,
                                &mut stats,
                            )?;
                        }
                    }
                    Ok::<_, BucketError>(stats.calls)
                }));
            }
            let mut calls = 0u64;
            for handle in handles {
                calls += handle.join().map_err(|_| BucketError::WorkerPanic)??;
            }
            Ok::<_, BucketError>(calls)
        })?;
        self.record_flush_stats(
            SharedBucketFlushStats { calls: flush_calls },
            started.elapsed(),
        );
        Ok(())
    }

    pub fn flush_colored_emitters(
        &self,
        emitters: Vec<SharedBucketEmitter>,
    ) -> Result<(), BucketError> {
        if !self.colored {
            return Err(BucketError::MalformedRecord);
        }

        let emitter_count = emitters.len();
        let mut pending_by_atlas = (0..self.atlases.len())
            .map(|_| Vec::with_capacity(emitter_count))
            .collect::<Vec<Vec<PendingColoredAtlas>>>();
        for emitter in emitters {
            if !std::ptr::eq(Arc::as_ptr(&emitter.sink), self) {
                return Err(BucketError::MalformedRecord);
            }
            if !emitter.pending.is_empty()
                || emitter.colored_pending.len() != pending_by_atlas.len()
            {
                return Err(BucketError::MalformedRecord);
            }
            for (atlas_id, pending) in emitter.colored_pending.into_iter().enumerate() {
                if !pending.bytes.is_empty() {
                    pending_by_atlas[atlas_id].push(pending);
                }
            }
        }

        let started = Instant::now();
        let workers = self.workers.min(pending_by_atlas.len().max(1));
        let chunk_size = pending_by_atlas.len().div_ceil(workers);
        let flush_calls = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(workers);
            for (chunk_id, atlas_work) in pending_by_atlas.chunks_mut(chunk_size).enumerate() {
                handles.push(scope.spawn(move || {
                    let mut stats = SharedBucketFlushStats::default();
                    let first_atlas = chunk_id * chunk_size;
                    for (offset, worker_buckets) in atlas_work.iter_mut().enumerate() {
                        let atlas_id = first_atlas + offset;
                        let mut atlas = self.atlases[atlas_id]
                            .lock()
                            .map_err(|_| BucketError::WorkerPanic)?;
                        let record_len = record_size(true, self.label_words);
                        for pending in worker_buckets.iter() {
                            if pending.graph_ids.len() * record_len != pending.bytes.len() {
                                return Err(BucketError::MalformedRecord);
                            }
                            for (&graph_id, record) in pending
                                .graph_ids
                                .iter()
                                .zip(pending.bytes.chunks_exact(record_len))
                            {
                                append_colored_atlas_record(&mut atlas, graph_id, record)?;
                            }
                        }
                        atlas.buffered_bytes += worker_buckets
                            .iter()
                            .map(|pending| pending.bytes.len())
                            .sum::<usize>();
                        for local_graph_id in 0..atlas.files.len() {
                            atlas.flush_subgraph(
                                local_graph_id,
                                &self.containers,
                                true,
                                self.label_words,
                                true,
                                &mut stats,
                            )?;
                            // This consumes the complete colored emitter set, so these
                            // buffers will not be reused. Releasing their capacity here
                            // avoids retaining every uncompressed subgraph bucket while
                            // the remaining worker atlas chunks are still resident.
                            atlas.files[local_graph_id].buffer = Vec::new();
                        }
                    }
                    Ok::<_, BucketError>(stats.calls)
                }));
            }
            let mut calls = 0u64;
            for handle in handles {
                calls += handle.join().map_err(|_| BucketError::WorkerPanic)??;
            }
            Ok::<_, BucketError>(calls)
        })?;
        self.record_flush_stats(
            SharedBucketFlushStats { calls: flush_calls },
            started.elapsed(),
        );
        Ok(())
    }

    pub fn finish(&self) -> Result<BucketEmitStats, BucketError> {
        let mut entries = Vec::new();
        let mut finish_flush_stats = SharedBucketFlushStats::default();
        let finish_started = Instant::now();
        let mut bytes_written = 0u64;
        for atlas in &self.atlases {
            let mut atlas = atlas.lock().map_err(|_| BucketError::WorkerPanic)?;
            atlas.flush_all(
                &self.containers,
                self.colored,
                self.label_words,
                self.compress_buckets,
                &mut finish_flush_stats,
            )?;
            let first_graph_id = atlas.first_graph_id;
            let container = (first_graph_id / ATLAS_GRAPH_COUNT) % self.containers.len();
            for (local_graph_id, meta) in atlas.files.iter_mut().enumerate() {
                let meta = std::mem::take(meta);
                if meta.total_records == 0 {
                    continue;
                }
                bytes_written += meta.bytes_written;
                entries.push(BucketManifestEntry {
                    graph_id: first_graph_id + local_graph_id,
                    records: meta.total_records,
                    location: BucketLocation::Container {
                        container,
                        segments: meta.segments,
                        bytes: meta.bytes_written,
                    },
                });
            }
        }
        self.record_flush_stats(finish_flush_stats, finish_started.elapsed());
        entries.sort_unstable_by_key(|entry| entry.graph_id);

        // Nothing to patch and nothing to reopen. The per-file layout finished
        // by reopening all 16,384 buckets to write the record count into each
        // header, parallelised across workers because it was slow enough to
        // matter; the count now lives in the manifest that has to be written
        // anyway.
        let header = ContainerManifestHeader {
            k: self.k,
            minimizer_len: self.minimizer_len,
            graph_count: self.graph_count,
            colored: self.colored,
            label_words: self.label_words,
            compressed: self.compress_buckets,
            interleaved_compression: self.compress_buckets
                && !self.colored
                && !force_split_compression(),
            segment_bytes: self.containers.segment_bytes(),
            // The number of containers actually created, which is not the
            // atlas count when the descriptor budget narrowed it.
            container_count: self.containers.len(),
        };
        write_container_manifest(&self.bucket_dir, &header, &entries)?;
        Ok(BucketEmitStats {
            bucket_dir: self.bucket_dir.clone(),
            bucket_files: entries.len(),
            bytes_written,
        })
    }

    fn record_flush_stats(&self, stats: SharedBucketFlushStats, elapsed: Duration) {
        if stats.calls == 0 {
            return;
        }
        self.flush_calls.fetch_add(stats.calls, Ordering::Relaxed);
        self.flush_nanos.fetch_add(
            u64::try_from(elapsed.as_nanos()).unwrap_or(u64::MAX),
            Ordering::Relaxed,
        );
    }
}

fn append_colored_atlas_record(
    atlas: &mut SharedBucketAtlas,
    graph_id: u16,
    record: &[u8],
) -> Result<(), BucketError> {
    let graph_id = usize::from(graph_id);
    let local_graph_id = graph_id
        .checked_sub(atlas.first_graph_id)
        .ok_or(BucketError::InvalidGraphId(graph_id))?;
    let file = atlas
        .files
        .get_mut(local_graph_id)
        .ok_or(BucketError::InvalidGraphId(graph_id))?;
    file.total_records = file
        .total_records
        .checked_add(1)
        .ok_or(BucketError::TooManyRecords)?;
    file.buffer_records = file
        .buffer_records
        .checked_add(1)
        .ok_or(BucketError::TooManyRecords)?;
    file.buffer.extend_from_slice(record);
    Ok(())
}

fn sort_colored_payload_by_source(
    payload: &mut Vec<u8>,
    record_len: usize,
    source_min: u32,
    source_max: u32,
) -> Result<(), BucketError> {
    if payload.is_empty() {
        return Ok(());
    }
    if payload.len() % record_len != 0 {
        return Err(BucketError::MalformedRecord);
    }
    let source_count =
        usize::try_from(source_max - source_min + 1).map_err(|_| BucketError::MalformedRecord)?;
    let mut offsets = vec![0usize; source_count];
    for record in payload.chunks_exact(record_len) {
        let attr = u32::from_le_bytes(record[..4].try_into().expect("colored attribute"));
        let source = attr >> 10;
        if source < source_min || source > source_max {
            return Err(BucketError::MalformedRecord);
        }
        offsets[(source - source_min) as usize] += 1;
    }
    let mut prefix = 0usize;
    for offset in &mut offsets {
        let count = *offset;
        *offset = prefix;
        prefix += count;
    }
    let mut sorted = vec![0u8; payload.len()];
    for record in payload.chunks_exact(record_len) {
        let attr = u32::from_le_bytes(record[..4].try_into().expect("colored attribute"));
        let source = ((attr >> 10) - source_min) as usize;
        let output = offsets[source] * record_len;
        sorted[output..output + record_len].copy_from_slice(record);
        offsets[source] += 1;
    }
    *payload = sorted;
    Ok(())
}

impl SharedBucketAtlas {
    fn append_bucket(
        &mut self,
        local_graph_id: usize,
        pending: PendingBucket,
    ) -> Result<(), BucketError> {
        let file = &mut self.files[local_graph_id];
        file.total_records = file
            .total_records
            .checked_add(pending.records)
            .ok_or(BucketError::TooManyRecords)?;
        file.buffer_records = file
            .buffer_records
            .checked_add(pending.records)
            .ok_or(BucketError::TooManyRecords)?;
        self.buffered_bytes += pending.bytes.len();
        file.buffer.extend_from_slice(&pending.bytes);
        Ok(())
    }

    fn flush_all(
        &mut self,
        containers: &BucketContainers,
        colored: bool,
        label_words: usize,
        compress_buckets: bool,
        stats: &mut SharedBucketFlushStats,
    ) -> Result<(), BucketError> {
        for local_graph_id in 0..self.files.len() {
            self.flush_subgraph(
                local_graph_id,
                containers,
                colored,
                label_words,
                compress_buckets,
                stats,
            )?;
        }
        Ok(())
    }

    /// Writes one bucket's staged records into its container.
    ///
    /// This is the whole syscall saving. The per-file path reached here with
    /// an `openat`, seven unbuffered reads to re-read the 42-byte header, a
    /// revalidation, an `lseek` to the end and a `close` -- around eleven
    /// syscalls for every 64 KiB flush, and a full-corpus build performs 14.4
    /// million of them. A container flush is the `pwrite` and nothing else:
    /// the header is in the manifest, the descriptor is already open, and the
    /// offset is known rather than sought.
    fn flush_subgraph(
        &mut self,
        local_graph_id: usize,
        containers: &BucketContainers,
        colored: bool,
        label_words: usize,
        compress_buckets: bool,
        stats: &mut SharedBucketFlushStats,
    ) -> Result<(), BucketError> {
        let container = (self.first_graph_id / ATLAS_GRAPH_COUNT) % containers.len();
        let Self { files, scratch, .. } = self;
        let file = &mut files[local_graph_id];
        if file.buffer.is_empty() {
            return Ok(());
        }
        let flushed_bytes = file.buffer.len();
        let record_size = record_size(colored, label_words) as u64;

        let written = if compress_buckets {
            let interleaved = !colored && !force_split_compression();
            let len = encode_compressed_block(
                &file.buffer,
                file.buffer_records,
                record_size,
                label_words,
                interleaved,
                scratch,
            )?;
            let block = std::mem::take(&mut scratch.block);
            let written = append_to_chain(containers, container, file, &block)?;
            scratch.block = block;
            debug_assert_eq!(written, len as u64);
            written
        } else {
            let buffer = std::mem::take(&mut file.buffer);
            let written = append_to_chain(containers, container, file, &buffer);
            file.buffer = buffer;
            written?
        };

        file.written_records += file.buffer_records;
        file.bytes_written += written;
        file.buffer.clear();
        file.buffer_records = 0;
        self.buffered_bytes -= flushed_bytes;
        stats.calls += 1;
        Ok(())
    }
}

/// Appends `bytes` to a bucket's segment chain, reserving as it goes.
///
/// A block may straddle a segment boundary rather than starting a fresh
/// segment. That costs a second `pwrite` for the rare block that spans one,
/// and saves refusing to fill the tail of every segment -- which at a 64 KiB
/// flush and a 256 KiB segment would waste up to a quarter of the directory.
/// The reader concatenates the chain in order, so a split block reassembles.
fn append_to_chain(
    containers: &BucketContainers,
    container: usize,
    file: &mut SharedBucketFileMeta,
    bytes: &[u8],
) -> Result<u64, BucketError> {
    let segment_bytes = containers.segment_bytes();
    let mut written = 0usize;
    while written < bytes.len() {
        if file.segments.is_empty() || file.segment_used == segment_bytes {
            let offset = containers.reserve_segment(container);
            let index =
                u32::try_from(offset / segment_bytes).map_err(|_| BucketError::TooManyRecords)?;
            file.segments.push(index);
            file.segment_used = 0;
        }
        let room = (segment_bytes - file.segment_used) as usize;
        let take = room.min(bytes.len() - written);
        let offset = u64::from(*file.segments.last().unwrap()) * segment_bytes + file.segment_used;
        containers.write_at(container, offset, &bytes[written..written + take])?;
        written += take;
        file.segment_used += take as u64;
    }
    Ok(written as u64)
}

impl SharedBucketEmitter {
    pub fn flush_colored_worker_if_required(&mut self) -> Result<(), BucketError> {
        if !self.sink.colored {
            return Err(BucketError::MalformedRecord);
        }
        for atlas_id in 0..self.colored_pending.len() {
            if self.colored_pending[atlas_id].bytes.len() >= SUBGRAPH_CHUNK_BYTES {
                self.flush_pending_colored_atlas(atlas_id)?;
            }
        }
        Ok(())
    }

    pub fn add(&mut self, superkmer: &WeakSuperKmer, seq: &[u8]) -> Result<(), BucketError> {
        self.add_impl(superkmer, seq, true)
    }

    pub fn add_valid(&mut self, superkmer: &WeakSuperKmer, seq: &[u8]) -> Result<(), BucketError> {
        debug_assert!(seq.iter().all(|&base| ascii_base_bits(base).is_some()));
        self.add_impl(superkmer, seq, false)
    }

    fn add_impl(
        &mut self,
        superkmer: &WeakSuperKmer,
        seq: &[u8],
        check_bases: bool,
    ) -> Result<(), BucketError> {
        if superkmer.graph_id >= self.sink.graph_count {
            return Err(BucketError::InvalidGraphId(superkmer.graph_id));
        }
        if seq.len() > u8::MAX as usize {
            return Err(BucketError::LabelTooLong(seq.len()));
        }

        let attr = if self.sink.colored {
            let source_id = superkmer.source_id.ok_or(BucketError::MissingSourceId)?;
            if source_id > MAX_SOURCE_ID {
                return Err(BucketError::SourceIdTooLarge(source_id));
            }
            pack_colored_attr(
                seq.len(),
                source_id,
                superkmer.left_discontinuous,
                superkmer.right_discontinuous,
            )
        } else {
            pack_uncolored_attr(
                seq.len(),
                superkmer.left_discontinuous,
                superkmer.right_discontinuous,
            )
        };
        let graph_id = superkmer.graph_id;
        let record_len = record_size(self.sink.colored, self.sink.label_words);
        if self.sink.colored {
            let atlas_id = graph_id / ATLAS_GRAPH_COUNT;
            let pending = &mut self.colored_pending[atlas_id];
            pending.graph_ids.push(graph_id as u16);
            if check_bases {
                append_record(
                    &mut pending.bytes,
                    attr,
                    graph_id,
                    seq,
                    self.sink.label_words,
                    true,
                )?;
            } else {
                append_record_valid(
                    &mut pending.bytes,
                    attr,
                    graph_id,
                    seq,
                    self.sink.label_words,
                    true,
                )?;
            }
        } else if self.deferred_uncolored {
            let atlas_id = graph_id / ATLAS_GRAPH_COUNT;
            let pending = &mut self.uncolored_pending[atlas_id];
            pending.graph_ids.push(graph_id as u16);
            if !check_bases {
                append_uncolored_record_valid(
                    &mut pending.bytes,
                    attr as u16,
                    graph_id as u16,
                    seq,
                    self.sink.label_words,
                )?;
            } else {
                append_record(
                    &mut pending.bytes,
                    attr,
                    graph_id,
                    seq,
                    self.sink.label_words,
                    false,
                )?;
            }
        } else {
            let pending = &mut self.pending[graph_id];
            pending.records = pending
                .records
                .checked_add(1)
                .ok_or(BucketError::TooManyRecords)?;
            if !check_bases {
                append_uncolored_record_valid(
                    &mut pending.bytes,
                    attr as u16,
                    graph_id as u16,
                    seq,
                    self.sink.label_words,
                )?;
            } else {
                append_record(
                    &mut pending.bytes,
                    attr,
                    graph_id,
                    seq,
                    self.sink.label_words,
                    false,
                )?;
            }
        }
        self.pending_bytes += record_len;

        // As in C++, colored worker-atlas chunks are checked and handed to the
        // shared atlas at the source boundary, not in the middle of a source.
        if self.deferred_uncolored {
            let atlas_id = graph_id / ATLAS_GRAPH_COUNT;
            if self.uncolored_pending[atlas_id].bytes.len() >= SUBGRAPH_CHUNK_BYTES {
                self.flush_pending_uncolored_atlas(atlas_id)?;
            }
        } else if !self.sink.colored
            && !self.deferred_uncolored
            && self.pending[graph_id].bytes.len() >= MAX_PENDING_BUCKET_BYTES
        {
            self.flush_pending_bucket(graph_id)?;
        } else if !self.sink.colored
            && !self.deferred_uncolored
            && self.pending_bytes >= MAX_TOTAL_PENDING_BYTES
        {
            self.flush_largest_pending_bucket()?;
        }
        Ok(())
    }

    fn flush_largest_pending_bucket(&mut self) -> Result<(), BucketError> {
        if self.sink.colored {
            let Some((atlas_id, _)) = self
                .colored_pending
                .iter()
                .enumerate()
                .max_by_key(|(_, pending)| pending.bytes.len())
            else {
                return Ok(());
            };
            return self.flush_pending_colored_atlas(atlas_id);
        }
        let Some((graph_id, _)) = self
            .pending
            .iter()
            .enumerate()
            .max_by_key(|(_, pending)| pending.bytes.len())
        else {
            return Ok(());
        };
        self.flush_pending_bucket(graph_id)
    }

    fn flush_pending_bucket(&mut self, graph_id: usize) -> Result<(), BucketError> {
        if self.pending[graph_id].bytes.is_empty() {
            return Ok(());
        }
        let pending = std::mem::take(&mut self.pending[graph_id]);
        self.pending_bytes -= pending.bytes.len();
        self.sink.append_bucket(graph_id, pending)
    }

    fn flush_pending_colored_atlas(&mut self, atlas_id: usize) -> Result<(), BucketError> {
        if self.colored_pending[atlas_id].bytes.is_empty() {
            return Ok(());
        }
        let pending = std::mem::take(&mut self.colored_pending[atlas_id]);
        self.pending_bytes -= pending.bytes.len();
        self.sink.append_colored_atlas(atlas_id, pending)
    }

    fn flush_pending_uncolored_atlas(&mut self, atlas_id: usize) -> Result<(), BucketError> {
        if self.uncolored_pending[atlas_id].bytes.is_empty() {
            return Ok(());
        }
        let pending = std::mem::take(&mut self.uncolored_pending[atlas_id]);
        self.pending_bytes -= pending.bytes.len();
        self.sink.append_uncolored_atlas(atlas_id, pending)
    }

    pub fn finish(mut self) -> Result<(), BucketError> {
        if self.sink.colored {
            for atlas_id in 0..self.colored_pending.len() {
                if !self.colored_pending[atlas_id].bytes.is_empty() {
                    self.flush_pending_colored_atlas(atlas_id)?;
                }
            }
        } else if self.deferred_uncolored {
            for atlas_id in 0..self.uncolored_pending.len() {
                if !self.uncolored_pending[atlas_id].bytes.is_empty() {
                    self.flush_pending_uncolored_atlas(atlas_id)?;
                }
            }
        } else {
            for graph_id in 0..self.pending.len() {
                if !self.pending[graph_id].bytes.is_empty() {
                    self.flush_pending_bucket(graph_id)?;
                }
            }
        }
        Ok(())
    }
}

pub fn bucket_dir(params: &BuildParams) -> PathBuf {
    let output_name = Path::new(&params.output_prefix)
        .file_name()
        .and_then(|s| s.to_str())
        .filter(|s| !s.is_empty())
        .unwrap_or("cuttlefish3");
    PathBuf::from(&params.work_dir).join(format!("{output_name}.cf3rs.wsk"))
}

pub fn write_manifest(
    bucket_dir: &Path,
    entries: &[(usize, u64, PathBuf)],
) -> Result<(), BucketError> {
    let path = bucket_dir.join("manifest.tsv");
    let mut out = File::create(&path).map_err(|source| BucketError::Io {
        path: path.clone(),
        source,
    })?;
    writeln!(out, "graph_id\trecords\tpath").map_err(|source| BucketError::Io {
        path: path.clone(),
        source,
    })?;
    for (graph_id, records, bucket_path) in entries {
        writeln!(out, "{graph_id}\t{records}\t{}", bucket_path.display()).map_err(|source| {
            BucketError::Io {
                path: path.clone(),
                source,
            }
        })?;
    }
    Ok(())
}

/// Encodes one compressed block into `scratch.block`, returning its length.
///
/// Assembling the whole block -- 12-byte header plus one or two LZ4 streams --
/// before it is written keeps a flush to a single `write` on a file that has
/// no buffering, and lets a container flush reuse the identical framing.
fn encode_compressed_block(
    bytes: &[u8],
    records: u64,
    record_size: u64,
    label_words: usize,
    interleaved: bool,
    scratch: &mut CompressionScratch,
) -> Result<usize, BucketError> {
    let records_u32 = u32::try_from(records).map_err(|_| BucketError::TooManyRecords)?;
    if records_u32 == 0 {
        scratch.block.clear();
        return Ok(0);
    }
    let (attr_len, label_len) = if interleaved {
        (CompressionScratch::encode(bytes, &mut scratch.encoded)?, 0)
    } else {
        let record_len = usize::try_from(record_size).unwrap();
        let fixed_len = record_len - label_words * 8;
        scratch.attrs.clear();
        scratch.labels.clear();
        scratch.attrs.reserve(records as usize * fixed_len);
        scratch.labels.reserve(records as usize * label_words * 8);
        for record in bytes.chunks_exact(record_len) {
            scratch.attrs.extend_from_slice(&record[..fixed_len]);
            scratch.labels.extend_from_slice(&record[fixed_len..]);
        }
        let attrs = std::mem::take(&mut scratch.attrs);
        let labels = std::mem::take(&mut scratch.labels);
        let attr_len = CompressionScratch::encode(&attrs, &mut scratch.encoded)?;
        let label_len = CompressionScratch::encode(&labels, &mut scratch.encoded_labels)?;
        scratch.attrs = attrs;
        scratch.labels = labels;
        (attr_len, label_len)
    };
    let attr_len_u32 = u32::try_from(attr_len).map_err(|_| BucketError::TooManyRecords)?;
    let label_len_u32 = u32::try_from(label_len).map_err(|_| BucketError::TooManyRecords)?;

    scratch.block.clear();
    scratch.block.extend_from_slice(&records_u32.to_le_bytes());
    scratch.block.extend_from_slice(&attr_len_u32.to_le_bytes());
    scratch
        .block
        .extend_from_slice(&label_len_u32.to_le_bytes());
    scratch
        .block
        .extend_from_slice(&scratch.encoded[..attr_len]);
    if label_len != 0 {
        scratch
            .block
            .extend_from_slice(&scratch.encoded_labels[..label_len]);
    }
    Ok(scratch.block.len())
}

fn append_record(
    out: &mut Vec<u8>,
    packed_attr: u32,
    _graph_id: usize,
    seq: &[u8],
    label_words: usize,
    colored: bool,
) -> Result<(), BucketError> {
    let mut words = [0u64; 4];
    if label_words > words.len() {
        return Err(BucketError::MalformedRecord);
    }
    for (idx, &ch) in seq.iter().enumerate() {
        let base_bits = ascii_base_bits(ch).ok_or(BucketError::InvalidBase(ch))?;
        let word_idx = idx / 32;
        let shift = 2 * (31 - (idx % 32));
        words[word_idx] |= (base_bits as u64) << shift;
    }

    if colored {
        out.extend_from_slice(&packed_attr.to_le_bytes());
    } else {
        out.extend_from_slice(&(packed_attr as u16).to_le_bytes());
    }
    for &word in &words[..label_words] {
        out.extend_from_slice(&word.to_le_bytes());
    }
    Ok(())
}

fn append_record_valid(
    out: &mut Vec<u8>,
    packed_attr: u32,
    _graph_id: usize,
    seq: &[u8],
    label_words: usize,
    colored: bool,
) -> Result<(), BucketError> {
    let words = pack_valid_label(seq, label_words)?;

    let mut record = [0u8; MAX_RECORD_BYTES];
    let attr_len = if colored {
        record[..4].copy_from_slice(&packed_attr.to_le_bytes());
        4
    } else {
        record[..2].copy_from_slice(&(packed_attr as u16).to_le_bytes());
        2
    };
    for (idx, &word) in words[..label_words].iter().enumerate() {
        let at = attr_len + idx * 8;
        record[at..at + 8].copy_from_slice(&word.to_le_bytes());
    }
    out.extend_from_slice(&record[..attr_len + label_words * 8]);
    Ok(())
}

#[inline]
fn append_uncolored_record_valid(
    out: &mut Vec<u8>,
    packed_attr: u16,
    _graph_id: u16,
    seq: &[u8],
    label_words: usize,
) -> Result<(), BucketError> {
    let words = pack_valid_label(seq, label_words)?;

    // Assemble the record on the stack and append it once. Reserving and then
    // extending per field costs a capacity check and a `memcpy` call for each of
    // the attribute and every label word; C++ writes its equivalent with plain
    // indexed stores into pre-reserved arrays.
    let mut record = [0u8; MAX_RECORD_BYTES];
    record[..2].copy_from_slice(&packed_attr.to_le_bytes());
    for (idx, &word) in words[..label_words].iter().enumerate() {
        let at = 2 + idx * 8;
        record[at..at + 8].copy_from_slice(&word.to_le_bytes());
    }
    out.extend_from_slice(&record[..2 + label_words * 8]);
    Ok(())
}

#[inline(always)]
/// Packs 32 valid ACGT bases into one 64-bit word, first base in the high bits.
///
/// The scalar form carries a loop-carried dependency (`word = (word << 2) | c`),
/// so it neither vectorizes nor pipelines: 32 serial iterations for one word,
/// which dominates record packing. Each 2-bit code is a pure bitwise function of
/// its byte, so eight bases can be reduced at a time inside a `u64` and the
/// codes gathered with a single `PEXT`.
fn pack_valid_word_32(seq: &[u8]) -> u64 {
    debug_assert!(seq.len() >= 32);
    #[cfg(all(target_arch = "x86_64", target_feature = "bmi2"))]
    {
        let mut word = 0u64;
        for chunk in 0..4 {
            let bytes = u64::from_le_bytes(
                seq[chunk * 8..chunk * 8 + 8]
                    .try_into()
                    .expect("eight bases"),
            );
            // Per byte: ((b >> 2) ^ (b >> 1)) & 0b11, evaluated eight at a time.
            let codes = ((bytes >> 2) ^ (bytes >> 1)) & 0x0303_0303_0303_0303;
            // Gather the eight 2-bit codes, lowest byte first, into 16 bits.
            let packed = unsafe { core::arch::x86_64::_pext_u64(codes, 0x0303_0303_0303_0303) };
            // The scalar form shifts the first base furthest left, and PEXT emits
            // the first (lowest-address) base in the least significant bits, so
            // reverse the 2-bit groups within this 16-bit lane.
            let reversed = reverse_2bit_groups_16(packed as u16);
            word = (word << 16) | u64::from(reversed);
        }
        word
    }
    #[cfg(not(all(target_arch = "x86_64", target_feature = "bmi2")))]
    {
        let mut word = 0u64;
        for &base in &seq[..32] {
            word = (word << 2) | u64::from(valid_ascii_base_bits(base));
        }
        word
    }
}

/// Reverses the order of the eight 2-bit groups in a 16-bit value.
#[cfg(all(target_arch = "x86_64", target_feature = "bmi2"))]
#[inline]
fn reverse_2bit_groups_16(value: u16) -> u16 {
    let v = value as u32;
    // Swap adjacent 2-bit pairs, then nibbles, then bytes.
    let v = ((v & 0x3333) << 2) | ((v >> 2) & 0x3333);
    let v = ((v & 0x0f0f) << 4) | ((v >> 4) & 0x0f0f);
    (((v & 0x00ff) << 8) | ((v >> 8) & 0x00ff)) as u16
}

#[inline]
fn pack_valid_label(seq: &[u8], label_words: usize) -> Result<[u64; 4], BucketError> {
    let mut words = [0u64; 4];
    if label_words > words.len() || seq.len() > label_words * 32 {
        return Err(BucketError::MalformedRecord);
    }

    let full_words = seq.len() / 32;
    for word_idx in 0..full_words {
        words[word_idx] = pack_valid_word_32(&seq[word_idx * 32..]);
    }
    let tail = &seq[full_words * 32..];
    if !tail.is_empty() {
        let mut word = 0u64;
        for &base in tail {
            word = (word << 2) | u64::from(valid_ascii_base_bits(base));
        }
        words[full_words] = word << (2 * (32 - tail.len()));
    }
    Ok(words)
}

struct BucketFile {
    path: PathBuf,
    file: File,
    records: u64,
    bytes_written: u64,
    record_size: u64,
    compressed: bool,
    interleaved_compression: bool,
    label_words: usize,
}

/// Reusable staging for compressed bucket writes.
///
/// A `BucketFile` is opened per flush, so these buffers belong to the emitter
/// that owns the flush loop; held there, a 64 KiB block costs no allocation.
#[derive(Default)]
pub(crate) struct CompressionScratch {
    /// Fixed-size attributes, de-interleaved out of the record stream.
    attrs: Vec<u8>,
    /// Label words, de-interleaved out of the record stream.
    labels: Vec<u8>,
    /// LZ4 output for the attribute stream, or for the whole block.
    encoded: Vec<u8>,
    /// LZ4 output for the label stream.
    encoded_labels: Vec<u8>,
    /// Header and payload assembled for a single `write_all`.
    block: Vec<u8>,
}

impl CompressionScratch {
    /// Compresses `input` into `out`, returning the encoded length.
    ///
    /// `out` only ever grows, so the zero-fill a resize implies is paid once
    /// per writer rather than once per block.
    fn encode(input: &[u8], out: &mut Vec<u8>) -> Result<usize, BucketError> {
        let bound = lz4_flex::block::get_maximum_output_size(input.len());
        if out.len() < bound {
            out.resize(bound, 0);
        }
        lz4_flex::block::compress_into(input, &mut out[..bound])
            .map_err(|_| BucketError::MalformedRecord)
    }
}

/// Whether uncolored buckets compress attributes and labels as separate
/// streams, the way the colored path and C++ both do.
///
/// The interleaved default compresses whole records in one stream. The header
/// records which was used, so readers stay correct under either setting.
fn force_split_compression() -> bool {
    static SPLIT: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *SPLIT.get_or_init(|| std::env::var_os("CF3_RS_SPLIT_COMPRESSION").is_some())
}

impl BucketFile {
    #[allow(clippy::too_many_arguments)]
    fn create(
        bucket_dir: &Path,
        k: u16,
        minimizer_len: u16,
        graph_count: usize,
        graph_id: usize,
        colored: bool,
        label_words: usize,
        compressed: bool,
    ) -> Result<Self, BucketError> {
        let path = bucket_dir.join(format!("{graph_id:05}.wsk"));
        let mut file = File::create(&path).map_err(|source| BucketError::Io {
            path: path.clone(),
            source,
        })?;

        file.write_all(MAGIC).map_err(|source| BucketError::Io {
            path: path.clone(),
            source,
        })?;
        write_u16(&mut file, &path, k)?;
        write_u16(&mut file, &path, minimizer_len)?;
        write_u64(&mut file, &path, graph_count as u64)?;
        write_u64(&mut file, &path, graph_id as u64)?;
        let interleaved_compression = compressed && !colored && !force_split_compression();
        file.write_all(&[
            u8::from(colored),
            label_words as u8,
            if interleaved_compression {
                2
            } else {
                u8::from(compressed)
            },
            0,
            0,
            0,
        ])
        .map_err(|source| BucketError::Io {
            path: path.clone(),
            source,
        })?;
        write_u64(&mut file, &path, 0)?;

        Ok(Self {
            path,
            file,
            records: 0,
            bytes_written: HEADER_LEN,
            record_size: record_size(colored, label_words) as u64,
            compressed,
            interleaved_compression,
            label_words,
        })
    }

    fn open_existing(path: &Path, records: u64, bytes_written: u64) -> Result<Self, BucketError> {
        let mut file = OpenOptions::new()
            .read(true)
            .write(true)
            .open(path)
            .map_err(|source| BucketError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        let header = read_header(&mut file, path)?;
        file.seek(SeekFrom::End(0))
            .map_err(|source| BucketError::Io {
                path: path.to_path_buf(),
                source,
            })?;

        Ok(Self {
            path: path.to_path_buf(),
            file,
            records,
            bytes_written,
            record_size: record_size(header.colored, header.label_words) as u64,
            compressed: header.compressed,
            interleaved_compression: header.interleaved_compression,
            label_words: header.label_words,
        })
    }

    fn write_records(
        &mut self,
        bytes: &[u8],
        records: u64,
        scratch: &mut CompressionScratch,
    ) -> Result<(u64, u64), BucketError> {
        debug_assert_eq!(bytes.len() as u64, records * self.record_size);
        let written = if self.compressed {
            self.write_compressed_block(bytes, records, scratch)?
        } else {
            self.file
                .write_all(bytes)
                .map_err(|source| BucketError::Io {
                    path: self.path.clone(),
                    source,
                })?;
            bytes.len() as u64
        };
        self.records = self
            .records
            .checked_add(records)
            .ok_or(BucketError::TooManyRecords)?;
        self.bytes_written += written;
        Ok((self.records, self.bytes_written))
    }

    fn write_compressed_block(
        &mut self,
        bytes: &[u8],
        records: u64,
        scratch: &mut CompressionScratch,
    ) -> Result<u64, BucketError> {
        let len = encode_compressed_block(
            bytes,
            records,
            self.record_size,
            self.label_words,
            self.interleaved_compression,
            scratch,
        )?;
        if len == 0 {
            return Ok(0);
        }
        self.file
            .write_all(&scratch.block)
            .map_err(|source| BucketError::Io {
                path: self.path.clone(),
                source,
            })?;
        Ok(len as u64)
    }

    fn flush(&mut self) -> Result<(), BucketError> {
        self.file.flush().map_err(|source| BucketError::Io {
            path: self.path.clone(),
            source,
        })
    }

    fn finish(&mut self) -> Result<(), BucketError> {
        write_record_count(&mut self.file, &self.path, self.records)?;
        self.file
            .seek(SeekFrom::End(0))
            .map_err(|source| BucketError::Io {
                path: self.path.clone(),
                source,
            })?;
        self.file.flush().map_err(|source| BucketError::Io {
            path: self.path.clone(),
            source,
        })
    }

    fn finish_closed(path: &Path, records: u64) -> Result<(), BucketError> {
        let mut file = OpenOptions::new()
            .read(true)
            .write(true)
            .open(path)
            .map_err(|source| BucketError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        write_record_count(&mut file, path, records)?;
        file.flush().map_err(|source| BucketError::Io {
            path: path.to_path_buf(),
            source,
        })
    }
}

fn write_record_count(file: &mut File, path: &Path, records: u64) -> Result<(), BucketError> {
    file.seek(SeekFrom::Start(RECORD_COUNT_OFFSET))
        .map_err(|source| BucketError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    write_u64(file, path, records)
}

fn label_word_count(k: u16, minimizer_len: u16) -> usize {
    let max_weak_superkmer_len = 2 * (usize::from(k) - 1) - usize::from(minimizer_len) + 2;
    max_weak_superkmer_len.div_ceil(32)
}

fn record_size(colored: bool, label_words: usize) -> usize {
    (if colored { 4 } else { 2 }) + label_words * 8
}

fn pack_uncolored_attr(len: usize, left_discontinuous: bool, right_discontinuous: bool) -> u32 {
    (len as u32) | ((left_discontinuous as u32) << 8) | ((right_discontinuous as u32) << 9)
}

fn pack_colored_attr(
    len: usize,
    source_id: u32,
    left_discontinuous: bool,
    right_discontinuous: bool,
) -> u32 {
    pack_uncolored_attr(len, left_discontinuous, right_discontinuous) | (source_id << 10)
}

fn decode_label_into(words: &[u64], len: usize, seq: &mut Vec<u8>) -> Result<(), BucketError> {
    if len > words.len() * 32 {
        return Err(BucketError::MalformedRecord);
    }

    seq.clear();
    seq.reserve(len);
    for idx in 0..len {
        let word_idx = idx / 32;
        let shift = 2 * (31 - (idx % 32));
        let base = match ((words[word_idx] >> shift) & 0b11) as u8 {
            0 => Base::A,
            1 => Base::C,
            2 => Base::G,
            3 => Base::T,
            _ => unreachable!(),
        };
        seq.push(base.to_ascii());
    }
    Ok(())
}

fn read_header(file: &mut impl Read, path: &Path) -> Result<BucketHeader, BucketError> {
    let mut magic = [0u8; 8];
    file.read_exact(&mut magic)
        .map_err(|source| BucketError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    if &magic != MAGIC {
        return Err(BucketError::BadMagic(path.to_path_buf()));
    }

    let k = read_u16(file, path)?;
    let minimizer_len = read_u16(file, path)?;
    let graph_count = usize::try_from(read_u64(file, path)?)
        .map_err(|_| BucketError::MalformedHeader(path.to_path_buf()))?;
    let graph_id = usize::try_from(read_u64(file, path)?)
        .map_err(|_| BucketError::MalformedHeader(path.to_path_buf()))?;

    let mut flags = [0u8; 6];
    file.read_exact(&mut flags)
        .map_err(|source| BucketError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    let colored = match flags[0] {
        0 => false,
        1 => true,
        _ => return Err(BucketError::MalformedHeader(path.to_path_buf())),
    };
    let label_words = flags[1] as usize;
    if label_words == 0 || label_words != label_word_count(k, minimizer_len) {
        return Err(BucketError::MalformedHeader(path.to_path_buf()));
    }
    let (compressed, interleaved_compression) = match flags[2] {
        0 => (false, false),
        1 => (true, false),
        2 if !colored => (true, true),
        _ => return Err(BucketError::MalformedHeader(path.to_path_buf())),
    };
    if flags[3..].iter().any(|&b| b != 0) {
        return Err(BucketError::MalformedHeader(path.to_path_buf()));
    }
    let records = read_u64(file, path)?;

    if graph_id >= graph_count {
        return Err(BucketError::MalformedHeader(path.to_path_buf()));
    }

    Ok(BucketHeader {
        k,
        minimizer_len,
        graph_count,
        graph_id,
        colored,
        compressed,
        interleaved_compression,
        label_words,
        records,
    })
}

fn read_u16(file: &mut impl Read, path: &Path) -> Result<u16, BucketError> {
    let mut bytes = [0u8; 2];
    file.read_exact(&mut bytes)
        .map_err(|source| BucketError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(u16::from_le_bytes(bytes))
}

fn read_u64(file: &mut impl Read, path: &Path) -> Result<u64, BucketError> {
    let mut bytes = [0u8; 8];
    file.read_exact(&mut bytes)
        .map_err(|source| BucketError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(u64::from_le_bytes(bytes))
}

fn write_u16(file: &mut File, path: &Path, value: u16) -> Result<(), BucketError> {
    file.write_all(&value.to_le_bytes())
        .map_err(|source| BucketError::Io {
            path: path.to_path_buf(),
            source,
        })
}

fn write_u64(file: &mut File, path: &Path, value: u64) -> Result<(), BucketError> {
    file.write_all(&value.to_le_bytes())
        .map_err(|source| BucketError::Io {
            path: path.to_path_buf(),
            source,
        })
}

#[derive(Debug)]
pub enum BucketError {
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    GraphCountTooLarge(usize),
    InvalidGraphId(usize),
    LabelTooLong(usize),
    MissingSourceId,
    SourceIdTooLarge(u32),
    InvalidBase(u8),
    TooManyRecords,
    BadMagic(PathBuf),
    MalformedHeader(PathBuf),
    MalformedManifest(PathBuf),
    MalformedRecord,
    RecordGraphMismatch {
        expected: usize,
        got: usize,
    },
    WorkerPanic,
}

impl std::fmt::Display for BucketError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io { path, source } => write!(f, "{}: {source}", path.display()),
            Self::GraphCountTooLarge(count) => write!(f, "graph count is too large: {count}"),
            Self::InvalidGraphId(graph_id) => write!(f, "invalid graph id: {graph_id}"),
            Self::LabelTooLong(len) => write!(f, "weak super-kmer label is too long: {len}"),
            Self::MissingSourceId => write!(f, "colored bucket record is missing source id"),
            Self::SourceIdTooLarge(source_id) => {
                write!(
                    f,
                    "source id exceeds 21-bit colored bucket limit: {source_id}"
                )
            }
            Self::InvalidBase(b) => write!(f, "invalid base in bucket label: '{}'", *b as char),
            Self::TooManyRecords => write!(f, "too many weak super-kmer records"),
            Self::BadMagic(path) => write!(
                f,
                "not a CF3 Rust weak-superkmer bucket: {}",
                path.display()
            ),
            Self::MalformedHeader(path) => {
                write!(
                    f,
                    "malformed weak-superkmer bucket header: {}",
                    path.display()
                )
            }
            Self::MalformedManifest(path) => {
                write!(
                    f,
                    "malformed weak-superkmer bucket manifest: {}",
                    path.display()
                )
            }
            Self::MalformedRecord => write!(f, "malformed weak-superkmer bucket record"),
            Self::RecordGraphMismatch { expected, got } => {
                write!(
                    f,
                    "bucket record graph id mismatch: expected {expected}, got {got}"
                )
            }
            Self::WorkerPanic => write!(f, "bucket worker thread panicked"),
        }
    }
}

impl std::error::Error for BucketError {}

#[cfg(test)]
mod tests {
    /// Reclaim must actually return blocks, not merely be called.
    ///
    /// A failed punch is ignored by design -- the container is unlinked
    /// wholesale at the end regardless -- which means a filesystem or platform
    /// where it silently does nothing costs peak disk with no other symptom.
    /// That is exactly what happened once already: deferring reclaim left the
    /// work directory 24.5 GB larger and presented as an unexplained number.
    /// This asserts the blocks come back, so the macOS `fcntl(F_PUNCHHOLE)`
    /// path is verified by CI rather than assumed from the fact that it
    /// compiles.
    #[test]
    fn releasing_segments_returns_blocks_to_the_filesystem() {
        use std::os::unix::fs::MetadataExt;

        let dir = std::env::temp_dir().join(format!(
            "cf3-punch-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();

        let containers = BucketContainers::create(&dir, 1).unwrap();
        let segment = containers.segment_bytes();
        let segments: Vec<u32> = (0..8).collect();
        let payload = vec![0xA5u8; segment as usize];
        for &index in &segments {
            containers
                .write_at(0, u64::from(index) * segment, &payload)
                .unwrap();
        }

        let path = dir.join("00000.wskc");
        let allocated = || fs::metadata(&path).unwrap().blocks() * 512;
        let before = allocated();
        assert!(
            before >= segments.len() as u64 * segment,
            "expected {} bytes of blocks before punching, saw {before}",
            segments.len() as u64 * segment
        );

        containers.release_segments(0, &segments);
        let after = allocated();

        assert!(
            after < before / 2,
            "punching {} segments freed {} of {before} bytes; this filesystem may \
             not support hole punching, in which case consumed bucket space is \
             held until the build ends and peak disk is higher than documented",
            segments.len(),
            before - after
        );
        // The file keeps its length; only the blocks behind it go away.
        assert_eq!(
            fs::metadata(&path).unwrap().len(),
            segments.len() as u64 * segment
        );

        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn packed_word_matches_scalar_reference() {
        fn scalar(seq: &[u8]) -> u64 {
            let mut word = 0u64;
            for &base in &seq[..32] {
                word = (word << 2) | u64::from(valid_ascii_base_bits(base));
            }
            word
        }
        let alphabet = b"ACGT";
        // Deterministic pseudo-random coverage plus structured edge cases.
        let mut state = 0x1234_5678_9abc_def0u64;
        for case in 0..2048 {
            let mut seq = [0u8; 32];
            for (i, slot) in seq.iter_mut().enumerate() {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                *slot = if case < 4 {
                    alphabet[(case + i) % 4]
                } else {
                    alphabet[(state >> 33) as usize % 4]
                };
            }
            assert_eq!(
                pack_valid_word_32(&seq),
                scalar(&seq),
                "mismatch for {:?}",
                std::str::from_utf8(&seq).unwrap()
            );
        }
    }

    use super::*;

    #[test]
    fn deferred_uncolored_atlas_chunks_preserve_graph_buckets() {
        let dir = std::env::temp_dir().join(format!(
            "cf3-uncolored-atlas-bucket-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        fs::create_dir_all(&dir).unwrap();
        let mut params = BuildParams::new(crate::GraphInput::References, "test".to_string());
        params.k = 31;
        params.minimizer_len = 15;
        params.threads = 2;
        params.work_dir = dir.to_string_lossy().into_owned();
        let sink = SharedBucketSink::create(&params, ATLAS_GRAPH_COUNT + 1).unwrap();
        let expected = [(0, b'A'), (7, b'C'), (ATLAS_GRAPH_COUNT, b'G')];
        let mut emitters = Vec::new();
        for repeat in 0..2 {
            let mut emitter = sink.deferred_uncolored_emitter();
            for &(graph_id, base) in &expected {
                for _ in 0..(repeat + 1) {
                    emitter
                        .add_valid(
                            &WeakSuperKmer {
                                graph_id,
                                offset: 0,
                                len: 31,
                                source_id: None,
                                left_discontinuous: false,
                                right_discontinuous: false,
                            },
                            &[base; 31],
                        )
                        .unwrap();
                }
            }
            emitters.push(emitter);
        }
        sink.flush_uncolored_emitters(emitters).unwrap();
        let stats = sink.finish().unwrap();

        let (store, entries) = BucketStore::open_dir(&stats.bucket_dir).unwrap();
        for &(graph_id, base) in &expected {
            let entry = entries
                .iter()
                .find(|entry| entry.graph_id == graph_id)
                .expect("bucket in manifest");
            let mut reader = store.reader(entry).unwrap();
            let mut record = BucketRecord::default();
            let mut count = 0;
            while reader.next_record_into(&mut record).unwrap() {
                assert_eq!(record.graph_id, graph_id);
                assert_eq!(record.label, vec![base; 31]);
                count += 1;
            }
            assert_eq!(count, 3);
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn colored_atlas_windows_preserve_global_source_order() {
        let dir = std::env::temp_dir().join(format!(
            "cf3-colored-bucket-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        fs::create_dir_all(&dir).unwrap();
        let mut params = BuildParams::new(crate::GraphInput::References, "test".to_string());
        params.color = true;
        params.k = 31;
        params.minimizer_len = 15;
        params.threads = 1;
        params.work_dir = dir.to_string_lossy().into_owned();
        let sink = SharedBucketSink::create(&params, 1).unwrap();

        for (source_min, source_max, sources) in
            [(1, 3, vec![3, 1, 3, 2, 1]), (4, 6, vec![6, 4, 5, 4])]
        {
            let mut emitter = sink.emitter();
            for source_id in sources {
                emitter
                    .add_valid(
                        &WeakSuperKmer {
                            graph_id: 0,
                            offset: 0,
                            len: 31,
                            source_id: Some(source_id),
                            left_discontinuous: false,
                            right_discontinuous: false,
                        },
                        &[b'A'; 31],
                    )
                    .unwrap();
            }
            emitter.finish().unwrap();
            sink.flush_colored_window(source_min, source_max).unwrap();
        }
        let stats = sink.finish().unwrap();

        let (store, entries) = BucketStore::open_dir(&stats.bucket_dir).unwrap();
        let entry = entries
            .iter()
            .find(|entry| entry.graph_id == 0)
            .expect("bucket in manifest");
        let mut reader = store.reader(entry).unwrap();
        assert!(reader.header().compressed);
        let mut sources = Vec::new();
        let mut record = BucketPackedRecord::default();
        while reader.next_packed_record_into(&mut record).unwrap() {
            sources.push(record.source_id.unwrap());
        }
        assert_eq!(sources, [1, 1, 2, 3, 3, 4, 4, 5, 6]);

        // Clipping the container leaves the manifest claiming a length the
        // chain can no longer supply, which is the container-shaped version of
        // a bucket file cut short by a killed run.
        let container_path = stats.bucket_dir.join("00000.wskc");
        let len = fs::metadata(&container_path).unwrap().len();
        OpenOptions::new()
            .write(true)
            .open(&container_path)
            .unwrap()
            .set_len(len - 1)
            .unwrap();
        let (store, entries) = BucketStore::open_dir(&stats.bucket_dir).unwrap();
        let entry = entries
            .iter()
            .find(|entry| entry.graph_id == 0)
            .expect("bucket in manifest");
        let mut truncated = store.reader(entry).unwrap();
        let mut record = BucketPackedRecord::default();
        let mut failed = false;
        loop {
            match truncated.next_packed_record_into(&mut record) {
                Ok(true) => {}
                Ok(false) => break,
                Err(_) => {
                    failed = true;
                    break;
                }
            }
        }
        assert!(failed, "truncated compressed block must be rejected");
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn colored_worker_tails_preserve_cpp_worker_order() {
        let dir = std::env::temp_dir().join(format!(
            "cf3-colored-worker-tails-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        fs::create_dir_all(&dir).unwrap();
        let mut params = BuildParams::new(crate::GraphInput::References, "test".to_string());
        params.color = true;
        params.k = 31;
        params.minimizer_len = 15;
        params.threads = 2;
        params.work_dir = dir.to_string_lossy().into_owned();
        let sink = SharedBucketSink::create(&params, 1).unwrap();

        let mut emitters = Vec::new();
        for sources in [[3, 1], [4, 2]] {
            let mut emitter = sink.emitter();
            for source_id in sources {
                emitter
                    .add_valid(
                        &WeakSuperKmer {
                            graph_id: 0,
                            offset: 0,
                            len: 31,
                            source_id: Some(source_id),
                            left_discontinuous: false,
                            right_discontinuous: false,
                        },
                        &[b'A'; 31],
                    )
                    .unwrap();
            }
            emitters.push(emitter);
        }
        sink.flush_colored_emitters(emitters).unwrap();
        let stats = sink.finish().unwrap();

        let (store, entries) = BucketStore::open_dir(&stats.bucket_dir).unwrap();
        let entry = entries
            .iter()
            .find(|entry| entry.graph_id == 0)
            .expect("bucket in manifest");
        let mut reader = store.reader(entry).unwrap();
        let mut sources = Vec::new();
        let mut record = BucketPackedRecord::default();
        while reader.next_packed_record_into(&mut record).unwrap() {
            sources.push(record.source_id.unwrap());
        }
        assert_eq!(sources, [3, 1, 4, 2]);
        fs::remove_dir_all(dir).unwrap();
    }
}

//! FASTA/FASTQ input discovery and zero-copy fragment parsing.
//!
//! Parsers split records at non-ACGT symbols. Borrowed callbacks are preferred
//! by the production partitioner so sequence data does not need to be copied.

use crate::dna::is_dna_ascii;
use crate::params::BuildParams;
use flate2::read::MultiGzDecoder;
use std::fs;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceFragment {
    pub source_id: u32,
    pub record_id: u64,
    pub offset: usize,
    pub seq: Vec<u8>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BorrowedSequenceFragment<'a> {
    pub source_id: u32,
    pub record_id: u64,
    pub offset: usize,
    pub seq: &'a [u8],
}

pub fn expand_input_paths(params: &BuildParams) -> Result<Vec<PathBuf>, InputError> {
    let mut paths = Vec::new();

    for path in &params.seqs {
        paths.push(PathBuf::from(path));
    }

    for list in &params.lists {
        let file = fs::File::open(list).map_err(|source| InputError::Io {
            path: PathBuf::from(list),
            source,
        })?;
        for line in BufReader::new(file).lines() {
            let line = line.map_err(|source| InputError::Io {
                path: PathBuf::from(list),
                source,
            })?;
            let trimmed = line.trim();
            if !trimmed.is_empty() {
                paths.push(PathBuf::from(trimmed));
            }
        }
    }

    for dir in &params.dirs {
        let mut entries = fs::read_dir(dir)
            .map_err(|source| InputError::Io {
                path: PathBuf::from(dir),
                source,
            })?
            .collect::<Result<Vec<_>, _>>()
            .map_err(|source| InputError::Io {
                path: PathBuf::from(dir),
                source,
            })?;
        entries.sort_by_key(|entry| entry.path());
        for entry in entries {
            let path = entry.path();
            if path.is_file() {
                paths.push(path);
            }
        }
    }

    if paths.is_empty() {
        return Err(InputError::NoInput);
    }

    Ok(paths)
}

pub fn parse_fragments<P, F>(
    path: P,
    source_id: u32,
    min_len: usize,
    mut on_fragment: F,
) -> Result<u64, InputError>
where
    P: AsRef<Path>,
    F: FnMut(SequenceFragment) -> Result<(), InputError>,
{
    parse_fragments_borrowed(path, source_id, min_len, |fragment| {
        on_fragment(SequenceFragment {
            source_id: fragment.source_id,
            record_id: fragment.record_id,
            offset: fragment.offset,
            seq: normalized_fragment_seq(fragment.seq),
        })
    })
}

pub fn parse_fragments_borrowed<P, F>(
    path: P,
    source_id: u32,
    min_len: usize,
    on_fragment: F,
) -> Result<u64, InputError>
where
    P: AsRef<Path>,
    F: for<'a> FnMut(BorrowedSequenceFragment<'a>) -> Result<(), InputError>,
{
    parse_fragments_borrowed_with(path, source_id, min_len, 1, on_fragment)
}

/// As [`parse_fragments_borrowed`], but permitted to use `inflate_workers`
/// threads to decompress gzip input.
///
/// The decompressed byte stream is identical either way; only the work of
/// producing it is shared out. BGZF blocks inflate on the dedicated
/// block-parallel reader; plain gzip members inflate through rapidgzip's
/// speculative parallel decoder, which finds block boundaries mid-stream.
pub fn parse_fragments_borrowed_with<P, F>(
    path: P,
    source_id: u32,
    min_len: usize,
    inflate_workers: usize,
    on_fragment: F,
) -> Result<u64, InputError>
where
    P: AsRef<Path>,
    F: for<'a> FnMut(BorrowedSequenceFragment<'a>) -> Result<(), InputError>,
{
    parse_fragments_borrowed_with_registry(
        path,
        source_id,
        min_len,
        inflate_workers,
        None,
        on_fragment,
    )
}

/// As [`parse_fragments_borrowed_with`], with an optional [`InflateRegistry`]
/// that gains control of every rapidgzip decoder opened for this file.
pub(crate) fn parse_fragments_borrowed_with_registry<P, F>(
    path: P,
    source_id: u32,
    min_len: usize,
    inflate_workers: usize,
    registry: Option<&std::sync::Arc<InflateRegistry>>,
    mut on_fragment: F,
) -> Result<u64, InputError>
where
    P: AsRef<Path>,
    F: for<'a> FnMut(BorrowedSequenceFragment<'a>) -> Result<(), InputError>,
{
    let path = path.as_ref();
    let file = fs::File::open(path).map_err(|source| InputError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let input: Box<dyn Read> = if path.extension().is_some_and(|ext| ext == "gz") {
        // rapidgzip covers the whole gzip family: BGZF's independent members
        // decode block-parallel, and plain single-member gzip decodes through
        // speculative mid-stream splits. Measured at parity with the retired
        // dedicated BGZF reader over three interleaved pairs.
        if inflate_workers > 1 && !sequential_inflate_diagnostic() {
            open_parallel_gzip_reader(path, inflate_workers, registry)?
        } else {
            Box::new(MultiGzDecoder::new(file))
        }
    } else {
        Box::new(file)
    };
    let mut reader = BufReader::with_capacity(1024 * 1024, input);
    let first_line = next_non_empty_line(&mut reader, path)?;
    match first_line.first().copied() {
        Some(b'>') => parse_fasta_reader(first_line, reader, source_id, min_len, &mut on_fragment),
        Some(b'@') => parse_fastq_reader(first_line, reader, source_id, min_len, &mut on_fragment),
        Some(_) if first_line.iter().copied().all(is_dna_ascii) => {
            parse_plain_sequence_reader(first_line, reader, source_id, min_len, &mut on_fragment)
        }
        Some(_) => Err(InputError::UnknownFormat(path.to_path_buf())),
        None => Err(InputError::EmptyFile(path.to_path_buf())),
    }
}

/// Whether to force single-threaded gzip inflation (measurement diagnostic).
fn sequential_inflate_diagnostic() -> bool {
    static SEQUENTIAL: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *SEQUENTIAL.get_or_init(|| std::env::var_os("CF3_RS_SEQUENTIAL_INFLATE").is_some())
}

/// Live rapidgzip decoders across every open input, as one controllable pool.
///
/// This is the seam a closed-loop thread broker drives: it aggregates busy
/// time and consumed bytes across decoders (including ones that have already
/// closed), and fans an aggregate worker limit out as an even per-decoder
/// share. Registration happens inside [`parse_fragments_borrowed_with_registry`]
/// when a rapidgzip reader is opened; a guard carried by the reader retires
/// the handle and banks its final counters on drop.
pub struct InflateRegistry {
    handles: std::sync::Mutex<Vec<Option<rapidgzip_core::DecoderHandle>>>,
    aggregate_limit: std::sync::atomic::AtomicUsize,
    finished_busy_nanos: std::sync::atomic::AtomicU64,
    finished_bytes: std::sync::atomic::AtomicU64,
}

impl InflateRegistry {
    pub fn new(aggregate_limit: usize) -> Self {
        Self {
            handles: std::sync::Mutex::new(Vec::new()),
            aggregate_limit: std::sync::atomic::AtomicUsize::new(aggregate_limit.max(1)),
            finished_busy_nanos: std::sync::atomic::AtomicU64::new(0),
            finished_bytes: std::sync::atomic::AtomicU64::new(0),
        }
    }

    /// Sets the aggregate decoder-worker ceiling and redistributes it.
    pub fn set_aggregate_limit(&self, limit: usize) {
        self.aggregate_limit
            .store(limit.max(1), std::sync::atomic::Ordering::Relaxed);
        self.redistribute(&self.handles.lock().unwrap_or_else(|p| p.into_inner()));
    }

    pub fn aggregate_limit(&self) -> usize {
        self.aggregate_limit
            .load(std::sync::atomic::Ordering::Relaxed)
    }

    /// Sums a stat over the live decoders.
    pub fn sum_stats<T: Default + std::ops::Add<Output = T>>(
        &self,
        mut field: impl FnMut(&rapidgzip_core::DecoderStats) -> T,
    ) -> T {
        let handles = self.handles.lock().unwrap_or_else(|p| p.into_inner());
        handles
            .iter()
            .flatten()
            .map(|handle| field(&handle.stats()))
            .fold(T::default(), |acc, v| acc + v)
    }

    /// Cumulative decoder busy nanoseconds and consumed bytes, closed
    /// decoders included.
    pub fn cumulative_work(&self) -> (u64, u64) {
        let handles = self.handles.lock().unwrap_or_else(|p| p.into_inner());
        let (mut nanos, mut bytes) = (
            self.finished_busy_nanos
                .load(std::sync::atomic::Ordering::Relaxed),
            self.finished_bytes
                .load(std::sync::atomic::Ordering::Relaxed),
        );
        for handle in handles.iter().flatten() {
            let stats = handle.stats();
            nanos += stats
                .accounted_busy_time
                .map_or(0, |t| u64::try_from(t.as_nanos()).unwrap_or(u64::MAX));
            bytes += stats.consumed_bytes;
        }
        (nanos, bytes)
    }

    /// Snapshot of every live decoder's pressure classification.
    pub fn pressures(&self) -> Vec<rapidgzip_core::DecoderPressure> {
        let handles = self.handles.lock().unwrap_or_else(|p| p.into_inner());
        handles
            .iter()
            .flatten()
            .map(|handle| handle.stats().pressure)
            .collect()
    }

    fn register(&self, handle: rapidgzip_core::DecoderHandle) -> usize {
        let mut handles = self.handles.lock().unwrap_or_else(|p| p.into_inner());
        let slot = handles.iter().position(Option::is_none).unwrap_or_else(|| {
            handles.push(None);
            handles.len() - 1
        });
        handles[slot] = Some(handle);
        self.redistribute(&handles);
        slot
    }

    fn unregister(&self, slot: usize) {
        let mut handles = self.handles.lock().unwrap_or_else(|p| p.into_inner());
        if let Some(handle) = handles[slot].take() {
            let stats = handle.stats();
            self.finished_busy_nanos.fetch_add(
                stats
                    .accounted_busy_time
                    .map_or(0, |t| u64::try_from(t.as_nanos()).unwrap_or(u64::MAX)),
                std::sync::atomic::Ordering::Relaxed,
            );
            self.finished_bytes
                .fetch_add(stats.consumed_bytes, std::sync::atomic::Ordering::Relaxed);
        }
        self.redistribute(&handles);
    }

    /// Fans the aggregate limit out as an even share per live decoder.
    fn redistribute(&self, handles: &[Option<rapidgzip_core::DecoderHandle>]) {
        let live = handles.iter().flatten().count();
        if live == 0 {
            return;
        }
        let share = (self.aggregate_limit() / live).max(1);
        for handle in handles.iter().flatten() {
            // A failed resize leaves the previous per-decoder limit standing,
            // which is safe; the broker re-issues limits every window.
            let _ = handle.set_worker_limit(share);
        }
    }
}

/// A rapidgzip reader that stays registered for as long as it is open.
struct RegisteredReader {
    inner: rapidgzip_core::DecoderReader,
    registry: std::sync::Arc<InflateRegistry>,
    slot: usize,
}

impl Read for RegisteredReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.inner.read(buf)
    }
}

impl Drop for RegisteredReader {
    fn drop(&mut self) {
        self.registry.unregister(self.slot);
    }
}

/// Opens a plain (non-BGZF) gzip file through rapidgzip's parallel decoder.
///
/// rapidgzip decodes ahead speculatively from mid-stream bit positions and
/// reconciles at member boundaries, so a single large gzip member -- which the
/// serial decoder is bound to walk alone at roughly 150 MB/s of compressed
/// input -- scales with the worker budget. The returned reader yields exactly
/// the bytes `MultiGzDecoder` would.
fn open_parallel_gzip_reader(
    path: &Path,
    inflate_workers: usize,
    registry: Option<&std::sync::Arc<InflateRegistry>>,
) -> Result<Box<dyn Read>, InputError> {
    let build_error = |message: String| InputError::Io {
        path: path.to_path_buf(),
        source: std::io::Error::other(message),
    };
    // Under a registry the per-decoder ceiling is owned by the broker, so the
    // builder gets the whole aggregate as immutable headroom and the live
    // limit is applied through the handle at registration.
    let headroom = registry.map_or(inflate_workers, |r| r.aggregate_limit());
    let reader = rapidgzip_core::Decoder::builder()
        .decoder_threads(headroom.max(inflate_workers))
        .build()
        .map_err(|error| build_error(format!("rapidgzip configuration: {error}")))?
        .open(path)
        .map_err(|error| build_error(format!("rapidgzip open: {error}")))?;
    Ok(match registry {
        None => Box::new(reader),
        Some(registry) => {
            let registry = std::sync::Arc::clone(registry);
            let slot = registry.register(reader.handle());
            Box::new(RegisteredReader {
                inner: reader,
                registry,
                slot,
            })
        }
    })
}

/// Appends a payload line, dropping embedded whitespace.
///
/// The filtering form is a per-byte predicate and push, which cannot vectorize.
/// Sequence lines essentially never contain interior whitespace, so probe first
/// with a scan the compiler can vectorize and fall back only when needed. The
/// slow path is kept so that interior whitespace still joins a record rather
/// than splitting it into fragments.
#[inline]
fn append_sequence_line(seq: &mut Vec<u8>, line: &[u8]) {
    if line.iter().any(u8::is_ascii_whitespace) {
        seq.extend(line.iter().copied().filter(|b| !b.is_ascii_whitespace()));
    } else {
        seq.extend_from_slice(line);
    }
}

fn parse_plain_sequence_reader<R, F>(
    first_line: Vec<u8>,
    mut reader: R,
    source_id: u32,
    min_len: usize,
    on_fragment: &mut F,
) -> Result<u64, InputError>
where
    R: BufRead,
    F: for<'a> FnMut(BorrowedSequenceFragment<'a>) -> Result<(), InputError>,
{
    let mut seq = first_line;
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader
            .read_until(b'\n', &mut line)
            .map_err(|source| InputError::Read { source })?
            == 0
        {
            break;
        }
        trim_ascii_line_in_place(&mut line);
        append_sequence_line(&mut seq, &line);
    }
    emit_actg_fragments(source_id, 1, &seq, min_len, on_fragment)?;
    Ok(1)
}

fn next_non_empty_line<R: BufRead>(reader: &mut R, path: &Path) -> Result<Vec<u8>, InputError> {
    let mut line = Vec::new();
    loop {
        line.clear();
        let n = reader
            .read_until(b'\n', &mut line)
            .map_err(|source| InputError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        if n == 0 {
            return Err(InputError::EmptyFile(path.to_path_buf()));
        }
        trim_ascii_line_in_place(&mut line);
        if !line.is_empty() {
            return Ok(line);
        }
    }
}

fn parse_fasta_reader<R, F>(
    first_header: Vec<u8>,
    mut reader: R,
    source_id: u32,
    min_len: usize,
    on_fragment: &mut F,
) -> Result<u64, InputError>
where
    R: BufRead,
    F: for<'a> FnMut(BorrowedSequenceFragment<'a>) -> Result<(), InputError>,
{
    debug_assert!(first_header.starts_with(b">"));
    let mut records = 0u64;
    let mut record_id = 1u64;
    let mut seq = Vec::new();
    let mut line = Vec::new();

    loop {
        line.clear();
        let n = reader
            .read_until(b'\n', &mut line)
            .map_err(|source| InputError::Read { source })?;
        if n == 0 {
            break;
        }
        trim_ascii_line_in_place(&mut line);
        if line.starts_with(b">") {
            records += 1;
            emit_actg_fragments(source_id, record_id, &seq, min_len, on_fragment)?;
            seq.clear();
            record_id += 1;
        } else if !line.is_empty() {
            append_sequence_line(&mut seq, &line);
        }
    }

    records += 1;
    emit_actg_fragments(source_id, record_id, &seq, min_len, on_fragment)?;

    Ok(records)
}

fn parse_fastq_reader<R, F>(
    first_header: Vec<u8>,
    mut reader: R,
    source_id: u32,
    min_len: usize,
    on_fragment: &mut F,
) -> Result<u64, InputError>
where
    R: BufRead,
    F: for<'a> FnMut(BorrowedSequenceFragment<'a>) -> Result<(), InputError>,
{
    debug_assert!(first_header.starts_with(b"@"));
    let mut record_id = 0u64;
    let mut header = first_header;
    let mut seq = Vec::new();
    let mut line = Vec::new();

    loop {
        record_id += 1;
        if !header.starts_with(b"@") {
            return Err(InputError::MalformedFastq(record_id));
        }

        seq.clear();
        loop {
            line.clear();
            if reader
                .read_until(b'\n', &mut line)
                .map_err(|source| InputError::Read { source })?
                == 0
            {
                return Err(InputError::MalformedFastq(record_id));
            }
            trim_ascii_line_in_place(&mut line);
            if line.starts_with(b"+") {
                break;
            }
            append_sequence_line(&mut seq, &line);
        }

        let mut qual_len = 0usize;
        while qual_len < seq.len() {
            line.clear();
            if reader
                .read_until(b'\n', &mut line)
                .map_err(|source| InputError::Read { source })?
                == 0
            {
                return Err(InputError::MalformedFastq(record_id));
            }
            trim_ascii_line_in_place(&mut line);
            qual_len += line.len();
        }

        emit_actg_fragments(source_id, record_id, &seq, min_len, on_fragment)?;

        header.clear();
        let bytes = reader
            .read_until(b'\n', &mut header)
            .map_err(|source| InputError::Read { source })?;
        if bytes == 0 {
            break;
        }
        trim_ascii_line_in_place(&mut header);
        if header.is_empty() {
            return Err(InputError::MalformedFastq(record_id + 1));
        }
    }

    Ok(record_id)
}

#[cfg(test)]
fn parse_fasta_bytes<F>(
    bytes: &[u8],
    source_id: u32,
    min_len: usize,
    on_fragment: &mut F,
) -> Result<u64, InputError>
where
    F: FnMut(SequenceFragment) -> Result<(), InputError>,
{
    let mut reader = BufReader::new(bytes);
    let first = next_non_empty_line(&mut reader, Path::new("<memory>"))?;
    let mut on_borrowed = |fragment: BorrowedSequenceFragment<'_>| {
        on_fragment(SequenceFragment {
            source_id: fragment.source_id,
            record_id: fragment.record_id,
            offset: fragment.offset,
            seq: normalized_fragment_seq(fragment.seq),
        })
    };
    match first.first().copied() {
        Some(b'>') => parse_fasta_reader(first, reader, source_id, min_len, &mut on_borrowed),
        Some(_) => Err(InputError::UnknownFormat(PathBuf::from("<memory>"))),
        None => Err(InputError::EmptyFile(PathBuf::from("<memory>"))),
    }
}

#[cfg(test)]
fn parse_fastq_bytes<F>(
    bytes: &[u8],
    source_id: u32,
    min_len: usize,
    on_fragment: &mut F,
) -> Result<u64, InputError>
where
    F: FnMut(SequenceFragment) -> Result<(), InputError>,
{
    let mut reader = BufReader::new(bytes);
    let first = next_non_empty_line(&mut reader, Path::new("<memory>"))?;
    let mut on_borrowed = |fragment: BorrowedSequenceFragment<'_>| {
        on_fragment(SequenceFragment {
            source_id: fragment.source_id,
            record_id: fragment.record_id,
            offset: fragment.offset,
            seq: normalized_fragment_seq(fragment.seq),
        })
    };
    match first.first().copied() {
        Some(b'@') => parse_fastq_reader(first, reader, source_id, min_len, &mut on_borrowed),
        Some(_) => Err(InputError::UnknownFormat(PathBuf::from("<memory>"))),
        None => Err(InputError::EmptyFile(PathBuf::from("<memory>"))),
    }
}

fn emit_actg_fragments<F>(
    source_id: u32,
    record_id: u64,
    seq: &[u8],
    min_len: usize,
    on_fragment: &mut F,
) -> Result<(), InputError>
where
    F: for<'a> FnMut(BorrowedSequenceFragment<'a>) -> Result<(), InputError>,
{
    let mut start = None;
    for (idx, &base) in seq.iter().enumerate() {
        if is_dna_ascii(base) {
            start.get_or_insert(idx);
        } else if let Some(beg) = start.take() {
            emit_fragment(
                source_id,
                record_id,
                beg,
                &seq[beg..idx],
                min_len,
                on_fragment,
            )?;
        }
    }

    if let Some(beg) = start {
        emit_fragment(source_id, record_id, beg, &seq[beg..], min_len, on_fragment)?;
    }

    Ok(())
}

fn emit_fragment<F>(
    source_id: u32,
    record_id: u64,
    offset: usize,
    seq: &[u8],
    min_len: usize,
    on_fragment: &mut F,
) -> Result<(), InputError>
where
    F: for<'a> FnMut(BorrowedSequenceFragment<'a>) -> Result<(), InputError>,
{
    if seq.len() >= min_len {
        on_fragment(BorrowedSequenceFragment {
            source_id,
            record_id,
            offset,
            seq,
        })?;
    }
    Ok(())
}

pub fn normalized_fragment_seq(seq: &[u8]) -> Vec<u8> {
    if seq.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
        seq.to_vec()
    } else {
        seq.iter().map(|b| b.to_ascii_uppercase()).collect()
    }
}

#[inline]
fn trim_ascii_line_in_place(line: &mut Vec<u8>) {
    while line.last().is_some_and(|b| b.is_ascii_whitespace()) {
        line.pop();
    }
}

#[derive(Debug)]
pub enum InputError {
    NoInput,
    EmptyFile(PathBuf),
    UnknownFormat(PathBuf),
    MalformedFastq(u64),
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    Read {
        source: std::io::Error,
    },
    Partition(crate::partition::PartitionError),
    Bucket(crate::buckets::BucketError),
}

impl std::fmt::Display for InputError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NoInput => write!(f, "no input files resolved"),
            Self::EmptyFile(path) => write!(f, "input file is empty: {}", path.display()),
            Self::UnknownFormat(path) => {
                write!(f, "unknown FASTA/FASTQ format: {}", path.display())
            }
            Self::MalformedFastq(record) => write!(f, "malformed FASTQ record {record}"),
            Self::Io { path, source } => write!(f, "{}: {source}", path.display()),
            Self::Read { source } => write!(f, "{source}"),
            Self::Partition(err) => write!(f, "{err}"),
            Self::Bucket(err) => write!(f, "{err}"),
        }
    }
}

impl std::error::Error for InputError {}

impl From<crate::buckets::BucketError> for InputError {
    fn from(value: crate::buckets::BucketError) -> Self {
        Self::Bucket(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    #[test]
    fn fasta_splits_on_non_actg() {
        let mut fragments = Vec::new();
        let records = parse_fasta_bytes(b">r1\nAACNNttg\n>r2\nCC\n", 1, 2, &mut |frag| {
            fragments.push(frag);
            Ok(())
        })
        .unwrap();

        assert_eq!(records, 2);
        assert_eq!(
            fragments
                .iter()
                .map(|f| f.seq.as_slice())
                .collect::<Vec<_>>(),
            vec![b"AAC".as_slice(), b"TTG".as_slice(), b"CC".as_slice()]
        );
        assert_eq!(fragments[1].offset, 5);
    }

    #[test]
    fn fastq_splits_and_counts_records() {
        let mut fragments = Vec::new();
        let records = parse_fastq_bytes(b"@r1\nACNTA\n+\nIIIII\n", 7, 2, &mut |frag| {
            fragments.push(frag);
            Ok(())
        })
        .unwrap();

        assert_eq!(records, 1);
        assert_eq!(fragments.len(), 2);
        assert_eq!(fragments[0].source_id, 7);
        assert_eq!(fragments[0].seq, b"AC");
        assert_eq!(fragments[1].seq, b"TA");
    }

    #[test]
    fn parses_gzipped_fastq() {
        let path =
            std::env::temp_dir().join(format!("cf3rs-input-{}.fastq.gz", std::process::id()));
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(b"@r1\nACGTNNTA\n+\nIIIIIIII\n").unwrap();
        fs::write(&path, encoder.finish().unwrap()).unwrap();

        let mut fragments = Vec::new();
        let records = parse_fragments(&path, 3, 2, |fragment| {
            fragments.push(fragment);
            Ok(())
        })
        .unwrap();

        assert_eq!(records, 1);
        assert_eq!(
            fragments
                .iter()
                .map(|fragment| fragment.seq.as_slice())
                .collect::<Vec<_>>(),
            vec![b"ACGT".as_slice(), b"TA".as_slice()]
        );

        let _ = fs::remove_file(path);
    }

    #[test]
    fn parallel_gzip_inflation_matches_serial() {
        // Large enough that rapidgzip's speculative workers actually engage,
        // and two members so the multi-member handling is covered too.
        let mut body = Vec::new();
        body.extend_from_slice(b">r1\n");
        for i in 0..200_000u32 {
            let base = match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            body.push(base);
            if i % 70 == 69 {
                body.push(b'\n');
            }
        }
        body.extend_from_slice(b"\n>r2\nACGTNNACGT\n");

        let mut first = GzEncoder::new(Vec::new(), Compression::default());
        first.write_all(&body[..body.len() / 2]).unwrap();
        let mut bytes = first.finish().unwrap();
        let mut second = GzEncoder::new(Vec::new(), Compression::default());
        second.write_all(&body[body.len() / 2..]).unwrap();
        bytes.extend_from_slice(&second.finish().unwrap());

        let path =
            std::env::temp_dir().join(format!("cf3rs-input-{}-parallel.fa.gz", std::process::id()));
        fs::write(&path, bytes).unwrap();

        let collect = |workers: usize| {
            let mut fragments = Vec::<(u32, usize, Vec<u8>)>::new();
            let records = parse_fragments_borrowed_with(&path, 5, 2, workers, |fragment| {
                fragments.push((fragment.source_id, fragment.offset, fragment.seq.to_vec()));
                Ok(())
            })
            .unwrap();
            (records, fragments)
        };

        let serial = collect(1);
        let parallel = collect(4);
        assert_eq!(serial, parallel);
        assert_eq!(serial.0, 2);

        let _ = fs::remove_file(path);
    }

    #[test]
    fn parses_headerless_wrapped_sequence() {
        let path =
            std::env::temp_dir().join(format!("cf3rs-input-{}-plain.fna.gz", std::process::id()));
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(b"\nAACGT\nTNNACGT\n").unwrap();
        fs::write(&path, encoder.finish().unwrap()).unwrap();

        let mut fragments = Vec::new();
        let records = parse_fragments(&path, 9, 3, |fragment| {
            fragments.push(fragment);
            Ok(())
        })
        .unwrap();
        assert_eq!(records, 1);
        assert_eq!(fragments[0].seq, b"AACGTT");
        assert_eq!(fragments[1].seq, b"ACGT");
        assert_eq!(fragments[1].offset, 8);

        let _ = fs::remove_file(path);
    }
}

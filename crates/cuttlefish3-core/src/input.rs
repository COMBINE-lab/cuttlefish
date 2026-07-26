//! FASTA/FASTQ input discovery and zero-copy fragment parsing.
//!
//! Parsers split records at non-ACGT symbols. Borrowed callbacks are preferred
//! by the production partitioner so sequence data does not need to be copied.

use crate::dna::is_dna_ascii;
use crate::params::BuildParams;
use flate2::read::MultiGzDecoder;
use std::fs;
use std::io::{BufRead, BufReader, Read, Seek};
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
/// threads to decompress block-structured (BGZF) gzip input.
///
/// The decompressed byte stream is identical either way; only the work of
/// producing it is shared out. Plain gzip ignores the budget.
pub fn parse_fragments_borrowed_with<P, F>(
    path: P,
    source_id: u32,
    min_len: usize,
    inflate_workers: usize,
    mut on_fragment: F,
) -> Result<u64, InputError>
where
    P: AsRef<Path>,
    F: for<'a> FnMut(BorrowedSequenceFragment<'a>) -> Result<(), InputError>,
{
    let path = path.as_ref();
    let mut file = fs::File::open(path).map_err(|source| InputError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let input: Box<dyn Read> = if path.extension().is_some_and(|ext| ext == "gz") {
        // BGZF concatenates independent gzip members, so its blocks can be
        // inflated concurrently. A plain member cannot be split, and falls back
        // to the serial decoder.
        let mut head = [0u8; crate::bgzf::PROBE_BYTES];
        let probed = read_probe(&mut file, &mut head).map_err(|source| InputError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        file.rewind().map_err(|source| InputError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        if inflate_workers > 1 && crate::bgzf::is_bgzf(&head[..probed]) {
            Box::new(crate::bgzf::ParallelBgzfReader::new(file, inflate_workers))
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

/// Fills as much of `head` as the file provides, tolerating short reads.
fn read_probe<R: Read>(source: &mut R, head: &mut [u8]) -> std::io::Result<usize> {
    let mut filled = 0;
    while filled < head.len() {
        match source.read(&mut head[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(error) => return Err(error),
        }
    }
    Ok(filled)
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
        seq.extend(line.iter().copied().filter(|b| !b.is_ascii_whitespace()));
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
            seq.extend(line.iter().copied().filter(|b| !b.is_ascii_whitespace()));
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
            seq.extend(line.iter().copied().filter(|b| !b.is_ascii_whitespace()));
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

//! `cuttlefish compare`: decide whether two unitig FASTA files describe the
//! same graph.
//!
//! A maximal unitig may be emitted in either orientation, and a cyclic one at
//! any rotation, so two correct builds of the same input can disagree
//! record-for-record. Comparison therefore canonicalizes each record before
//! matching, and compares multisets rather than files.
//!
//! Graphs outgrow memory, so the two sides are chunk-sorted to disk in
//! parallel and then merged, which bounds resident set by the chunk size
//! rather than by the graph.

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::error::Error;
use std::fs::{self, File};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::mem;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, mpsc};

type Result<T> = std::result::Result<T, Box<dyn Error + Send + Sync>>;

struct Params {
    a: PathBuf,
    b: PathBuf,
    work_dir: PathBuf,
    k: usize,
    chunk_records: usize,
    threads: usize,
    reuse_chunks: bool,
    dump_mismatch: Option<PathBuf>,
    a_chunks: Option<PathBuf>,
    b_chunks: Option<PathBuf>,
    full_diff: bool,
}

struct ChunkTask {
    index: usize,
    records: Vec<Vec<u8>>,
}

/// Entry point for the `compare` subcommand. `args` is the argument list with
/// the subcommand name already consumed.
pub fn run<I>(args: I) -> Result<i32>
where
    I: Iterator<Item = String>,
{
    let Some(params) = parse_args(args)? else {
        return Ok(0); // --help
    };
    fs::create_dir_all(&params.work_dir)?;
    let a_dir = params.work_dir.join("a");
    let b_dir = params.work_dir.join("b");
    if !params.reuse_chunks
        && ((params.b_chunks.is_none() && b_dir.exists())
            || (params.a_chunks.is_none() && a_dir.exists()))
    {
        return Err(format!(
            "comparison directories already exist under {}",
            params.work_dir.display()
        )
        .into());
    }

    let (a_chunks, b_chunks) = if params.reuse_chunks {
        (existing_chunks(&a_dir)?, existing_chunks(&b_dir)?)
    } else {
        let a_chunks = if let Some(directory) = &params.a_chunks {
            existing_chunks(directory)?
        } else {
            eprintln!("cuttlefish compare: sorting A");
            build_chunks(
                &params.a,
                &a_dir,
                params.k,
                params.chunk_records,
                params.threads,
            )?
        };
        let b_chunks = if let Some(directory) = &params.b_chunks {
            existing_chunks(directory)?
        } else {
            eprintln!("cuttlefish compare: sorting B");
            build_chunks(
                &params.b,
                &b_dir,
                params.k,
                params.chunk_records,
                params.threads,
            )?
        };
        (a_chunks, b_chunks)
    };

    let mut a = ExternalMerge::new(&a_chunks)?;
    let mut b = ExternalMerge::new(&b_chunks)?;
    if params.full_diff {
        return compare_full_diff(&mut a, &mut b);
    }
    let mut count = 0u64;
    loop {
        match (a.next_record()?, b.next_record()?) {
            (Some(left), Some(right)) if left == right => count += 1,
            (Some(left), Some(right)) => {
                if let Some(directory) = &params.dump_mismatch {
                    fs::create_dir_all(directory)?;
                    fs::write(directory.join("a.seq"), &left)?;
                    fs::write(directory.join("b.seq"), &right)?;
                }
                return Err(format!(
                    "topology differs at sorted record {}:\n  A {}\n  B {}\n  common prefix={} common suffix={}",
                    count + 1,
                    describe_record(&left),
                    describe_record(&right),
                    common_prefix(&left, &right),
                    common_suffix(&left, &right),
                )
                .into());
            }
            (None, None) => break,
            (left, right) => {
                return Err(format!(
                    "topology counts differ after {count} matching records (A remaining={}, B remaining={})",
                    left.is_some(),
                    right.is_some()
                )
                .into());
            }
        }
    }
    println!("OK: {count} strand-normalized unitigs match");
    Ok(0)
}

/// Returns `None` when `--help` was requested and already printed.
fn parse_args<I>(args: I) -> Result<Option<Params>>
where
    I: Iterator<Item = String>,
{
    let mut a = None;
    let mut b = None;
    let mut work_dir = None;
    let mut k = None;
    let mut chunk_records = 1_000_000usize;
    let mut threads = 8usize;
    let mut reuse_chunks = false;
    let mut dump_mismatch = None;
    let mut a_chunks = None;
    let mut b_chunks = None;
    let mut full_diff = false;
    let mut args = args;
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-a" | "--a" => a = Some(PathBuf::from(argument_value(&mut args, &arg)?)),
            "-b" | "--b" => b = Some(PathBuf::from(argument_value(&mut args, &arg)?)),
            "-w" | "--work-dir" => work_dir = Some(PathBuf::from(argument_value(&mut args, &arg)?)),
            "-k" | "--kmer-len" => k = Some(argument_value(&mut args, &arg)?.parse()?),
            "--chunk-records" => chunk_records = argument_value(&mut args, &arg)?.parse()?,
            "-t" | "--threads" => threads = argument_value(&mut args, &arg)?.parse()?,
            "--reuse-chunks" => reuse_chunks = true,
            "--dump-mismatch" => {
                dump_mismatch = Some(PathBuf::from(argument_value(&mut args, &arg)?))
            }
            "--a-chunks" => a_chunks = Some(PathBuf::from(argument_value(&mut args, &arg)?)),
            "--b-chunks" => b_chunks = Some(PathBuf::from(argument_value(&mut args, &arg)?)),
            "--full-diff" => full_diff = true,
            "-h" | "--help" => {
                print_help();
                return Ok(None);
            }
            _ => return Err(format!("unknown argument {arg}").into()),
        }
    }
    if chunk_records == 0 || threads == 0 {
        return Err("--chunk-records and --threads must be positive".into());
    }
    Ok(Some(Params {
        a: a.ok_or("--a is required")?,
        b: b.ok_or("--b is required")?,
        work_dir: work_dir.ok_or("--work-dir is required")?,
        k: k.ok_or("-k is required")?,
        chunk_records,
        threads,
        reuse_chunks,
        dump_mismatch,
        a_chunks,
        b_chunks,
        full_diff,
    }))
}

pub fn print_help() {
    println!(
        "Compare two unitig FASTA files for graph equivalence, up to strand and cyclic rotation."
    );
    println!("Usage:");
    println!("  cuttlefish compare -a FILE -b FILE -w DIR -k K [OPTION...]");
    println!();
    println!(" common options:");
    println!("  -a, --a <arg>          first unitig FASTA");
    println!("  -b, --b <arg>          second unitig FASTA");
    println!("  -k, --kmer-len <arg>   k-mer length the graphs were built at");
    println!("  -w, --work-dir <arg>   scratch directory for the sorted chunks");
    println!("  -t, --threads <arg>    sorting threads (default: 8)");
    println!("      --full-diff        report all differences instead of the first");
    println!();
    println!(" diagnostic options:");
    println!("      --chunk-records <arg> records per sorted chunk (default: 1000000)");
    println!("      --reuse-chunks     reuse the chunks under --work-dir from a prior run");
    println!("      --a-chunks <arg>   reuse pre-sorted chunks for A from this directory");
    println!("      --b-chunks <arg>   reuse pre-sorted chunks for B from this directory");
    println!("      --dump-mismatch <arg> write the first differing pair to this directory");
    println!("  -h, --help             print usage");
}

fn compare_full_diff(a: &mut ExternalMerge, b: &mut ExternalMerge) -> Result<i32> {
    const MAX_SAMPLES: usize = 5;
    let mut left = a.next_record()?;
    let mut right = b.next_record()?;
    let mut matching = 0u64;
    let mut a_only = 0u64;
    let mut b_only = 0u64;
    let mut a_samples = Vec::new();
    let mut b_samples = Vec::new();
    while left.is_some() || right.is_some() {
        match (&left, &right) {
            (Some(a_record), Some(b_record)) if a_record == b_record => {
                matching += 1;
                left = a.next_record()?;
                right = b.next_record()?;
            }
            (Some(a_record), Some(b_record)) if a_record < b_record => {
                a_only += 1;
                if a_samples.len() < MAX_SAMPLES {
                    a_samples.push(describe_record(a_record));
                }
                left = a.next_record()?;
            }
            (Some(_), Some(b_record)) => {
                b_only += 1;
                if b_samples.len() < MAX_SAMPLES {
                    b_samples.push(describe_record(b_record));
                }
                right = b.next_record()?;
            }
            (Some(a_record), None) => {
                a_only += 1;
                if a_samples.len() < MAX_SAMPLES {
                    a_samples.push(describe_record(a_record));
                }
                left = a.next_record()?;
            }
            (None, Some(b_record)) => {
                b_only += 1;
                if b_samples.len() < MAX_SAMPLES {
                    b_samples.push(describe_record(b_record));
                }
                right = b.next_record()?;
            }
            (None, None) => break,
        }
    }
    eprintln!("cuttlefish compare: matching={matching} A-only={a_only} B-only={b_only}");
    for sample in a_samples {
        eprintln!("cuttlefish compare: A-only {sample}");
    }
    for sample in b_samples {
        eprintln!("cuttlefish compare: B-only {sample}");
    }
    if a_only == 0 && b_only == 0 {
        println!("OK: {matching} strand-normalized unitigs match");
        Ok(0)
    } else {
        Err("topology multisets differ".into())
    }
}

fn existing_chunks(directory: &Path) -> Result<Vec<PathBuf>> {
    let mut chunks = fs::read_dir(directory)?
        .map(|entry| entry.map(|entry| entry.path()))
        .collect::<io::Result<Vec<_>>>()?;
    chunks.retain(|path| path.is_file());
    chunks.sort_unstable();
    if chunks.is_empty() {
        return Err(format!("no chunks found in {}", directory.display()).into());
    }
    eprintln!(
        "cuttlefish compare: reusing {} chunk(s) from {}",
        chunks.len(),
        directory.display()
    );
    Ok(chunks)
}

fn describe_record(record: &[u8]) -> String {
    const EXCERPT: usize = 48;
    let mut hasher = DefaultHasher::new();
    record.hash(&mut hasher);
    let head = &record[..record.len().min(EXCERPT)];
    let tail = &record[record.len().saturating_sub(EXCERPT)..];
    format!(
        "len={} hash={:016x} head={} tail={}",
        record.len(),
        hasher.finish(),
        String::from_utf8_lossy(head),
        String::from_utf8_lossy(tail)
    )
}

fn common_prefix(left: &[u8], right: &[u8]) -> usize {
    left.iter()
        .zip(right)
        .take_while(|(left, right)| left == right)
        .count()
}

fn common_suffix(left: &[u8], right: &[u8]) -> usize {
    left.iter()
        .rev()
        .zip(right.iter().rev())
        .take_while(|(left, right)| left == right)
        .count()
}

fn argument_value<I>(args: &mut I, flag: &str) -> Result<String>
where
    I: Iterator<Item = String>,
{
    args.next()
        .ok_or_else(|| format!("missing value after {flag}").into())
}

fn build_chunks(
    fasta: &Path,
    directory: &Path,
    k: usize,
    chunk_records: usize,
    threads: usize,
) -> Result<Vec<PathBuf>> {
    fs::create_dir(directory)?;
    let (task_tx, task_rx) = mpsc::sync_channel::<ChunkTask>(threads * 2);
    let task_rx = Arc::new(Mutex::new(task_rx));
    let (result_tx, result_rx) = mpsc::channel::<Result<(usize, PathBuf, usize)>>();
    let directory = Arc::new(directory.to_path_buf());
    let mut task_count = 0usize;

    std::thread::scope(|scope| -> Result<()> {
        for _ in 0..threads {
            let task_rx = Arc::clone(&task_rx);
            let result_tx = result_tx.clone();
            let directory = Arc::clone(&directory);
            scope.spawn(move || {
                loop {
                    let task = match task_rx.lock().expect("chunk receiver poisoned").recv() {
                        Ok(task) => task,
                        Err(_) => break,
                    };
                    let result = write_chunk(task, &directory, k);
                    if result_tx.send(result).is_err() {
                        break;
                    }
                }
            });
        }
        drop(result_tx);
        read_fasta_chunks(fasta, chunk_records, |records| {
            let index = task_count;
            task_count += 1;
            task_tx
                .send(ChunkTask { index, records })
                .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "chunk workers stopped"))
        })?;
        drop(task_tx);
        Ok(())
    })?;

    let mut chunks = vec![PathBuf::new(); task_count];
    let mut record_count = 0u64;
    for result in result_rx {
        let (index, path, records) = result?;
        chunks[index] = path;
        record_count += records as u64;
    }
    if chunks.iter().any(|path| path.as_os_str().is_empty()) {
        return Err("a chunk worker exited without producing output".into());
    }
    eprintln!(
        "cuttlefish compare: wrote {} sorted chunk(s), {} record(s)",
        chunks.len(),
        record_count
    );
    Ok(chunks)
}

fn read_fasta_chunks<F>(path: &Path, chunk_records: usize, mut send: F) -> Result<()>
where
    F: FnMut(Vec<Vec<u8>>) -> io::Result<()>,
{
    let mut reader = BufReader::with_capacity(8 * 1024 * 1024, File::open(path)?);
    let mut line = Vec::new();
    let mut sequence = Vec::new();
    let mut records = Vec::with_capacity(chunk_records);
    let mut saw_header = false;
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        if line.first() == Some(&b'>') {
            if saw_header {
                records.push(mem::take(&mut sequence));
                if records.len() == chunk_records {
                    send(mem::replace(
                        &mut records,
                        Vec::with_capacity(chunk_records),
                    ))?;
                }
            }
            saw_header = true;
        } else if saw_header && !line.is_empty() {
            sequence.extend(line.iter().map(u8::to_ascii_uppercase));
        }
    }
    if saw_header {
        records.push(sequence);
    }
    if !records.is_empty() {
        send(records)?;
    }
    Ok(())
}

fn write_chunk(mut task: ChunkTask, directory: &Path, k: usize) -> Result<(usize, PathBuf, usize)> {
    for record in &mut task.records {
        *record = canonical_unitig(record, k);
    }
    task.records.sort_unstable();
    let path = directory.join(format!("topology-{:06}.txt", task.index));
    let mut output = BufWriter::with_capacity(8 * 1024 * 1024, File::create(&path)?);
    for record in &task.records {
        output.write_all(record)?;
        output.write_all(b"\n")?;
    }
    output.flush()?;
    Ok((task.index, path, task.records.len()))
}

fn canonical_unitig(sequence: &[u8], k: usize) -> Vec<u8> {
    if sequence.len() > k && sequence[..k - 1] == sequence[sequence.len() - k + 1..] {
        let body_len = sequence.len() - k + 1;
        let forward = minimal_rotation(&sequence[..body_len]);
        let reverse = minimal_rotation(&reverse_complement(&sequence[..body_len]));
        let mut body = if forward <= reverse { forward } else { reverse };
        for index in 0..k - 1 {
            body.push(body[index % body_len]);
        }
        body
    } else {
        let reverse = reverse_complement(sequence);
        if sequence <= reverse.as_slice() {
            sequence.to_vec()
        } else {
            reverse
        }
    }
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => *base,
        })
        .collect()
}

fn minimal_rotation(sequence: &[u8]) -> Vec<u8> {
    if sequence.len() <= 1 {
        return sequence.to_vec();
    }
    let n = sequence.len();
    let (mut i, mut j, mut matched) = (0usize, 1usize, 0usize);
    while i < n && j < n && matched < n {
        match sequence[(i + matched) % n].cmp(&sequence[(j + matched) % n]) {
            std::cmp::Ordering::Equal => matched += 1,
            std::cmp::Ordering::Greater => {
                i += matched + 1;
                if i <= j {
                    i = j + 1;
                }
                matched = 0;
            }
            std::cmp::Ordering::Less => {
                j += matched + 1;
                if j <= i {
                    j = i + 1;
                }
                matched = 0;
            }
        }
    }
    let start = i.min(j);
    sequence[start..]
        .iter()
        .chain(&sequence[..start])
        .copied()
        .collect()
}

struct ExternalMerge {
    readers: Vec<BufReader<File>>,
    heap: BinaryHeap<Reverse<(Vec<u8>, usize)>>,
}

impl ExternalMerge {
    fn new(paths: &[PathBuf]) -> io::Result<Self> {
        let mut readers = Vec::with_capacity(paths.len());
        let mut heap = BinaryHeap::new();
        for (index, path) in paths.iter().enumerate() {
            let mut reader = BufReader::with_capacity(1024 * 1024, File::open(path)?);
            if let Some(record) = read_chunk_record(&mut reader)? {
                heap.push(Reverse((record, index)));
            }
            readers.push(reader);
        }
        Ok(Self { readers, heap })
    }

    fn next_record(&mut self) -> io::Result<Option<Vec<u8>>> {
        let Some(Reverse((record, index))) = self.heap.pop() else {
            return Ok(None);
        };
        if let Some(next) = read_chunk_record(&mut self.readers[index])? {
            self.heap.push(Reverse((next, index)));
        }
        Ok(Some(record))
    }
}

fn read_chunk_record(reader: &mut BufReader<File>) -> io::Result<Option<Vec<u8>>> {
    let mut record = Vec::new();
    if reader.read_until(b'\n', &mut record)? == 0 {
        return Ok(None);
    }
    if record.last() == Some(&b'\n') {
        record.pop();
    }
    Ok(Some(record))
}

#[cfg(test)]
mod tests {
    use super::{canonical_unitig, reverse_complement};

    #[test]
    fn canonicalizes_linear_strands() {
        let sequence = b"AACCGGTTACG";
        assert_eq!(
            canonical_unitig(sequence, 5),
            canonical_unitig(&reverse_complement(sequence), 5)
        );
    }

    #[test]
    fn canonicalizes_short_period_cycles_without_truncation() {
        let body = b"ACGTC";
        let mut sequence = body.to_vec();
        for index in 0..30 {
            sequence.push(body[index % body.len()]);
        }
        let canonical = canonical_unitig(&sequence, 31);
        assert_eq!(canonical.len(), sequence.len());

        let rotated_body = b"GTCAC";
        let mut rotated = rotated_body.to_vec();
        for index in 0..30 {
            rotated.push(rotated_body[index % rotated_body.len()]);
        }
        assert_eq!(canonical, canonical_unitig(&rotated, 31));
    }
}

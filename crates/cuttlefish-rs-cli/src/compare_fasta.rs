use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::env;
use std::error::Error;
use std::fs::{self, File};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::mem;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, mpsc};

type Result<T> = std::result::Result<T, Box<dyn Error + Send + Sync>>;

struct Params {
    cpp: PathBuf,
    rust: PathBuf,
    work_dir: PathBuf,
    k: usize,
    chunk_records: usize,
    threads: usize,
    reuse_chunks: bool,
    dump_mismatch: Option<PathBuf>,
    cpp_chunks: Option<PathBuf>,
    rust_chunks: Option<PathBuf>,
    full_diff: bool,
}

struct ChunkTask {
    index: usize,
    records: Vec<Vec<u8>>,
}

fn main() {
    if let Err(error) = run() {
        eprintln!("cf3-compare-fasta: {error}");
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let params = parse_args()?;
    fs::create_dir_all(&params.work_dir)?;
    let cpp_dir = params.work_dir.join("cpp");
    let rust_dir = params.work_dir.join("rust");
    if !params.reuse_chunks
        && ((params.rust_chunks.is_none() && rust_dir.exists())
            || (params.cpp_chunks.is_none() && cpp_dir.exists()))
    {
        return Err(format!(
            "comparison directories already exist under {}",
            params.work_dir.display()
        )
        .into());
    }

    let (cpp_chunks, rust_chunks) = if params.reuse_chunks {
        (existing_chunks(&cpp_dir)?, existing_chunks(&rust_dir)?)
    } else {
        let cpp_chunks = if let Some(directory) = &params.cpp_chunks {
            existing_chunks(directory)?
        } else {
            eprintln!("cf3-compare-fasta: sorting C++ topology");
            build_chunks(
                &params.cpp,
                &cpp_dir,
                params.k,
                params.chunk_records,
                params.threads,
            )?
        };
        let rust_chunks = if let Some(directory) = &params.rust_chunks {
            existing_chunks(directory)?
        } else {
            eprintln!("cf3-compare-fasta: sorting Rust topology");
            build_chunks(
                &params.rust,
                &rust_dir,
                params.k,
                params.chunk_records,
                params.threads,
            )?
        };
        (cpp_chunks, rust_chunks)
    };

    let mut cpp = ExternalMerge::new(&cpp_chunks)?;
    let mut rust = ExternalMerge::new(&rust_chunks)?;
    if params.full_diff {
        return compare_full_diff(&mut cpp, &mut rust);
    }
    let mut count = 0u64;
    loop {
        match (cpp.next_record()?, rust.next_record()?) {
            (Some(left), Some(right)) if left == right => count += 1,
            (Some(left), Some(right)) => {
                if let Some(directory) = &params.dump_mismatch {
                    fs::create_dir_all(directory)?;
                    fs::write(directory.join("cpp.seq"), &left)?;
                    fs::write(directory.join("rust.seq"), &right)?;
                }
                return Err(format!(
                    "topology differs at sorted record {}:\n  C++ {}\n  Rust {}\n  common prefix={} common suffix={}",
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
                    "topology counts differ after {count} matching records (C++ remaining={}, Rust remaining={})",
                    left.is_some(),
                    right.is_some()
                )
                .into());
            }
        }
    }
    println!("OK: {count} strand-normalized unitigs match C++ topology");
    Ok(())
}

fn parse_args() -> Result<Params> {
    let mut cpp = None;
    let mut rust = None;
    let mut work_dir = None;
    let mut k = None;
    let mut chunk_records = 1_000_000usize;
    let mut threads = 8usize;
    let mut reuse_chunks = false;
    let mut dump_mismatch = None;
    let mut cpp_chunks = None;
    let mut rust_chunks = None;
    let mut full_diff = false;
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--cpp-fasta" => cpp = Some(PathBuf::from(argument_value(&mut args, &arg)?)),
            "--rust-fasta" => rust = Some(PathBuf::from(argument_value(&mut args, &arg)?)),
            "--work-dir" => work_dir = Some(PathBuf::from(argument_value(&mut args, &arg)?)),
            "-k" | "--kmer-len" => k = Some(argument_value(&mut args, &arg)?.parse()?),
            "--chunk-records" => chunk_records = argument_value(&mut args, &arg)?.parse()?,
            "--threads" => threads = argument_value(&mut args, &arg)?.parse()?,
            "--reuse-chunks" => reuse_chunks = true,
            "--dump-mismatch" => {
                dump_mismatch = Some(PathBuf::from(argument_value(&mut args, &arg)?))
            }
            "--cpp-chunks" => cpp_chunks = Some(PathBuf::from(argument_value(&mut args, &arg)?)),
            "--rust-chunks" => rust_chunks = Some(PathBuf::from(argument_value(&mut args, &arg)?)),
            "--full-diff" => full_diff = true,
            "-h" | "--help" => {
                println!(
                    "usage: cf3-compare-fasta --cpp-fasta FILE --rust-fasta FILE --work-dir DIR -k K [--chunk-records N] [--threads N]"
                );
                std::process::exit(0);
            }
            _ => return Err(format!("unknown argument {arg}").into()),
        }
    }
    if chunk_records == 0 || threads == 0 {
        return Err("--chunk-records and --threads must be positive".into());
    }
    Ok(Params {
        cpp: cpp.ok_or("--cpp-fasta is required")?,
        rust: rust.ok_or("--rust-fasta is required")?,
        work_dir: work_dir.ok_or("--work-dir is required")?,
        k: k.ok_or("-k is required")?,
        chunk_records,
        threads,
        reuse_chunks,
        dump_mismatch,
        cpp_chunks,
        rust_chunks,
        full_diff,
    })
}

fn compare_full_diff(cpp: &mut ExternalMerge, rust: &mut ExternalMerge) -> Result<()> {
    const MAX_SAMPLES: usize = 5;
    let mut left = cpp.next_record()?;
    let mut right = rust.next_record()?;
    let mut matching = 0u64;
    let mut cpp_only = 0u64;
    let mut rust_only = 0u64;
    let mut cpp_samples = Vec::new();
    let mut rust_samples = Vec::new();
    while left.is_some() || right.is_some() {
        match (&left, &right) {
            (Some(cpp_record), Some(rust_record)) if cpp_record == rust_record => {
                matching += 1;
                left = cpp.next_record()?;
                right = rust.next_record()?;
            }
            (Some(cpp_record), Some(rust_record)) if cpp_record < rust_record => {
                cpp_only += 1;
                if cpp_samples.len() < MAX_SAMPLES {
                    cpp_samples.push(describe_record(cpp_record));
                }
                left = cpp.next_record()?;
            }
            (Some(_), Some(rust_record)) => {
                rust_only += 1;
                if rust_samples.len() < MAX_SAMPLES {
                    rust_samples.push(describe_record(rust_record));
                }
                right = rust.next_record()?;
            }
            (Some(cpp_record), None) => {
                cpp_only += 1;
                if cpp_samples.len() < MAX_SAMPLES {
                    cpp_samples.push(describe_record(cpp_record));
                }
                left = cpp.next_record()?;
            }
            (None, Some(rust_record)) => {
                rust_only += 1;
                if rust_samples.len() < MAX_SAMPLES {
                    rust_samples.push(describe_record(rust_record));
                }
                right = rust.next_record()?;
            }
            (None, None) => break,
        }
    }
    eprintln!("cf3-compare-fasta: matching={matching} C++-only={cpp_only} Rust-only={rust_only}");
    for sample in cpp_samples {
        eprintln!("cf3-compare-fasta: C++-only {sample}");
    }
    for sample in rust_samples {
        eprintln!("cf3-compare-fasta: Rust-only {sample}");
    }
    if cpp_only == 0 && rust_only == 0 {
        println!("OK: {matching} strand-normalized unitigs match C++ topology");
        Ok(())
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
        "cf3-compare-fasta: reusing {} chunk(s) from {}",
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
        "cf3-compare-fasta: wrote {} sorted chunk(s), {} record(s)",
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

use crate::hash::{FastBuildHasher, hash_u64};
use crate::state::{ColorCoordinate, UnitigColor};
use scc::{HashMap as SccHashMap, hash_map::Entry as SccEntry};
use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Seek, Write};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::{
    Mutex,
    atomic::{AtomicBool, AtomicU64, Ordering},
};

const COLOR_BUFFER_BYTES: usize = 128 * 1024;

pub fn reverse_color_runs(runs: &[UnitigColor], vertex_count: u32) -> Vec<UnitigColor> {
    let mut reversed = Vec::with_capacity(runs.len());
    for index in (0..runs.len()).rev() {
        let end = runs
            .get(index + 1)
            .map_or(vertex_count, |next| next.offset());
        reversed.push(UnitigColor::new(
            vertex_count - end,
            ColorCoordinate::from_u40(runs[index].coordinate()),
        ));
    }
    reversed
}

pub fn reverse_color_runs_in_place(runs: &mut [UnitigColor], vertex_count: u32) {
    if runs.is_empty() {
        return;
    }
    runs.reverse();
    for index in (1..runs.len()).rev() {
        runs[index] = UnitigColor::new(
            vertex_count - runs[index - 1].offset(),
            ColorCoordinate::from_u40(runs[index].coordinate()),
        );
    }
    runs[0] = UnitigColor::new(0, ColorCoordinate::from_u40(runs[0].coordinate()));
}

pub fn append_color_runs(
    output: &mut Vec<UnitigColor>,
    output_vertex_count: u32,
    runs: &[UnitigColor],
    unitig_vertex_count: u32,
    reverse: bool,
) {
    if reverse {
        output.reserve(runs.len());
        for index in (0..runs.len()).rev() {
            if index == runs.len() - 1
                && output
                    .last()
                    .is_some_and(|left| left.coordinate() == runs[index].coordinate())
            {
                continue;
            }
            let end = runs
                .get(index + 1)
                .map_or(unitig_vertex_count, |next| next.offset());
            output.push(UnitigColor::new(
                output_vertex_count + unitig_vertex_count - end - 1,
                ColorCoordinate::from_u40(runs[index].coordinate()),
            ));
        }
    } else {
        let start = usize::from(
            output
                .last()
                .zip(runs.first())
                .is_some_and(|(left, right)| left.coordinate() == right.coordinate()),
        );
        output.reserve(runs.len().saturating_sub(start));
        for run in &runs[start..] {
            output.push(UnitigColor::new(
                output_vertex_count + run.offset() - 1,
                ColorCoordinate::from_u40(run.coordinate()),
            ));
        }
    }
}

pub fn rotate_cycle_color_runs(
    runs: &[UnitigColor],
    vertex_count: u32,
    pivot: u32,
    reverse: bool,
) -> Vec<UnitigColor> {
    if runs.is_empty() || vertex_count == 0 {
        return Vec::new();
    }
    let oriented = if reverse {
        reverse_color_runs(runs, vertex_count)
    } else {
        runs.to_vec()
    };
    let mut colors = vec![0u64; vertex_count as usize];
    for (index, run) in oriented.iter().enumerate() {
        let end = oriented
            .get(index + 1)
            .map_or(vertex_count, |next| next.offset());
        colors[run.offset() as usize..end as usize].fill(run.coordinate());
    }
    colors.rotate_left((pivot % vertex_count) as usize);
    compress_color_coordinates(&colors)
}

fn compress_color_coordinates(colors: &[u64]) -> Vec<UnitigColor> {
    let mut runs = Vec::new();
    for (offset, &coordinate) in colors.iter().enumerate() {
        if runs
            .last()
            .is_none_or(|previous: &UnitigColor| previous.coordinate() != coordinate)
        {
            runs.push(UnitigColor::new(
                offset as u32,
                ColorCoordinate::from_u40(coordinate),
            ));
        }
    }
    runs
}

pub struct ColorRunSidecarWriter {
    run_path: PathBuf,
    runs: BufWriter<File>,
    run_bytes: u64,
    run_count: u64,
    unitigs: u64,
}

impl ColorRunSidecarWriter {
    pub fn create(path_prefix: impl AsRef<Path>) -> Result<Self, ColorError> {
        let run_path = path_prefix.as_ref().with_extension("color-runs");
        let run_file = OpenOptions::new()
            .create(true)
            .truncate(true)
            .read(true)
            .write(true)
            .open(&run_path)
            .map_err(|source| ColorError::Io {
                path: run_path.clone(),
                source,
            })?;
        Ok(Self {
            run_path,
            runs: BufWriter::with_capacity(COLOR_BUFFER_BYTES, run_file),
            run_bytes: 0,
            run_count: 0,
            unitigs: 0,
        })
    }

    pub fn position(&self) -> u64 {
        self.run_bytes
    }

    pub fn write_unitig(&mut self, colors: &[UnitigColor]) -> Result<(), ColorError> {
        let count = u32::try_from(colors.len()).map_err(|_| ColorError::TooManyColorRuns)?;
        let mut encoded_count = Vec::with_capacity(5);
        append_varint_u32(&mut encoded_count, count);
        self.runs
            .write_all(&encoded_count)
            .map_err(|source| ColorError::Io {
                path: self.run_path.clone(),
                source,
            })?;
        for color in colors {
            self.runs
                .write_all(&color.raw().to_le_bytes())
                .map_err(|source| ColorError::Io {
                    path: self.run_path.clone(),
                    source,
                })?;
        }
        self.run_bytes += encoded_count.len() as u64 + u64::from(count) * 8;
        self.run_count += u64::from(count);
        self.unitigs += 1;
        Ok(())
    }

    pub fn finish(mut self) -> Result<ColorRunSidecar, ColorError> {
        self.runs.flush().map_err(|source| ColorError::Io {
            path: self.run_path.clone(),
            source,
        })?;
        self.runs.into_inner().map_err(|error| ColorError::Io {
            path: self.run_path.clone(),
            source: error.into_error(),
        })?;
        Ok(ColorRunSidecar {
            run_path: self.run_path,
            unitigs: self.unitigs,
            runs: self.run_count,
        })
    }
}

pub struct ConcurrentColorRunSidecarWriter {
    run_path: PathBuf,
    runs: File,
    run_bytes: AtomicU64,
    run_count: AtomicU64,
    unitigs: AtomicU64,
}

impl ConcurrentColorRunSidecarWriter {
    pub fn create(path_prefix: impl AsRef<Path>) -> Result<Self, ColorError> {
        let run_path = path_prefix.as_ref().with_extension("color-runs");
        let runs = OpenOptions::new()
            .create(true)
            .truncate(true)
            .read(true)
            .write(true)
            .open(&run_path)
            .map_err(|source| ColorError::Io {
                path: run_path.clone(),
                source,
            })?;
        Ok(Self {
            run_path,
            runs,
            run_bytes: AtomicU64::new(0),
            run_count: AtomicU64::new(0),
            unitigs: AtomicU64::new(0),
        })
    }

    pub fn write_unitigs(&self, unitigs: &[Vec<UnitigColor>]) -> Result<u64, ColorError> {
        let run_count = unitigs.iter().map(Vec::len).sum::<usize>();
        let byte_len = unitigs
            .len()
            .checked_mul(5)
            .and_then(|bytes| bytes.checked_add(run_count.checked_mul(8)?))
            .ok_or(ColorError::TooManyColorRuns)?;
        let mut encoded = Vec::with_capacity(byte_len);
        for colors in unitigs {
            let count = u32::try_from(colors.len()).map_err(|_| ColorError::TooManyColorRuns)?;
            append_varint_u32(&mut encoded, count);
            for color in colors {
                encoded.extend_from_slice(&color.raw().to_le_bytes());
            }
        }
        let offset = self
            .run_bytes
            .fetch_add(encoded.len() as u64, Ordering::Relaxed);
        self.runs
            .write_all_at(&encoded, offset)
            .map_err(|source| ColorError::Io {
                path: self.run_path.clone(),
                source,
            })?;
        self.run_count
            .fetch_add(run_count as u64, Ordering::Relaxed);
        self.unitigs
            .fetch_add(unitigs.len() as u64, Ordering::Relaxed);
        Ok(offset)
    }

    pub fn finish(self) -> Result<ColorRunSidecar, ColorError> {
        self.runs.sync_data().map_err(|source| ColorError::Io {
            path: self.run_path.clone(),
            source,
        })?;
        Ok(ColorRunSidecar {
            run_path: self.run_path,
            unitigs: self.unitigs.load(Ordering::Relaxed),
            runs: self.run_count.load(Ordering::Relaxed),
        })
    }
}

#[derive(Debug, Clone)]
pub struct ColorRunSidecar {
    pub run_path: PathBuf,
    pub unitigs: u64,
    pub runs: u64,
}

impl ColorRunSidecar {
    pub fn reader_at(&self, byte_offset: u64) -> Result<ColorRunStreamReader, ColorError> {
        let mut file = File::open(&self.run_path).map_err(|source| ColorError::Io {
            path: self.run_path.clone(),
            source,
        })?;
        file.seek(std::io::SeekFrom::Start(byte_offset))
            .map_err(|source| ColorError::Io {
                path: self.run_path.clone(),
                source,
            })?;
        Ok(ColorRunStreamReader {
            path: self.run_path.clone(),
            input: BufReader::with_capacity(COLOR_BUFFER_BYTES, file),
        })
    }

    pub fn read_unitig(&self, unitig: u64) -> Result<Vec<UnitigColor>, ColorError> {
        if unitig >= self.unitigs {
            return Err(ColorError::MalformedUnitigIndex(unitig));
        }
        let mut reader = self.reader_at(0)?;
        for index in 0..=unitig {
            let colors = reader.read_next()?;
            if index == unitig {
                return Ok(colors);
            }
        }
        Err(ColorError::MalformedUnitigIndex(unitig))
    }
}

pub struct ColorRunStreamReader {
    path: PathBuf,
    input: BufReader<File>,
}

impl ColorRunStreamReader {
    pub fn read_next(&mut self) -> Result<Vec<UnitigColor>, ColorError> {
        let mut colors = Vec::new();
        self.read_next_into(&mut colors)?;
        Ok(colors)
    }

    pub fn read_next_into(&mut self, colors: &mut Vec<UnitigColor>) -> Result<(), ColorError> {
        let count = read_varint(&mut self.input).map_err(|source| ColorError::Io {
            path: self.path.clone(),
            source,
        })?;
        colors.clear();
        colors.reserve(count as usize);
        for _ in 0..count {
            let mut raw = [0u8; 8];
            self.input
                .read_exact(&mut raw)
                .map_err(|source| ColorError::Io {
                    path: self.path.clone(),
                    source,
                })?;
            let raw = u64::from_le_bytes(raw);
            colors.push(UnitigColor::new(
                (raw & 0xff_ffff) as u32,
                ColorCoordinate::from_u40(raw >> 24),
            ));
        }
        Ok(())
    }
}

pub struct ConcurrentColorRepository {
    dir: PathBuf,
    table: AtomicColorTable,
    overflow: SccHashMap<u64, ColorCoordinate, FastBuildHasher>,
    workers: Vec<Mutex<ColorWorkerWriter>>,
}

struct AtomicColorSlot {
    state: AtomicU64,
    key: AtomicU64,
    value: AtomicU64,
}

struct AtomicColorTable {
    slots: Vec<AtomicColorSlot>,
    mask: usize,
    saturated: AtomicBool,
}

enum AtomicColorEntry<'a> {
    Occupied(ColorCoordinate),
    Vacant(&'a AtomicColorSlot),
    Full,
}

impl AtomicColorTable {
    const EMPTY: u64 = 0;
    const BUSY: u64 = 1;

    fn with_expected_entries(expected: usize) -> Self {
        let capacity = expected
            .max(8)
            .saturating_mul(4)
            .div_ceil(3)
            .next_power_of_two()
            .min(64 * 1024 * 1024);
        Self {
            slots: (0..capacity)
                .map(|_| AtomicColorSlot {
                    state: AtomicU64::new(Self::EMPTY),
                    key: AtomicU64::new(0),
                    value: AtomicU64::new(0),
                })
                .collect(),
            mask: capacity - 1,
            saturated: AtomicBool::new(false),
        }
    }

    #[inline]
    fn entry(&self, key: u64) -> AtomicColorEntry<'_> {
        let hash = hash_u64(key, 0);
        let tag = hash.wrapping_add(2).max(2);
        let mut index = hash as usize & self.mask;
        for _ in 0..self.slots.len() {
            let slot = &self.slots[index];
            let mut state = slot.state.load(Ordering::Acquire);
            while state == Self::BUSY {
                std::hint::spin_loop();
                state = slot.state.load(Ordering::Acquire);
            }
            if state == Self::EMPTY {
                if slot
                    .state
                    .compare_exchange(
                        Self::EMPTY,
                        Self::BUSY,
                        Ordering::Acquire,
                        Ordering::Relaxed,
                    )
                    .is_ok()
                {
                    slot.key.store(key, Ordering::Relaxed);
                    return AtomicColorEntry::Vacant(slot);
                }
                continue;
            }
            if state == tag && slot.key.load(Ordering::Relaxed) == key {
                return AtomicColorEntry::Occupied(ColorCoordinate::from_u40(
                    slot.value.load(Ordering::Relaxed),
                ));
            }
            index = (index + 1) & self.mask;
        }
        self.saturated.store(true, Ordering::Release);
        AtomicColorEntry::Full
    }

    #[inline]
    fn get(&self, key: u64) -> Option<ColorCoordinate> {
        let hash = hash_u64(key, 0);
        let tag = hash.wrapping_add(2).max(2);
        let mut index = hash as usize & self.mask;
        for _ in 0..self.slots.len() {
            let slot = &self.slots[index];
            let mut state = slot.state.load(Ordering::Acquire);
            while state == Self::BUSY {
                std::hint::spin_loop();
                state = slot.state.load(Ordering::Acquire);
            }
            if state == Self::EMPTY {
                return None;
            }
            if state == tag && slot.key.load(Ordering::Relaxed) == key {
                return Some(ColorCoordinate::from_u40(
                    slot.value.load(Ordering::Relaxed),
                ));
            }
            index = (index + 1) & self.mask;
        }
        None
    }

    fn publish(slot: &AtomicColorSlot, coordinate: ColorCoordinate) {
        let key = slot.key.load(Ordering::Relaxed);
        let tag = hash_u64(key, 0).wrapping_add(2).max(2);
        slot.value.store(coordinate.as_u40(), Ordering::Relaxed);
        slot.state.store(tag, Ordering::Release);
    }

    fn abort(slot: &AtomicColorSlot) {
        slot.state.store(Self::EMPTY, Ordering::Release);
    }

    #[inline]
    fn is_saturated(&self) -> bool {
        self.saturated.load(Ordering::Acquire)
    }
}

struct ColorWorkerWriter {
    records: u32,
    output: BufWriter<File>,
}

impl ConcurrentColorRepository {
    pub(crate) fn worker_count(&self) -> usize {
        self.workers.len()
    }

    pub fn create(
        dir: impl AsRef<Path>,
        workers: usize,
        expected_colors: usize,
    ) -> Result<Self, ColorError> {
        if workers == 0 || workers > 256 {
            return Err(ColorError::InvalidWorkerCount(workers));
        }
        let dir = dir.as_ref().to_path_buf();
        if dir.exists() {
            fs::remove_dir_all(&dir).map_err(|source| ColorError::Io {
                path: dir.clone(),
                source,
            })?;
        }
        fs::create_dir_all(&dir).map_err(|source| ColorError::Io {
            path: dir.clone(),
            source,
        })?;
        let mut worker_writers = Vec::with_capacity(workers);
        for worker in 0..workers {
            let path = color_worker_path(&dir, worker);
            let file = File::create(&path).map_err(|source| ColorError::Io {
                path: path.clone(),
                source,
            })?;
            worker_writers.push(Mutex::new(ColorWorkerWriter {
                records: 0,
                output: BufWriter::with_capacity(COLOR_BUFFER_BYTES, file),
            }));
        }
        Ok(Self {
            dir,
            table: AtomicColorTable::with_expected_entries(expected_colors),
            overflow: SccHashMap::with_capacity_and_hasher(8, FastBuildHasher::default()),
            workers: worker_writers,
        })
    }

    pub fn resolve_or_insert(
        &self,
        color_hash: u64,
        sources: &[u32],
        worker: usize,
    ) -> Result<ColorCoordinate, ColorError> {
        if sources.is_empty() {
            return self.get(color_hash).ok_or(ColorError::MalformedSourceSet);
        }
        if !sources.windows(2).all(|pair| pair[0] < pair[1]) {
            return Err(ColorError::MalformedSourceSet);
        }
        if worker >= self.workers.len() {
            return Err(ColorError::InvalidWorkerCount(worker + 1));
        }
        match self.table.entry(color_hash) {
            AtomicColorEntry::Occupied(coordinate) => Ok(coordinate),
            AtomicColorEntry::Vacant(slot) => {
                let mut writer = self.workers[worker]
                    .lock()
                    .map_err(|_| ColorError::PoisonedWriter)?;
                let index = writer.records;
                if let Err(source) = write_source_set(&mut writer.output, sources) {
                    AtomicColorTable::abort(slot);
                    return Err(ColorError::Io {
                        path: color_worker_path(&self.dir, worker),
                        source,
                    });
                }
                let Some(next_records) = writer.records.checked_add(1) else {
                    AtomicColorTable::abort(slot);
                    return Err(ColorError::TooManyColors);
                };
                writer.records = next_records;
                let coordinate = ColorCoordinate::discovered(worker as u64, index as u64);
                AtomicColorTable::publish(slot, coordinate);
                Ok(coordinate)
            }
            AtomicColorEntry::Full => match self.overflow.entry_sync(color_hash) {
                SccEntry::Occupied(entry) => Ok(*entry.get()),
                SccEntry::Vacant(entry) => {
                    let mut writer = self.workers[worker]
                        .lock()
                        .map_err(|_| ColorError::PoisonedWriter)?;
                    let index = writer.records;
                    write_source_set(&mut writer.output, sources).map_err(|source| {
                        ColorError::Io {
                            path: color_worker_path(&self.dir, worker),
                            source,
                        }
                    })?;
                    writer.records = writer
                        .records
                        .checked_add(1)
                        .ok_or(ColorError::TooManyColors)?;
                    let coordinate = ColorCoordinate::discovered(worker as u64, index as u64);
                    entry.insert_entry(coordinate);
                    Ok(coordinate)
                }
            },
        }
    }

    pub fn get(&self, color_hash: u64) -> Option<ColorCoordinate> {
        self.table.get(color_hash).or_else(|| {
            self.table.is_saturated().then(|| {
                self.overflow
                    .read_sync(&color_hash, |_, coordinate| *coordinate)
            })?
        })
    }

    pub fn finish(&self) -> Result<ColorRepositoryManifest, ColorError> {
        eprintln!(
            "cuttlefish3-rs: color repository table capacity {}, overflow entries {}",
            self.table.slots.len(),
            self.overflow.len()
        );
        let mut records = Vec::with_capacity(self.workers.len());
        for (worker, writer) in self.workers.iter().enumerate() {
            let mut writer = writer.lock().map_err(|_| ColorError::PoisonedWriter)?;
            writer.output.flush().map_err(|source| ColorError::Io {
                path: color_worker_path(&self.dir, worker),
                source,
            })?;
            records.push(writer.records);
        }
        let manifest = ColorRepositoryManifest {
            dir: self.dir.clone(),
            records,
        };
        manifest.write()?;
        Ok(manifest)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ColorRepositoryManifest {
    pub dir: PathBuf,
    pub records: Vec<u32>,
}

impl ColorRepositoryManifest {
    fn write(&self) -> Result<(), ColorError> {
        let path = self.dir.join("manifest.tsv");
        let mut output = BufWriter::new(File::create(&path).map_err(|source| ColorError::Io {
            path: path.clone(),
            source,
        })?);
        writeln!(output, "worker\trecords\tpath").map_err(|source| ColorError::Io {
            path: path.clone(),
            source,
        })?;
        for (worker, records) in self.records.iter().enumerate() {
            writeln!(output, "{worker}\t{records}\t{:03}.colors", worker,).map_err(|source| {
                ColorError::Io {
                    path: path.clone(),
                    source,
                }
            })?;
        }
        output
            .flush()
            .map_err(|source| ColorError::Io { path, source })
    }

    pub fn write_metadata(
        &self,
        k: u16,
        fasta_path: &Path,
        sources: &[PathBuf],
    ) -> Result<(), ColorError> {
        let metadata_path = self.dir.join("metadata.tsv");
        let mut metadata =
            BufWriter::new(
                File::create(&metadata_path).map_err(|source| ColorError::Io {
                    path: metadata_path.clone(),
                    source,
                })?,
            );
        writeln!(metadata, "format\tcf3rs-color-repository-v1")
            .and_then(|_| writeln!(metadata, "k\t{k}"))
            .and_then(|_| writeln!(metadata, "fasta\t{}", fasta_path.display()))
            .and_then(|_| writeln!(metadata, "coordinate\tworker:u8,index:u32"))
            .and_then(|_| writeln!(metadata, "encoding\tcount-varint,delta-source-varints"))
            .and_then(|_| writeln!(metadata, "source_count\t{}", sources.len()))
            .map_err(|source| ColorError::Io {
                path: metadata_path.clone(),
                source,
            })?;
        for (index, source_path) in sources.iter().enumerate() {
            writeln!(metadata, "source\t{}\t{}", index + 1, source_path.display()).map_err(
                |source| ColorError::Io {
                    path: metadata_path.clone(),
                    source,
                },
            )?;
        }
        metadata.flush().map_err(|source| ColorError::Io {
            path: metadata_path,
            source,
        })
    }

    pub fn read_color(&self, coordinate: ColorCoordinate) -> Result<Vec<u32>, ColorError> {
        let worker = coordinate.worker();
        let target = coordinate.index();
        if worker >= self.records.len() || target >= self.records[worker] {
            return Err(ColorError::MalformedCoordinate(coordinate.as_u40()));
        }
        let path = color_worker_path(&self.dir, worker);
        let mut input = BufReader::with_capacity(
            COLOR_BUFFER_BYTES,
            File::open(&path).map_err(|source| ColorError::Io {
                path: path.clone(),
                source,
            })?,
        );
        for index in 0..=target {
            let sources = read_source_set(&mut input).map_err(|source| ColorError::Io {
                path: path.clone(),
                source,
            })?;
            if index == target {
                return Ok(sources);
            }
        }
        Err(ColorError::MalformedCoordinate(coordinate.as_u40()))
    }
}

fn color_worker_path(dir: &Path, worker: usize) -> PathBuf {
    dir.join(format!("{worker:03}.colors"))
}

fn write_source_set(output: &mut impl Write, sources: &[u32]) -> std::io::Result<()> {
    write_varint(output, sources.len() as u32)?;
    let mut previous = 0u32;
    for &source in sources {
        write_varint(output, source - previous)?;
        previous = source;
    }
    Ok(())
}

fn read_source_set(input: &mut impl Read) -> std::io::Result<Vec<u32>> {
    let len = read_varint(input)? as usize;
    let mut sources = Vec::with_capacity(len);
    let mut previous = 0u32;
    for _ in 0..len {
        previous = previous
            .checked_add(read_varint(input)?)
            .ok_or_else(|| std::io::Error::from(std::io::ErrorKind::InvalidData))?;
        sources.push(previous);
    }
    Ok(sources)
}

fn write_varint(output: &mut impl Write, mut value: u32) -> std::io::Result<()> {
    while value >= 0x80 {
        output.write_all(&[((value as u8) & 0x7f) | 0x80])?;
        value >>= 7;
    }
    output.write_all(&[value as u8])
}

fn append_varint_u32(output: &mut Vec<u8>, mut value: u32) {
    while value >= 0x80 {
        output.push(((value as u8) & 0x7f) | 0x80);
        value >>= 7;
    }
    output.push(value as u8);
}

fn read_varint(input: &mut impl Read) -> std::io::Result<u32> {
    let mut value = 0u32;
    for shift in (0..35).step_by(7) {
        let mut byte = [0u8; 1];
        input.read_exact(&mut byte)?;
        value |= u32::from(byte[0] & 0x7f) << shift;
        if byte[0] & 0x80 == 0 {
            return Ok(value);
        }
    }
    Err(std::io::Error::from(std::io::ErrorKind::InvalidData))
}

#[derive(Debug)]
pub enum ColorError {
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    InvalidWorkerCount(usize),
    MalformedSourceSet,
    MalformedCoordinate(u64),
    TooManyColors,
    TooManyColorRuns,
    MalformedUnitigIndex(u64),
    PoisonedWriter,
}

impl std::fmt::Display for ColorError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io { path, source } => write!(f, "{}: {source}", path.display()),
            Self::InvalidWorkerCount(count) => write!(f, "invalid color worker count: {count}"),
            Self::MalformedSourceSet => write!(f, "color source set is empty or not sorted"),
            Self::MalformedCoordinate(coord) => write!(f, "invalid color coordinate: {coord}"),
            Self::TooManyColors => write!(f, "worker color repository exceeds 2^32 entries"),
            Self::TooManyColorRuns => write!(f, "a local unitig exceeds 2^32 color runs"),
            Self::MalformedUnitigIndex(unitig) => {
                write!(f, "invalid local-unitig color index: {unitig}")
            }
            Self::PoisonedWriter => write!(f, "color repository writer lock is poisoned"),
        }
    }
}

impl std::error::Error for ColorError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repository_deduplicates_hashes_and_round_trips_sparse_sets() {
        let dir = std::env::temp_dir().join(format!(
            "cf3-color-repo-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        let repository = ConcurrentColorRepository::create(&dir, 2, 8).unwrap();
        let first = repository.resolve_or_insert(17, &[1, 2, 150], 0).unwrap();
        let duplicate = repository.resolve_or_insert(17, &[1, 2, 150], 1).unwrap();
        assert_eq!(first, duplicate);
        let second = repository.resolve_or_insert(23, &[3, 1000], 1).unwrap();
        let manifest = repository.finish().unwrap();
        assert_eq!(manifest.read_color(first).unwrap(), [1, 2, 150]);
        assert_eq!(manifest.read_color(second).unwrap(), [3, 1000]);
        manifest
            .write_metadata(31, Path::new("graph.fa"), &[PathBuf::from("source.fa")])
            .unwrap();
        let metadata = fs::read_to_string(dir.join("metadata.tsv")).unwrap();
        assert!(metadata.contains("format\tcf3rs-color-repository-v1"));
        assert!(metadata.contains("source\t1\tsource.fa"));
        let manifest_text = fs::read_to_string(dir.join("manifest.tsv")).unwrap();
        assert!(manifest_text.contains("000.colors"));
        assert!(!manifest_text.contains(&dir.display().to_string()));
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn color_runs_reverse_append_and_rotate_like_vertex_colors() {
        let a = ColorCoordinate::discovered(0, 1);
        let b = ColorCoordinate::discovered(0, 2);
        let c = ColorCoordinate::discovered(0, 3);
        let runs = vec![UnitigColor::new(0, a), UnitigColor::new(2, b)];
        let reversed = reverse_color_runs(&runs, 5);
        assert_eq!(
            reversed
                .iter()
                .map(|run| (run.offset(), run.coordinate()))
                .collect::<Vec<_>>(),
            [(0, b.as_u40()), (3, a.as_u40())]
        );

        let mut joined = vec![UnitigColor::new(0, c), UnitigColor::new(3, a)];
        append_color_runs(&mut joined, 4, &runs, 5, false);
        assert_eq!(
            joined
                .iter()
                .map(|run| (run.offset(), run.coordinate()))
                .collect::<Vec<_>>(),
            [(0, c.as_u40()), (3, a.as_u40()), (5, b.as_u40())]
        );

        let rotated = rotate_cycle_color_runs(&runs, 5, 2, false);
        assert_eq!(
            rotated
                .iter()
                .map(|run| (run.offset(), run.coordinate()))
                .collect::<Vec<_>>(),
            [(0, b.as_u40()), (3, a.as_u40())]
        );
    }

    #[test]
    fn color_run_sidecar_supports_range_streams() {
        let dir = std::env::temp_dir().join(format!(
            "cf3-color-runs-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        fs::create_dir_all(&dir).unwrap();
        let mut writer = ColorRunSidecarWriter::create(dir.join("local")).unwrap();
        let coordinate = ColorCoordinate::discovered(1, 9);
        writer.write_unitig(&[]).unwrap();
        let second_offset = writer.position();
        writer
            .write_unitig(&[
                UnitigColor::new(0, coordinate),
                UnitigColor::new(7, coordinate),
            ])
            .unwrap();
        let sidecar = writer.finish().unwrap();
        assert!(sidecar.read_unitig(0).unwrap().is_empty());
        assert_eq!(
            sidecar
                .reader_at(second_offset)
                .unwrap()
                .read_next()
                .unwrap()
                .len(),
            2
        );
        assert_eq!(
            sidecar
                .read_unitig(1)
                .unwrap()
                .iter()
                .map(|color| color.raw())
                .collect::<Vec<_>>(),
            [
                UnitigColor::new(0, coordinate).raw(),
                UnitigColor::new(7, coordinate).raw()
            ]
        );
        fs::remove_dir_all(dir).unwrap();
    }
}

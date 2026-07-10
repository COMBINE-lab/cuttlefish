use crate::hash::FastBuildHasher;
use crate::state::ColorCoordinate;
use scc::{HashMap as SccHashMap, hash_map::Entry as SccEntry};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

const COLOR_BUFFER_BYTES: usize = 128 * 1024;

pub struct ConcurrentColorRepository {
    dir: PathBuf,
    table: SccHashMap<u64, ColorCoordinate, FastBuildHasher>,
    workers: Vec<Mutex<ColorWorkerWriter>>,
}

struct ColorWorkerWriter {
    records: u32,
    output: BufWriter<File>,
}

impl ConcurrentColorRepository {
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
            table: SccHashMap::with_capacity_and_hasher(
                expected_colors.max(8),
                FastBuildHasher::default(),
            ),
            workers: worker_writers,
        })
    }

    pub fn resolve_or_insert(
        &self,
        color_hash: u64,
        sources: &[u32],
        worker: usize,
    ) -> Result<ColorCoordinate, ColorError> {
        if sources.is_empty() || !sources.windows(2).all(|pair| pair[0] < pair[1]) {
            return Err(ColorError::MalformedSourceSet);
        }
        if worker >= self.workers.len() {
            return Err(ColorError::InvalidWorkerCount(worker + 1));
        }
        match self.table.entry_sync(color_hash) {
            SccEntry::Occupied(entry) => Ok(*entry.get()),
            SccEntry::Vacant(entry) => {
                let mut writer = self.workers[worker]
                    .lock()
                    .map_err(|_| ColorError::PoisonedWriter)?;
                let index = writer.records;
                write_source_set(&mut writer.output, sources).map_err(|source| ColorError::Io {
                    path: color_worker_path(&self.dir, worker),
                    source,
                })?;
                writer.records = writer
                    .records
                    .checked_add(1)
                    .ok_or(ColorError::TooManyColors)?;
                let coordinate = ColorCoordinate::discovered(worker as u64, index as u64);
                entry.insert_entry(coordinate);
                Ok(coordinate)
            }
        }
    }

    pub fn finish(&self) -> Result<ColorRepositoryManifest, ColorError> {
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
            writeln!(
                output,
                "{worker}\t{records}\t{}",
                color_worker_path(&self.dir, worker).display()
            )
            .map_err(|source| ColorError::Io {
                path: path.clone(),
                source,
            })?;
        }
        output
            .flush()
            .map_err(|source| ColorError::Io { path, source })
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
        fs::remove_dir_all(dir).unwrap();
    }
}

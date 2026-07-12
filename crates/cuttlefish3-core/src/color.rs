use crate::hash::FastBuildHasher;
use crate::state::{ColorCoordinate, UnitigColor};
use scc::{HashMap as SccHashMap, hash_map::Entry as SccEntry};
use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Seek, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

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

pub fn append_color_runs(
    output: &mut Vec<UnitigColor>,
    output_vertex_count: u32,
    runs: &[UnitigColor],
    unitig_vertex_count: u32,
    reverse: bool,
) {
    let oriented = if reverse {
        reverse_color_runs(runs, unitig_vertex_count)
    } else {
        runs.to_vec()
    };
    let mut start = 0usize;
    if output
        .last()
        .zip(oriented.first())
        .is_some_and(|(left, right)| left.coordinate() == right.coordinate())
    {
        start = 1;
    }
    for run in &oriented[start..] {
        output.push(UnitigColor::new(
            output_vertex_count + run.offset() - 1,
            ColorCoordinate::from_u40(run.coordinate()),
        ));
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
        self.runs
            .write_all(&count.to_le_bytes())
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
        self.run_bytes += 4 + u64::from(count) * 8;
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
        let mut count = [0u8; 4];
        self.input
            .read_exact(&mut count)
            .map_err(|source| ColorError::Io {
                path: self.path.clone(),
                source,
            })?;
        let count = u32::from_le_bytes(count);
        let mut colors = Vec::with_capacity(count as usize);
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
        Ok(colors)
    }
}

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

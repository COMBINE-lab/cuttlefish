//! `cuttlefish colors`: read the color repository a colored build leaves
//! behind.
//!
//! A colored graph is written as two artifacts: the unitig FASTA, whose header
//! carries one packed run per color change, and the repository directory,
//! which holds the deduplicated source sets those runs point at. Neither is
//! meant to be read by hand -- the runs are packed `u64`s and the sets are
//! Elias-delta bit streams -- so these subcommands are the supported way to
//! get at them.
//!
//! Everything here streams, and the access pattern is chosen for the scale
//! these graphs actually reach: a 150000-genome colored build holds 87 million
//! distinct color sets across a 41 GB repository. `dump` emits a line per
//! color run per unitig -- more output than the graph itself -- through an
//! optional gzip encoder, resolving coordinates through a bounded cache since
//! consecutive runs repeat them heavily. `grep` never resolves per run at all:
//! it decides the predicate once per set in one sequential pass and carries
//! the answer in a bitset. `sets` is a pure stream and never builds an index.

use cuttlefish_rs::color::ColorRepositoryReader;
use cuttlefish_rs::state::ColorCoordinate;
use flate2::Compression;
use flate2::write::GzEncoder;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

type Result<T> = std::result::Result<T, Box<dyn Error + Send + Sync>>;

/// How many decoded sets to keep resident.
///
/// Sized so the cache stays small next to the graph even when every entry is a
/// wide set: at 64 Ki entries and a few hundred sources each it is tens of
/// megabytes, against the tens of gigabytes a dump can emit.
const SET_CACHE_ENTRIES: usize = 64 * 1024;

enum Mode {
    /// One line per color run per unitig.
    Dump,
    /// One line per distinct set in the repository.
    Sets,
    /// Unitigs whose colors satisfy a source predicate.
    Grep,
}

struct Params {
    mode: Mode,
    repository: PathBuf,
    fasta: Option<PathBuf>,
    output: Option<PathBuf>,
    gzip: Option<bool>,
    names: bool,
    any_of: Vec<u32>,
    all_of: Vec<u32>,
    none_of: Vec<u32>,
    count_only: bool,
}

/// Entry point for the `colors` subcommand. `args` is the argument list with
/// the subcommand name already consumed.
pub fn run<I>(args: I) -> Result<i32>
where
    I: Iterator<Item = String>,
{
    let Some(params) = parse_args(args)? else {
        return Ok(0); // --help
    };
    let mut reader = ColorRepositoryReader::open(&params.repository)?;
    let mut output = open_output(params.output.as_deref(), params.gzip)?;
    let result = match params.mode {
        Mode::Dump => dump(&mut reader, &params, &mut output),
        Mode::Sets => sets(&reader, &params, &mut output),
        Mode::Grep => grep(&mut reader, &params, &mut output),
    };
    // `colors sets | head` is how anyone looks at a repository for the first
    // time, and a reader that stops early is not an error.
    if is_broken_pipe(&result) {
        return Ok(0);
    }
    result?;
    let flushed = output.flush().map_err(Box::from).and_then(|()| {
        output.finish()?;
        Ok(())
    });
    if is_broken_pipe(&flushed) {
        return Ok(0);
    }
    flushed?;
    Ok(0)
}

fn is_broken_pipe(result: &Result<()>) -> bool {
    let Err(error) = result else {
        return false;
    };
    let mut source: Option<&(dyn Error + 'static)> = Some(error.as_ref());
    while let Some(current) = source {
        if current
            .downcast_ref::<std::io::Error>()
            .is_some_and(|io| io.kind() == std::io::ErrorKind::BrokenPipe)
        {
            return true;
        }
        source = current.source();
    }
    false
}

/// One line per color run: unitig, the vertex range it covers, and its sources.
fn dump(reader: &mut ColorRepositoryReader, params: &Params, out: &mut Output) -> Result<()> {
    let fasta = params.fasta.as_ref().ok_or("-i/--input is required")?;
    let mut cache = SetCache::default();
    let mut line = String::new();
    writeln!(out, "#unitig\tvertex_start\tvertex_end\tsources")?;
    let mut unitig = 0u64;
    let mut input = open_fasta(fasta)?;
    let mut label = String::new();
    while read_record(&mut input, &mut line, &mut label)? {
        let runs = parse_runs(&line)?;
        let vertices = vertex_count(label.len(), reader.k())?;
        for (index, &(offset, coordinate)) in runs.iter().enumerate() {
            let end = runs.get(index + 1).map_or(vertices, |next| next.0);
            let sources = cache.get(reader, coordinate)?;
            write!(out, "{unitig}\t{offset}\t{end}\t")?;
            write_sources(out, sources, reader, params.names)?;
            out.write_all(b"\n")?;
        }
        unitig += 1;
    }
    Ok(())
}

/// One line per distinct set, in repository order.
fn sets(reader: &ColorRepositoryReader, params: &Params, out: &mut Output) -> Result<()> {
    writeln!(out, "#coordinate\tworker\tindex\tsize\tsources")?;
    let mut error = None;
    for worker in 0..reader.records_per_worker().len() {
        reader.for_each_record(worker, |coordinate, sources| {
            if error.is_some() {
                return Ok(());
            }
            let mut emit = || -> Result<()> {
                write!(
                    out,
                    "{}\t{}\t{}\t{}\t",
                    coordinate.as_u40(),
                    coordinate.worker(),
                    coordinate.index(),
                    sources.len()
                )?;
                write_sources(out, sources, reader, params.names)?;
                out.write_all(b"\n")?;
                Ok(())
            };
            if let Err(failed) = emit() {
                error = Some(failed);
            }
            Ok(())
        })?;
        if let Some(failed) = error.take() {
            return Err(failed);
        }
    }
    Ok(())
}

/// Unitigs with at least one run whose set satisfies the source predicate.
///
/// The predicate is decided once per *set*, not once per run, by streaming the
/// repository and recording the answer in a bitset over coordinates. That is
/// what keeps this linear on real data: a 150000-genome build has 87 million
/// distinct sets in 41 GB, which is one sequential pass and an 11 MB bitset,
/// against 41 GB of random reads if each run were resolved on demand. The
/// FASTA pass that follows then never touches the repository at all --
/// except, when the matched sets are to be printed, for the ones that matched.
fn grep(reader: &mut ColorRepositoryReader, params: &Params, out: &mut Output) -> Result<()> {
    let fasta = params.fasta.as_ref().ok_or("-i/--input is required")?;
    if params.any_of.is_empty() && params.all_of.is_empty() && params.none_of.is_empty() {
        return Err("grep needs at least one of --any-of, --all-of, --none-of".into());
    }
    let selected = select_coordinates(reader, params)?;
    let mut cache = SetCache::default();
    let mut line = String::new();
    let mut label = String::new();
    let mut input = open_fasta(fasta)?;
    let mut unitig = 0u64;
    let mut matched = 0u64;
    if !params.count_only {
        writeln!(out, "#unitig\tvertex_start\tvertex_end\tlabel\tsources")?;
    }
    while read_record(&mut input, &mut line, &mut label)? {
        let runs = parse_runs(&line)?;
        let vertices = vertex_count(label.len(), reader.k())?;
        let mut hit = false;
        for (index, &(offset, coordinate)) in runs.iter().enumerate() {
            let end = runs.get(index + 1).map_or(vertices, |next| next.0);
            if !selected.contains(coordinate) {
                continue;
            }
            hit = true;
            if params.count_only {
                break;
            }
            // The run's own substring, not the whole unitig: a colored unitig
            // can be megabases long and the run is what actually matched.
            let start = offset as usize;
            let stop = (end as usize + reader.k() as usize - 1).min(label.len());
            write!(out, "{unitig}\t{offset}\t{end}\t")?;
            out.write_all(label.as_bytes().get(start..stop).unwrap_or_default())?;
            out.write_all(b"\t")?;
            let sources = cache.get(reader, coordinate)?;
            write_sources(out, sources, reader, params.names)?;
            out.write_all(b"\n")?;
        }
        matched += u64::from(hit);
        unitig += 1;
    }
    if params.count_only {
        writeln!(out, "{matched}")?;
    } else {
        eprintln!("cuttlefish colors: {matched} of {unitig} unitig(s) matched");
    }
    Ok(())
}

/// One bit per repository record: does this set satisfy the predicate?
struct CoordinateSet {
    words: Vec<Vec<u64>>,
}

impl CoordinateSet {
    fn new(records_per_worker: &[u32]) -> Self {
        Self {
            words: records_per_worker
                .iter()
                .map(|&count| vec![0u64; (count as usize).div_ceil(64)])
                .collect(),
        }
    }

    fn insert(&mut self, coordinate: ColorCoordinate) {
        let index = coordinate.index() as usize;
        self.words[coordinate.worker()][index / 64] |= 1u64 << (index % 64);
    }

    #[inline]
    fn contains(&self, coordinate: u64) -> bool {
        let coordinate = ColorCoordinate::from_u40(coordinate);
        let index = coordinate.index() as usize;
        self.words
            .get(coordinate.worker())
            .and_then(|words| words.get(index / 64))
            .is_some_and(|word| word & (1u64 << (index % 64)) != 0)
    }
}

/// Streams every set once, marking the coordinates the predicate accepts.
fn select_coordinates(reader: &ColorRepositoryReader, params: &Params) -> Result<CoordinateSet> {
    let mut selected = CoordinateSet::new(reader.records_per_worker());
    for worker in 0..reader.records_per_worker().len() {
        reader.for_each_record(worker, |coordinate, sources| {
            if matches(sources, params) {
                selected.insert(coordinate);
            }
            Ok(())
        })?;
    }
    Ok(selected)
}

#[inline]
fn matches(sources: &[u32], params: &Params) -> bool {
    if !params.any_of.is_empty() && !params.any_of.iter().any(|s| sources.contains(s)) {
        return false;
    }
    if !params.all_of.is_empty() && !params.all_of.iter().all(|s| sources.contains(s)) {
        return false;
    }
    if params.none_of.iter().any(|s| sources.contains(s)) {
        return false;
    }
    true
}

fn write_sources(
    out: &mut Output,
    sources: &[u32],
    reader: &ColorRepositoryReader,
    names: bool,
) -> Result<()> {
    for (index, &source) in sources.iter().enumerate() {
        if index > 0 {
            out.write_all(b",")?;
        }
        if names {
            // Source ids are 1-based into the metadata's source list.
            match reader.sources().get(source as usize - 1) {
                Some(path) => write!(out, "{}", path.display())?,
                None => write!(out, "?{source}")?,
            }
        } else {
            write!(out, "{source}")?;
        }
    }
    Ok(())
}

/// Bounded map from coordinate to decoded set.
///
/// Runs repeat coordinates heavily -- a set shared by many unitigs is the
/// common case in a pangenome -- so this turns most lookups into a hash hit
/// rather than a seek and a bit-stream decode. It is cleared wholesale when
/// full rather than evicting by age: the access pattern has no reuse distance
/// worth tracking, and a clear costs nothing next to the decodes it saves.
#[derive(Default)]
struct SetCache {
    entries: HashMap<u64, Vec<u32>>,
}

impl SetCache {
    fn get(&mut self, reader: &mut ColorRepositoryReader, coordinate: u64) -> Result<&[u32]> {
        if !self.entries.contains_key(&coordinate) {
            if self.entries.len() >= SET_CACHE_ENTRIES {
                self.entries.clear();
            }
            let sources = reader.read_color(ColorCoordinate::from_u40(coordinate))?;
            self.entries.insert(coordinate, sources);
        }
        Ok(self
            .entries
            .get(&coordinate)
            .expect("entry was just inserted"))
    }
}

/// Reads one FASTA record, returning false at end of input.
///
/// Labels are single-line in cuttlefish output, which is what lets this hand
/// back a borrowed label instead of assembling one.
fn read_record(
    input: &mut BufReader<File>,
    header: &mut String,
    label: &mut String,
) -> Result<bool> {
    header.clear();
    label.clear();
    if input.read_line(header)? == 0 {
        return Ok(false);
    }
    if !header.starts_with('>') {
        return Err(format!("expected a FASTA header, found {header:?}").into());
    }
    if input.read_line(label)? == 0 {
        return Err("FASTA header with no sequence".into());
    }
    while label.ends_with('\n') || label.ends_with('\r') {
        label.pop();
    }
    Ok(true)
}

/// Splits a colored header into its `(vertex offset, coordinate)` runs.
fn parse_runs(header: &str) -> Result<Vec<(u32, u64)>> {
    let mut runs = Vec::new();
    for field in header.split_whitespace().skip(1) {
        let packed: u64 = field
            .parse()
            .map_err(|_| format!("colored header field {field:?} is not a packed color run"))?;
        // Low 24 bits are the vertex offset; the rest is the coordinate.
        runs.push(((packed & 0xFF_FFFF) as u32, packed >> 24));
    }
    if runs.is_empty() {
        return Err(
            "no color runs in FASTA header; is this the FASTA from a --color build?".into(),
        );
    }
    Ok(runs)
}

fn vertex_count(label_len: usize, k: u16) -> Result<u32> {
    let k = k as usize;
    if label_len < k {
        return Err(format!("unitig shorter than k = {k}").into());
    }
    Ok((label_len - k + 1) as u32)
}

fn open_fasta(path: &Path) -> Result<BufReader<File>> {
    Ok(BufReader::with_capacity(1 << 20, File::open(path)?))
}

/// Output sink: stdout or a file, optionally gzipped.
enum Output {
    Plain(Box<dyn Write + Send>),
    Gzip(Box<GzEncoder<Box<dyn Write + Send>>>),
}

impl Output {
    fn finish(self) -> Result<()> {
        match self {
            Self::Plain(mut out) => out.flush()?,
            Self::Gzip(encoder) => {
                encoder.finish()?.flush()?;
            }
        }
        Ok(())
    }
}

impl Write for Output {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Plain(out) => out.write(buf),
            Self::Gzip(out) => out.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Plain(out) => out.flush(),
            Self::Gzip(out) => out.flush(),
        }
    }
}

/// Opens the sink, compressing when asked or when the name says `.gz`.
fn open_output(path: Option<&Path>, gzip: Option<bool>) -> Result<Output> {
    let named_gz = path.is_some_and(|path| {
        path.extension()
            .is_some_and(|extension| extension.eq_ignore_ascii_case("gz"))
    });
    let compress = gzip.unwrap_or(named_gz);
    let sink: Box<dyn Write + Send> = match path {
        Some(path) => Box::new(BufWriter::with_capacity(1 << 20, File::create(path)?)),
        None => Box::new(BufWriter::with_capacity(1 << 20, std::io::stdout())),
    };
    Ok(if compress {
        // Level 1: these are wide, highly repetitive TSV lines, where the
        // fast level already gets most of the ratio and the slow one would
        // make compression, not decoding, the bottleneck.
        Output::Gzip(Box::new(GzEncoder::new(sink, Compression::new(1))))
    } else {
        Output::Plain(sink)
    })
}

fn parse_args<I>(args: I) -> Result<Option<Params>>
where
    I: Iterator<Item = String>,
{
    let mut args = args;
    let mode = match args.next().as_deref() {
        Some("dump") => Mode::Dump,
        Some("sets") => Mode::Sets,
        Some("grep") => Mode::Grep,
        Some("-h") | Some("--help") | None => {
            print_help();
            return Ok(None);
        }
        Some(other) => return Err(format!("unknown colors mode {other}").into()),
    };
    let mut repository = None;
    let mut fasta = None;
    let mut output = None;
    let mut gzip = None;
    let mut names = false;
    let mut any_of = Vec::new();
    let mut all_of = Vec::new();
    let mut none_of = Vec::new();
    let mut count_only = false;
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-r" | "--repository" => {
                repository = Some(PathBuf::from(argument_value(&mut args, &arg)?))
            }
            "-i" | "--input" => fasta = Some(PathBuf::from(argument_value(&mut args, &arg)?)),
            "-o" | "--output" => output = Some(PathBuf::from(argument_value(&mut args, &arg)?)),
            "-z" | "--gzip" => gzip = Some(true),
            "--no-gzip" => gzip = Some(false),
            "--names" => names = true,
            "-c" | "--count" => count_only = true,
            "--any-of" => any_of = parse_sources(&argument_value(&mut args, &arg)?)?,
            "--all-of" => all_of = parse_sources(&argument_value(&mut args, &arg)?)?,
            "--none-of" => none_of = parse_sources(&argument_value(&mut args, &arg)?)?,
            "-h" | "--help" => {
                print_help();
                return Ok(None);
            }
            _ => return Err(format!("unknown argument {arg}").into()),
        }
    }
    Ok(Some(Params {
        mode,
        repository: repository.ok_or("-r/--repository is required")?,
        fasta,
        output,
        gzip,
        names,
        any_of,
        all_of,
        none_of,
        count_only,
    }))
}

fn parse_sources(value: &str) -> Result<Vec<u32>> {
    let mut sources = Vec::new();
    for field in value.split(',').filter(|field| !field.is_empty()) {
        let source: u32 = field
            .parse()
            .map_err(|_| format!("source id {field:?} is not a number"))?;
        if source == 0 {
            return Err("source ids are 1-based; 0 is not a source".into());
        }
        sources.push(source);
    }
    if sources.is_empty() {
        return Err("empty source list".into());
    }
    Ok(sources)
}

fn argument_value<I>(args: &mut I, flag: &str) -> Result<String>
where
    I: Iterator<Item = String>,
{
    args.next()
        .ok_or_else(|| format!("{flag} needs a value").into())
}

fn print_help() {
    println!(
        "Read the color repository written by a colored build.

Usage:
  cuttlefish colors dump -r DIR -i FASTA [OPTION...]
  cuttlefish colors sets -r DIR [OPTION...]
  cuttlefish colors grep -r DIR -i FASTA (--any-of|--all-of|--none-of LIST) [OPTION...]

 modes:
  dump                   every color run of every unitig, as
                         unitig, vertex range, sources
  sets                   every distinct source set in the repository
  grep                   unitigs with a run matching a source predicate

 common options:
  -r, --repository <arg> the .cf3rs.color-repository directory
  -i, --input <arg>      the colored unitig FASTA (dump, grep)
  -o, --output <arg>     write here instead of stdout
  -z, --gzip             gzip the output (implied by an .gz output name)
      --no-gzip          never gzip, whatever the output is named
      --names            print source paths instead of 1-based ids

 grep options:
      --any-of <arg>     comma-separated source ids; keep runs containing any
      --all-of <arg>     keep runs containing all of these
      --none-of <arg>    drop runs containing any of these
  -c, --count            print only how many unitigs matched

Output is a TSV with a `#` header line, written as a stream: a dump of a real
corpus is larger than the graph, so send it somewhere with room, or gzip it."
    );
}

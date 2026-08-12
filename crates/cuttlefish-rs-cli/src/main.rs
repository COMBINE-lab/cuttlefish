use cuttlefish_rs::colored::{ColoredBuildError, build_colored_from_buckets};
use cuttlefish_rs::discontinuity::{raise_open_file_limit, report_process_memory};
use cuttlefish_rs::input::InputError;
use cuttlefish_rs::params::{BuildParams, ParamError};
use cuttlefish_rs::partition::{
    PartitionEmissionStats, PartitionRunError, emit_weak_superkmer_buckets,
};
use cuttlefish_rs::uncolored::{UncoloredBuildError, build_uncolored_from_buckets};
use cuttlefish_rs::{
    DEFAULT_SUBGRAPH_COUNT, GraphInput, MAX_K, configure_global_parallelism, default_threads,
    default_work_dir,
};
use std::time::Instant;

mod colors;
mod compare;

#[cfg(all(feature = "jemalloc", feature = "mimalloc"))]
compile_error!("allocator features `jemalloc` and `mimalloc` are mutually exclusive");

#[cfg(feature = "jemalloc")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() {
    std::process::exit(match run(std::env::args().skip(1)) {
        Ok(code) => code,
        Err(err) => {
            eprintln!("{err}");
            1
        }
    });
}

fn run<I>(mut args: I) -> Result<i32, CliError>
where
    I: Iterator<Item = String>,
{
    let Some(command) = args.next() else {
        print_top_help();
        return Ok(0);
    };

    match command.to_ascii_lowercase().as_str() {
        "build" => {
            let params = parse_build(args)?;
            params.validate()?;
            // The compact k <= 31 path uses phase-local bounded pools. Avoid
            // keeping an otherwise idle global Rayon pool contending with
            // those workers for the entire build.
            configure_global_parallelism(if params.k <= 31 { 1 } else { params.threads })
                .map_err(|error| CliError::Parallelism(error.to_string()))?;
            // Nobody should have to discover `ulimit -n` to build a graph.
            // Raising the soft limit to the hard limit needs no privileges,
            // and the phases that fan out adapt to whatever budget they see,
            // so a generous one simply lets them stay wide.
            let (fd_before, fd_after) = raise_open_file_limit();
            if fd_after > fd_before {
                eprintln!("cuttlefish: raised open-file limit from {fd_before} to {fd_after}");
            }
            eprintln!(
                "cuttlefish: parsed build request for k={}, l={}, cutoff={}, color={}",
                params.k,
                params.minimizer_len,
                params.cutoff(),
                params.color
            );
            let partition_start = Instant::now();
            report_process_memory("before partition");
            let emission = dispatch_partition(&params)?;
            let partition_elapsed = partition_start.elapsed();
            report_process_memory("after partition");
            let partition_stats = &emission.partition;
            eprintln!(
                "cuttlefish: resolved {} input file(s)",
                partition_stats.input_files
            );
            eprintln!(
                "cuttlefish: partitioned {} record(s), {} ACTG fragment(s), {} fragment base(s), {} weak super-kmer(s), {} non-empty subgraph(s), max subgraph bucket {}",
                partition_stats.records,
                partition_stats.fragments,
                partition_stats.fragment_bases,
                partition_stats.weak_superkmers,
                partition_stats.non_empty_graphs(),
                partition_stats.max_graph_superkmers()
            );
            eprintln!(
                "cuttlefish: wrote {} weak-superkmer bucket(s), {} byte(s), in {}",
                emission.buckets.bucket_files,
                emission.buckets.bytes_written,
                emission.buckets.bucket_dir.display()
            );
            eprintln!(
                "cuttlefish: partition detail: parse/send {:.3}s, worker scan+pack {:.3}s, bucket flush {} call(s) {:.3}s, finish {:.3}s",
                emission.parse_elapsed.as_secs_f64(),
                emission.worker_elapsed.as_secs_f64(),
                emission.bucket_flushes,
                emission.bucket_flush_elapsed.as_secs_f64(),
                emission.bucket_finish_elapsed.as_secs_f64()
            );
            eprintln!(
                "cuttlefish: partition and bucket emission completed in {:.3}s",
                partition_elapsed.as_secs_f64()
            );
            // Diagnostic: stop after partitioning so a profiler can attribute
            // work to that phase alone.
            if std::env::var_os("CF3_RS_EXIT_AFTER_PARTITION").is_some() {
                eprintln!("cuttlefish: exiting after partition (diagnostic)");
                std::process::exit(0);
            }
            let build_start = Instant::now();
            if params.color {
                report_process_memory("before colored build");
                let build = dispatch_colored_build(&params, &emission)?;
                report_process_memory("after colored build");
                eprintln!(
                    "cuttlefish: wrote {} colored unitig(s), {} base(s), FASTA at {}; color repository at {}",
                    build.unitigs,
                    build.unitig_bases,
                    build.output_path.display(),
                    build.color_repository.display(),
                );
            } else {
                report_process_memory("before uncolored build");
                let build = dispatch_uncolored_build(&params, &emission)?;
                report_process_memory("after uncolored build");
                eprintln!(
                    "cuttlefish: graph handoff recorded {} retained/active edge(s) from {} observed exit/edge event(s)",
                    build.retained_edges, build.observed_edges
                );
                eprintln!(
                    "cuttlefish: wrote {} unitig(s), {} base(s), FASTA at {}",
                    build.unitigs,
                    build.unitig_bases,
                    build.output_path.display()
                );
            }
            eprintln!(
                "cuttlefish: graph build completed in {:.3}s",
                build_start.elapsed().as_secs_f64()
            );
            // The build path has already flushed its output files. Exiting here
            // lets the OS reclaim large graph buffers instead of spending wall
            // time recursively dropping them on the hot benchmark path.
            std::process::exit(0);
        }
        "compare" => compare::run(args).map_err(|err| CliError::Compare(err.to_string())),
        "colors" => colors::run(args).map_err(|err| CliError::Colors(err.to_string())),
        "help" | "--help" | "-h" => {
            print_top_help();
            Ok(0)
        }
        "version" | "--version" | "-V" => {
            println!("cuttlefish {}", env!("CARGO_PKG_VERSION"));
            Ok(0)
        }
        _ => {
            print_top_help();
            Ok(1)
        }
    }
}

fn dispatch_partition(params: &BuildParams) -> Result<PartitionEmissionStats, CliError> {
    macro_rules! dispatch {
        ($($k:literal),* $(,)?) => {
            match params.k {
                $($k => emit_weak_superkmer_buckets::<$k>(params, DEFAULT_SUBGRAPH_COUNT),)*
                other => Err(PartitionRunError::Partition(
                    cuttlefish_rs::partition::PartitionError::InvalidK(other as usize),
                )),
            }
        };
    }

    Ok(dispatch!(
        3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49,
        51, 53, 55, 57, 59, 61, 63,
    )?)
}

fn dispatch_uncolored_build(
    params: &BuildParams,
    emission: &PartitionEmissionStats,
) -> Result<cuttlefish_rs::uncolored::UncoloredBuildStats, CliError> {
    macro_rules! dispatch {
        ($($k:literal),* $(,)?) => {
            match params.k {
                $($k => build_uncolored_from_buckets::<$k>(params, &emission.buckets.bucket_dir),)*
                other => Err(UncoloredBuildError::Kmer(
                    cuttlefish_rs::kmer::KmerError::UnsupportedK(other as usize),
                )),
            }
        };
    }

    Ok(dispatch!(
        3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49,
        51, 53, 55, 57, 59, 61, 63,
    )?)
}

fn dispatch_colored_build(
    params: &BuildParams,
    emission: &PartitionEmissionStats,
) -> Result<cuttlefish_rs::colored::ColoredBuildStats, CliError> {
    macro_rules! dispatch {
        ($($k:literal),* $(,)?) => {
            match params.k {
                $($k => build_colored_from_buckets::<$k>(params, &emission.buckets.bucket_dir),)*
                other => Err(ColoredBuildError::Kmer(
                    cuttlefish_rs::kmer::KmerError::UnsupportedK(other as usize),
                )),
            }
        };
    }
    Ok(dispatch!(
        3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49,
        51, 53, 55, 57, 59, 61, 63,
    )?)
}

fn parse_build<I>(args: I) -> Result<BuildParams, CliError>
where
    I: Iterator<Item = String>,
{
    let mut raw = RawBuild::default();
    let mut args = args.peekable();
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                print_build_help();
                return Err(CliError::Help);
            }
            "-s" | "--seq" => split_values(take_value("--seq", &mut args)?, &mut raw.seqs),
            "-l" | "--list" => split_values(take_value("--list", &mut args)?, &mut raw.lists),
            "-d" | "--dir" => split_values(take_value("--dir", &mut args)?, &mut raw.dirs),
            "-k" | "--kmer-len" => raw.k = Some(parse_value("--kmer-len", &mut args)?),
            "--min-len" => raw.minimizer_len = Some(parse_value("--min-len", &mut args)?),
            "-c" | "--cutoff" => raw.cutoff = Some(parse_value("--cutoff", &mut args)?),
            "-t" | "--threads" => raw.threads = Some(parse_value("--threads", &mut args)?),
            "-m" | "--max-memory" => {
                raw.max_memory_gb = Some(parse_value("--max-memory", &mut args)?)
            }
            "-o" | "--output" => raw.output = Some(take_value("--output", &mut args)?),
            "-w" | "--work-dir" => raw.work_dir = Some(take_value("--work-dir", &mut args)?),
            "--read" => raw.read = true,
            "--ref" => raw.reference = true,
            "--color" => raw.color = true,
            "--compress-buckets" => raw.compress_buckets = Some(true),
            "--no-compress-buckets" => raw.compress_buckets = Some(false),
            "--skip-unreadable" => raw.skip_unreadable = true,
            _ if arg.starts_with("--seq=") => split_values(arg[6..].to_string(), &mut raw.seqs),
            _ if arg.starts_with("--list=") => split_values(arg[7..].to_string(), &mut raw.lists),
            _ if arg.starts_with("--dir=") => split_values(arg[6..].to_string(), &mut raw.dirs),
            _ if arg.starts_with("--kmer-len=") => raw.k = Some(parse_inline(&arg, "--kmer-len=")?),
            _ if arg.starts_with("--min-len=") => {
                raw.minimizer_len = Some(parse_inline(&arg, "--min-len=")?)
            }
            _ if arg.starts_with("--cutoff=") => {
                raw.cutoff = Some(parse_inline(&arg, "--cutoff=")?)
            }
            _ if arg.starts_with("--threads=") => {
                raw.threads = Some(parse_inline(&arg, "--threads=")?)
            }
            _ if arg.starts_with("--max-memory=") => {
                raw.max_memory_gb = Some(parse_inline(&arg, "--max-memory=")?)
            }
            _ if arg.starts_with("--output=") => raw.output = Some(arg[9..].to_string()),
            _ if arg.starts_with("--work-dir=") => raw.work_dir = Some(arg[11..].to_string()),
            _ => return Err(CliError::UnknownArg(arg)),
        }
    }

    if raw.read == raw.reference {
        return Err(CliError::InputMode);
    }

    let output = raw.output.ok_or(CliError::MissingArg("--output"))?;
    let mut params = BuildParams::new(
        if raw.read {
            GraphInput::Reads
        } else {
            GraphInput::References
        },
        output,
    );

    params.seqs = raw.seqs;
    params.lists = raw.lists;
    params.dirs = raw.dirs;
    if let Some(k) = raw.k {
        params.k = k;
    }
    if let Some(l) = raw.minimizer_len {
        params.minimizer_len = l;
    }
    params.cutoff = raw.cutoff;
    params.color = raw.color;
    if let Some(compress) = raw.compress_buckets {
        params.compress_buckets = compress;
    }
    params.skip_unreadable = raw.skip_unreadable;
    params.max_memory_gb = raw.max_memory_gb;
    if let Some(threads) = raw.threads {
        params.threads = threads;
    }
    if let Some(w) = raw.work_dir {
        params.work_dir = w;
    }

    Ok(params)
}

#[derive(Default)]
struct RawBuild {
    seqs: Vec<String>,
    lists: Vec<String>,
    dirs: Vec<String>,
    k: Option<u16>,
    minimizer_len: Option<u16>,
    cutoff: Option<u32>,
    threads: Option<usize>,
    max_memory_gb: Option<usize>,
    output: Option<String>,
    work_dir: Option<String>,
    read: bool,
    reference: bool,
    color: bool,
    compress_buckets: Option<bool>,
    skip_unreadable: bool,
}

fn split_values(value: String, out: &mut Vec<String>) {
    out.extend(
        value
            .split(',')
            .filter(|s| !s.is_empty())
            .map(ToString::to_string),
    );
}

fn take_value<I>(flag: &'static str, args: &mut std::iter::Peekable<I>) -> Result<String, CliError>
where
    I: Iterator<Item = String>,
{
    args.next().ok_or(CliError::MissingArg(flag))
}

fn parse_value<T, I>(flag: &'static str, args: &mut std::iter::Peekable<I>) -> Result<T, CliError>
where
    T: std::str::FromStr,
    I: Iterator<Item = String>,
{
    take_value(flag, args)?
        .parse()
        .map_err(|_| CliError::InvalidValue(flag))
}

fn parse_inline<T>(arg: &str, prefix: &'static str) -> Result<T, CliError>
where
    T: std::str::FromStr,
{
    arg[prefix.len()..]
        .parse()
        .map_err(|_| CliError::InvalidValue(prefix))
}

fn print_top_help() {
    println!("cuttlefish {}", env!("CARGO_PKG_VERSION"));
    println!("Supported commands: `build`, `compare`, `colors`, `help`, `version`.");
    println!("Usage:");
    println!("    cuttlefish build [options]");
    println!("    cuttlefish compare [options]");
    println!("    cuttlefish colors dump|sets|grep [options]");
}

fn print_build_help() {
    let work_dir = default_work_dir();
    let threads = default_threads();
    println!(
        "Efficiently construct the (colored) compacted de Bruijn graph from reference sequences or sequencing reads."
    );
    println!("Usage:");
    println!("  cuttlefish build [OPTION...]");
    println!();
    println!(" common options:");
    println!("  -s, --seq <arg>       input files");
    println!("  -l, --list <arg>      input file lists");
    println!("  -d, --dir <arg>       input file directories");
    println!("  -k, --kmer-len <arg>  k-mer length (default: 31; odd <= {MAX_K})");
    println!("      --min-len <arg>   minimizer length (default: 12)");
    println!("  -o, --output <arg>    output file prefix");
    println!("  -w, --work-dir <arg>  working directory (default: {work_dir})");
    println!("  -t, --threads <arg>   worker threads for parallel phases (default: {threads})");
    println!("  -m, --max-memory <arg> soft maximum memory budget in GiB");
    println!("      --read            construct a compacted read de Bruijn graph");
    println!("      --ref             construct a compacted reference de Bruijn graph");
    println!("  -c, --cutoff <arg>    frequency cutoff for (k + 1)-mers");
    println!("      --color           color the compacted graph");
    println!("      --compress-buckets compress uncolored temporary buckets (default)");
    println!("      --no-compress-buckets store uncolored temporary buckets uncompressed");
    println!("      --skip-unreadable skip inputs that fail to parse instead of aborting");
    println!("  -h, --help            print usage");
}

#[derive(Debug)]
enum CliError {
    Help,
    UnknownArg(String),
    MissingArg(&'static str),
    InvalidValue(&'static str),
    InputMode,
    Parallelism(String),
    Compare(String),
    Colors(String),
    Param(ParamError),
    Input(InputError),
    Partition(PartitionRunError),
    Uncolored(UncoloredBuildError),
    Colored(ColoredBuildError),
}

impl From<ParamError> for CliError {
    fn from(value: ParamError) -> Self {
        Self::Param(value)
    }
}

impl From<InputError> for CliError {
    fn from(value: InputError) -> Self {
        Self::Input(value)
    }
}

impl From<PartitionRunError> for CliError {
    fn from(value: PartitionRunError) -> Self {
        Self::Partition(value)
    }
}

impl From<UncoloredBuildError> for CliError {
    fn from(value: UncoloredBuildError) -> Self {
        Self::Uncolored(value)
    }
}

impl From<ColoredBuildError> for CliError {
    fn from(value: ColoredBuildError) -> Self {
        Self::Colored(value)
    }
}

impl std::fmt::Display for CliError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Help => write!(f, "help requested"),
            Self::UnknownArg(arg) => write!(f, "unknown argument: {arg}"),
            Self::MissingArg(arg) => write!(f, "missing value for {arg}"),
            Self::InvalidValue(arg) => write!(f, "invalid value for {arg}"),
            Self::InputMode => write!(f, "select exactly one of --read or --ref"),
            Self::Parallelism(err) => write!(f, "failed to configure worker threads: {err}"),
            Self::Compare(err) => write!(f, "cuttlefish compare: {err}"),
            Self::Colors(err) => write!(f, "cuttlefish colors: {err}"),
            Self::Param(err) => write!(f, "{err}"),
            Self::Input(err) => write!(f, "{err}"),
            Self::Partition(err) => write!(f, "{err}"),
            Self::Uncolored(err) => write!(f, "{err}"),
            Self::Colored(err) => write!(f, "{err}"),
        }
    }
}

impl std::error::Error for CliError {}

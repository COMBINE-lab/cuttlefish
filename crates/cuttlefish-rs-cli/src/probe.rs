//! `cuttlefish probe`: what a build would learn about its work directory.
//!
//! `build --compress-intermediates auto` (the default) classifies the work
//! directory's storage and times a short write and read-back before deciding
//! whether to compress intermediate streams. This runs the same measurement
//! on its own and prints both the numbers and the decision, so a choice that
//! looks wrong can be checked, and a slow or shared scratch path can be
//! sized up before a long build.

use std::error::Error;
use std::path::PathBuf;
use std::time::Instant;

use cuttlefish_rs::intermediates;
use cuttlefish_rs::params::IntermediateCompression;
use cuttlefish_rs::{default_threads, default_work_dir};

type Result<T> = std::result::Result<T, Box<dyn Error + Send + Sync>>;

/// Entry point for the `probe` subcommand. `args` is the argument list with
/// the subcommand name already consumed.
pub fn run<I>(args: I) -> Result<i32>
where
    I: Iterator<Item = String>,
{
    let mut work_dir = PathBuf::from(default_work_dir());
    let mut threads = default_threads();
    let mut args = args.peekable();
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                print_help();
                return Ok(0);
            }
            "-w" | "--work-dir" => {
                work_dir = PathBuf::from(args.next().ok_or("missing value for --work-dir")?);
            }
            "-t" | "--threads" => {
                threads = args
                    .next()
                    .ok_or("missing value for --threads")?
                    .parse()
                    .map_err(|_| "invalid value for --threads")?;
            }
            _ if arg.starts_with("--work-dir=") => work_dir = PathBuf::from(&arg[11..]),
            _ if arg.starts_with("--threads=") => {
                threads = arg[10..]
                    .parse()
                    .map_err(|_| "invalid value for --threads")?;
            }
            _ => return Err(format!("unknown argument: {arg}").into()),
        }
    }
    if !work_dir.is_dir() {
        return Err(format!("{} is not a directory", work_dir.display()).into());
    }

    let storage = scratch_probe::classify(&work_dir);
    println!("work directory   {}", work_dir.display());
    println!(
        "filesystem       {} ({:?}{})",
        storage.filesystem,
        storage.kind,
        match storage.rotational {
            Some(true) => ", rotational",
            Some(false) => ", solid state",
            None => "",
        }
    );
    if let Some(bytes) = scratch_probe::mem_available() {
        println!("memory available {:.1} GB", bytes as f64 / 1e9);
    }
    let started = Instant::now();
    let choice = intermediates::choose(IntermediateCompression::Auto, &work_dir, threads);
    // `auto` skips the timing on network and rotational storage; time it
    // anyway, since that is what this command is for.
    let probe = match choice.probe.clone() {
        Some(probe) => Ok(probe),
        None => scratch_probe::probe(&work_dir, &scratch_probe::Limits::default()),
    };
    match probe {
        Ok(probe) => {
            let io = if probe.direct_io {
                "direct"
            } else {
                "buffered"
            };
            println!(
                "write            {:.2} GB/s ({:.0} MB in {:.3}s, {io}, with fdatasync)",
                probe.write.bytes_per_second() / 1e9,
                probe.write.bytes as f64 / 1e6,
                probe.write.elapsed.as_secs_f64(),
            );
            println!(
                "read             {:.2} GB/s ({:.0} MB in {:.3}s, {io})",
                probe.read.bytes_per_second() / 1e9,
                probe.read.bytes as f64 / 1e6,
                probe.read.elapsed.as_secs_f64(),
            );
        }
        Err(error) => println!("probe failed     {error}"),
    }
    println!("probe took       {:.2}s", started.elapsed().as_secs_f64());
    println!(
        "build produces   about {:.2} GB/s of intermediates at {threads} thread(s)",
        intermediates::intermediate_bytes_per_second(threads) / 1e9
    );
    println!(
        "auto chooses     compression {} ({})",
        if choice.compress { "on" } else { "off" },
        choice.reason
    );
    Ok(0)
}

fn print_help() {
    println!("Measure a work directory's storage and show whether `build` would compress");
    println!("intermediate streams there (`--compress-intermediates auto`).");
    println!("Usage:");
    println!("  cuttlefish probe [OPTION...]");
    println!();
    println!(
        "  -w, --work-dir <arg>  directory to measure (default: {})",
        default_work_dir()
    );
    println!(
        "  -t, --threads <arg>   worker threads the build would use (default: {})",
        default_threads()
    );
    println!("  -h, --help            print usage");
}

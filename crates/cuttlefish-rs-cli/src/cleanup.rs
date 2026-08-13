//! `cuttlefish cleanup`: remove what a bailed run left in a working directory.
//!
//! A completed build unlinks its own intermediates, so a working directory is
//! empty when it finishes. A build that is killed, runs out of disk, or hits an
//! error partway leaves them, and they are large -- the external-memory
//! intermediates for a 150000-genome graph peak in the hundreds of gigabytes.
//!
//! The awkward part is that the working directory is usually a shared scratch
//! path, so "delete everything under it" is not an option. This removes only
//! names that cuttlefish itself can produce, matched exactly, and reports
//! everything else it finds so a stale directory is legible rather than merely
//! smaller.

use std::error::Error;
use std::fs;
use std::io::ErrorKind;
use std::path::{Path, PathBuf};

type Result<T> = std::result::Result<T, Box<dyn Error + Send + Sync>>;

/// The intermediate suffixes a build can create, all under `<name>.cf3rs.`.
///
/// Kept as one table so the docs, this command, and the code that creates them
/// can be checked against each other. `stitch-coords.cpp-expansion` and
/// `lmtig-unitigs.edge-matrix` are derived from their bases by the code that
/// creates them, so they are listed rather than inferred.
const INTERMEDIATE_SUFFIXES: &[&str] = &[
    "wsk",
    "lmtig-labels",
    "lmtig-unitigs",
    "lmtig-unitigs.edge-matrix",
    "local-unitig-buckets",
    "stitch-coords",
    "stitch-coords.cpp-expansion",
    "final-unitigs",
    "colors",
    "color-runs",
    "trivial.fa",
];

/// Colored builds write this beside the *output*, and a completed one needs it
/// to interpret the FASTA, so it is never removed without being asked for.
const COLOR_REPOSITORY_SUFFIX: &str = "color-repository";

struct Params {
    work_dir: PathBuf,
    prefix: Option<String>,
    dry_run: bool,
    include_repository: bool,
}

struct Found {
    path: PathBuf,
    bytes: u64,
    entries: u64,
    is_dir: bool,
}

/// Entry point for the `cleanup` subcommand. `args` is the argument list with
/// the subcommand name already consumed.
pub fn run<I>(args: I) -> Result<i32>
where
    I: Iterator<Item = String>,
{
    let Some(params) = parse_args(args)? else {
        return Ok(0); // --help
    };
    if !params.work_dir.is_dir() {
        return Err(format!("{} is not a directory", params.work_dir.display()).into());
    }

    let mut intermediates = Vec::new();
    let mut repositories = Vec::new();
    let mut unrecognized = Vec::new();
    for entry in fs::read_dir(&params.work_dir)? {
        let entry = entry?;
        let path = entry.path();
        let Some(name) = path.file_name().and_then(|name| name.to_str()) else {
            continue;
        };
        match classify(name, params.prefix.as_deref()) {
            Some(Kind::Intermediate) => intermediates.push(measure(&path)?),
            Some(Kind::Repository) => repositories.push(measure(&path)?),
            None => unrecognized.push(name.to_string()),
        }
    }
    // Largest first: on a bailed run the interesting line is what is actually
    // holding the disk.
    intermediates.sort_by_key(|found| std::cmp::Reverse(found.bytes));
    repositories.sort_by_key(|found| std::cmp::Reverse(found.bytes));

    let mut removed_bytes = 0u64;
    let mut removed = 0usize;
    let targets = intermediates.iter().chain(
        params
            .include_repository
            .then_some(&repositories)
            .into_iter()
            .flatten(),
    );
    for found in targets {
        let kind = if found.is_dir { "dir " } else { "file" };
        println!(
            "{} {kind} {} ({}, {} entr{})",
            if params.dry_run {
                "would remove"
            } else {
                "removing"
            },
            found.path.display(),
            human_bytes(found.bytes),
            found.entries,
            if found.entries == 1 { "y" } else { "ies" },
        );
        if !params.dry_run {
            remove(&found.path)?;
        }
        removed_bytes += found.bytes;
        removed += 1;
    }

    if !params.include_repository {
        for found in &repositories {
            println!(
                "keeping dir  {} ({}) -- a colored graph's FASTA cannot be read \
                 without it; pass --include-repository to remove it too",
                found.path.display(),
                human_bytes(found.bytes),
            );
        }
    }
    if !unrecognized.is_empty() {
        println!(
            "leaving {} unrecognized entr{} alone: {}",
            unrecognized.len(),
            if unrecognized.len() == 1 { "y" } else { "ies" },
            preview(&unrecognized),
        );
    }
    println!(
        "{} {removed} intermediate(s), {}",
        if params.dry_run {
            "would remove"
        } else {
            "removed"
        },
        human_bytes(removed_bytes),
    );
    if removed == 0 && repositories.is_empty() {
        println!(
            "nothing to clean in {}; a finished build removes its own \
             intermediates",
            params.work_dir.display()
        );
    }
    Ok(0)
}

enum Kind {
    Intermediate,
    Repository,
}

/// Decides whether a directory entry is cuttlefish's to delete.
///
/// The match is exact on the whole `<prefix>.cf3rs.<suffix>` shape rather than
/// a substring test: this runs against shared scratch directories, so a name
/// that merely mentions `cf3rs` is somebody else's.
fn classify(name: &str, prefix: Option<&str>) -> Option<Kind> {
    let (found_prefix, rest) = name.split_once(".cf3rs.")?;
    if let Some(prefix) = prefix
        && found_prefix != prefix
    {
        return None;
    }
    if INTERMEDIATE_SUFFIXES.contains(&rest) {
        return Some(Kind::Intermediate);
    }
    if rest == COLOR_REPOSITORY_SUFFIX {
        return Some(Kind::Repository);
    }
    None
}

fn measure(path: &Path) -> Result<Found> {
    let metadata = fs::symlink_metadata(path)?;
    if !metadata.is_dir() {
        return Ok(Found {
            path: path.to_path_buf(),
            bytes: metadata.len(),
            entries: 1,
            is_dir: false,
        });
    }
    let mut bytes = 0u64;
    let mut entries = 0u64;
    let mut stack = vec![path.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let listing = match fs::read_dir(&dir) {
            Ok(listing) => listing,
            // A directory racing with its own creator is still ours to remove.
            Err(error) if error.kind() == ErrorKind::NotFound => continue,
            Err(error) => return Err(error.into()),
        };
        for entry in listing {
            let entry = entry?;
            let metadata = entry.metadata()?;
            if metadata.is_dir() {
                stack.push(entry.path());
            } else {
                bytes += metadata.len();
                entries += 1;
            }
        }
    }
    Ok(Found {
        path: path.to_path_buf(),
        bytes,
        entries,
        is_dir: true,
    })
}

fn remove(path: &Path) -> Result<()> {
    let metadata = fs::symlink_metadata(path)?;
    if metadata.is_dir() {
        fs::remove_dir_all(path)?;
    } else {
        fs::remove_file(path)?;
    }
    Ok(())
}

fn preview(names: &[String]) -> String {
    const SHOWN: usize = 4;
    let mut shown = names
        .iter()
        .take(SHOWN)
        .cloned()
        .collect::<Vec<_>>()
        .join(", ");
    if names.len() > SHOWN {
        shown.push_str(&format!(", and {} more", names.len() - SHOWN));
    }
    shown
}

fn human_bytes(bytes: u64) -> String {
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];
    let mut value = bytes as f64;
    let mut unit = 0;
    while value >= 1024.0 && unit + 1 < UNITS.len() {
        value /= 1024.0;
        unit += 1;
    }
    if unit == 0 {
        format!("{bytes} B")
    } else {
        format!("{value:.1} {}", UNITS[unit])
    }
}

fn parse_args<I>(args: I) -> Result<Option<Params>>
where
    I: Iterator<Item = String>,
{
    let mut work_dir = None;
    let mut prefix = None;
    let mut dry_run = false;
    let mut include_repository = false;
    let mut args = args;
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-w" | "--work-dir" => work_dir = Some(PathBuf::from(argument_value(&mut args, &arg)?)),
            "-p" | "--prefix" => prefix = Some(argument_value(&mut args, &arg)?),
            "-n" | "--dry-run" => dry_run = true,
            "--include-repository" => include_repository = true,
            "-h" | "--help" => {
                print_help();
                return Ok(None);
            }
            _ => return Err(format!("unknown argument {arg}").into()),
        }
    }
    Ok(Some(Params {
        work_dir: work_dir.ok_or("-w/--work-dir is required")?,
        prefix,
        dry_run,
        include_repository,
    }))
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
        "Remove the intermediates a bailed build left in a working directory.

Usage:
  cuttlefish cleanup -w DIR [OPTION...]

 options:
  -w, --work-dir <arg>   the --work-dir the build was given
  -p, --prefix <arg>     only this build's artifacts, by output name
  -n, --dry-run          report what would be removed, remove nothing
      --include-repository
                         also remove <name>.cf3rs.color-repository, which a
                         completed colored build needs to interpret its FASTA
  -h, --help             print this help

A finished build removes its own intermediates, so a working directory is empty
when it succeeds. Only names cuttlefish can produce are touched, matched exactly
as <name>.cf3rs.<suffix>; anything else in the directory is reported and left
alone. The output FASTA is never touched -- after a bailed run it is partial,
and it is yours to inspect or delete."
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classifies_only_cuttlefish_names() {
        assert!(matches!(
            classify("graph.cf3rs.wsk", None),
            Some(Kind::Intermediate)
        ));
        assert!(matches!(
            classify("graph.cf3rs.stitch-coords.cpp-expansion", None),
            Some(Kind::Intermediate)
        ));
        assert!(matches!(
            classify("graph.cf3rs.color-repository", None),
            Some(Kind::Repository)
        ));
        // Neighbours in a shared scratch directory.
        assert!(classify("graph.fa", None).is_none());
        assert!(classify("cf3rs-notes.txt", None).is_none());
        assert!(classify("graph.cf3rs.something-else", None).is_none());
        assert!(classify("my.cf3rs.wsk.backup", None).is_none());
        // Prefix filtering picks one build out of a shared directory.
        assert!(classify("graph.cf3rs.wsk", Some("graph")).is_some());
        assert!(classify("other.cf3rs.wsk", Some("graph")).is_none());
    }

    #[test]
    fn reports_sizes_in_units() {
        assert_eq!(human_bytes(512), "512 B");
        assert_eq!(human_bytes(2048), "2.0 KiB");
        assert_eq!(human_bytes(3 << 30), "3.0 GiB");
    }
}

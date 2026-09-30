//! Whether to compress intermediate record streams, chosen from the storage
//! under the work directory.
//!
//! lz4-blocking the local-unitig, coordinate and colour streams costs CPU
//! and saves I/O. On a host whose storage outruns the build it is a small
//! loss; on slower storage it wins. Measured on 10k Salmonella at 16 threads,
//! over two striped NVMe drives with the intermediates held in page cache,
//! compression cost about 0.26 s of wall time per GB of writes saved:
//! colored builds wrote 6.5 GB less and ran 2.4% slower. It pays wherever
//! writing and reading those bytes costs more than that.
//!
//! [`choose`] makes the call with [`scratch_probe`].
//! * A network filesystem or a rotational disk: compress, without timing
//!   anything.
//! * A FUSE filesystem: compress without timing, with a warning. It may
//!   front a fast local disk or a remote object store, and a short probe
//!   cannot tell a cache from the store behind it.
//! * Otherwise, a probe of at most 256 MiB and about half a second each way
//!   times direct writes and reads. Storage that sustains twice the rate the build
//!   produces intermediates ([`intermediate_bytes_per_second`]) absorbs them
//!   without stalling the build, so compression is left off; slower storage
//!   gets it.
//! * A probe that fails compresses, the safer choice on storage it cannot
//!   measure.

use std::path::Path;

use crate::params::IntermediateCompression;

/// Intermediate bytes per second a build produces at 16 threads: 10k
/// Salmonella wrote 63 GB in 70 s colored and 49 GB in 49 s uncolored. Larger
/// inputs run at lower average rates (150k: 0.50 GB/s colored, 0.55 GB/s
/// uncolored), so the smaller corpus gives the demanding estimate.
pub const BYTES_PER_SECOND_AT_16_THREADS: f64 = 0.96e9;

/// How the rate grows with threads. The phases do not scale linearly: 150k
/// wrote 1.33 GB/s colored and 1.24 GB/s uncolored at 64 threads against
/// 0.50 and 0.55 at 16, a growth of threads^0.65.
pub const THREAD_SCALING_EXPONENT: f64 = 0.65;

/// The rate, in bytes per second, at which a build of `threads` workers
/// produces intermediates.
pub fn intermediate_bytes_per_second(threads: usize) -> f64 {
    BYTES_PER_SECOND_AT_16_THREADS * (threads.max(1) as f64 / 16.0).powf(THREAD_SCALING_EXPONENT)
}

/// How far storage must outrun that rate to be left uncompressed.
pub const STORAGE_MARGIN: f64 = 2.0;

/// The decision, and why.
#[derive(Debug, Clone)]
pub struct CompressionChoice {
    pub compress: bool,
    pub reason: String,
    /// The storage measurement, when one was made.
    pub probe: Option<scratch_probe::Probe>,
    /// Something the user should know about the choice, such as a guess
    /// made without measuring.
    pub warning: Option<String>,
}

/// Decides whether to compress intermediates for a build of `threads` workers
/// writing into `work_dir`, which must exist.
pub fn choose(
    setting: IntermediateCompression,
    work_dir: &Path,
    threads: usize,
) -> CompressionChoice {
    match setting {
        IntermediateCompression::On => CompressionChoice {
            compress: true,
            reason: "requested".to_string(),
            probe: None,
            warning: None,
        },
        IntermediateCompression::Off => CompressionChoice {
            compress: false,
            reason: "requested".to_string(),
            probe: None,
            warning: None,
        },
        IntermediateCompression::Auto => choose_auto(work_dir, threads),
    }
}

fn choose_auto(work_dir: &Path, threads: usize) -> CompressionChoice {
    if let Some(choice) = choose_without_timing(&scratch_probe::classify(work_dir)) {
        return choice;
    }
    let probe = match scratch_probe::probe(work_dir, &scratch_probe::Limits::default()) {
        Ok(probe) => probe,
        Err(error) => {
            return CompressionChoice {
                compress: true,
                reason: format!("could not measure the work directory's storage ({error})"),
                probe: None,
                warning: None,
            };
        }
    };
    decide(probe, intermediate_bytes_per_second(threads))
}

/// The storage that is compressed for without timing, from metadata alone:
/// network filesystems and rotational disks, which are slow as a rule, and
/// FUSE mounts, whose speed says little about what they front.
fn choose_without_timing(storage: &scratch_probe::Storage) -> Option<CompressionChoice> {
    if storage.kind == scratch_probe::StorageKind::Network {
        return Some(CompressionChoice {
            compress: true,
            reason: format!(
                "work directory is on a network filesystem ({})",
                storage.filesystem
            ),
            probe: None,
            warning: None,
        });
    }
    if storage.kind == scratch_probe::StorageKind::Fuse {
        return Some(CompressionChoice {
            compress: true,
            reason: format!(
                "work directory is on a FUSE filesystem ({})",
                storage.filesystem
            ),
            probe: None,
            warning: Some(
                "the work directory is on a FUSE filesystem, which may be a local disk or a \
                 remote store; intermediates are compressed without timing it. Pass \
                 --compress-intermediates off if it is fast local storage."
                    .to_string(),
            ),
        });
    }
    if storage.rotational == Some(true) {
        return Some(CompressionChoice {
            compress: true,
            reason: format!(
                "work directory is on a rotational disk ({})",
                storage.filesystem
            ),
            probe: None,
            warning: None,
        });
    }
    None
}

/// The timed half of the decision, on a finished probe.
fn decide(probe: scratch_probe::Probe, need: f64) -> CompressionChoice {
    let compress = !probe.sustains(need, STORAGE_MARGIN);
    let reason = format!(
        "{} storage ({}{}) writes {:.2} GB/s and reads {:.2} GB/s; the build produces about {:.2} GB/s, so storage must sustain {:.2} GB/s",
        probe.storage.filesystem,
        if probe.direct_io {
            "direct I/O"
        } else {
            "buffered I/O"
        },
        match probe.storage.rotational {
            Some(false) => ", solid state",
            _ => "",
        },
        probe.write.bytes_per_second() / 1e9,
        probe.read.bytes_per_second() / 1e9,
        need / 1e9,
        need * STORAGE_MARGIN / 1e9,
    );
    CompressionChoice {
        compress,
        reason,
        probe: Some(probe),
        warning: None,
    }
}

/// Makes `compress` the setting for every intermediate stream opened after
/// this call.
pub fn apply(compress: bool) {
    crate::block_io::set_compress_blocks(compress);
}

/// Applies a build's own setting before it opens any intermediate stream:
/// `on` and `off` as given, and `auto` by probing the work directory, which
/// is reported. The library's build entry points call this, so a caller that
/// wants the decision made earlier (the CLI makes it before partitioning)
/// resolves `auto` itself and passes `on` or `off` down.
///
/// The setting is process-wide. Builds running concurrently in one process
/// with different settings may each write some streams under the other's;
/// every block records how it is stored, so that changes speed, never the
/// output.
pub fn apply_params(params: &crate::params::BuildParams) {
    let compress = match params.compress_intermediates {
        IntermediateCompression::On => true,
        IntermediateCompression::Off => false,
        IntermediateCompression::Auto => {
            let work_dir = Path::new(&params.work_dir);
            let _ = std::fs::create_dir_all(work_dir);
            let choice = choose(IntermediateCompression::Auto, work_dir, params.threads);
            if let Some(warning) = &choice.warning {
                eprintln!("cuttlefish: warning: {warning}");
            }
            eprintln!(
                "cuttlefish: intermediate compression {} (auto: {})",
                if choice.compress { "on" } else { "off" },
                choice.reason
            );
            choice.compress
        }
    };
    apply(compress);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    fn probe_at(write_gb_s: f64, read_gb_s: f64) -> scratch_probe::Probe {
        let throughput = |gb_s: f64| scratch_probe::Throughput {
            bytes: (gb_s * 1e9) as u64,
            elapsed: Duration::from_secs(1),
        };
        scratch_probe::Probe {
            dir: std::path::PathBuf::from("/work"),
            storage: scratch_probe::Storage {
                filesystem: "xfs",
                kind: scratch_probe::StorageKind::Local,
                rotational: Some(false),
            },
            write: throughput(write_gb_s),
            read: throughput(read_gb_s),
            direct_io: true,
            mem_available: None,
        }
    }

    #[test]
    fn fast_storage_skips_compression_and_slow_storage_gets_it() {
        // 16 threads produce about 0.96 GB/s; with the margin, 1.92 GB/s.
        let need = intermediate_bytes_per_second(16);
        assert!((need - 0.96e9).abs() < 1.0);
        assert!(!decide(probe_at(5.7, 6.8), need).compress, "striped NVMe");
        assert!(decide(probe_at(0.5, 0.55), need).compress, "SATA SSD");
        assert!(decide(probe_at(5.7, 1.5), need).compress, "slow reads");
        // The same NVMe still keeps up at 64 threads (about 2.4 GB/s), but
        // not at 256 (about 5.8 GB/s).
        let nvme = probe_at(5.7, 6.8);
        assert!(!decide(nvme.clone(), intermediate_bytes_per_second(64)).compress);
        assert!(decide(nvme, intermediate_bytes_per_second(256)).compress);
    }

    #[test]
    fn explicit_settings_do_not_probe() {
        let missing = Path::new("/nonexistent/cuttlefish-work");
        assert!(choose(IntermediateCompression::On, missing, 16).compress);
        let off = choose(IntermediateCompression::Off, missing, 16);
        assert!(!off.compress && off.probe.is_none());
    }

    /// A build's own setting reaches the writers it opens. The one test that
    /// changes the process-wide setting; others that care set theirs per
    /// writer, and readers take either.
    #[test]
    fn a_builds_setting_reaches_its_writers() {
        use std::io::Write;
        let stored = |setting: IntermediateCompression| {
            let mut params = crate::params::BuildParams::new(
                crate::GraphInput::References,
                "unused".to_string(),
            );
            params.compress_intermediates = setting;
            apply_params(&params);
            let mut writer = crate::block_io::Lz4BlockWriter::new(Vec::new());
            writer.write_all(&[7u8; 64 * 1024]).unwrap();
            writer.into_inner().len()
        };
        assert!(
            stored(IntermediateCompression::Off) > 64 * 1024,
            "stored raw"
        );
        assert!(stored(IntermediateCompression::On) < 4096, "compressed");
    }

    #[test]
    fn some_storage_compresses_without_timing() {
        let storage = |kind, rotational| scratch_probe::Storage {
            filesystem: "test",
            kind,
            rotational,
        };
        use scratch_probe::StorageKind::{Fuse, Local, Memory, Network, Unknown};
        let fuse = choose_without_timing(&storage(Fuse, None)).unwrap();
        assert!(fuse.compress && fuse.probe.is_none() && fuse.warning.is_some());
        for (kind, rotational) in [(Network, None), (Local, Some(true))] {
            let choice = choose_without_timing(&storage(kind, rotational)).unwrap();
            assert!(choice.compress && choice.warning.is_none());
        }
        for (kind, rotational) in [
            (Local, Some(false)),
            (Local, None),
            (Memory, None),
            (Unknown, None),
        ] {
            assert!(choose_without_timing(&storage(kind, rotational)).is_none());
        }
    }

    #[test]
    fn unmeasurable_storage_compresses() {
        let choice = choose(
            IntermediateCompression::Auto,
            Path::new("/nonexistent/cuttlefish-work"),
            16,
        );
        assert!(choice.compress && choice.probe.is_none());
    }
}

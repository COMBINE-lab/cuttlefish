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
//! * Otherwise, a probe of at most 256 MiB and about half a second times
//!   direct writes and reads. Storage that sustains twice the rate the build
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
        },
        IntermediateCompression::Off => CompressionChoice {
            compress: false,
            reason: "requested".to_string(),
            probe: None,
        },
        IntermediateCompression::Auto => choose_auto(work_dir, threads),
    }
}

fn choose_auto(work_dir: &Path, threads: usize) -> CompressionChoice {
    let storage = scratch_probe::classify(work_dir);
    if storage.kind == scratch_probe::StorageKind::Network {
        return CompressionChoice {
            compress: true,
            reason: format!(
                "work directory is on a network filesystem ({})",
                storage.filesystem
            ),
            probe: None,
        };
    }
    if storage.rotational == Some(true) {
        return CompressionChoice {
            compress: true,
            reason: format!(
                "work directory is on a rotational disk ({})",
                storage.filesystem
            ),
            probe: None,
        };
    }
    let probe = match scratch_probe::probe(work_dir, &scratch_probe::Limits::default()) {
        Ok(probe) => probe,
        Err(error) => {
            return CompressionChoice {
                compress: true,
                reason: format!("could not measure the work directory's storage ({error})"),
                probe: None,
            };
        }
    };
    decide(probe, intermediate_bytes_per_second(threads))
}

/// The timed half of the decision, on a finished probe.
fn decide(probe: scratch_probe::Probe, need: f64) -> CompressionChoice {
    let compress = !probe.sustains(need, STORAGE_MARGIN);
    let reason = format!(
        "{} storage ({}{}) writes {:.2} GB/s and reads {:.2} GB/s; the build produces about {:.2} GB/s",
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
    );
    CompressionChoice {
        compress,
        reason,
        probe: Some(probe),
    }
}

/// Makes `compress` the setting for every intermediate stream opened after
/// this call.
pub fn apply(compress: bool) {
    crate::block_io::set_compress_blocks(compress);
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

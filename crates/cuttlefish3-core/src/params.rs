//! Build configuration and phase-local resource planning.
//!
//! [`BuildParams::threads`] is an upper bound. When `max_memory_gb` is set,
//! partition, local contraction, and post-local processing independently lower
//! concurrency according to replicated state estimates. The limit is soft:
//! workload-sized shared tables can impose a higher minimum RSS.

use crate::{
    DEFAULT_CUTOFF_READS, DEFAULT_CUTOFF_REFS, DEFAULT_GMTIG_BUCKETS, DEFAULT_K,
    DEFAULT_LMTIG_BUCKETS, DEFAULT_MINIMIZER_LEN, DEFAULT_VERTEX_PARTITIONS, GraphInput, MAX_K,
    MAX_MINIMIZER_LEN, default_threads, default_work_dir,
};

const GIB: usize = 1024 * 1024 * 1024;
// This is the soft-memory policy for the Rust-only `--max-memory` extension,
// not a C++ algorithm parameter. These estimates cover phase-local replicated
// state; workload-sized graph tables can impose a higher minimum.
const PARTITION_FIXED_MEMORY: usize = 2 * GIB;
const PARTITION_MEMORY_PER_WORKER: usize = 32 * 1024 * 1024;
const LOCAL_FIXED_MEMORY: usize = 4 * GIB;
const LOCAL_MEMORY_PER_WORKER: usize = 160 * 1024 * 1024;
const POST_LOCAL_FIXED_MEMORY: usize = 4 * GIB;
const POST_LOCAL_MEMORY_PER_WORKER: usize = 64 * 1024 * 1024;

/// Configuration shared by partitioning and graph construction phases.
///
/// Bucket and partition counts must be powers of two. The public fields support
/// programmatic construction, while [`BuildParams::validate`] enforces the
/// cross-field invariants expected by the production pipeline.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BuildParams {
    /// Whether inputs are sequencing reads or references.
    pub input: GraphInput,
    /// Explicit sequence file paths.
    pub seqs: Vec<String>,
    /// Files containing one input path per line.
    pub lists: Vec<String>,
    /// Directories whose regular files are used as inputs.
    pub dirs: Vec<String>,
    /// K-mer length. Must be odd and no greater than [`crate::MAX_K`].
    pub k: u16,
    /// Minimizer length used for weak-super-k-mer partitioning.
    pub minimizer_len: u16,
    /// Optional edge-frequency cutoff; input-class defaults apply when absent.
    pub cutoff: Option<u32>,
    /// Whether to construct and emit positional colors.
    pub color: bool,
    /// Whether weak-super-k-mer buckets use LZ4 block compression.
    pub compress_buckets: bool,
    /// Whether an input that fails to parse is skipped instead of aborting.
    ///
    /// A skipped source retains its position in the input list, so colored
    /// source assignments are unaffected. Records read before the failure are
    /// retained, so a file that fails part-way through contributes its prefix.
    pub skip_unreadable: bool,
    /// Prefix for final graph files.
    pub output_prefix: String,
    /// Directory for external-memory intermediates.
    pub work_dir: String,
    /// Number of blocked discontinuity-graph vertex partitions.
    pub vertex_partitions: usize,
    /// Logical local-unitig bucket count.
    pub lmtig_buckets: usize,
    /// Logical maximal-unitig coordinate bucket count.
    pub gmtig_buckets: usize,
    /// Maximum worker count requested for each phase.
    pub threads: usize,
    /// Optional soft memory budget in GiB.
    ///
    /// The budget controls replicated phase state and is not a hard RSS limit.
    pub max_memory_gb: Option<usize>,
}

impl BuildParams {
    /// Constructs parameters with Cuttlefish-compatible defaults.
    pub fn new(input: GraphInput, output_prefix: String) -> Self {
        Self {
            input,
            seqs: Vec::new(),
            lists: Vec::new(),
            dirs: Vec::new(),
            k: DEFAULT_K,
            minimizer_len: DEFAULT_MINIMIZER_LEN,
            cutoff: None,
            color: false,
            compress_buckets: false,
            skip_unreadable: false,
            output_prefix,
            work_dir: default_work_dir(),
            vertex_partitions: DEFAULT_VERTEX_PARTITIONS,
            lmtig_buckets: DEFAULT_LMTIG_BUCKETS,
            gmtig_buckets: DEFAULT_GMTIG_BUCKETS,
            threads: default_threads(),
            max_memory_gb: None,
        }
    }

    /// Returns the explicit cutoff or the input-class default.
    #[inline]
    pub fn cutoff(&self) -> u32 {
        self.cutoff.unwrap_or_else(|| self.input.default_cutoff())
    }

    /// Returns the soft memory budget in bytes, if configured.
    pub fn max_memory_bytes(&self) -> Option<usize> {
        self.max_memory_gb.and_then(|gb| gb.checked_mul(GIB))
    }

    fn memory_bounded_workers(&self, fixed: usize, per_worker: usize) -> usize {
        let requested = self.threads.max(1);
        let Some(budget) = self.max_memory_bytes() else {
            return requested;
        };
        let bounded = requested
            .min(budget.saturating_sub(fixed) / per_worker)
            .max(1);
        if bounded == requested {
            return requested;
        }
        // When the estimate requires a reduction, round down to leave headroom
        // for workload-sized shared tables outside this replicated-state estimate.
        1usize << bounded.ilog2()
    }

    /// Returns the partition worker count after input and memory bounds.
    pub fn partition_workers(&self, input_files: usize) -> usize {
        self.memory_bounded_workers(PARTITION_FIXED_MEMORY, PARTITION_MEMORY_PER_WORKER)
            .min(input_files.max(1))
    }

    /// Returns the memory-bounded local-contraction worker count.
    pub fn local_workers(&self) -> usize {
        self.memory_bounded_workers(LOCAL_FIXED_MEMORY, LOCAL_MEMORY_PER_WORKER)
    }

    /// Returns the memory-bounded discontinuity/collation worker count.
    pub fn post_local_workers(&self) -> usize {
        self.memory_bounded_workers(POST_LOCAL_FIXED_MEMORY, POST_LOCAL_MEMORY_PER_WORKER)
    }

    /// Validates input, graph dimensions, and resource parameters.
    pub fn validate(&self) -> Result<(), ParamError> {
        if self.seqs.is_empty() && self.lists.is_empty() && self.dirs.is_empty() {
            return Err(ParamError::NoInput);
        }

        if self.output_prefix.is_empty() {
            return Err(ParamError::MissingOutput);
        }

        if self.k <= 1 || self.k > MAX_K || self.k % 2 == 0 {
            return Err(ParamError::InvalidK(self.k));
        }

        if self.minimizer_len == 0
            || self.minimizer_len >= self.k
            || self.minimizer_len > MAX_MINIMIZER_LEN
        {
            return Err(ParamError::InvalidMinimizerLen {
                k: self.k,
                l: self.minimizer_len,
            });
        }

        if self.cutoff() == 0 {
            return Err(ParamError::InvalidCutoff);
        }

        if !self.vertex_partitions.is_power_of_two()
            || !self.lmtig_buckets.is_power_of_two()
            || !self.gmtig_buckets.is_power_of_two()
        {
            return Err(ParamError::BucketCountsMustBePowersOfTwo);
        }
        if self.threads == 0 {
            return Err(ParamError::InvalidThreadCount);
        }
        if self.max_memory_gb == Some(0)
            || self.max_memory_gb.is_some_and(|gb| gb > usize::MAX / GIB)
        {
            return Err(ParamError::InvalidMemoryLimit);
        }

        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParamError {
    NoInput,
    MissingOutput,
    InvalidK(u16),
    InvalidMinimizerLen { k: u16, l: u16 },
    InvalidCutoff,
    BucketCountsMustBePowersOfTwo,
    InvalidThreadCount,
    InvalidMemoryLimit,
}

impl std::fmt::Display for ParamError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NoInput => write!(f, "no sequence input provided"),
            Self::MissingOutput => write!(f, "missing output prefix"),
            Self::InvalidK(k) => write!(f, "k-mer length {k} is invalid; expected odd 3..={MAX_K}"),
            Self::InvalidMinimizerLen { k, l } => {
                write!(f, "minimizer length {l} is invalid for k={k}")
            }
            Self::InvalidCutoff => write!(f, "cutoff frequency must be at least 1"),
            Self::BucketCountsMustBePowersOfTwo => write!(f, "bucket counts must be powers of two"),
            Self::InvalidThreadCount => write!(f, "thread count must be at least 1"),
            Self::InvalidMemoryLimit => write!(f, "maximum memory must be at least 1 GiB"),
        }
    }
}

impl std::error::Error for ParamError {}

pub const _DEFAULT_CUTOFFS_FOR_DOCS: (u32, u32) = (DEFAULT_CUTOFF_REFS, DEFAULT_CUTOFF_READS);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validates_current_cpp_defaults() {
        let mut p = BuildParams::new(GraphInput::References, "out".to_string());
        p.seqs.push("data/refs1.fa".to_string());
        assert_eq!(p.cutoff(), 1);
        assert!(p.validate().is_ok());
    }

    #[test]
    fn rejects_even_k_and_zero_cutoff() {
        let mut p = BuildParams::new(GraphInput::Reads, "out".to_string());
        p.seqs.push("data/reads.fq".to_string());
        p.k = 32;
        assert_eq!(p.validate(), Err(ParamError::InvalidK(32)));
        p.k = 31;
        p.cutoff = Some(0);
        assert_eq!(p.validate(), Err(ParamError::InvalidCutoff));
    }

    #[test]
    fn explicit_memory_budget_bounds_phase_concurrency() {
        let mut p = BuildParams::new(GraphInput::References, "out".to_string());
        p.threads = 256;
        p.max_memory_gb = Some(16);

        assert_eq!(p.partition_workers(1_000), 256);
        assert_eq!(p.local_workers(), 64);
        assert_eq!(p.post_local_workers(), 128);
        assert_eq!(p.partition_workers(17), 17);
    }

    #[test]
    fn unrestricted_concurrency_follows_the_user_request() {
        let mut p = BuildParams::new(GraphInput::References, "out".to_string());
        p.threads = 96;

        assert_eq!(p.partition_workers(1_000), 96);
        assert_eq!(p.local_workers(), 96);
        assert_eq!(p.post_local_workers(), 96);
    }

    #[test]
    fn ample_memory_budget_preserves_non_power_of_two_request() {
        let mut p = BuildParams::new(GraphInput::References, "out".to_string());
        p.threads = 96;
        p.max_memory_gb = Some(24);

        assert_eq!(p.partition_workers(1_000), 96);
        assert_eq!(p.local_workers(), 96);
        assert_eq!(p.post_local_workers(), 96);
    }

    #[test]
    fn binding_memory_budget_only_reduces_constrained_phases() {
        let mut p = BuildParams::new(GraphInput::References, "out".to_string());
        p.threads = 96;
        p.max_memory_gb = Some(12);

        assert_eq!(p.partition_workers(1_000), 96);
        assert_eq!(p.local_workers(), 32);
        assert_eq!(p.post_local_workers(), 96);
    }

    #[test]
    fn rejects_zero_memory_budget() {
        let mut p = BuildParams::new(GraphInput::References, "out".to_string());
        p.seqs.push("data/refs1.fa".to_string());
        p.max_memory_gb = Some(0);
        assert_eq!(p.validate(), Err(ParamError::InvalidMemoryLimit));
    }
}

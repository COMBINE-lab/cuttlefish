//! External-memory compacted de Bruijn graph construction.
//!
//! This crate implements the Rust Cuttlefish 3 pipeline. The production path
//! intentionally follows the phase ordering of the C++ implementation:
//!
//! 1. [`partition`] parses sequences and emits weak super-k-mers into atlas
//!    buckets.
//! 2. [`subgraph`] constructs and contracts each local de Bruijn subgraph.
//! 3. [`discontinuity`] contracts and expands the external discontinuity graph,
//!    then collates local unitigs into maximal unitigs.
//! 4. [`color`] stores deduplicated source sets and positional color runs.
//!
//! Intermediate formats are private to this implementation. Their compact
//! layouts, ordering, and bucket fanouts are performance-sensitive; changing
//! them requires both topology tests and matched scale benchmarks.
#![warn(rustdoc::broken_intra_doc_links)]

pub mod buckets;
pub mod color;
pub mod colored;
pub mod discontinuity;
pub mod dna;
pub mod hash;
pub mod input;
pub mod kmer;
pub mod minimizer;
pub mod params;
pub mod partition;
pub mod state;
pub mod subgraph;
pub mod uncolored;

/// Default k-mer length used by the CLI.
pub const DEFAULT_K: u16 = 31;
/// Default minimizer length used to partition k-mers.
pub const DEFAULT_MINIMIZER_LEN: u16 = 12;
/// Default edge-frequency cutoff for sequencing-read input.
pub const DEFAULT_CUTOFF_READS: u32 = 2;
/// Default edge-frequency cutoff for reference-sequence input.
pub const DEFAULT_CUTOFF_REFS: u32 = 1;
/// Number of atlases in the C++-compatible default partition layout.
pub const DEFAULT_ATLAS_COUNT: usize = 128;
/// Number of subgraphs buffered by each atlas.
pub const DEFAULT_GRAPHS_PER_ATLAS: usize = 128;
/// Total number of local subgraphs in the default partition layout.
pub const DEFAULT_SUBGRAPH_COUNT: usize = DEFAULT_ATLAS_COUNT * DEFAULT_GRAPHS_PER_ATLAS;
/// Number of vertex partitions in the external discontinuity graph.
pub const DEFAULT_VERTEX_PARTITIONS: usize = 128;
/// Default number of local-unitig buckets.
pub const DEFAULT_LMTIG_BUCKETS: usize = 1024;
/// Default number of maximal-unitig coordinate buckets.
pub const DEFAULT_GMTIG_BUCKETS: usize = 1024;
/// Largest k-mer length supported by the packed two-word representation.
pub const MAX_K: u16 = 63;
/// Largest supported minimizer length.
pub const MAX_MINIMIZER_LEN: u16 = 32;

/// Returns the platform temporary directory used for intermediate files.
pub fn default_work_dir() -> String {
    std::env::temp_dir().to_string_lossy().into_owned()
}

/// Returns the default worker count, matching Cuttlefish's quarter-machine policy.
pub fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(|threads| (threads.get() / 4).max(1))
        .unwrap_or(8)
}

/// Configures the process-wide Rayon pool for code paths that use it.
///
/// The compact production pipeline primarily uses phase-local pools. Call this
/// at most once, before performing parallel work.
pub fn configure_global_parallelism(threads: usize) -> Result<(), rayon::ThreadPoolBuildError> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build_global()
}

/// Side of a canonical k-mer vertex.
///
/// `Front` and `Back` describe incidence in the canonical representation, not
/// the strand chosen for a final unitig label.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Side {
    /// The prefix side of the canonical vertex.
    Front = 0,
    /// The suffix side of the canonical vertex.
    Back = 1,
}

impl Side {
    /// Returns the opposite side of the vertex.
    #[inline]
    pub const fn inverse(self) -> Self {
        match self {
            Self::Front => Self::Back,
            Self::Back => Self::Front,
        }
    }
}

/// Semantic class of sequence input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GraphInput {
    /// Sequencing reads, for which low-frequency edges are normally filtered.
    Reads,
    /// Reference sequences, whose edges are retained by default.
    References,
}

impl GraphInput {
    /// Returns the default edge-frequency cutoff for this input class.
    #[inline]
    pub const fn default_cutoff(self) -> u32 {
        match self {
            Self::Reads => DEFAULT_CUTOFF_READS,
            Self::References => DEFAULT_CUTOFF_REFS,
        }
    }
}

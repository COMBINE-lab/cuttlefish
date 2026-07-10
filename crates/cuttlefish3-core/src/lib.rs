pub mod buckets;
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

pub const DEFAULT_K: u16 = 31;
pub const DEFAULT_MINIMIZER_LEN: u16 = 12;
pub const DEFAULT_CUTOFF_READS: u32 = 2;
pub const DEFAULT_CUTOFF_REFS: u32 = 1;
pub const DEFAULT_WORK_DIR: &str = "/scratch3/tmp";
pub const DEFAULT_ATLAS_COUNT: usize = 128;
pub const DEFAULT_GRAPHS_PER_ATLAS: usize = 128;
pub const DEFAULT_SUBGRAPH_COUNT: usize = DEFAULT_ATLAS_COUNT * DEFAULT_GRAPHS_PER_ATLAS;
pub const DEFAULT_VERTEX_PARTITIONS: usize = 128;
pub const DEFAULT_LMTIG_BUCKETS: usize = 1024;
pub const DEFAULT_GMTIG_BUCKETS: usize = 1024;
pub const DEFAULT_THREADS: usize = 1;
pub const MAX_K: u16 = 63;
pub const MAX_MINIMIZER_LEN: u16 = 32;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Side {
    Front = 0,
    Back = 1,
}

impl Side {
    #[inline]
    pub const fn inverse(self) -> Self {
        match self {
            Self::Front => Self::Back,
            Self::Back => Self::Front,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GraphInput {
    Reads,
    References,
}

impl GraphInput {
    #[inline]
    pub const fn default_cutoff(self) -> u32 {
        match self {
            Self::Reads => DEFAULT_CUTOFF_READS,
            Self::References => DEFAULT_CUTOFF_REFS,
        }
    }
}

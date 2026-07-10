use crate::{
    DEFAULT_CUTOFF_READS, DEFAULT_CUTOFF_REFS, DEFAULT_GMTIG_BUCKETS, DEFAULT_K,
    DEFAULT_LMTIG_BUCKETS, DEFAULT_MINIMIZER_LEN, DEFAULT_THREADS, DEFAULT_VERTEX_PARTITIONS,
    DEFAULT_WORK_DIR, GraphInput, MAX_K, MAX_MINIMIZER_LEN,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BuildParams {
    pub input: GraphInput,
    pub seqs: Vec<String>,
    pub lists: Vec<String>,
    pub dirs: Vec<String>,
    pub k: u16,
    pub minimizer_len: u16,
    pub cutoff: Option<u32>,
    pub color: bool,
    pub output_prefix: String,
    pub work_dir: String,
    pub vertex_partitions: usize,
    pub lmtig_buckets: usize,
    pub gmtig_buckets: usize,
    pub threads: usize,
}

impl BuildParams {
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
            output_prefix,
            work_dir: DEFAULT_WORK_DIR.to_string(),
            vertex_partitions: DEFAULT_VERTEX_PARTITIONS,
            lmtig_buckets: DEFAULT_LMTIG_BUCKETS,
            gmtig_buckets: DEFAULT_GMTIG_BUCKETS,
            threads: DEFAULT_THREADS,
        }
    }

    #[inline]
    pub fn cutoff(&self) -> u32 {
        self.cutoff.unwrap_or_else(|| self.input.default_cutoff())
    }

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
}

use crate::color::ColorError;
use crate::discontinuity::{
    DiscontinuityInputError, SerialCollationError, SerialUncoloredCollator,
    emit_colored_external_discontinuity_inputs_with_threads_in_dir, report_process_memory,
    trim_process_allocations,
};
use crate::input::{InputError, expand_input_paths};
use crate::kmer::KmerError;
use crate::params::BuildParams;
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ColoredBuildStats {
    pub input_buckets: usize,
    pub bucket_records: u64,
    pub unitigs: u64,
    pub unitig_bases: u64,
    pub output_path: PathBuf,
    pub color_repository: PathBuf,
}

pub fn build_colored_from_buckets<const K: usize>(
    params: &BuildParams,
    bucket_dir: impl AsRef<Path>,
) -> Result<ColoredBuildStats, ColoredBuildError> {
    if !params.color {
        return Err(ColoredBuildError::ColorRequired);
    }
    let bucket_dir = bucket_dir.as_ref().to_path_buf();
    let output_name = Path::new(&params.output_prefix)
        .file_name()
        .and_then(|name| name.to_str())
        .filter(|name| !name.is_empty())
        .unwrap_or("cuttlefish3");
    let work_dir = PathBuf::from(&params.work_dir);
    let label_path = work_dir.join(format!("{output_name}.cf3rs.lmtig-labels"));
    let color_path = work_dir.join(format!("{output_name}.cf3rs.colors"));

    let local_started = Instant::now();
    report_process_memory("before colored local contraction");
    let mut inputs = emit_colored_external_discontinuity_inputs_with_threads_in_dir::<K>(
        &bucket_dir,
        params.cutoff(),
        params.threads,
        &label_path,
        &color_path,
    )?;
    eprintln!(
        "cuttlefish3-rs: colored local contraction completed in {:.3}s; {} local unitig(s)",
        local_started.elapsed().as_secs_f64(),
        inputs.stats.local_unitigs,
    );
    if std::env::var_os("CF3_RS_KEEP_INTERMEDIATES").is_none() {
        std::fs::remove_dir_all(&bucket_dir).map_err(|source| ColoredBuildError::Cleanup {
            path: bucket_dir.clone(),
            source,
        })?;
    }
    trim_process_allocations();

    let output_path = PathBuf::from(format!("{}.fa", params.output_prefix));
    let coord_dir = work_dir.join(format!("{output_name}.cf3rs.stitch-coords"));
    let final_dir = work_dir.join(format!("{output_name}.cf3rs.final-unitigs"));
    let stats = SerialUncoloredCollator::collate_external_stitched_to_fasta_with_threads_in_dir(
        &mut inputs,
        params.threads.min(64),
        &coord_dir,
        &final_dir,
        &output_path,
    )?;
    report_process_memory("after colored collation");
    let color_repository = inputs
        .color_repository()
        .ok_or(ColoredBuildError::MissingColorArtifacts)?
        .clone();
    color_repository.write_metadata(params.k, &output_path, &expand_input_paths(params)?)?;
    Ok(ColoredBuildStats {
        input_buckets: inputs.stats.input_buckets,
        bucket_records: inputs.stats.weak_superkmers,
        unitigs: stats.emitted_unitigs,
        unitig_bases: stats.emitted_bases,
        output_path,
        color_repository: color_repository.dir,
    })
}

#[derive(Debug)]
pub enum ColoredBuildError {
    ColorRequired,
    MissingColorArtifacts,
    SourceInput(InputError),
    Input(DiscontinuityInputError),
    Collation(SerialCollationError),
    Kmer(KmerError),
    Color(ColorError),
    Cleanup {
        path: PathBuf,
        source: std::io::Error,
    },
}

impl From<DiscontinuityInputError> for ColoredBuildError {
    fn from(value: DiscontinuityInputError) -> Self {
        Self::Input(value)
    }
}

impl From<SerialCollationError> for ColoredBuildError {
    fn from(value: SerialCollationError) -> Self {
        Self::Collation(value)
    }
}

impl From<InputError> for ColoredBuildError {
    fn from(value: InputError) -> Self {
        Self::SourceInput(value)
    }
}

impl From<ColorError> for ColoredBuildError {
    fn from(value: ColorError) -> Self {
        Self::Color(value)
    }
}

impl std::fmt::Display for ColoredBuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ColorRequired => write!(f, "colored build requires --color"),
            Self::MissingColorArtifacts => {
                write!(f, "colored build did not produce color artifacts")
            }
            Self::SourceInput(err) => write!(f, "{err}"),
            Self::Input(err) => write!(f, "{err}"),
            Self::Collation(err) => write!(f, "{err}"),
            Self::Kmer(err) => write!(f, "{err}"),
            Self::Color(err) => write!(f, "{err}"),
            Self::Cleanup { path, source } => write!(f, "{}: {source}", path.display()),
        }
    }
}

impl std::error::Error for ColoredBuildError {}

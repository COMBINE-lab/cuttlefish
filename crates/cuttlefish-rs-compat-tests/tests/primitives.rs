use cuttlefish_rs::Side;
use cuttlefish_rs::buckets::{BucketLocation, BucketStore, read_container_manifest};
use cuttlefish_rs::discontinuity::{
    DISCONTINUITY_PARALLELIZATION_OPPORTUNITIES, DiscontinuityEdge, DiscontinuityEndpoint,
    DiscontinuityInputStats, DiscontinuityInputs, EdgePathInfo, FullSerialContractionStats,
    FullSerialDiscontinuityContraction, MatrixEndpoint, OwnedDiscontinuityUnitig, PathInfo,
    SerialDiscontinuityContractor, SerialDiscontinuityExpander, SerialEdgeMatrix, SerialExpansion,
    SerialExpansionStats, SerialMetaVertex, SerialUncoloredCollator,
    emit_colored_external_discontinuity_inputs_with_threads_in_dir,
    emit_uncolored_discontinuity_inputs,
};
use cuttlefish_rs::dna::Base;
use cuttlefish_rs::input::{expand_input_paths, parse_fragments};
use cuttlefish_rs::kmer::Kmer;
use cuttlefish_rs::minimizer::canonical_minimizer;
use cuttlefish_rs::params::BuildParams;
use cuttlefish_rs::partition::{emit_weak_superkmer_buckets, partition_fragment, partition_inputs};
use cuttlefish_rs::state::VertexState;
use cuttlefish_rs::subgraph::LocalSubgraph;
use cuttlefish_rs::uncolored::{
    build_uncolored_from_buckets, build_uncolored_with_serial_discontinuity_pipeline,
};
use cuttlefish_rs::{DEFAULT_SUBGRAPH_COUNT, GraphInput};
use std::fs;
use std::path::PathBuf;

fn fixture(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(path)
}

fn scratch_prefix(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("cf3-rs-{}-{name}", std::process::id()))
}

/// Expected unitigs for `data/compat-small.fa` at k=7, shared by the
/// default-threads test and the low-worker sweep.
const COMPAT_SMALL_UNITIGS: &[&str] = &[
    "AACGGTAA",
    "AACGGTCAA",
    "ACCGGTAA",
    "ACCGGTCCGATGGA",
    "ACCGGTTA",
    "ACCGTTAACCG",
    "ATCGATCG",
    "CGGTAACC",
    "CGGTTAC",
    "GAATTCGAAAGCTGGATCC",
    "GTAACCTGGA",
    "TACTGACATGCAACGTACTGA",
    "TAGCCTTAGCCTACTGA",
];

fn normalized_uncolored_fixture(
    fixture_path: &str,
    name: &str,
    expected_unitigs: &[&str],
    expected_bases: u64,
) {
    normalized_uncolored_fixture_with_threads(
        fixture_path,
        name,
        expected_unitigs,
        expected_bases,
        None,
    );
}

/// Builds a fixture and checks its unitigs, optionally pinning the worker count.
///
/// Pinning matters because the local-contraction sink has two implementations
/// and the thread count picks between them: one worker takes a serial path that
/// a developer machine, sizing workers from its core count, effectively never
/// runs. That path shipped broken -- it derived the unitig file's record index
/// from a counter that also included trivial unitigs, which are written to the
/// FASTA and never to that file, so every read seeked past the records that
/// existed. CI on a two-core runner caught it; nothing else did.
fn normalized_uncolored_fixture_with_threads(
    fixture_path: &str,
    name: &str,
    expected_unitigs: &[&str],
    expected_bases: u64,
    threads: Option<usize>,
) {
    normalized_uncolored_fixture_tuned(
        fixture_path,
        name,
        expected_unitigs,
        expected_bases,
        threads,
        DEFAULT_SUBGRAPH_COUNT,
        true,
    );
}

/// As above, but also choosing how many buckets the partitioner may use.
///
/// The bucket count is the data-driven lever on `workers`: local contraction
/// takes `threads.min(groups.len())`, and a group is a non-empty bucket, so
/// one bucket pins it to a single worker however many threads are asked for.
fn normalized_uncolored_fixture_tuned(
    fixture_path: &str,
    name: &str,
    expected_unitigs: &[&str],
    expected_bases: u64,
    threads: Option<usize>,
    subgraph_count: usize,
    expect_discontinuity_edges: bool,
) {
    let output_prefix = scratch_prefix(name);
    let bucket_dir = output_prefix
        .parent()
        .unwrap()
        .join(format!("{name}.cf3rs.wsk"));
    let fasta_path = output_prefix.with_extension("fa");
    let _ = fs::remove_dir_all(&bucket_dir);
    let _ = fs::remove_file(&fasta_path);

    let mut params = BuildParams::new(
        GraphInput::References,
        output_prefix.with_extension("").display().to_string(),
    );
    params.k = 7;
    params.minimizer_len = 3;
    params
        .seqs
        .push(fixture(fixture_path).display().to_string());
    params.work_dir = output_prefix.parent().unwrap().display().to_string();
    if let Some(threads) = threads {
        params.threads = threads;
    }

    let emitted = emit_weak_superkmer_buckets::<7>(&params, subgraph_count).unwrap();
    let built = build_uncolored_from_buckets::<7>(&params, &emitted.buckets.bucket_dir).unwrap();
    let actual = normalized_fasta_labels(&built.output_path);
    let mut expected = expected_unitigs
        .iter()
        .map(|label| canonical_label(label.as_bytes()))
        .collect::<Vec<_>>();
    expected.sort_unstable();

    assert_eq!(actual, expected);
    assert_eq!(built.unitigs, expected_unitigs.len() as u64);
    assert_eq!(built.unitig_bases, expected_bases);
    assert_eq!(built.bucket_records, emitted.partition.weak_superkmers);
    assert_eq!(built.observed_edges > 0, expect_discontinuity_edges);
    assert_eq!(built.retained_edges > 0, expect_discontinuity_edges);

    let _ = fs::remove_dir_all(emitted.buckets.bucket_dir);
    let _ = fs::remove_file(built.output_path);
}

/// As `normalized_uncolored_fixture`, but at an arbitrary K.
///
/// The k <= 31 pipeline runs on one-word k-mers; K > 31 flips the local
/// vertex map from the flat table to the two-word hash map and widens every
/// k-mer container, so the large-K arms need their own end-to-end coverage.
/// Expectations come from the C++ implementation run on the same fixture.
/// As `normalized_uncolored_fixture_k`, but for a read graph.
///
/// Read mode is a distinct pipeline configuration -- `GraphInput::Reads`
/// carries its own `(k + 1)`-mer cutoff, which discards edges seen fewer than
/// `cutoff` times -- and until this existed nothing in CI exercised it at any
/// k. Expectations come from the C++ implementation, run at two threads and
/// confirmed identical over three repeats, since it is only trustworthy as an
/// oracle at low thread counts.
fn normalized_read_fixture_k<const K: usize>(
    fixture_path: &str,
    name: &str,
    minimizer_len: u16,
    expected_unitigs: &[&str],
    expected_bases: u64,
) {
    let output_prefix = scratch_prefix(name);
    let bucket_dir = output_prefix
        .parent()
        .unwrap()
        .join(format!("{name}.cf3rs.wsk"));
    let fasta_path = output_prefix.with_extension("fa");
    let _ = fs::remove_dir_all(&bucket_dir);
    let _ = fs::remove_file(&fasta_path);

    let mut params = BuildParams::new(
        GraphInput::Reads,
        output_prefix.with_extension("").display().to_string(),
    );
    params.k = K as u16;
    params.minimizer_len = minimizer_len;
    params
        .seqs
        .push(fixture(fixture_path).display().to_string());
    params.work_dir = output_prefix.parent().unwrap().display().to_string();
    // The read default, stated rather than inherited: this fixture exists to
    // check that a depth-one read is discarded, which only happens above 1.
    assert_eq!(params.cutoff(), 2);

    let emitted = emit_weak_superkmer_buckets::<K>(&params, DEFAULT_SUBGRAPH_COUNT).unwrap();
    let built = build_uncolored_from_buckets::<K>(&params, &emitted.buckets.bucket_dir).unwrap();
    let actual = normalized_fasta_labels(&built.output_path);
    let mut expected = expected_unitigs
        .iter()
        .map(|label| canonical_label(label.as_bytes()))
        .collect::<Vec<_>>();
    expected.sort_unstable();

    assert_eq!(actual, expected);
    assert_eq!(built.unitigs, expected_unitigs.len() as u64);
    assert_eq!(built.unitig_bases, expected_bases);

    let _ = fs::remove_dir_all(emitted.buckets.bucket_dir);
    let _ = fs::remove_file(built.output_path);
}

fn normalized_uncolored_fixture_k<const K: usize>(
    fixture_path: &str,
    name: &str,
    minimizer_len: u16,
    expected_unitigs: &[&str],
    expected_bases: u64,
) {
    let output_prefix = scratch_prefix(name);
    let bucket_dir = output_prefix
        .parent()
        .unwrap()
        .join(format!("{name}.cf3rs.wsk"));
    let fasta_path = output_prefix.with_extension("fa");
    let _ = fs::remove_dir_all(&bucket_dir);
    let _ = fs::remove_file(&fasta_path);

    let mut params = BuildParams::new(
        GraphInput::References,
        output_prefix.with_extension("").display().to_string(),
    );
    params.k = K as u16;
    params.minimizer_len = minimizer_len;
    params
        .seqs
        .push(fixture(fixture_path).display().to_string());
    params.work_dir = output_prefix.parent().unwrap().display().to_string();

    let emitted = emit_weak_superkmer_buckets::<K>(&params, DEFAULT_SUBGRAPH_COUNT).unwrap();
    let built = build_uncolored_from_buckets::<K>(&params, &emitted.buckets.bucket_dir).unwrap();
    let actual = normalized_fasta_labels(&built.output_path);
    let mut expected = expected_unitigs
        .iter()
        .map(|label| canonical_label(label.as_bytes()))
        .collect::<Vec<_>>();
    expected.sort_unstable();

    assert_eq!(actual, expected);
    assert_eq!(built.unitigs, expected_unitigs.len() as u64);
    assert_eq!(built.unitig_bases, expected_bases);
    assert_eq!(built.bucket_records, emitted.partition.weak_superkmers);

    let _ = fs::remove_dir_all(emitted.buckets.bucket_dir);
    let _ = fs::remove_file(built.output_path);
}

/// Colored local runs resolve, across the low end of the worker range.
///
/// Colored has its own local-contraction sink and its own single-worker
/// branches, and it reads every bucket a second time to collect source
/// sets, so it is worth sweeping for the same reason the uncolored path is.
#[test]
fn colored_bucket_local_runs_resolve_after_external_collation() {
    for threads in [1, 2, 3] {
        colored_bucket_local_runs_at(threads);
    }
}

fn colored_bucket_local_runs_at(threads: usize) {
    let root = scratch_prefix(&format!("colored-external-local-t{threads}"));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let first = root.join("first.fa");
    let second = root.join("second.fa");
    fs::write(&first, b">first\nAACCGGTTAACCGGTT\n").unwrap();
    fs::write(&second, b">second\nAACCGGTTAACCTTAA\n").unwrap();

    let mut params = BuildParams::new(
        GraphInput::References,
        root.join("colored").display().to_string(),
    );
    params.k = 7;
    params.minimizer_len = 3;
    params.color = true;
    params.threads = threads;
    params.work_dir = root.display().to_string();
    params.seqs = vec![first.display().to_string(), second.display().to_string()];
    let emitted = emit_weak_superkmer_buckets::<7>(&params, DEFAULT_SUBGRAPH_COUNT).unwrap();
    let mut inputs = emit_colored_external_discontinuity_inputs_with_threads_in_dir::<7>(
        &emitted.buckets.bucket_dir,
        1,
        threads,
        root.join("local.labels"),
        root.join("local.colors"),
        root.join("local.color-repository"),
        3,
    )
    .unwrap();
    let repository = inputs.color_repository().unwrap().clone();
    let output_path = root.join("colored.fa");
    SerialUncoloredCollator::collate_external_stitched_to_fasta_with_threads_in_dir(
        &mut inputs,
        threads,
        &root.join("stitch-coords"),
        &root.join("final-unitigs"),
        &output_path,
    )
    .unwrap();

    let fasta = fs::read_to_string(output_path).unwrap();
    let mut runs = 0usize;
    for header in fasta.lines().filter(|line| line.starts_with('>')) {
        assert_eq!(header.split_whitespace().next(), Some(">0"));
        for encoded_run in header.split_whitespace().skip(1) {
            let raw = encoded_run.parse::<u64>().unwrap();
            let coordinate = cuttlefish_rs::state::ColorCoordinate::from_u40(raw >> 24);
            let sources = repository.read_color(coordinate).unwrap();
            assert!(!sources.is_empty());
            assert!(sources.windows(2).all(|pair| pair[0] < pair[1]));
            assert!(sources.iter().all(|&source| source == 1 || source == 2));
            runs += 1;
        }
    }
    assert!(runs > 0);
    fs::remove_dir_all(root).unwrap();
}

/// Every vertex's colors, read back through the repository reader, must equal
/// the set of sources that actually contain that k-mer.
///
/// This is the check `colors dump` rests on, at the level below the CLI: the
/// reader has to agree with the encoder about coordinates, and the FASTA
/// header's packed runs have to name the right vertex ranges. Ground truth
/// here is the input itself, not another implementation.
#[test]
fn colored_runs_resolve_to_the_sources_that_contain_each_kmer() {
    colored_run_sources_at::<7>(
        "colored-run-sources-k7",
        3,
        &[
            b"AACCGGTTAACCGGTTACGT",
            b"AACCGGTTAACCTTAAGGCC",
            b"TTAACCGGTTAACCGGTTAA",
        ],
    );
}

/// The same guarantee above the one-word k-mer boundary.
///
/// Colored and K > 31 are separate configurations, and nothing in CI crossed
/// them: the wide arms' three fatal defects survived exactly because every
/// automated test was either colored at k = 7 or uncolored at k = 33/63. The
/// sources share a long core so classes of more than one source actually
/// arise at k = 63, which needs sequences far longer than the k = 7 fixture's.
#[test]
fn colored_runs_resolve_to_the_sources_at_k33() {
    colored_run_sources_at::<33>("colored-run-sources-k33", 3, &COLORED_WIDE_SOURCES);
}

#[test]
fn colored_runs_resolve_to_the_sources_at_k63() {
    colored_run_sources_at::<63>("colored-run-sources-k63", 3, &COLORED_WIDE_SOURCES);
}

/// Three sources over a shared 150-base core with distinct flanks, so every
/// class -- core (all three) and flank (one each) -- is exercised at k = 63.
const COLORED_WIDE_SOURCES: [&[u8]; 3] = [
    b"GCTTAAGGCACCTTGAACGTGCATCGGTAACCTTGGCAAGTCCATGAACGGTTACCAGGTCATGCTAAGGCCATTGAACCGGTTAAGCCTTGGAACCTTAGGCCAATTGGCCTTAAGGTTCCAAGGTTCCAAGGTTCCAATTGGCCAAGGTTAACCGGTTAACCGGAATTCCGGAATTCC",
    b"TTGGCCAATTGGCCAATTGGCCAAGGTTCCAAGCTTAAGGCACCTTGAACGTGCATCGGTAACCTTGGCAAGTCCATGAACGGTTACCAGGTCATGCTAAGGCCATTGAACCGGTTAAGCCTTGGAACCTTAGGCCAATTGGCCTTAAGGTTCCAAGGTTCCAAGGTTCCAA",
    b"AACCGGTTAACCGGTTAACCGGGCTTAAGGCACCTTGAACGTGCATCGGTAACCTTGGCAAGTCCATGAACGGTTACCAGGTCATGCTAAGGCCATTGAACCGGTTAAGCCTTGGAACCTTAGGCCAATTGGCCTTAAGGTTCCAAGGTTCCAAGGTTCCAATTAAGGCCTTAA",
];

/// Builds a colored graph from `sequences` and checks every vertex's colours
/// against the sources that actually contain that k-mer.
fn colored_run_sources_at<const K: usize>(name: &str, threads: usize, sequences: &[&[u8]]) {
    let root = scratch_prefix(name);
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let mut seqs = Vec::new();
    for (index, sequence) in sequences.iter().enumerate() {
        let path = root.join(format!("source{index}.fa"));
        let mut record = format!(">source{index}\n").into_bytes();
        record.extend_from_slice(sequence);
        record.push(b'\n');
        fs::write(&path, record).unwrap();
        seqs.push(path.display().to_string());
    }

    let mut params = BuildParams::new(
        GraphInput::References,
        root.join("colored").display().to_string(),
    );
    params.k = K as u16;
    // Small enough for the k = 7 case, and the wide cases keep their sources
    // long so the minimizer still has room.
    params.minimizer_len = 3;
    params.color = true;
    params.threads = threads;
    params.work_dir = root.display().to_string();
    params.seqs = seqs;
    let emitted = emit_weak_superkmer_buckets::<K>(&params, DEFAULT_SUBGRAPH_COUNT).unwrap();
    let mut inputs = emit_colored_external_discontinuity_inputs_with_threads_in_dir::<K>(
        &emitted.buckets.bucket_dir,
        1,
        params.threads,
        root.join("local.labels"),
        root.join("local.colors"),
        root.join("local.color-repository"),
        // Source IDs are one-based, so the colour alphabet is one wider than
        // the source count -- what `colored::build` passes. Getting this wrong
        // silently mis-encodes sets containing the highest source id.
        u32::try_from(sequences.len() + 1).unwrap(),
    )
    .unwrap();
    let repository = inputs.color_repository().unwrap().clone();
    let output_path = root.join("colored.fa");
    SerialUncoloredCollator::collate_external_stitched_to_fasta_with_threads_in_dir(
        &mut inputs,
        params.threads,
        &root.join("stitch-coords"),
        &root.join("final-unitigs"),
        &output_path,
    )
    .unwrap();
    repository
        .write_metadata(
            K as u16,
            &output_path,
            &sequences
                .iter()
                .enumerate()
                .map(|(index, _)| root.join(format!("source{index}.fa")))
                .collect::<Vec<_>>(),
        )
        .unwrap();

    // Which sources contain each canonical k-mer, straight from the inputs.
    let mut truth = std::collections::HashMap::<Vec<u8>, std::collections::BTreeSet<u32>>::new();
    for (index, sequence) in sequences.iter().enumerate() {
        for window in sequence.windows(K) {
            truth
                .entry(canonical_label(window))
                .or_default()
                .insert(index as u32 + 1);
        }
    }

    let mut reader =
        cuttlefish_rs::color::ColorRepositoryReader::open(root.join("local.color-repository"))
            .unwrap();
    assert_eq!(reader.k(), K as u16);
    assert_eq!(reader.sources().len(), sequences.len());

    let fasta = fs::read_to_string(&output_path).unwrap();
    let mut lines = fasta.lines();
    let mut vertices_checked = 0usize;
    while let (Some(header), Some(label)) = (lines.next(), lines.next()) {
        let runs = header
            .split_whitespace()
            .skip(1)
            .map(|field| {
                let packed = field.parse::<u64>().unwrap();
                ((packed & 0xFF_FFFF) as usize, packed >> 24)
            })
            .collect::<Vec<_>>();
        assert!(!runs.is_empty());
        let vertex_count = label.len() - K + 1;
        for (index, &(start, coordinate)) in runs.iter().enumerate() {
            let end = runs.get(index + 1).map_or(vertex_count, |next| next.0);
            let sources = reader
                .read_color(cuttlefish_rs::state::ColorCoordinate::from_u40(coordinate))
                .unwrap();
            for vertex in start..end {
                let kmer = canonical_label(&label.as_bytes()[vertex..vertex + K]);
                let expected = truth
                    .get(&kmer)
                    .unwrap_or_else(|| panic!("emitted a k-mer no source contains"));
                assert_eq!(
                    sources,
                    expected.iter().copied().collect::<Vec<_>>(),
                    "vertex {vertex} of {label}",
                );
                vertices_checked += 1;
            }
        }
    }
    assert_eq!(vertices_checked, truth.len());

    fs::remove_dir_all(root).unwrap();
}

fn normalized_fasta_labels(path: &std::path::Path) -> Vec<Vec<u8>> {
    let mut labels = Vec::new();
    let mut current = Vec::new();
    for line in fs::read_to_string(path).unwrap().lines() {
        if line.starts_with('>') {
            if !current.is_empty() {
                labels.push(canonical_label(&current));
                current.clear();
            }
        } else if !line.is_empty() {
            current.extend_from_slice(line.as_bytes());
        }
    }
    if !current.is_empty() {
        labels.push(canonical_label(&current));
    }
    labels.sort_unstable();
    labels
}

fn canonical_label(label: &[u8]) -> Vec<u8> {
    let rc = label
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => panic!("invalid DNA base in FASTA"),
        })
        .collect::<Vec<_>>();
    if rc.as_slice() < label {
        rc
    } else {
        label.to_vec()
    }
}

fn endpoint<const K: usize>(label: &[u8], side: Side) -> DiscontinuityEndpoint<K> {
    DiscontinuityEndpoint {
        vertex: Kmer::<K>::from_ascii(label).unwrap(),
        side,
    }
}

fn endpoint_in_partition<const K: usize>(
    matrix: &SerialEdgeMatrix<K>,
    partition: usize,
    ordinal: usize,
    side: Side,
) -> DiscontinuityEndpoint<K> {
    let mut found = 0usize;
    for encoded in 0..(1usize << (2 * K)) {
        let mut label = vec![b'A'; K];
        let mut value = encoded;
        for base in label.iter_mut().rev() {
            *base = b"ACGT"[value & 3];
            value >>= 2;
        }
        let candidate = endpoint::<K>(&label, side);
        if matrix.partition(MatrixEndpoint::Vertex(candidate)) == partition {
            if found == ordinal {
                return candidate;
            }
            found += 1;
        }
    }
    panic!("no k-mer found for partition {partition} ordinal {ordinal}")
}

#[test]
fn uncolored_refs2_matches_cpp_validated_unitigs() {
    normalized_uncolored_fixture(
        "data/refs2.fa",
        "uncolored-refs2-expected",
        &["AAGTATAATAGACATG", "ACATGTC", "CATGTCG"],
        30,
    );
}

#[test]
fn uncolored_compat_small_matches_cpp_validated_unitigs() {
    normalized_uncolored_fixture(
        "data/compat-small.fa",
        "uncolored-compat-small-expected",
        COMPAT_SMALL_UNITIGS,
        148,
    );
}

/// One bucket: `workers` collapses to 1 by input shape, not by `--threads`.
///
/// This is the route a thread floor could not have closed, since `workers` is
/// `threads.min(groups.len())` and a group is a non-empty bucket. Forcing the
/// partitioner to a single bucket reproduces it deterministically at a thread
/// count that would otherwise take the concurrent path.
///
/// It does *not* reach the code the single-worker bug lived in, and cannot:
/// a discontinuity vertex is one whose neighbour is absent from its subgraph,
/// so once every weak super-k-mer shares a bucket the junctions between them
/// are ordinary vertices and no unitig has an exit. Everything is then trivial
/// and goes straight to the FASTA, never to the unitig file whose index was
/// wrong. One group and a non-empty unitig file are mutually exclusive, which
/// is why `--threads 1` remains the only way to exercise that path and why the
/// sweep below is the test that guards it.
///
/// What this does cover is that the degenerate shape -- one subgraph, no
/// discontinuity edges, an empty unitig file and an empty edge matrix -- still
/// produces the right FASTA rather than dividing by zero somewhere.
#[test]
fn uncolored_single_bucket_collapses_workers_without_discontinuity_edges() {
    normalized_uncolored_fixture_tuned(
        "data/compat-small.fa",
        "uncolored-compat-small-one-bucket",
        COMPAT_SMALL_UNITIGS,
        148,
        Some(8),
        1,
        false,
    );
}

/// The same fixture across the low end of the worker range.
///
/// The worker count selects between two local-contraction sinks, and the
/// serial one is unreachable on a machine with cores to spare, so it went
/// wrong unnoticed until a two-core CI runner tripped over it: it indexed the
/// unitig file with a counter that also included trivial unitigs, which are
/// written to the FASTA and never to that file.
///
/// Sweeping rather than pinning one value because `workers` is derived, not
/// given -- `threads.min(groups.len())` here, and elsewhere from block counts,
/// partition counts and bucket counts -- so the interesting boundaries sit at
/// the bottom of the range wherever the data is small.
#[test]
fn uncolored_compat_small_matches_across_low_worker_counts() {
    for threads in [1, 2, 3, 4] {
        normalized_uncolored_fixture_with_threads(
            "data/compat-small.fa",
            &format!("uncolored-compat-small-t{threads}"),
            COMPAT_SMALL_UNITIGS,
            148,
            Some(threads),
        );
    }
}

#[test]
fn uncolored_generated_synthetic_matches_cpp_validated_unitigs() {
    normalized_uncolored_fixture(
        "data/generated/uncolored-synthetic.fa",
        "uncolored-generated-synthetic-expected",
        &[
            "AACGTACTGACATGCAACGTA",
            "ACCGGTAA",
            "ACCGGTTAA",
            "ACGGTAATCTGTAATC",
            "ATCGATCG",
            "GGTTAAC",
            "GGTTAAGTTAAC",
            "TACCGTAAGGCTA",
            "TACCGTCCGGCTA",
        ],
        107,
    );
}

/// The fixture for the two large-K tests below: four ~200-300 bp records
/// built around a shared 120-base core so the graph branches even at k = 63,
/// plus an embedded `N` (fragment split) and a reverse-complement repeat.
/// Expectations were generated by the C++ implementation:
/// `cuttlefish build --ref -s data/generated/k-large-synthetic.fa -k K
/// --min-len 12`.
const K_LARGE_FIXTURE: &str = "data/generated/k-large-synthetic.fa";

/// Overlapping reads over a 200-base core, plus one depth-one read the read
/// cutoff must discard. Expectations were generated by the C++ implementation:
/// `cuttlefish build --read -s data/generated/reads-synthetic.fq -k K
/// --min-len 5`, at two threads and identical over three repeats.
const READ_FIXTURE: &str = "data/generated/reads-synthetic.fq";

/// Reads at k = 7: the mode's own pipeline, and its cutoff.
///
/// The fixture tiles a 200-base core with overlapping reads so its k-mers are
/// seen repeatedly, and appends one depth-one read. At the read default cutoff
/// of 2 that read contributes nothing, which is what distinguishes this from
/// the reference pipeline over the same bases.
#[test]
fn read_mode_matches_cpp_validated_unitigs() {
    normalized_read_fixture_k::<7>(
        READ_FIXTURE,
        "read-mode-k7-expected",
        5,
        &[
            "CTTATAA",
            "CACCAGCAA",
            "ACTTATA",
            "CTGGTGGCAAGCTCGGAGCGAGATGATCTACCTGTCTGATTACGCGTCGCAATTAAGTCCAAGATAACACATAGCGCGTCATACGTCTGGTGATTGAAGCTTTAGGCCACTCTTTATGATTGTGAACTTTGCTTTTTTA",
            "CTGGTGCAGGTTGATGCGTAGTCTCTGAATTGTTCTTCGGGCCTTATA",
            "CAGCAAACCATTGCTG",
        ],
        226,
    );
}

/// Reads above the one-word k-mer boundary.
///
/// Read mode and K > 31 are separate configurations that had never been
/// exercised together; the wide arms carried three independently fatal defects
/// precisely because no test crossed a pipeline mode with a width.
#[test]
fn read_mode_k33_matches_cpp_validated_unitigs() {
    normalized_read_fixture_k::<33>(
        READ_FIXTURE,
        "read-mode-k33-expected",
        5,
        &[
            "ACTTATAAGGCCCGAAGAACAATTCAGAGACTACGCATCAACCTGCACCAGCAATGGTTTGCTGGTGGCAAGCTCGGAGCGAGATGATCTACCTGTCTGATTACGCGTCGCAATTAAGTCCAAGATAACACATAGCGCGTCATACGTCTGGTGATTGAAGCTTTAGGCCACTCTTTATGATTGTGAACTTTGCTTTTTTA",
        ],
        200,
    );
}

#[test]
fn read_mode_k63_matches_cpp_validated_unitigs() {
    normalized_read_fixture_k::<63>(
        READ_FIXTURE,
        "read-mode-k63-expected",
        5,
        &[
            "ACTTATAAGGCCCGAAGAACAATTCAGAGACTACGCATCAACCTGCACCAGCAATGGTTTGCTGGTGGCAAGCTCGGAGCGAGATGATCTACCTGTCTGATTACGCGTCGCAATTAAGTCCAAGATAACACATAGCGCGTCATACGTCTGGTGATTGAAGCTTTAGGCCACTCTTTATGATTGTGAACTTTGCTTTTTTA",
        ],
        200,
    );
}

#[test]
fn uncolored_k33_matches_cpp_validated_unitigs() {
    normalized_uncolored_fixture_k::<33>(
        K_LARGE_FIXTURE,
        "uncolored-k33-expected",
        12,
        &[
            "GCCTAAACCTATTTGAAGGAGTCTAGCAGCCGACTATGTGTTCCGCAAGAATCAACAACTACAA",
            "CTTCTGTGAACTTCGTGCCAGTAGGGTGGTTTGCTTATCACACCAGCACTGAGATTATGAGACATAACGA",
            "ACTATGTGTTCCGCAAGAATCAACAACTACAATGGCGCGTCGTGAATAACGCGACGGCTGAGACGAACGGCGCGTGAATGAAGCGCTTAA",
            "CGTGCCCACCACACCGAGCTTAAAGAAGTGTCCCCCGAGCATCAAGACTGCCAAACGAGACATAATCTCA",
            "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTATTTGAAGGAGTCTAGCAGCCG",
            "TAGCGTACCGTCAAATGTGGGTCGGTTCACAAAGATGTCGGTTAATTACTGCAAACTTAGCTTCGTCGCGCGTCCAGGGTACTATGTGTTCCGCAAGAATCAACAACTACAA",
            "CAGTAAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTTTAAATGACCCTCTCGTCATAAAACCTTTCTAAGCCCAATAAACCACTCTGACTGGCCGAATA",
            "GAACCTAATAAGTAACAATGGAGGGTCCCGCTACGACGCAAGAAACAGATGTAGGATTGGAAAGAATAGTTGGGCGTCCTGAGAAACAGCCGGCTGCTAGACTCCTTCAAATAGGTTTAGGC",
            "ACAGCTCAGGAGCCAGTCCCCTACGTCGCATATCCTGGCCACTGGAGGTGAAGCGAATGGTATCGATACGTAGGAGGTGTGCCTTCGTAGAAGCCCAATAAACCACTCTGACTGGCCGAATA",
        ],
        892,
    );
}

#[test]
fn uncolored_k63_matches_cpp_validated_unitigs() {
    normalized_uncolored_fixture_k::<63>(
        K_LARGE_FIXTURE,
        "uncolored-k63-expected",
        12,
        &[
            "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTATTTGAAGGAGTCTAGCAGCCG",
            "CTCAGCCGTCGCGTTATTCACGACGCGCCATTGTAGTTGTTGATTCTTGCGGAACACATAGTACCCTGGACGCGCGACGAAGCTAAGTTTGCAGTAATTAACCGACATCTTTGTGAACCGACCCACATTTGACGGTACGCTA",
            "GAACCTAATAAGTAACAATGGAGGGTCCCGCTACGACGCAAGAAACAGATGTAGGATTGGAAAGAATAGTTGGGCGTCCTGAGAAACAGCCGGCTGCTAGACTCCTTCAAATAGGTTTAGGCAACGGCGAAAGCGTCACTGTCGCAAGGGTC",
            "ACAGCTCAGGAGCCAGTCCCCTACGTCGCATATCCTGGCCACTGGAGGTGAAGCGAATGGTATCGATACGTAGGAGGTGTGCCTTCGTAGAAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACC",
            "CGTGCCCACCACACCGAGCTTAAAGAAGTGTCCCCCGAGCATCAAGACTGCCAAACGAGACATAATCTCA",
            "CTTCTGTGAACTTCGTGCCAGTAGGGTGGTTTGCTTATCACACCAGCACTGAGATTATGAGACATAACGA",
            "ACTATGTGTTCCGCAAGAATCAACAACTACAATGGCGCGTCGTGAATAACGCGACGGCTGAGACGAACGGCGCGTGAATGAAGCGCTTAA",
            "CTCAGCCGTCGCGTTATTCACGACGCGCCATTGTAGTTGTTGATTCTTGCGGAACACATAGTCGGCTGCTAGACTCCTTCAAATAGGTTTAGGCAACGGCGAAAGCGTCACTGTCGCAAGGGTC",
            "CAGTAAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTTTAAATGACCCTCTCGTCATAAAACCTTTCTAAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACC",
        ],
        1072,
    );
}

#[test]
fn canonical_minimizer_is_strand_invariant_for_fixture_sequence() {
    let kmer = Kmer::<31>::from_ascii(b"ACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
    assert_eq!(
        canonical_minimizer(kmer, 12, 0),
        canonical_minimizer(kmer.reverse_complement(), 12, 0)
    );
}

#[test]
fn vertex_state_matches_cutoff_model() {
    let mut state = VertexState::<()>::default();
    state.update_edges(Base::E, Base::A);
    state.update_edges(Base::E, Base::A);
    state.update_edges(Base::E, Base::C);

    assert_eq!(state.edge_at(Side::Back, 1), Base::N);
    assert_eq!(state.edge_at(Side::Back, 2), Base::A);
    assert!(!state.is_branching_side(Side::Back, 2));

    state.mark_discontinuous(Side::Front);
    assert!(state.is_discontinuous(Side::Front));
    assert!(!state.is_discontinuous(Side::Back));
}

#[test]
fn expands_input_and_parses_reference_fragments() {
    let mut params = BuildParams::new(GraphInput::References, "out".to_string());
    params
        .seqs
        .push(fixture("data/refs1.fa").display().to_string());

    let paths = expand_input_paths(&params).unwrap();
    assert_eq!(paths.len(), 1);

    let mut fragments = Vec::new();
    let records = parse_fragments(&paths[0], 1, 3, |fragment| {
        fragments.push(fragment);
        Ok(())
    })
    .unwrap();

    assert!(records >= 1);
    assert!(fragments.iter().all(|fragment| fragment.source_id == 1));
    assert!(fragments.iter().all(|fragment| fragment.seq.len() >= 3));
}

#[test]
fn partitions_reference_fragments_into_weak_superkmers() {
    let mut fragments = Vec::new();
    parse_fragments(fixture("data/refs2.fa"), 1, 8, |fragment| {
        fragments.push(fragment);
        Ok(())
    })
    .unwrap();

    let mut superkmers = Vec::new();
    for fragment in fragments {
        superkmers.extend(
            partition_fragment::<7>(&fragment.seq, 3, DEFAULT_SUBGRAPH_COUNT, Some(1)).unwrap(),
        );
    }

    assert!(!superkmers.is_empty());
    assert!(superkmers.iter().all(|sk| sk.len >= 7));
    assert!(
        superkmers
            .iter()
            .all(|sk| sk.graph_id < DEFAULT_SUBGRAPH_COUNT)
    );
}

#[test]
fn partition_phase_reports_fixture_stats() {
    let mut params = BuildParams::new(GraphInput::References, "out".to_string());
    params.k = 7;
    params.minimizer_len = 3;
    params
        .seqs
        .push(fixture("data/refs2.fa").display().to_string());

    let stats = partition_inputs::<7>(&params, DEFAULT_SUBGRAPH_COUNT).unwrap();
    assert_eq!(stats.input_files, 1);
    assert!(stats.records >= 1);
    assert!(stats.fragments >= 1);
    assert!(stats.weak_superkmers >= 1);
    assert!(stats.non_empty_graphs() >= 1);
    assert_eq!(
        stats.graph_histogram.iter().sum::<u64>(),
        stats.weak_superkmers
    );
}

#[test]
fn emits_external_memory_weak_superkmer_buckets() {
    let output_prefix = scratch_prefix("bucket-refs2");
    let bucket_dir = output_prefix
        .parent()
        .unwrap()
        .join("bucket-refs2.cf3rs.wsk");
    let _ = fs::remove_dir_all(&bucket_dir);

    let mut params = BuildParams::new(
        GraphInput::References,
        output_prefix.with_extension("").display().to_string(),
    );
    params.k = 7;
    params.minimizer_len = 3;
    params.compress_buckets = true;
    params
        .seqs
        .push(fixture("data/refs2.fa").display().to_string());
    params.work_dir = output_prefix.parent().unwrap().display().to_string();

    let emitted = emit_weak_superkmer_buckets::<7>(&params, DEFAULT_SUBGRAPH_COUNT).unwrap();
    assert_eq!(emitted.partition.input_files, 1);
    assert!(emitted.partition.weak_superkmers > 0);
    assert!(emitted.buckets.bucket_files > 0);
    assert!(emitted.buckets.bytes_written > 0);

    // Buckets share container files, so the directory holds one file per atlas
    // plus the manifest -- not one file per bucket -- and the per-bucket header
    // that a whole file carried inline now lives in that manifest.
    let (header, manifest) = read_container_manifest(&emitted.buckets.bucket_dir)
        .unwrap()
        .expect("container manifest");
    assert_eq!(header.k, 7);
    assert_eq!(header.minimizer_len, 3);
    assert!(!header.colored);
    assert!(header.compressed);
    assert_eq!(manifest.len(), emitted.buckets.bucket_files);
    // The file count tracks the atlas count, not the bucket count: one
    // container per atlas plus the manifest, whatever DEFAULT_SUBGRAPH_COUNT
    // is. Containers are created eagerly, so a fixture this small writes more
    // files than it fills -- at real bucket counts every one of them is used,
    // and creating them lazily would put a check on the flush path to save
    // 128 empty inodes.
    let physical = fs::read_dir(&emitted.buckets.bucket_dir).unwrap().count();
    assert_eq!(physical, header.container_count + 1);
    assert_eq!(
        header.container_count,
        DEFAULT_SUBGRAPH_COUNT.div_ceil(128),
        "one container per atlas"
    );
    for entry in &manifest {
        let BucketLocation::Container {
            container,
            segments,
            bytes,
        } = &entry.location
        else {
            panic!("expected a container location");
        };
        assert!(*container < header.container_count);
        assert!(!segments.is_empty() && *bytes > 0);
        // The chain holds the payload with less than one segment to spare.
        let capacity = segments.len() as u64 * header.segment_bytes;
        assert!(*bytes <= capacity && capacity - *bytes < header.segment_bytes);
    }

    let (store, entries) = BucketStore::open_dir(&emitted.buckets.bucket_dir).unwrap();
    assert_eq!(entries.len(), emitted.buckets.bucket_files);
    assert_eq!(
        entries.iter().map(|entry| entry.records).sum::<u64>(),
        emitted.partition.weak_superkmers
    );

    let mut decoded_records = 0u64;
    for entry in entries {
        let mut reader = store.reader(&entry).unwrap();
        assert_eq!(reader.header().k, 7);
        assert_eq!(reader.header().minimizer_len, 3);
        assert_eq!(reader.header().graph_id, entry.graph_id);
        assert_eq!(reader.header().records, entry.records);
        assert!(!reader.header().colored);
        assert!(reader.header().compressed);

        while let Some(record) = reader.next_record().unwrap() {
            assert_eq!(record.graph_id, entry.graph_id);
            assert_eq!(record.len, record.label.len());
            assert!(record.len >= 7);
            assert!(
                record
                    .label
                    .iter()
                    .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
            );
            assert_eq!(record.source_id, None);
            decoded_records += 1;
        }
    }
    assert_eq!(decoded_records, emitted.partition.weak_superkmers);

    let _ = fs::remove_dir_all(emitted.buckets.bucket_dir);
}

#[test]
fn constructs_local_subgraphs_from_emitted_buckets() {
    let output_prefix = scratch_prefix("local-subgraph-refs2");
    let bucket_dir = output_prefix
        .parent()
        .unwrap()
        .join("local-subgraph-refs2.cf3rs.wsk");
    let _ = fs::remove_dir_all(&bucket_dir);

    let mut params = BuildParams::new(
        GraphInput::References,
        output_prefix.with_extension("").display().to_string(),
    );
    params.k = 7;
    params.minimizer_len = 3;
    params
        .seqs
        .push(fixture("data/refs2.fa").display().to_string());
    params.work_dir = output_prefix.parent().unwrap().display().to_string();

    let emitted = emit_weak_superkmer_buckets::<7>(&params, DEFAULT_SUBGRAPH_COUNT).unwrap();
    let (store, entries) = BucketStore::open_dir(&emitted.buckets.bucket_dir).unwrap();

    let mut weak_superkmers = 0u64;
    let mut observed_vertices = 0u64;
    let mut observed_edges = 0u64;
    let mut unique_vertices = 0u64;
    for entry in entries {
        let mut subgraph = LocalSubgraph::<7>::from_manifest_entries(
            &store,
            std::slice::from_ref(&entry),
            params.cutoff(),
        )
        .unwrap();
        assert_eq!(subgraph.graph_id, entry.graph_id);
        assert_eq!(subgraph.stats.weak_superkmers, entry.records);
        assert_eq!(
            subgraph.stats.unique_vertices,
            subgraph.vertices.len() as u64
        );
        assert_eq!(subgraph.stats.unique_edges, subgraph.edges.len() as u64);
        assert!(subgraph.stats.observed_vertices >= subgraph.stats.unique_vertices);
        assert!(subgraph.stats.observed_edges >= subgraph.stats.unique_edges);

        weak_superkmers += subgraph.stats.weak_superkmers;
        observed_vertices += subgraph.stats.observed_vertices;
        observed_edges += subgraph.stats.observed_edges;
        unique_vertices += subgraph.stats.unique_vertices;

        let unitigs = subgraph.contract().unwrap();
        assert_eq!(subgraph.stats.unitigs, unitigs.len() as u64);
        assert!(
            unitigs
                .iter()
                .all(|unitig| unitig.label.len() >= params.k as usize)
        );
        assert_eq!(
            subgraph.stats.unitig_bases,
            unitigs
                .iter()
                .map(|unitig| unitig.label.len() as u64)
                .sum::<u64>()
        );
    }

    assert_eq!(weak_superkmers, emitted.partition.weak_superkmers);
    assert!(observed_vertices >= weak_superkmers);
    assert!(observed_edges > 0);
    assert!(unique_vertices > 0);

    let _ = fs::remove_dir_all(emitted.buckets.bucket_dir);
}

#[test]
fn emits_discontinuity_inputs_from_local_subgraphs() {
    let output_prefix = scratch_prefix("discontinuity-inputs-refs2");
    let bucket_dir = output_prefix
        .parent()
        .unwrap()
        .join("discontinuity-inputs-refs2.cf3rs.wsk");
    let _ = fs::remove_dir_all(&bucket_dir);

    let mut params = BuildParams::new(
        GraphInput::References,
        output_prefix.with_extension("").display().to_string(),
    );
    params.k = 7;
    params.minimizer_len = 3;
    params
        .seqs
        .push(fixture("data/refs2.fa").display().to_string());
    params.work_dir = output_prefix.parent().unwrap().display().to_string();

    let emitted = emit_weak_superkmer_buckets::<7>(&params, DEFAULT_SUBGRAPH_COUNT).unwrap();
    let discontinuities =
        emit_uncolored_discontinuity_inputs::<7>(&emitted.buckets.bucket_dir, params.cutoff())
            .unwrap();

    assert_eq!(
        discontinuities.stats.input_buckets,
        emitted.buckets.bucket_files
    );
    assert_eq!(
        discontinuities.stats.weak_superkmers,
        emitted.partition.weak_superkmers
    );
    assert_eq!(
        discontinuities.stats.local_unitigs,
        discontinuities.unitigs.len() as u64
    );
    assert!(discontinuities.stats.local_unitigs > 0);
    assert!(discontinuities.stats.discontinuity_exits > 0);
    assert!(
        discontinuities
            .unitigs
            .iter()
            .all(|unitig| unitig.label(&discontinuities).len() >= params.k as usize)
    );

    let _ = fs::remove_dir_all(emitted.buckets.bucket_dir);
}

#[test]
fn serial_edge_matrix_uses_phi_and_upper_triangular_blocks() {
    let mut matrix = SerialEdgeMatrix::<7>::new(4).unwrap();
    let low = endpoint_in_partition(&matrix, 1, 0, Side::Front);
    let high = endpoint_in_partition(&matrix, 4, 0, Side::Back);

    assert_eq!(matrix.partition(MatrixEndpoint::Phi), 0);
    assert_eq!(matrix.partition(MatrixEndpoint::Vertex(low)), 1);
    assert_eq!(matrix.partition(MatrixEndpoint::Vertex(high)), 4);

    matrix.add_edge(
        MatrixEndpoint::Vertex(high),
        MatrixEndpoint::Vertex(low),
        11,
        3,
    );
    matrix.add_edge(MatrixEndpoint::Phi, MatrixEndpoint::Vertex(high), 5, 4);

    let internal_block = matrix.block(1, 4);
    assert_eq!(internal_block.len(), 1);
    assert!(internal_block[0].swapped);
    assert_eq!(internal_block[0].first, MatrixEndpoint::Vertex(low));
    assert_eq!(internal_block[0].second, MatrixEndpoint::Vertex(high));
    assert_eq!(internal_block[0].weight, 11);

    let phi_block = matrix.block(0, 4);
    assert_eq!(phi_block.len(), 1);
    assert_eq!(phi_block[0].first, MatrixEndpoint::Phi);
    assert_eq!(phi_block[0].second, MatrixEndpoint::Vertex(high));

    let stats = matrix.stats();
    assert_eq!(stats.edges, 2);
    assert_eq!(stats.phi_edges, 1);
    assert_eq!(stats.diagonal_edges, 0);
}

#[test]
fn serial_edge_matrix_uses_unit_weight_for_original_local_unitigs() {
    let left = endpoint::<7>(b"AAAAAAA", Side::Front);
    let right = endpoint::<7>(b"AAAAAAT", Side::Back);
    let inputs = DiscontinuityInputs::from_unitigs(
        vec![OwnedDiscontinuityUnitig {
            graph_id: 0,
            label: b"AAAAAAATTT".to_vec(),
            left_exit: Some(left),
            right_exit: Some(right),
            is_cycle: false,
        }],
        DiscontinuityInputStats::default(),
    );

    let matrix = SerialEdgeMatrix::from_inputs(&inputs, 4).unwrap();
    let edge = matrix.edges().next().unwrap();

    assert_eq!(matrix.stats().edges, 1);
    assert_eq!(edge.weight, 1);
    assert_eq!(edge.unitig_index, 0);
}

#[test]
fn serial_discontinuity_contractor_groups_path_and_cycle_components() {
    let a = endpoint::<7>(b"AAAAAAA", Side::Front);
    let b = endpoint::<7>(b"AAAAAAC", Side::Back);
    let c = endpoint::<7>(b"AAAAAAG", Side::Front);
    let x = endpoint::<7>(b"CCCCCCC", Side::Back);
    let y = endpoint::<7>(b"CCCCCCA", Side::Front);
    let z = endpoint::<7>(b"CCCCCCG", Side::Back);
    let mut matrix = SerialEdgeMatrix::<7>::new(8).unwrap();

    matrix.add_edge(MatrixEndpoint::Vertex(a), MatrixEndpoint::Vertex(b), 2, 0);
    matrix.add_edge(MatrixEndpoint::Vertex(b), MatrixEndpoint::Vertex(c), 3, 1);
    matrix.add_edge(MatrixEndpoint::Phi, MatrixEndpoint::Vertex(a), 5, 2);
    matrix.add_edge(MatrixEndpoint::Vertex(x), MatrixEndpoint::Vertex(y), 7, 3);
    matrix.add_edge(MatrixEndpoint::Vertex(y), MatrixEndpoint::Vertex(z), 11, 4);
    matrix.add_edge(MatrixEndpoint::Vertex(z), MatrixEndpoint::Vertex(x), 13, 5);

    let contracted = SerialDiscontinuityContractor::contract(&matrix);

    assert_eq!(contracted.stats.input_edges, 6);
    assert_eq!(contracted.stats.components, 2);
    assert_eq!(contracted.stats.phi_edges, 1);
    assert_eq!(contracted.stats.cyclic_components, 1);
    assert!(
        contracted.components.iter().any(|component| {
            component.edges == 3 && component.weight == 31 && component.cyclic
        })
    );
    assert!(contracted.components.iter().any(|component| {
        component.edges == 3
            && component.weight == 10
            && component.phi_edges == 1
            && !component.cyclic
    }));
}

#[test]
fn diagonal_compression_collapses_same_partition_chains() {
    let mut matrix = SerialEdgeMatrix::<7>::new(4).unwrap();
    let partition = 2;
    let a = endpoint_in_partition(&matrix, partition, 0, Side::Front);
    let b = endpoint_in_partition(&matrix, partition, 1, Side::Back);
    let c = endpoint_in_partition(&matrix, partition, 2, Side::Front);

    assert_eq!(partition, matrix.partition(MatrixEndpoint::Vertex(b)));
    assert_eq!(partition, matrix.partition(MatrixEndpoint::Vertex(c)));

    matrix.add_edge(MatrixEndpoint::Vertex(a), MatrixEndpoint::Vertex(b), 2, 0);
    matrix.add_edge(MatrixEndpoint::Vertex(b), MatrixEndpoint::Vertex(c), 3, 1);

    let compressed = SerialDiscontinuityContractor::compress_diagonal_block(&matrix, partition);

    assert_eq!(compressed.partition, partition);
    assert_eq!(compressed.stats.input_edges, 2);
    assert_eq!(compressed.stats.compressed_edges, 1);
    assert_eq!(compressed.stats.meta_vertices, 0);
    assert_eq!(compressed.edges.len(), 1);
    assert_eq!(compressed.edges[0].weight, 5);
    assert_eq!(compressed.edges[0].first, MatrixEndpoint::Vertex(a));
    assert_eq!(compressed.edges[0].second, MatrixEndpoint::Vertex(c));
}

#[test]
fn diagonal_compression_records_isolated_cycle_meta_vertices() {
    let mut matrix = SerialEdgeMatrix::<7>::new(4).unwrap();
    let partition = 2;
    let a = endpoint_in_partition(&matrix, partition, 0, Side::Front);
    let b = endpoint_in_partition(&matrix, partition, 1, Side::Back);
    let c = endpoint_in_partition(&matrix, partition, 2, Side::Back);

    matrix.add_edge(MatrixEndpoint::Vertex(a), MatrixEndpoint::Vertex(b), 2, 0);
    matrix.add_edge(MatrixEndpoint::Vertex(b), MatrixEndpoint::Vertex(c), 3, 1);
    matrix.add_edge(MatrixEndpoint::Vertex(c), MatrixEndpoint::Vertex(a), 5, 2);

    let compressed = SerialDiscontinuityContractor::compress_diagonal_block(&matrix, partition);

    assert_eq!(compressed.stats.input_edges, 3);
    assert_eq!(compressed.stats.compressed_edges, 0);
    assert_eq!(compressed.stats.meta_vertices, 1);
    assert_eq!(compressed.stats.isolated_cordless_cycles, 1);
    assert_eq!(compressed.meta_vertices.len(), 1);
    assert_eq!(compressed.meta_vertices[0].vertex, a.vertex);
    assert_eq!(compressed.meta_vertices[0].partition, partition);
    assert!(compressed.meta_vertices[0].is_cycle);
}

#[test]
fn partition_contraction_joins_non_diagonal_edges_through_partition_vertices() {
    let mut matrix = SerialEdgeMatrix::<7>::new(4).unwrap();
    let partition = 4;
    let left = endpoint_in_partition(&matrix, 1, 0, Side::Front);
    let current = endpoint_in_partition(&matrix, partition, 0, Side::Back);
    let right = endpoint_in_partition(&matrix, 2, 0, Side::Front);

    assert!(matrix.partition(MatrixEndpoint::Vertex(left)) < partition);
    assert!(matrix.partition(MatrixEndpoint::Vertex(right)) < partition);

    matrix.add_edge(
        MatrixEndpoint::Vertex(left),
        MatrixEndpoint::Vertex(current),
        7,
        0,
    );
    matrix.add_edge(
        MatrixEndpoint::Vertex(right),
        MatrixEndpoint::Vertex(current),
        11,
        1,
    );

    let contracted = SerialDiscontinuityContractor::contract_partition(&matrix, partition);

    assert_eq!(contracted.partition, partition);
    assert_eq!(contracted.stats.input_non_diagonal_edges, 2);
    assert_eq!(contracted.stats.compressed_diagonal_edges, 0);
    assert_eq!(contracted.stats.phantom_edges, 0);
    assert_eq!(contracted.stats.output_edges, 1);
    assert_eq!(contracted.edges[0].weight, 18);
    assert_eq!(contracted.edges[0].first, MatrixEndpoint::Vertex(right));
    assert_eq!(contracted.edges[0].second, MatrixEndpoint::Vertex(left));
}

#[test]
fn partition_contraction_emits_false_phantoms_for_unpaired_endpoints() {
    let mut matrix = SerialEdgeMatrix::<7>::new(4).unwrap();
    let partition = 4;
    let lower = endpoint_in_partition(&matrix, 1, 0, Side::Front);
    let current = endpoint_in_partition(&matrix, partition, 0, Side::Back);

    matrix.add_edge(
        MatrixEndpoint::Vertex(lower),
        MatrixEndpoint::Vertex(current),
        7,
        0,
    );

    let contracted = SerialDiscontinuityContractor::contract_partition(&matrix, partition);

    assert_eq!(contracted.stats.input_non_diagonal_edges, 1);
    assert_eq!(contracted.stats.phantom_edges, 1);
    assert_eq!(contracted.stats.output_edges, 1);
    assert!(contracted.expansion_edges.iter().any(|edge| {
        let phantom = DiscontinuityEndpoint {
            vertex: current.vertex,
            side: current.side.inverse(),
        };
        edge.weight == 1
            && edge.first == MatrixEndpoint::Phi
            && edge.second == MatrixEndpoint::Vertex(phantom)
            && edge.phantom_unitig == Some(phantom)
    }));
    assert!(contracted.edges.iter().any(|edge| {
        edge.weight == 8
            && edge.first == MatrixEndpoint::Phi
            && edge.second == MatrixEndpoint::Vertex(lower)
    }));
}

#[test]
fn full_serial_contraction_reinserts_edges_for_lower_partitions() {
    let mut matrix = SerialEdgeMatrix::<7>::new(4).unwrap();
    let low = endpoint_in_partition(&matrix, 1, 0, Side::Front);
    let mid = endpoint_in_partition(&matrix, 2, 0, Side::Back);
    let high = endpoint_in_partition(&matrix, 4, 0, Side::Front);

    assert!(
        matrix.partition(MatrixEndpoint::Vertex(low))
            < matrix.partition(MatrixEndpoint::Vertex(high))
    );
    assert!(
        matrix.partition(MatrixEndpoint::Vertex(mid))
            < matrix.partition(MatrixEndpoint::Vertex(high))
    );
    assert_ne!(
        matrix.partition(MatrixEndpoint::Vertex(low)),
        matrix.partition(MatrixEndpoint::Vertex(mid))
    );

    matrix.add_edge(
        MatrixEndpoint::Vertex(low),
        MatrixEndpoint::Vertex(high),
        2,
        0,
    );
    matrix.add_edge(
        MatrixEndpoint::Vertex(mid),
        MatrixEndpoint::Vertex(high),
        3,
        1,
    );

    let contracted = SerialDiscontinuityContractor::contract_all_partitions(&matrix);

    assert_eq!(contracted.stats.input_edges, 2);
    assert_eq!(contracted.stats.reinserted_edges, 2);
    assert_eq!(contracted.stats.phantom_edges, 2);
    assert_eq!(contracted.stats.final_edges, 0);
    assert_eq!(contracted.stats.meta_vertices, 1);
    assert!(contracted.partitions.iter().any(|partition| {
        partition.stats.output_edges == 1 && partition.stats.phantom_edges == 0
    }));
    assert!(contracted.expansion_edges.iter().any(|edge| {
        edge.weight == 1
            && (edge.first == MatrixEndpoint::Phi || edge.second == MatrixEndpoint::Phi)
            && edge.phantom_unitig.is_some()
    }));
    assert!(
        contracted
            .meta_vertices
            .iter()
            .any(|meta| { meta.vertex == low.vertex && meta.weight > 0 && !meta.is_cycle })
    );
}

#[test]
fn serial_expander_infers_path_info_like_cpp_rule() {
    let path_id = Kmer::<7>::from_ascii(b"AAAAAAA").unwrap();
    let source = PathInfo {
        path_id,
        rank: 5,
        exit_side: Side::Back,
        is_cycle: false,
    };

    let forward = SerialDiscontinuityExpander::infer(source, Side::Back, Side::Front, 3);
    assert_eq!(forward.path_id, path_id);
    assert_eq!(forward.rank, 8);
    assert_eq!(forward.exit_side, Side::Back);

    let backward = SerialDiscontinuityExpander::infer(source, Side::Front, Side::Back, 3);
    assert_eq!(backward.path_id, path_id);
    assert_eq!(backward.rank, 2);
    assert_eq!(backward.exit_side, Side::Back);
}

#[test]
fn serial_expansion_seeds_meta_vertices_and_propagates_final_edges() {
    let low = endpoint::<7>(b"AAAAAAA", Side::Front);
    let mid = endpoint::<7>(b"AAAAAAC", Side::Back);
    let edge = DiscontinuityEdge {
        first: MatrixEndpoint::Vertex(low),
        second: MatrixEndpoint::Vertex(mid),
        weight: 1,
        unitig_bucket: 0,
        unitig_index: 42,
        unitig_exit_side: Side::Back,
        phantom_unitig: None,
        swapped: false,
    };
    let mut expansion_matrix = SerialEdgeMatrix::new(128).unwrap();
    expansion_matrix.add_edge_with_orientation(
        edge.first,
        edge.second,
        edge.weight,
        edge.unitig_index,
        edge.unitig_exit_side,
    );
    let contracted = FullSerialDiscontinuityContraction {
        vertex_partitions: 128,
        final_edges: vec![edge.clone()],
        expansion_edges: vec![edge],
        expansion_matrix,
        compressed_diagonal_edges: vec![Vec::new(); 128],
        meta_vertices: vec![SerialMetaVertex {
            vertex: low.vertex,
            partition: 1,
            entry_side: Side::Back,
            weight: 3,
            is_cycle: false,
        }],
        partitions: Vec::new(),
        stats: FullSerialContractionStats {
            final_edges: 1,
            meta_vertices: 1,
            ..FullSerialContractionStats::default()
        },
    };
    let expanded = SerialDiscontinuityExpander::expand(&contracted);

    assert_eq!(expanded.stats.seed_vertices, contracted.stats.meta_vertices);
    assert!(expanded.stats.inferred_vertices > 0);
    assert_eq!(expanded.stats.unresolved_edges, 0);
    assert!(
        expanded
            .vertices
            .iter()
            .any(|entry| entry.vertex == low.vertex)
    );
    assert!(
        expanded
            .edges
            .iter()
            .any(|entry| entry.unitig_index == 42 && entry.info.path_id == low.vertex)
    );
}

#[test]
fn serial_uncolored_collator_joins_ranked_local_unitigs() {
    let path_id = Kmer::<7>::from_ascii(b"AAAAAAA").unwrap();
    let inputs = DiscontinuityInputs::from_unitigs(
        vec![
            OwnedDiscontinuityUnitig {
                graph_id: 0,
                label: b"AAAAAAAC".to_vec(),
                left_exit: None,
                right_exit: None,
                is_cycle: false,
            },
            OwnedDiscontinuityUnitig {
                graph_id: 0,
                label: b"AAAAACGG".to_vec(),
                left_exit: None,
                right_exit: None,
                is_cycle: false,
            },
        ],
        DiscontinuityInputStats::default(),
    );
    let expansion = SerialExpansion {
        vertices: Vec::new(),
        edges: vec![
            EdgePathInfo {
                unitig_index: 0,
                phantom_unitig: None,
                info: PathInfo {
                    path_id,
                    rank: 1,
                    exit_side: Side::Back,
                    is_cycle: false,
                },
            },
            EdgePathInfo {
                unitig_index: 1,
                phantom_unitig: None,
                info: PathInfo {
                    path_id,
                    rank: 2,
                    exit_side: Side::Back,
                    is_cycle: false,
                },
            },
        ],
        stats: SerialExpansionStats::default(),
    };

    let collated = SerialUncoloredCollator::collate(&inputs, &expansion);

    assert_eq!(collated.stats.input_path_infos, 2);
    assert_eq!(collated.stats.missing_unitig_labels, 0);
    assert_eq!(collated.unitigs, vec![b"AAAAAAACG".to_vec()]);
}

#[test]
fn serial_uncolored_collator_emits_completed_local_unitigs() {
    let inputs = DiscontinuityInputs::from_unitigs(
        vec![
            OwnedDiscontinuityUnitig {
                graph_id: 0,
                label: b"CCCCAAA".to_vec(),
                left_exit: None,
                right_exit: None,
                is_cycle: false,
            },
            OwnedDiscontinuityUnitig {
                graph_id: 0,
                label: b"AAAAAAC".to_vec(),
                left_exit: Some(DiscontinuityEndpoint {
                    vertex: Kmer::<7>::from_ascii(b"AAAAAAC").unwrap(),
                    side: Side::Back,
                }),
                right_exit: None,
                is_cycle: false,
            },
        ],
        DiscontinuityInputStats::default(),
    );
    let expansion = SerialExpansion {
        vertices: Vec::new(),
        edges: Vec::new(),
        stats: SerialExpansionStats::default(),
    };

    let collated = SerialUncoloredCollator::collate(&inputs, &expansion);

    assert_eq!(collated.stats.input_path_infos, 0);
    assert_eq!(collated.stats.direct_local_unitigs, 1);
    assert_eq!(collated.stats.missing_unitig_labels, 0);
    assert_eq!(
        collated.unitigs,
        vec![b"AAAAAAC".to_vec(), b"CCCCAAA".to_vec()]
    );
}

#[test]
fn discontinuity_parallelization_opportunities_are_recorded() {
    assert!(
        DISCONTINUITY_PARALLELIZATION_OPPORTUNITIES
            .iter()
            .any(|opportunity| opportunity.contains("diagonal block compression"))
    );
    assert!(
        DISCONTINUITY_PARALLELIZATION_OPPORTUNITIES
            .iter()
            .any(|opportunity| opportunity.contains("false-phantom"))
    );
}

#[test]
fn builds_serial_edge_matrix_from_fixture_discontinuity_inputs() {
    let output_prefix = scratch_prefix("serial-matrix-refs2");
    let bucket_dir = output_prefix
        .parent()
        .unwrap()
        .join("serial-matrix-refs2.cf3rs.wsk");
    let _ = fs::remove_dir_all(&bucket_dir);

    let mut params = BuildParams::new(
        GraphInput::References,
        output_prefix.with_extension("").display().to_string(),
    );
    params.k = 7;
    params.minimizer_len = 3;
    params
        .seqs
        .push(fixture("data/refs2.fa").display().to_string());
    params.work_dir = output_prefix.parent().unwrap().display().to_string();

    let emitted = emit_weak_superkmer_buckets::<7>(&params, DEFAULT_SUBGRAPH_COUNT).unwrap();
    let discontinuities =
        emit_uncolored_discontinuity_inputs::<7>(&emitted.buckets.bucket_dir, params.cutoff())
            .unwrap();
    let expected_edges = discontinuities
        .unitigs
        .iter()
        .filter(|unitig| unitig.left_exit().is_some() || unitig.right_exit().is_some())
        .count() as u64;

    let matrix = SerialEdgeMatrix::from_inputs(&discontinuities, 128).unwrap();
    let contracted = SerialDiscontinuityContractor::contract(&matrix);
    let full_contracted = SerialDiscontinuityContractor::contract_all_partitions(&matrix);
    let expanded = SerialDiscontinuityExpander::expand(&full_contracted);
    let collated = SerialUncoloredCollator::collate(&discontinuities, &expanded);
    let expected_collated = vec![
        b"AAGTATAATAGACATG".to_vec(),
        b"ACATGTC".to_vec(),
        b"CATGTCG".to_vec(),
    ];

    assert_eq!(matrix.vertex_partitions(), 128);
    assert_eq!(matrix.partition_count(), 129);
    assert_eq!(matrix.stats().edges, expected_edges);
    assert!(matrix.stats().edges > 0);
    assert_eq!(contracted.stats.input_edges, matrix.stats().edges);
    assert_eq!(contracted.stats.phi_edges, matrix.stats().phi_edges);
    assert!(contracted.stats.components > 0);
    assert_eq!(full_contracted.stats.partitions, 128);
    assert_eq!(full_contracted.stats.input_edges, matrix.stats().edges);
    assert_eq!(
        expanded.stats.seed_vertices,
        full_contracted.stats.meta_vertices
    );
    assert_eq!(
        collated.stats.input_path_infos,
        expanded.stats.edge_path_infos
    );
    assert_eq!(collated.stats.direct_local_unitigs, 1);
    assert_eq!(collated.stats.stitched_discontinuity_unitigs, 2);
    assert_eq!(collated.unitigs, expected_collated);

    let _ = fs::remove_dir_all(emitted.buckets.bucket_dir);
}

#[test]
fn serial_discontinuity_pipeline_builds_fasta_from_emitted_buckets() {
    let output_prefix = scratch_prefix("serial-discontinuity-refs2");
    let bucket_dir = output_prefix
        .parent()
        .unwrap()
        .join("serial-discontinuity-refs2.cf3rs.wsk");
    let fasta_path = output_prefix.with_extension("fa");
    let _ = fs::remove_dir_all(&bucket_dir);
    let _ = fs::remove_file(&fasta_path);

    let mut params = BuildParams::new(
        GraphInput::References,
        output_prefix.with_extension("").display().to_string(),
    );
    params.k = 7;
    params.minimizer_len = 3;
    params
        .seqs
        .push(fixture("data/refs2.fa").display().to_string());
    params.work_dir = output_prefix.parent().unwrap().display().to_string();

    let emitted = emit_weak_superkmer_buckets::<7>(&params, DEFAULT_SUBGRAPH_COUNT).unwrap();
    let built = build_uncolored_with_serial_discontinuity_pipeline::<7>(
        &params,
        &emitted.buckets.bucket_dir,
    )
    .unwrap();

    assert_eq!(built.input_buckets, emitted.buckets.bucket_files);
    assert_eq!(built.bucket_records, emitted.partition.weak_superkmers);
    assert_eq!(built.output_path, fasta_path);

    let fasta = fs::read_to_string(&built.output_path).unwrap();
    assert_eq!(
        fasta.lines().filter(|line| line.starts_with('>')).count() as u64,
        built.unitigs
    );

    let _ = fs::remove_dir_all(emitted.buckets.bucket_dir);
    let _ = fs::remove_file(built.output_path);
}

#[test]
fn builds_uncolored_fasta_from_emitted_buckets() {
    let output_prefix = scratch_prefix("uncolored-refs2");
    let bucket_dir = output_prefix
        .parent()
        .unwrap()
        .join("uncolored-refs2.cf3rs.wsk");
    let fasta_path = output_prefix.with_extension("fa");
    let _ = fs::remove_dir_all(&bucket_dir);
    let _ = fs::remove_file(&fasta_path);

    let mut params = BuildParams::new(
        GraphInput::References,
        output_prefix.with_extension("").display().to_string(),
    );
    params.k = 7;
    params.minimizer_len = 3;
    params
        .seqs
        .push(fixture("data/refs2.fa").display().to_string());
    params.work_dir = output_prefix.parent().unwrap().display().to_string();

    let emitted = emit_weak_superkmer_buckets::<7>(&params, DEFAULT_SUBGRAPH_COUNT).unwrap();
    let built = build_uncolored_from_buckets::<7>(&params, &emitted.buckets.bucket_dir).unwrap();

    assert_eq!(built.input_buckets, emitted.buckets.bucket_files);
    assert_eq!(built.bucket_records, emitted.partition.weak_superkmers);
    assert!(built.observed_edges > 0);
    assert!(built.retained_edges > 0);
    assert!(built.unitigs > 0);
    assert!(built.unitig_bases >= built.unitigs * params.k as u64);
    assert_eq!(built.output_path, fasta_path);

    let fasta = fs::read_to_string(&built.output_path).unwrap();
    assert!(fasta.starts_with(">0\n"));
    let labels: Vec<_> = fasta
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect();
    assert_eq!(labels.len() as u64, built.unitigs);
    assert!(labels.iter().all(|label| label.len() >= params.k as usize));
    assert!(labels.iter().all(|label| {
        label
            .as_bytes()
            .iter()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
    }));

    let _ = fs::remove_dir_all(emitted.buckets.bucket_dir);
    let _ = fs::remove_file(built.output_path);
}

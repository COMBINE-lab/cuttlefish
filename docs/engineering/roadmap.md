# Rust rewrite roadmap

This is the working backlog for bringing the Rust implementation to Cuttlefish 3 parity.

## Current checkpoint

- The production colored and uncolored paths are end-to-end, parallel, and
  external-memory.
- On the 1000-genome HumGut2 colored benchmark, Rust completes in 19.92 seconds
  with 13,840,356 KiB peak RSS versus C++ at 20.60 seconds and 17,234,540 KiB.
  Counts and all 40,370,131 strand-normalized unitigs match exactly.
- Colored partitioning now uses C++-style source barriers and atlas collation,
  with in-place worker-stream partitioning and bounded long-lived consumers.
- Local contraction emits directly into a blocked discontinuity edge matrix;
  contraction, expansion, mapping, and reduction follow C++ phase ordering.
- See `performance-record.md` for architecture decisions, commands,
  measurements, and the scale-validation protocol.

- Uncolored Rust output matches C++ normalized FASTA on:
  - `data/refs2.fa`, `k=7`, `l=3`
  - `data/compat-small.fa`, `k=7`, `l=3`
  - `data/generated/uncolored-synthetic.fa`, `k=7`, `l=3`
- Weak-superkmer bucket emission is consumable.
- Local subgraph construction and contraction exist.
- Serial discontinuity edge-matrix construction exists.
- Serial discontinuity contraction now has diagonal-block compression, partition-level non-diagonal contraction, false-phantom handling, meta-vertex records, and a high-to-low full contraction driver with edge reinsertion.
- Serial discontinuity expansion now has path-info types and an initial path-info propagation pass over contracted edges and meta-vertices.
- The production uncolored Rust path now defaults to the discontinuity pipeline. The global in-memory canonical contractor remains available behind `CF3_RS_DEBUG_GLOBAL_CONTRACTOR=1` as a correctness oracle.
- The first discontinuity-graph handoff exists as `emit_uncolored_discontinuity_inputs`.
- Weak-superkmer emission and stats-only partitioning are parallel across input files. The current implementation preserves one consumable bucket manifest through a mutex-protected bucket emitter; replacing that with worker-local bucket shards plus manifest merge is the next scaling step for high-thread-count runs.

## Test scaling

- Keep fast Rust-only regression tests for C++-validated expected FASTA.
- Keep the C++ comparator harness for integration checks and benchmark gates.
- Add fixture classes incrementally:
  - linear paths,
  - reverse-complement duplicate paths,
  - branches,
  - bubbles,
  - tips and dead ends,
  - detached cycles,
  - palindromic and near-palindromic k-mers,
  - repeated k-mers within one reference,
  - multi-record shared paths,
  - cutoff-sensitive duplicate observations,
  - `N`-split fragments.
- Compare stats in tiers:
  - must match: unitig count, unitig bases, normalized FASTA;
  - should match after architecture parity: local unitigs, discontinuity exits, edge matrix counts;
  - informational: timing, bucket compression, peak memory.

## Uncolored parity path

1. Preserve the debug global contractor as an oracle while porting the Cuttlefish 3 pipeline.
2. Emit local maximal unitigs and discontinuity exits from local subgraphs.
3. Add serial edge-matrix construction.
4. Add serial discontinuity graph contraction.
5. Expand contracted discontinuity graph into final unitigs.
6. Collate final unitig buckets into FASTA/GFA.
7. Replace the debug global contractor in the default uncolored path. Done; keep it as an opt-in oracle until larger benchmark parity is stable.
8. Parallelize and externalize each phase after serial parity is established.

## External-memory work

- Add binary manifests and round-trip tests for:
  - local unitig buckets,
  - edge matrix buckets,
  - contracted path-info buckets,
  - final unitig buckets.
- Track file descriptor use explicitly.
- Keep deterministic manifests and sorted output for reproducible tests.

## Parallelization map

- Diagonal block compression can run concurrently with the non-diagonal column scan for the same partition.
- Non-diagonal blocks in one matrix column can be scanned in parallel against a shared partition table.
- Compressed diagonal-chain edges can be processed independently after the column table is built.
- False-phantom filtering is a parallel scan over unprocessed partition-table entries.
- Output edge and meta-vertex emission should use worker-local buffers, followed by batched external-memory bucket flushes.
- Partition processing still follows the C++ high-to-low partition dependency order until the expansion algorithm proves a wider scheduling window is safe.
- Weak-superkmer emission currently parallelizes by input-file chunks. The shared emitter is correct but can serialize writes under load; use worker-local bucket files and a deterministic manifest merge before serious throughput benchmarking.
- Local subgraph contraction is parallel across weak-superkmer bucket manifest entries and preserves manifest order in the handoff to the discontinuity pipeline.

## Benchmark ladder

1. Committed tiny fixtures.
2. Generated synthetic fixtures.
3. Small license-safe real references.
4. Medium generated references.
5. Public microbial references.
6. Larger datasets from `data/Datasets.md`.

For each rung, collect wall time, peak RSS, temporary disk bytes, unitig count, unitig bases, and comparator status.

## Colored parity path

Start colored work after uncolored discontinuity graph parity.

1. Preserve source IDs through bucket emission.
2. Finish color accumulation in vertex state and local unitigs.
3. Implement color-equivalence or color-coordinate buckets matching C++ semantics.
4. Carry colors through discontinuity contraction and expansion.
5. Emit colored output formats.
6. Add a colored comparator for sequence plus color-set parity.

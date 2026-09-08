# Issue #63: compact wide label offsets — 3.0.2 validation

The fix passes the correctness and performance checks below and is prepared
for 3.0.2. No end-to-end performance regression was measured on the tested
workloads. Coordinate records remain **24 bytes** in memory and on disk.

## Change

[Issue #63](https://github.com/COMBINE-lab/cuttlefish/issues/63) exposed a
32-bit **label byte offset**, rather than a path-ID limit. Release builds
truncated offsets while writing materialized coordinates; merging retained
tails subsequently failed when the label buffer exceeded `u32::MAX`.

`LoadedMaterializedStitchedCoordRecord` now stores the low 32 offset bits in
`label_offset_low` and the high 14 bits in the previously unused upper bits of
`flags`. Orientation and cycle remain in bits 0–1. This addresses a 64 TiB
label buffer per bucket without widening the record. A direct replacement
with `u64`, retaining the other fields, would grow this record to 32 bytes.
Path IDs, rank, lengths, color indices, and color counts retain their widths.

Writers and loaders use checked wide offset arithmetic. Both individual
records and batch writes preserve high bits. Shard and retained-tail merging
share the same loader, and disk label ranges are checked before assembly.
The loader removes coordinate files only after all shards and tails have
loaded successfully. Oversized offsets report the bucket and supported range.

The batch writer rebases its already packed records in place and copies the
native records once per batch. This avoids expanding each record and encoding
it again. Zero-base rebasing leaves the first shard's records untouched. These
details prevent the normal small-bucket path from paying for unnecessary work.

The private coordinate magic advances from `CF3MCB2` to `CF3MCB3`, so old
readers cannot silently ignore high bits. Old failed-run files are not a
supported resume checkpoint. This change does not add resumability or extend
the separate color-index limits. Final FASTA/color-repository formats and
public interfaces are unchanged.

The patch also applies equivalent `as_chunks` iterator forms and removes a
redundant cast where current/MSRV Clippy flagged pre-existing code. These
changes preserve the consumed chunks and remainders and are included in the
end-to-end candidate measurements.

## Correctness

- The complete workspace suite passes in debug and release: **73 library,
  4 CLI, and 38 compatibility tests**. The compatibility tests include
  C++-validated reference/read outputs, colored source truth, low worker
  counts, and k=33/63. Full tests also pass with the minimum supported
  compiler, Rust 1.91.0; the main validation compiler was Rust 1.98.0.
- New ordinary regressions cover serialization and native decoding at the
  32- and 46-bit boundaries, all orientation/cycle combinations, color-index
  preservation, checked rebasing and rejection, obsolete formats, label
  bounds, and retaining files when a tail is invalid. They do not require
  multi-gigabyte allocations.
- An explicitly run scalability regression loads a real label buffer above
  4 GiB. It exercises the individual and batch writers, native bulk loading,
  multiple shards, retained tails, and assembly in uncolored and colored
  modes. Four overlapping fragments carry distinct color transitions, with
  agreement at shared endpoint k-mers. The resulting label and all five
  expected color runs match. This test is ignored by default because it
  allocates approximately 4 GiB.
- Full external topology comparisons match **every unitig** in all six
  timed workloads, including strand/cycle normalization and multiplicity:

| Workload | Matching unitigs | Matching bases |
| --- | ---: | ---: |
| Salmonella 1,000, uncolored, k=31 | 18,910,541 | 1,783,360,006 |
| Salmonella 1,000, colored, k=31 | 18,910,541 | 1,783,360,006 |
| Salmonella 1,000, colored, k=55 | 10,373,413 | 1,887,406,978 |
| SRR105788 first 1M reads, uncolored, k=31 | 161,346 | 7,415,825 |
| Salmonella 10,000, colored, k=31 | 51,644,203 | 4,130,164,227 |
| HumGut 1,000, colored, k=31 | 40,370,131 | 2,726,597,540 |

Decoded source colors additionally match at every k-mer of **17,265 sampled
linear unitigs at k=31** and **9,821 at k=55**. Sampling selected canonical
labels beginning `AAAACC`; coordinates were resolved independently through
each build's color repository. This is supplementary color sampling; the
topology comparisons above cover the entire graphs.

## Real data above 4 GiB

A second colored 10,000-Salmonella build used
`CF3_RS_MCOORD_BUCKETS=1 CF3_RS_KEEP_INTERMEDIATES=1`. This forced the normal
pipeline through a single large coordinate bucket and the single-worker
reducer. Its persisted materialized bucket contained:

| Quantity | Value |
| --- | ---: |
| Coordinate records | 293,285,032 |
| Label bytes | **11,573,040,368** |
| Last persisted label offset | **11,573,040,322** |
| Last label length | 46 |
| Record width | 24 bytes |

The header, file sizes, and last encoded offset were checked directly. The
last offset crosses both the 4 and 8 GiB boundaries. The build completed and
all **51,644,203 unitigs** matched the 3.0.1 graph built with the normal fanout.
Decoded colors also matched at every k-mer of **46,693 sampled linear unitigs**
against that baseline. This verifies wide offsets using real sequences, in
addition to the sparse-file regression. It is a correctness experiment;
forced-fanout timings are not part of the performance comparison.

## Performance

Baseline: unmodified production code from
`f9d91d6f1d780cb0486f6029470438329241801b` (3.0.1). Candidate: the completed
fix, measured before its version-only update to 3.0.2. Both were rebuilt with
Rust 1.98.0, `cargo build --release`, the repository's x86-64-v3/AVX2 flags,
and default jemalloc. Host: AMD EPYC 9575F, Linux x86-64. Builds requested
32 threads and were pinned to CPUs 0–31.

Each workload warmed both binaries, then ran three interleaved pairs in
alternating order. Validation/compilation did not run concurrently with timed
builds. The table reports medians; negative wall-time changes indicate a
faster candidate. Peak RSS comes from `/usr/bin/time`. All these ordinary
builds used the default coordinate fanout, where 32-bit offsets suffice.

| Workload | Baseline seconds | Candidate seconds | Change | Baseline GiB RSS | Candidate GiB RSS |
| --- | ---: | ---: | ---: | ---: | ---: |
| Salmonella 1,000, uncolored, k=31 | 15.72 | 15.60 | −0.76% | 5.333 | 5.331 |
| Salmonella 1,000, colored, k=31 | 19.51 | 19.39 | −0.62% | 12.123 | 12.104 |
| Salmonella 1,000, colored, k=55 | 20.76 | 20.64 | −0.58% | 10.264 | 10.290 |
| SRR105788 first 1M reads, uncolored, k=31 | 0.81 | 0.81 | 0.00% | 0.464 | 0.466 |
| Salmonella 10,000, colored, k=31 | 65.06 | 64.65 | −0.63% | 19.188 | 19.230 |
| HumGut 1,000, colored, k=31 | 22.34 | 22.11 | −1.03% | 8.626 | 8.632 |

Individual phase timers fluctuate more than whole-run times. In the 10,000
case, median map time was 2.376→2.319 seconds and reduction 1.513→1.464.
In HumGut, mapping was 1.599→1.638 while reduction was 1.042→0.887.
The 1M-read map timer is only 12–17 ms, so percentage changes there are not
a useful regression measure. Complete samples, including warmups, are in
[issue-63-benchmarks.json](issue-63-benchmarks.json).

A separate production-path benchmark writes, loads, sorts, and assembles
2,000,000 records with small label buffers. The exact same benchmark function
was added to the baseline snapshot. Five interleaved pairs, after warmups,
were pinned to CPU 0. Median seconds:

| Mode | Phase | Baseline | Candidate |
| --- | --- | ---: | ---: |
| Uncolored | Write | 0.033256 | 0.032100 |
| Uncolored | Load | 0.044629 | 0.045171 |
| Uncolored | Sort | 0.046377 | 0.046231 |
| Uncolored | Assemble | 0.312820 | 0.312770 |
| Colored | Write | 0.040258 | 0.038960 |
| Colored | Load | 0.048762 | 0.048992 |
| Colored | Sort | 0.046385 | 0.046194 |
| Colored | Assemble | 0.322205 | 0.322533 |

The packed layout has zero additional record bytes; the shifts/checks still
have a CPU cost. These measurements show no measurable end-to-end regression
on this host and workload set, rather than guaranteeing identical performance
on every architecture, compiler, or input. The original multi-terabyte plant
dataset was not available for replay.

## Reproduction and remaining release steps

The checked-in regression and focused benchmark can be run with:

```sh
cargo test --workspace
cargo test --release --workspace
cargo test --release -p cuttlefish-rs --lib materialized_large_bucket_round_trip -- --ignored --test-threads=1
cargo test --release -p cuttlefish-rs --lib benchmark_materialized_coordinate_round_trip -- --ignored --nocapture --test-threads=1
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo doc --workspace --no-deps
```

On a many-core host, give parallel compatibility tests enough descriptors
(`ulimit -n 65536`) or run them with fewer test threads. The initial default
run exhausted the shell's descriptor limit; with the adjusted test environment
all tests passed. Format, strict Clippy, and documentation checks pass.
The final 3.0.2 CLI also produced the expected fixture graph with hard
descriptor limits of both 256 and 384.

The local campaign scripts, binaries, logs, dataset lists, and raw results are
under `/scratch3/rob/cuttlefish-rewrite/issue-63-validation`. Timed reference
commands use `--ref --list INPUT -k K --min-len 12 -t 32`, adding `--color`
where indicated. The read command uses `--read --min-len 15`. Saved baseline
and candidate outputs were compared using `cuttlefish compare --full-diff`.
The real wide experiment adds the two environment settings given above.
Only the final comparison outputs were retained; temporary work, sort chunks,
and superseded outputs were removed after verification. The wide bucket's
header and final record are preserved in the benchmark JSON alongside its
sizes and decoded offset.

The validation campaign built and checked a local **3.0.2** binary. For release
preparation, the Cargo manifests and lockfile remain at 3.0.1 so that
`scripts/bump_and_publish.sh` can create its separate 3.0.2 version-bump commit.
The changelog and displayed documentation version are already prepared.

Submit the fix to `develop`, as required by `CONTRIBUTING.md`, and integrate
the reviewed changes into `main` without losing its newer installation docs.
From a clean release checkout, use the repository script:

```sh
./scripts/bump_and_publish.sh 3.0.2 --publish --dry-run
# After review and approval to publish:
./scripts/bump_and_publish.sh 3.0.2 --publish
```

The dry run checks and tests the workspace, then validates both crate packages
at the current version; it does not bump, commit, tag, push, or upload. The
live command bumps both workspace and internal dependency versions, commits,
creates the annotated `v3.0.2` tag, pushes the branch and tag, and publishes
`cuttlefish-rs` followed by `cuttlefish-rs-cli`. The pushed tag triggers the
four-platform GitHub release workflow, using this version's changelog section
as release notes. Watch that workflow and verify both registry versions before
considering the release complete.

Use `--publish` on the live invocation when publishing to crates.io is intended:
omitting it still pushes the release tag, and a later invocation with the same
version is rejected by the script's existing-version/tag checks. No tagging or
publishing was performed by the validation campaign.

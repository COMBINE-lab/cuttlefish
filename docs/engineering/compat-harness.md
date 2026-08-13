# Rust Rewrite Compatibility Harness

This note tracks the comparison strategy for the Rust rewrite against the C++
`cuttlefish3` implementation.

## Goals

The harness should catch semantic drift before performance work starts:

- run the same fixture inputs through C++ and Rust;
- compare normalized graph artifacts rather than logs or incidental ordering;
- keep small deterministic fixtures in CI and leave large benchmarks as explicit
  manual jobs;
- record command lines and tool versions with each comparison result.

## Fixtures

Start with the existing checked-in fixtures:

- `data/refs1.fa`
- `data/refs2.fa`
- `data/rob-palindrome.fa`
- `data/reads.fq`
- `data/compat-small.fa`

Use small odd `k` values first (`3`, `5`, `7`) to make expected differences
easy to inspect. Add `k=31` fixture smoke tests once Rust reaches full unitig
emission.

## Current Rust Stage

The Rust implementation currently emits weak-superkmer bucket files and builds
local subgraph state from them. C++ does not expose an equivalent stable bucket
artifact, so direct C++ parity starts at later algorithm boundaries.

The default Rust uncolored path now uses the discontinuity pipeline. The older
global in-memory contractor remains as an opt-in oracle with
`CF3_RS_DEBUG_GLOBAL_CONTRACTOR=1`.

Rust-only invariants should cover:

- manifest record counts match partition stats;
- bucket headers match command parameters;
- decoded weak-superkmer labels are valid ACTG strings;
- local subgraph vertex and edge counts are internally consistent;
- discontinuity flags survive bucket round trips.
- threaded runs with `-t > 1` produce the same normalized FASTA multiset as
  single-threaded runs and C++.

## Full Graph Comparison

Once Rust emits maximal unitigs, compare C++ and Rust FASTA output as normalized
multisets:

1. Build the C++ binary, then run:

   ```bash
   module load conda/202211 cmake/3.29.2 gcc/14.2.0 boost/1.80 nasm/2.14.02
   export LZ4_PKG=/fs/cbcb-software/RedHat-8-x86_64/local/conda/202211/pkgs/lz4-c-1.9.4-h6a678d5_1
   export PKG_CONFIG_PATH="$LZ4_PKG/lib/pkgconfig:$PKG_CONFIG_PATH"
   export CPATH="$LZ4_PKG/include:$CPATH"
   export LIBRARY_PATH="$LZ4_PKG/lib:$LIBRARY_PATH"
   export LD_LIBRARY_PATH="$LZ4_PKG/lib:$LD_LIBRARY_PATH"

   cmake -S . -B build-cpp-compare -DCMAKE_BUILD_TYPE=Release
   cmake --build build-cpp-compare --target cuttlefish -j 8

   scripts/compare_uncolored_fasta.py \
     --cpp-bin build-cpp-compare/src/cuttlefish \
     -s data/compat-small.fa \
     -k 7 \
     --min-len 3 \
     --metrics-tsv /scratch2/rob/tmp/cf3-metrics.tsv \
     --keep
   ```

2. The script runs both implementations, then normalizes each FASTA:

   - parse records;
   - uppercase sequence labels;
   - canonicalize each unitig by taking `min(seq, reverse_complement(seq))`;
   - sort;
   - compare sequence multisets and total bases.

This avoids false failures from record IDs, output order, and orientation.

## Read Graphs And Cutoffs

For FASTQ/read graphs, include cutoff-sensitive fixtures:

- cutoff below all repeated edges;
- cutoff above singleton edges;
- a branch with two outgoing edges where only one reaches cutoff.

The normalized comparison remains the same, but the fixture should also record
expected retained `(k + 1)`-mer counts so failures can be localized.

## Colored Graphs

Colored output should be a separate harness phase after uncolored graph parity:

- compare uncolored unitig sequences first;
- parse color-coded FASTA headers;
- normalize source/color coordinates only after the Rust color repository format
  is finalized.

## Performance Harness

Performance comparisons should run only after semantic parity:

- pin thread count through `PARLAY_NUM_THREADS` for C++;
- pass `--threads` to the comparator so it forwards the same count to Rust
  through `-t`;
- use `/scratch2/tmp` or `/scratch3/tmp` for working directories;
- record wall time, peak RSS, bytes written, and temporary directory size;
- use `--metrics-tsv` on `scripts/compare_uncolored_fasta.py` to append elapsed
  time, max RSS, output size, and command lines for both implementations;
- run at least one small fixture, one medium reference, and one read dataset.

## The C++ implementation is not a reliable oracle

Expectations in the compat tests were generated from the C++ build, and that
requires care: it drops unitigs nondeterministically, at every thread count
tested, in both modes. Roughly a third of 64-thread runs on 10,000 genomes
disagree with their own siblings, and the lost records are always exactly
one k-mer long. Full evidence in [cpp-nondeterminism.md](cpp-nondeterminism.md).

Consequences for anyone regenerating expectations:

- **Pin threads low.** Every expectation currently checked in was produced at
  two threads and confirmed identical over three repeats. Say so in the test's
  doc comment, as the read-mode and k = 33/63 fixtures do.
- **Repeat before trusting.** A single run is not evidence. If repeats disagree,
  take the majority and record that they disagreed.
- **Prefer input-derived truth where the property allows it.** The colored tests
  check each vertex against the sources that actually contain that k-mer,
  computed from the fixture itself, which needs no oracle at all. That is the
  stronger pattern; use it when the expected value is derivable.
- **CI does not build or run C++.** Expectations are baked in as literals, so a
  regression shows up as a test failure rather than as a comparison against a
  moving target.

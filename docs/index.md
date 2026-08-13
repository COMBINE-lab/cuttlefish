# Cuttlefish 3

Cuttlefish 3 is a parallel, external-memory compacted de Bruijn graph
constructor written in Rust. It builds uncolored and positional colored graphs
from reference sequences or sequencing reads, including collections that are
too large to keep entirely in RAM.

The production pipeline partitions input into weak super-k-mer buckets,
contracts independent local subgraphs, resolves their discontinuities through
a blocked external graph, and emits maximal unitigs directly to FASTA.

> Cuttlefish 3 is under active development. The command-line interface and
> private intermediate formats may still change before a stable release.

## Features

- Parallel uncolored and colored compacted graph construction
- FASTA, FASTQ, and gzip-compressed input
- Reference graphs and cutoff-filtered read graphs
- External-memory intermediates with adaptive I/O fanout
- User-controlled worker count and soft memory budget
- Optional LZ4 compression for uncolored partition buckets
- Odd `k` values from 3 through 63
- Color queries over a finished graph, and cleanup after an interrupted run

## Build

Cuttlefish requires Rust 1.85 or newer and the standard linker and native build
tools for the selected Rust target.

```bash
git clone https://github.com/COMBINE-lab/cuttlefish.git
cd cuttlefish
cargo build --release
```

The executable is `target/release/cuttlefish`. It can also be installed
under the Cargo binary prefix:

```bash
cargo install --path crates/cuttlefish-rs-cli
```

## Quick Start

Build an uncolored reference graph:

```bash
mkdir -p work
target/release/cuttlefish build \
  --ref \
  --seq data/refs1.fa \
  --kmer-len 7 \
  --min-len 3 \
  --threads 16 \
  --work-dir work \
  --output graph
```

Build a colored graph in which each listed reference is a source color:

```bash
target/release/cuttlefish build \
  --ref \
  --list genomes.list \
  --kmer-len 31 \
  --threads 32 \
  --max-memory 128 \
  --work-dir work \
  --output colored-graph \
  --color
```

Build a graph from compressed sequencing reads:

```bash
target/release/cuttlefish build \
  --read \
  --seq reads.fastq.gz \
  --cutoff 2 \
  --work-dir work \
  --output reads-graph
```

## Command Line

```text
cuttlefish build [OPTIONS]

  -s, --seq <PATH>          input sequence file; may be repeated
  -l, --list <PATH>         file containing one input path per line
  -d, --dir <PATH>          directory containing input files
  -k, --kmer-len <INT>      odd k-mer length, 3..=63 (default: 31)
      --min-len <INT>       minimizer length (default: 12)
  -c, --cutoff <INT>        (k + 1)-mer frequency cutoff
  -t, --threads <INT>       maximum workers for parallel phases
  -m, --max-memory <GiB>    soft memory budget
  -o, --output <PREFIX>     output prefix
  -w, --work-dir <PATH>     external-memory working directory
      --read                build from sequencing reads
      --ref                 build from references
      --color               emit positional colors
      --compress-buckets    LZ4-compress uncolored partition buckets
      --skip-unreadable     skip inputs that fail to parse
  -h, --help                print build help
```

Exactly one of `--read` and `--ref` is required. Input options may be mixed and
repeated, and they accept comma-separated values. `--dir` includes regular
files directly within a directory and is not recursive.

FASTA and FASTQ are detected from file contents. A `.gz` suffix enables gzip
decompression. Characters outside `A`, `C`, `G`, and `T` split records into
independent graph fragments.

An input that cannot be parsed aborts the build, so a truncated or empty file
never silently yields an incomplete graph. `--skip-unreadable` reports such
inputs and continues, which is useful for large downloaded corpora where a
single bad file would otherwise cost a full restart. A skipped source keeps its
position in the input list, so colored source assignments are unaffected; a
source that fails part-way through contributes the records read before the
failure.

Reference builds use a default `(k + 1)`-mer cutoff of 1; read builds use 2.
Use `--cutoff` to override either default.

Besides `build`, the binary carries `compare`, which decides whether two unitig
FASTA files describe the same graph up to strand and rotation; `colors`, which
reads a colored build's repository back (see [Output](#output)); and `cleanup`,
which removes the intermediates a bailed run left behind (see [Working
Directory and Cleanup](#working-directory-and-cleanup)). Each prints its own
`--help`.

## Resource Control

`--threads` is a maximum rather than a promise that every phase uses the same
worker count. Cuttlefish independently adapts partitioning, local contraction,
and final collation to available work and memory. The default is one quarter of
the available hardware threads.

`--max-memory` supplies a soft budget in GiB for replicated phase state. Shared
graph tables can impose a higher workload-dependent minimum, so this is not a
hard operating-system memory limit.

Inputs presenting fewer files than the requested worker count — typically a
handful of large FASTQ files — are streamed: one reader thread per file handles
decompression and record assembly while every requested worker performs the
minimizer scan and bucket packing.

A plain gzip member cannot be split, so a single very large `.gz` is bounded by
single-threaded decompression no matter how many workers are requested. BGZF
input is detected automatically and its blocks are inflated in parallel, which
removes that bound; compressing input with `bgzip` is worthwhile for files that
will be read more than once. Splitting input into several files also exposes
more reader parallelism.

External-memory intermediates are placed in `--work-dir`. Fast local storage is
recommended. Cuttlefish adapts its bucket fanout to the process descriptor
limit, and narrows the maximal-unitig fanout at high worker counts, so a large
`ulimit -n` is not required. Use `--compress-buckets` to trade CPU for lower
uncolored temporary-disk usage.

## Output

Every build writes `<PREFIX>.fa`, with one maximal unitig per FASTA record.
Headers currently begin with `>0`; order and record identifiers are not stable
identifiers. Unitigs may be emitted in either forward or reverse-complement
orientation. Cyclic-unitig comparisons must account for rotation as well as
strand.

Colored FASTA headers additionally contain packed positional color runs. Each
run combines a unitig offset with a coordinate into
`<WORK_DIR>/<OUTPUT_NAME>.cf3rs.color-repository/`. The repository contains:

- `metadata.tsv` with graph parameters and one-based source assignments
- `manifest.tsv` with shard sizes and paths
- `NNN.colors` files containing delta-coded source sets

The current color repository format is `cf3rs-color-repository-v2`: source
sets are stored in a hybrid Elias-delta encoding that picks per record between
sparse gap codes, a plain bitmap, and gap codes over the complement.

### Reading colors

The packed runs and the coded sets are not meant to be read by hand.
`cuttlefish colors` is the supported way in:

```bash
# every color run of every unitig, gzipped because a dump outgrows the graph
cuttlefish colors dump -r graph.cf3rs.color-repository -i graph.fa -o dump.tsv.gz

# the distinct source sets the repository holds
cuttlefish colors sets -r graph.cf3rs.color-repository --names

# unitigs carrying source 3 but not source 7
cuttlefish colors grep -r graph.cf3rs.color-repository -i graph.fa \
    --all-of 3 --none-of 7
```

All three stream and accept `-o`/`-z`; `--names` prints source paths instead of
one-based ids. Queries by *sequence* rather than by color would need a k-mer
locator, which is sketched in [color query index](color-query-index.md) and not
built.

## Working Directory and Cleanup

A build stages its external-memory intermediates under `--work-dir` and unlinks
them as it consumes them, so a successful run leaves the working directory
empty. A run that is killed, fills the disk, or fails partway leaves them
behind, and they are large: the intermediates for a 1000-genome graph peak
around 17 GiB, and a 150000-genome graph reaches hundreds.

Everything a build creates there is named `<OUTPUT_NAME>.cf3rs.<SUFFIX>`:

| name | what it holds | phase |
|---|---|---|
| `.cf3rs.wsk` | weak super-k-mer buckets | partition |
| `.cf3rs.lmtig-labels` | concatenated local unitig labels | local contraction |
| `.cf3rs.lmtig-unitigs` | local unitig records into the labels | local contraction |
| `.cf3rs.lmtig-unitigs.edge-matrix` | blocked discontinuity edges | local contraction |
| `.cf3rs.local-unitig-buckets` | local unitigs, bucketed for collation | local contraction |
| `.cf3rs.colors`, `.cf3rs.color-runs` | positional color runs (colored builds) | local contraction |
| `.cf3rs.stitch-coords`, `.cf3rs.stitch-coords.cpp-expansion` | discontinuity path coordinates | expansion |
| `.cf3rs.final-unitigs` | final unitig buckets before the FASTA | collation |
| `.cf3rs.trivial.fa` | unitigs with no discontinuity exits, when not written straight to the output | local contraction |

`CF3_RS_KEEP_INTERMEDIATES` retains all of them, which is what to set when
inspecting a run rather than cleaning up after one.

`<OUTPUT_PREFIX>.cf3rs.color-repository` is different: it sits beside the
*output*, not in the working directory, and a completed colored build needs it
to interpret its FASTA. It is durable output, not an intermediate.

To reclaim a bailed run's space:

```bash
cuttlefish cleanup -w /path/to/work-dir --dry-run   # report, remove nothing
cuttlefish cleanup -w /path/to/work-dir             # remove them
```

Cleanup removes only names matching `<name>.cf3rs.<suffix>` exactly, since a
working directory is usually shared scratch; anything else is reported and left
alone. Use `-p/--prefix` to restrict it to one build when several share a
directory. It never touches the output FASTA -- after a bailed run that file is
partial, and whether to keep it is yours to decide -- and it never removes a
color repository without `--include-repository`.

## Pipeline

1. Parse records and emit weak super-k-mers through atlas-buffered writers.
2. Construct and contract local de Bruijn subgraphs in parallel.
3. Contract and expand the blocked external discontinuity graph.
4. Materialize maximal-unitig coordinates and reduce them directly to FASTA
   and positional colors.

This organization bounds resident state while retaining parallel work across
input parsing, local subgraphs, discontinuity partitions, and final buckets.

## Rust API

Generate the crate documentation locally:

```bash
cargo doc --workspace --no-deps --open
```

Implementation notes include the [module decomposition plan](rust-rewrite-modules.md),
[performance record](rust-rewrite-performance.md),
[compatibility harness](rust-rewrite-harness.md), and the
[color query index sketch](color-query-index.md).

## Testing

```bash
cargo test --workspace
cargo build --release
```

Topology tests canonicalize forward and reverse-complement orientations and
normalize cyclic rotations.

## Citation

> Jamshed Khan, Laxman Dhulipala, and Rob Patro. Fast and Scalable Parallel
> External-Memory Construction of Colored Compacted de Bruijn Graphs with
> Cuttlefish 3. bioRxiv (2025).
> <https://doi.org/10.1101/2025.02.02.636161>

## License

Cuttlefish is distributed under the BSD 3-Clause License.

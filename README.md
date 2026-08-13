# Cuttlefish 3

Cuttlefish 3 constructs uncolored and colored compacted de Bruijn graphs from
sequencing reads or reference sequences. It is a parallel, external-memory
implementation written in Rust and designed for collections that are too large
to keep entirely in RAM.

The implementation partitions input into weak super-k-mer buckets, contracts
local subgraphs in parallel, and stitches the resulting paths through a blocked
external discontinuity graph. Both uncolored and colored builds use the same
production pipeline.

> **Release candidate 1.** Cuttlefish 3 is feature-complete and validated
> against the C++ implementation on reference and read inputs, uncolored and
> colored, for k from 3 to 63. The command-line interface and the private
> intermediate formats may still change before the stable release.

## Highlights

- FASTA, FASTQ, and gzip-compressed input
- Reference and sequencing-read graph construction
- Uncolored and positional colored compacted graphs
- Parallel partitioning, local contraction, discontinuity contraction, and
  final collation
- External-memory intermediates with adaptive file-descriptor fanout
- User-controlled worker count and soft memory budget
- Optional LZ4 compression for uncolored weak super-k-mer buckets
- Odd `k` values from 3 through 63
- Color queries over a finished graph, and cleanup after an interrupted run

## Requirements

- A 64-bit Linux or macOS system
- Rust 1.91 or newer
- Standard linker and native build tools for the selected Rust target
- Sufficient temporary disk space for external-memory intermediates

## Build

```bash
git clone https://github.com/COMBINE-lab/cuttlefish.git
cd cuttlefish
cargo build --release
```

The executable is written to `target/release/cuttlefish`. To install it
under a Cargo-style prefix:

```bash
cargo install --path crates/cuttlefish-rs-cli
```

The default build uses jemalloc. The system allocator and mimalloc are also
available (but generally will not perform as well):

```bash
cargo build --release -p cuttlefish-rs-cli --no-default-features
cargo build --release -p cuttlefish-rs-cli --no-default-features --features mimalloc
```

## Quick Start

Construct an uncolored compacted graph from references:

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

The maximal unitigs are written to `graph.fa`.

Construct a colored graph from a collection of references:

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

Each path in `genomes.list` becomes one source color. The build writes
`colored-graph.fa` and a color repository under the working directory.

For sequencing reads, select `--read`; the default `(k + 1)`-mer frequency
cutoff then changes from 1 to 2:

```bash
target/release/cuttlefish build \
  --read \
  --seq reads.fastq.gz \
  --cutoff 2 \
  --work-dir work \
  --output reads-graph
```

## Command Line

`cuttlefish` has four subcommands: `build` constructs a graph, `compare`
decides whether two unitig FASTA files describe the same one, `colors` reads a
colored build's repository back, and `cleanup` removes the intermediates an
interrupted build left behind.

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
      --compress-buckets    LZ4-compress uncolored partition buckets (default)
      --no-compress-buckets store uncolored partition buckets uncompressed
      --skip-unreadable     report and skip inputs that fail to parse
  -h, --help                print build help
```

Exactly one of `--read` and `--ref` is required. Input options may be mixed and
repeated. Comma-separated values are accepted by `--seq`, `--list`, and
`--dir`. Directory inputs include regular files directly within the directory;
directory traversal is not recursive.

FASTA and FASTQ formats are detected from their contents. Files ending in
`.gz` are decompressed automatically. Characters outside `A`, `C`, `G`, and
`T` split a record into independent graph fragments.

## Resource Control

`--threads` is an upper bound for every parallel phase. By default, Cuttlefish
uses one quarter of the machine's available hardware threads. Input count,
work availability, and `--max-memory` may reduce phase-local concurrency.

`--max-memory` is a soft budget in GiB. It bounds replicated worker state, but
workload-sized shared tables can impose a higher minimum resident set. It is
not an operating-system-enforced memory limit.

The working directory holds partition buckets, blocked edge matrices, path
information, and collation buckets. Required disk space depends on input
redundancy and graph structure. Place it on fast local storage when possible.
The implementation adapts bucket fanout to the process file-descriptor limit;
raising `ulimit -n` can expose more I/O parallelism on large runs.

Uncolored partition buckets are LZ4-compressed by default, which reduces
temporary-disk usage; whether it also improves wall time depends on storage
bandwidth against CPU capacity, so `--no-compress-buckets` turns it off.
Colored buckets are always compressed regardless of either flag.

## Output

### Uncolored graphs

`<PREFIX>.fa` contains one maximal unitig per FASTA record. Headers currently
use `>0`; record order and identifiers are not stable identifiers.

A unitig may be emitted in either its forward or reverse-complement orientation.
Correct graph comparisons must canonicalize strand, and cyclic unitigs must
also account for rotation. `cuttlefish compare` does both, so two runs that
differ only in thread count or bucket layout still compare equal:

```bash
cuttlefish compare -a graph-a.fa -b graph-b.fa --kmer-len 31 --work-dir cmp
```

It sorts each side to disk and merges, so the comparison is bounded by
`--chunk-records` rather than by the size of the graph. `--full-diff` reports
every difference instead of stopping at the first.

### Colored graphs

Colored builds write the same unitig FASTA plus positional color runs in each
FASTA header. A run is stored as a packed decimal integer containing its unitig
offset and a coordinate into the color repository.

The repository is written as `<PREFIX>.cf3rs.color-repository/`, beside the
output FASTA rather than in the working directory, because the FASTA cannot be
interpreted without it. It contains:

- `metadata.tsv`: graph parameters, source numbering, and source paths
- `manifest.tsv`: repository shard metadata
- `NNN.colors`: delta-coded source sets

Source IDs are one-based and follow resolved input order. The repository format
is currently versioned as `cf3rs-color-repository-v2`.

`cuttlefish colors` reads all of this back without unpacking headers by hand:

```bash
cuttlefish colors dump -r colored-graph.cf3rs.color-repository \
    -i colored-graph.fa -o dump.tsv.gz          # every run of every unitig
cuttlefish colors sets -r colored-graph.cf3rs.color-repository --names
cuttlefish colors grep -r colored-graph.cf3rs.color-repository \
    -i colored-graph.fa --all-of 3 --none-of 7  # unitigs by source predicate
```

A dump is larger than the graph it describes, so it streams and can gzip on the
way out.

### Interrupted runs

A build unlinks its intermediates as it consumes them, so a successful run
leaves `--work-dir` empty. One that is killed or fails partway does not, and
the leftovers are large. `cuttlefish cleanup -w DIR` removes them, matching only
names cuttlefish itself produces and reporting anything else it finds; add
`--dry-run` to look first. The output FASTA is never touched.

## Architecture

The production build has four external-memory phases:

1. Parse input and emit weak super-k-mers into atlas-buffered partition buckets.
2. Build and contract independent local de Bruijn subgraphs in parallel.
3. Contract and expand the blocked discontinuity graph connecting local paths.
4. Map local labels into maximal-unitig coordinate buckets and reduce them
   directly to FASTA and color runs.

Rust API documentation can be generated with:

```bash
cargo doc --workspace --no-deps --open
```

User documentation is in [`docs/index.md`](docs/index.md). The development
record -- measurements, attempts, and the reasoning behind each decision,
including the reverted ones -- is catalogued in
[`docs/engineering/`](docs/engineering/README.md).

## Testing

```bash
cargo test --workspace
cargo build --release
```

The compatibility test crate canonicalizes unitig strand and cyclic rotation
when comparing graph topology.

## Citation

If you use Cuttlefish 3, please cite:

> Jamshed Khan, Laxman Dhulipala, and Rob Patro. Fast and Scalable Parallel
> External-Memory Construction of Colored Compacted de Bruijn Graphs with
> Cuttlefish 3. bioRxiv (2025).
> <https://doi.org/10.1101/2025.02.02.636161>

## License

Cuttlefish is distributed under the BSD 3-Clause License. See [`LICENSE`](LICENSE).

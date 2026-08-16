# Changelog

## 3.0.1

Performance release; outputs and interfaces are unchanged.

- The weak-super-k-mer label packer's `PEXT` fast path now dispatches on
  runtime CPU-feature detection instead of a compile-time `bmi2` gate. Every
  distributed 3.0.0 artifact (GitHub binaries, `cargo install`, bioconda) was
  silently running the scalar fallback — forfeiting a measured ~6.5% of the
  partition phase that only `target-cpu=native` builds kept. Any BMI2-capable
  CPU (Haswell, 2013+) now gets the fast path regardless of build flags.
- Builds now carry portable per-target CPU baselines, mirroring piscem:
  x86-64-v3 with AVX2 on x86-64, Neoverse-N1 / Apple-A14 on arm64.
  `.cargo/config.toml` is tracked with these baselines and the release
  workflow applies them to the prebuilt binaries; local `target-cpu=native`
  tuning moves to an explicit, gitignored `.cargo/config.local.toml`.
  The prebuilt x86-64 binaries consequently assume x86-64-v3; on older
  machines, `cargo install cuttlefish-rs-cli` still produces a baseline
  binary whose runtime dispatch keeps the packer correct everywhere.

## 3.0.0

First release of the Rust implementation of Cuttlefish 3, and the release in
which it becomes the canonical Cuttlefish: `main` is now the Rust
implementation, and this version succeeds the C++ Cuttlefish 2 (2.2.0) as the
`cuttlefish` everyone installs. It constructs uncolored and colored compacted
de Bruijn graphs from reference sequences or sequencing reads, in external
memory, for odd `k` from 3 through 63.

The Cuttlefish 3 algorithm was first carefully implemented in C++; that
initial implementation is preserved on the `cuttlefish3-cpp` branch. The C++
Cuttlefish 1 & 2 line lives on the `cuttlefish-1-2` branch (integration
history on `develop-legacy`), and its `v1.x`/`v2.x` tags are unaffected.

Versioning, stated once: the major version tracks the product generation, so
a backward-incompatible change to what a user depends on — the output FASTA,
the color repository format, or the command line — bumps the *minor* version
and is called out here. Rust library dependents should pin an exact version.

### Highlights

- **Reference and read graphs**, uncolored and colored, from FASTA, FASTQ, and
  gzip-compressed input.
- **Exact positional colors.** A vertex's color set is the set of inputs that
  contain it; the implementation does not introduce approximation beyond a
  genuine hash collision between distinct observed color sets.
- **Parallel decompression within the thread budget.** The whole gzip family is
  handled by one decoder: BGZF members decode block-parallel, and a plain
  single-member `.gz` decodes through speculative mid-stream splits, so a large
  compressed input is not bounded by one decompressing thread. Decode workers
  come out of `--threads`, never in addition to it.
- **`cuttlefish colors`** reads a colored build's repository back: `dump` for
  every color run of every unitig, `sets` for the distinct source sets, and
  `grep` for unitigs matching a source predicate. All three stream and can gzip
  their output, since a dump is larger than the graph it describes.
- **`cuttlefish cleanup`** removes the intermediates an interrupted build left
  behind, matching only names cuttlefish itself produces and reporting anything
  else it finds.
- **`cuttlefish compare`** decides whether two unitig FASTA files describe the
  same graph, canonicalizing strand and cyclic rotation, sorting to disk so the
  comparison is bounded by chunk size rather than by the graph.

### Performance

Measured against the C++ implementation on the same host and inputs, 64
threads, 149,998 Salmonella assemblies at k = 31:

| | Rust | C++ |
| --- | ---: | ---: |
| uncolored, wall | 4:27.6 | 4:57.6 |
| uncolored, peak RSS | 10.0 GB | 23.3 GB |
| colored, partition phase | 112.6 s | 119.5 s |

Roughly half the peak memory and half the intermediate disk of the C++
implementation, with identical unitig and base counts. At k = 55 on 10,000
assemblies the uncolored build is 40.8 s against 42.2 s, and colored is at
parity.

### Notes

- Odd `k` only, up to 63. The k-mer representation switches from one word to
  two above k = 31; both widths are covered by the test suite.
- Uncolored partition buckets are LZ4-compressed by default;
  `--no-compress-buckets` turns that off. Colored buckets are always compressed.
- A successful build leaves `--work-dir` empty. An interrupted one does not; see
  `cuttlefish cleanup`.
- The color repository format is `cf3rs-color-repository-v2` and sits beside the
  output FASTA, which cannot be interpreted without it.

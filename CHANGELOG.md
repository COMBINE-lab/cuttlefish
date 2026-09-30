# Changelog

## 3.1.0

Faster builds that write far less to the working directory, with the same
graphs. On 149,998 Salmonella assemblies (k = 31), measured back to back on
one host against 3.0.3, with intermediate compression off (3.0.3 has none):

| | 16 threads | 64 threads | written (t16) |
| --- | ---: | ---: | ---: |
| uncolored | 12:12 to 6:05 (-50%) | 4:45 to 2:45 (-42%) | 472 to 204 GB |
| colored | 17:46 to 11:56 (-33%) | 6:30 to 4:20 (-33%) | 616 to 347 GB |

Peak memory stays close to 3.0.3: within 0.2 GB in three of the four
configurations, and about 0.4 GB higher colored at 64 threads (18.5 against
18.1 GB). For comparison, C++ Cuttlefish 3 took 14:43 uncolored and 25:06
colored at 16 threads, writing 644 and 970 GB. Output matches 3.0.3: the same
unitigs, and matching color digests (a cyclic unitig may start at a different
position).

### Speed

- The partition's minimizer scan is vectorized: AVX2 on x86-64, chosen at
  run time with an identical scalar fallback, and NEON on aarch64. It follows
  simd-minimizers (Groot Koerkamp and Martayan, SEA 2025). Partitioning is
  23-28% faster.
- Reference FASTA records made only of upper-case ACGT skip the per-byte
  fragment scan (SSE2 on x86-64, NEON on aarch64).
- The 32-base label packer gains a NEON form on aarch64.
- At cutoff 1, the reference default, local contraction keeps edges as
  presence bits in the vertex state. Vertex-table slots shrink from 24 to
  20 bytes colored and from 16 to 12 uncolored, making uncolored local
  contraction 11% faster at 150k.
- For k > 31 the vertex table no longer pads every slot to 32 bytes.
  Uncolored local contraction is 14% faster at k = 55, and read mode's peak
  memory falls from 7.1 to 6.1 GB.

### Intermediate I/O

- Much less is written to the working directory:
  - partition buckets store each label once and reference it after that;
  - local contraction replays its color pass from memory instead of
    re-reading the bucket;
  - intermediate labels are stored at 2 bits per base;
  - the path-info, coordinate and edge records are narrower.
- `--compress-intermediates auto|on|off` lz4-compresses the local-unitig,
  coordinate and color streams. `auto`, the default, decides at startup from
  the storage under `--work-dir`:
  - network filesystems, rotational disks and FUSE mounts (with a warning)
    are compressed without timing;
  - other storage is timed for about a second, and compressed unless it
    keeps well ahead of the build.
  On a spinning disk, compression made colored builds 24% faster. On fast
  NVMe, forcing it on costs about 2%, and `auto` leaves it off at moderate
  thread counts. The bar rises with threads (up to about 4.7 GB/s from 64
  threads on), so near it the choice can go either way; it never changes the
  graph.
- `cuttlefish probe -w DIR -t N` shows that measurement and what `auto`
  would choose. It is built on a new crate, `scratch-probe`, supported on
  Linux and macOS.
- `cuttlefish cleanup` also removes files left by an interrupted probe, but
  skips any modified in the last minute.

### Compatibility

- Final outputs are unchanged. The private intermediate formats change, so
  intermediates left by an interrupted 3.0.x build are not readable by 3.1.0.
- Super-k-mers are assigned to subgraphs by a new hash, so subgraph contents
  differ from earlier releases and from C++ Cuttlefish 3. The unitigs do not.
- Library users: `BuildParams` has a new public field,
  `compress_intermediates`, so code that builds it with a struct literal must
  set it (`BuildParams::new` does). Its default, `Auto`, makes the build entry
  points time the work directory once and print their choice to stderr; set
  `On` or `Off` to skip that.

### Fixes

- A colored build whose super-k-mers all fall in one subgraph, and so has no
  discontinuity edges, wrote an empty FASTA. This affected only very small
  inputs.

## 3.0.3

- Colored builds use far less memory at low and moderate thread counts. On
  149,998 Salmonella assemblies at 16 threads, peak RSS falls from 22.4 GB to
  8.9 GB (C++ Cuttlefish 3: 11.5 GB), with wall time unchanged; at 64 threads
  it falls from 25.5 GB to 18.0 GB. Colored partition buckets are now written
  as they fill instead of being staged for a 12 GiB source window. The window
  is kept only as the fallback for a source too large to stage whole, which it
  still regroups exactly. Output is unchanged.
- The maximal-unitig coordinate map holds smaller per-bucket writer buffers
  (128 KiB instead of 1 MiB), and colored worker batches flush when any of
  their streams fills. This removes about 3 GB from the colored map-phase peak
  and about 1.2 GB from uncolored builds.
- Local contraction workers reuse their per-subgraph output buffers instead of
  allocating and freeing them for every one of the 16,384 subgraphs, which
  removes about a third of the phase's page faults, and the unitig walk no
  longer builds error values it then has to drop. Together these keep colored
  local contraction at 3.0.2's speed; end to end, 150k-assembly colored builds
  are as fast as 3.0.2 at 16 threads and 8.5% faster at 256.

## 3.0.2

- Fix [#63](https://github.com/COMBINE-lab/cuttlefish/issues/63): materialized
  coordinate buckets can address label buffers larger than 4 GiB. Label
  offsets now use 46 bits by storing their high bits in the record's unused
  flag bits, preserving the 24-byte record size and addressing up to 64 TiB
  per bucket. Path IDs and color fields retain their existing widths.
- Check label-offset limits in release builds, including disk serialization,
  shard concatenation, and retained-tail rebasing. Oversized offsets now report
  the affected bucket and supported range instead of wrapping silently.
- Rebase and copy packed records in batches, and avoid rewriting offsets when
  loading the first shard. These keep the ordinary small-bucket path efficient.
- Validate the complete materialized bucket and its retained tails before
  removing its coordinate files. This does not add checkpoint/resume support.
  The private coordinate format advances to `CF3MCB3`; failed 3.0.1 runs cannot
  be resumed using their old coordinate files. Final output formats and CLI
  interfaces are unchanged.

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
- `cuttlefish build` now reports its CPU-capability tier at startup: silent
  when the fast paths are compiled in statically, a one-line note on portable
  builds (dispatch active, with the `RUSTFLAGS` recipe for a native rebuild),
  and a loud warning on pre-BMI2 CPUs where the scalar fallback runs.

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

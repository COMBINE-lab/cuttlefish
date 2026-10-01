---
title: Command line
description: The cuttlefish commands — build, compare, colors, cleanup, probe — and how input is resolved.
---

The `cuttlefish` binary carries seven commands:

```bash
cuttlefish build [OPTION...]      # construct a (colored) compacted graph
cuttlefish compare [OPTION...]    # decide whether two unitig FASTAs match
cuttlefish colors dump|sets|grep [OPTION...]   # query a colored build
cuttlefish cleanup [OPTION...]    # remove an interrupted build's intermediates
cuttlefish probe [OPTION...]      # measure a work directory's storage
cuttlefish help                   # print the command summary
cuttlefish version                # print the release version
```

`build`, `compare` and `probe` are detailed below. `colors` has [its own
page](../colors/), as does [`cleanup`](../cleanup/). Every command prints its
own `--help`.

## `cuttlefish build`

```text
 common options:
  -s, --seq <arg>       input files
  -l, --list <arg>      input file lists
  -d, --dir <arg>       input file directories
  -k, --kmer-len <arg>  k-mer length (default: 31; odd <= 63)
      --min-len <arg>   minimizer length (default: 12)
  -o, --output <arg>    output file prefix
  -w, --work-dir <arg>  working directory (default: the system temp directory)
  -t, --threads <arg>   worker threads for parallel phases
  -m, --max-memory <arg> soft maximum memory budget in GiB
      --read            construct a compacted read de Bruijn graph
      --ref             construct a compacted reference de Bruijn graph
  -c, --cutoff <arg>    frequency cutoff for (k + 1)-mers
      --color           color the compacted graph
      --compress-buckets compress uncolored temporary buckets (default)
      --no-compress-buckets store uncolored temporary buckets uncompressed
      --compress-intermediates <auto|on|off>
                        lz4-compress intermediate record streams (default: auto,
                        which times the work directory's storage at startup;
                        see `cuttlefish probe`)
      --skip-unreadable skip inputs that fail to parse instead of aborting
  -h, --help            print usage
```

Exactly one of `--read` and `--ref` is required. `--threads` defaults to one
quarter of the machine's available hardware threads.

### Input

Input may be given as individual files (`--seq`), as newline-separated list
files (`--list`), or as directories (`--dir`). The options may be mixed and
repeated, and each accepts comma-separated values:

```bash
cuttlefish build --ref --seq a.fa,b.fa --seq c.fa --list more.txt --dir refs/
```

Directory inputs include the regular files directly within the directory;
traversal is **not** recursive.

FASTA and FASTQ are detected from file contents rather than from the extension.
Files ending in `.gz` are decompressed automatically. Characters outside `A`,
`C`, `G`, and `T` split a record into independent graph fragments.

`--skip-unreadable` turns a file that fails to parse into a warning instead of
an aborted run, which is what you want when sweeping a large heterogeneous
collection. A skipped source keeps its position in the input list, so colored
source numbering is unaffected.

### *k* and the minimizer length

`-k` must be **odd** and at most 63. The default is 31. Odd values only is a
consequence of canonical *k*-mers: an even *k* admits a *k*-mer that is its own
reverse complement.

`--min-len` sets the minimizer length used to partition the input into weak
super-*k*-mers. The default is 12.

### Frequency cutoff

`-c` drops (*k* + 1)-mers occurring fewer than the given number of times. It
defaults to **1 for references** and **2 for reads** — the read default exists
to discard singleton *k*-mers produced by sequencing error.

### Coloring

`--color` records, for each stretch of each unitig, which input sources cover
it. Each input file is one source, numbered from 1 in resolved input order.
The colored build writes the unitig FASTA plus a color repository directory
beside it; see [Output formats](../output/) for the layout and
[Querying colors](../colors/) for reading it back.

### Bucket compression

`--compress-buckets` is the default and LZ4-compresses uncolored partition
buckets, reducing working-directory size. `--no-compress-buckets` turns it off.
Whether compression improves wall time depends on your storage bandwidth
relative to spare CPU. Colored buckets are always compressed, regardless of
either flag.

### Intermediate compression

`--compress-intermediates` controls LZ4 compression of three later
intermediate streams: local unitigs, unitig coordinates, and color runs.
Compressing them costs a little CPU (about 2% of wall time where storage is
fast) and cuts what the build writes and reads back, which pays on slower
storage.

- `auto` (the default) decides at startup from the storage under
  `--work-dir`:
  - network filesystems and rotational disks are compressed without timing
    (rotational disks are detected on Linux; on macOS they are timed);
  - FUSE mounts are compressed too, with a warning, because they may be a
    fast local disk or a remote store;
  - anything else is timed for about a second, and left uncompressed only if
    it keeps well ahead of the rate the build writes;
  - if the measurement fails, intermediates are compressed.

  The build log states the choice and why:

  ```text
  cuttlefish: intermediate compression off (auto: xfs storage (direct I/O, solid state) writes 6.09 GB/s and reads 13.22 GB/s; the build produces about 0.96 GB/s, so storage must sustain 1.92 GB/s; decided in 0.07s)
  ```
- `on` and `off` skip the measurement. Pass `off` for the last few percent on
  fast local storage that `auto` classified cautiously (a FUSE mount, say),
  and `on` when scratch space is tight.

The setting never changes the graph, only how the intermediates are stored.

## `cuttlefish compare`

```text
 common options:
  -a, --a <arg>          first unitig FASTA
  -b, --b <arg>          second unitig FASTA
  -k, --kmer-len <arg>   k-mer length the graphs were built at
  -w, --work-dir <arg>   scratch directory for the sorted chunks
  -t, --threads <arg>    sorting threads (default: 8)
      --full-diff        report all differences instead of the first

 diagnostic options:
      --chunk-records <arg> records per sorted chunk (default: 1000000)
      --reuse-chunks     reuse the chunks under --work-dir from a prior run
      --a-chunks <arg>   reuse pre-sorted chunks for A from this directory
      --b-chunks <arg>   reuse pre-sorted chunks for B from this directory
      --dump-mismatch <arg> write the first differing pair to this directory
  -h, --help             print usage
```

See [Comparing graphs](../comparing-graphs/) for what it does and why a plain
`diff` will not do.

## `cuttlefish probe`

```text
  -w, --work-dir <arg>  directory to measure (default: the system temp directory)
  -t, --threads <arg>   worker threads the build would use (default: as for build)
  -h, --help            print usage
```

`probe` runs the measurement behind `--compress-intermediates auto` and says
what a build with those threads would choose, without building anything. An
excerpt of its output:

```text
write            6.09 GB/s (268 MB in 0.044s, direct, with fdatasync)
read             13.22 GB/s (268 MB in 0.020s, direct)
probe took       0.08s
build produces   about 0.96 GB/s of intermediates at 16 thread(s)
auto chooses     compression off (xfs storage (direct I/O, solid state) ...)
```

It writes, then reads back, at most 256 MiB of random data in a hidden
`.scratch-probe-*` file, bypassing the page cache where the filesystem allows,
and removes the file before it exits. Unlike `auto`, it times network, FUSE and
rotational storage too, since measuring is what it is for.

## Diagnostics

A handful of `CF3_RS_*` environment variables exist for profiling and ablation
during development (for example, stopping after the partition phase). They are
unsupported and may change without notice; the one worth knowing is
`CF3_RS_KEEP_INTERMEDIATES`, which retains the working-directory intermediates
a build would otherwise delete — see
[Cleaning up](../cleanup/#inspecting-a-run-instead).

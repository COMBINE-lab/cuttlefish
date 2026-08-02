---
title: Command line
description: The cuttlefish build and cuttlefish compare subcommands, and how input is resolved.
---

`cuttlefish` has two subcommands: `build` constructs a graph, and `compare`
decides whether two unitig FASTA files describe the same one.

```bash
cuttlefish build [OPTION...]
cuttlefish compare [OPTION...]
cuttlefish version
```

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
      --skip-unreadable skip inputs that fail to parse instead of aborting
  -h, --help            print usage
```

Exactly one of `--read` and `--ref` is required.

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
collection.

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

### Bucket compression

`--compress-buckets` is the default and LZ4-compresses uncolored partition
buckets, reducing working-directory size. `--no-compress-buckets` turns it off.
Whether compression improves wall time depends on your storage bandwidth
relative to spare CPU.

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

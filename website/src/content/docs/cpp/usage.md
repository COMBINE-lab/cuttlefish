---
title: Usage
description: The cuttlefish build command line for the C++ Cuttlefish 1 and 2.
---

`cuttlefish build --help` prints the following (the default `threads` value is
machine specific):

```txt
Efficiently construct the compacted de Bruijn graph from sequencing reads or reference sequences
Usage:
  cuttlefish build [OPTION...]

 common options:
  -s, --seq arg            input files
  -l, --list arg           input file lists
  -d, --dir arg            input file directories
  -k, --kmer-len arg       k-mer length (default: 27)
  -t, --threads arg        number of threads to use (default: 22)
  -o, --output arg         output file
  -w, --work-dir arg       working directory (default: .)
  -m, --max-memory arg     soft maximum memory limit in GB (default: 3)
      --unrestrict-memory  do not impose memory usage restriction
  -h, --help               print usage

 cuttlefish_1 options:
  -f, --format arg  output format (0: FASTA, 1: GFA 1.0, 2: GFA 2.0, 3:
                    GFA-reduced)

 cuttlefish_2 options:
      --read        construct a compacted read de Bruijn graph (for FASTQ
                    input)
      --ref         construct a compacted reference de Bruijn graph (for
                    FASTA input)
  -c, --cutoff arg  frequency cutoff for (k + 1)-mers (default: refs: 1,
                    reads: 2)
      --path-cover  extract a maximal path cover of the de Bruijn graph

 debug options:
      --vertex-set arg  set of vertices, i.e. k-mers (KMC database) prefix
                        (default: "")
      --edge-set arg    set of edges, i.e. (k + 1)-mers (KMC database) prefix
                        (default: "")

 specialized options:
      --save-mph       save the minimal perfect hash (BBHash) over the vertex
                       set
      --save-buckets   save the DFA-states collection of the vertices
      --save-vertices  save the vertex set of the graph
```

It supports GNU style arguments: `--` for long options, `-` for short. A long
option taking a parameter may be written `--opt=parameter` or `--opt
parameter`; a short option as `-o parameter`.

## Common arguments

These apply to both Cuttlefish 1 and 2.

### Input

Input files can be passed in any of these ways, and the options may be mixed:

- `-s <data files>`
- `-l <newline-separated list files of data files>`
- `-d <directories containing only the data files>`

Multiple values per option are accepted as `--seq=s1,s2,...`, `--seq s1 --seq
s2 ...`, `-s s1,s2 ...`, or `-s s1 -s s2` — and likewise for `list` and `dir`.

Sequencing reads must be **FASTQ**; reference sequences must be **FASTA**.
Input files may be gzipped.

### *k*-mer length

`k` must be **odd** and at most **127** — or **63** if installed from source
without raising the limit; see [Larger *k*-mer sizes](../large-k/). The default
is 27.

### Threads

`t` defaults to a quarter of the supported concurrent threads. Using
high-enough values is recommended.

### Output

Cuttlefish generates **two** output files:

- The maximal unitigs of the graph, in FASTA / GFA1 / GFA2 (extension `.fa` /
  `.gfa1` / `.gfa2`). The GFA formats are exclusive to Cuttlefish 1.
- A metadata file with structural characteristics of the graph and its
  compacted form (extension `.json`).

### Working directory

`w` is used for temporary files. It is **not created by Cuttlefish and must
exist beforehand**. The current directory is the default.

### Memory

A soft maximum memory limit `m` (in GB) trades RAM against execution time. It
is only adhered to if the limit is at least the minimum Cuttlefish determines
it needs internally.

`--unrestrict-memory` lifts the restriction entirely, trading extra RAM for
speed.

## Cuttlefish 1 arguments

The output format `f` is one of:

- `0` — only the maximal unitig (non-branching path) fragments, in FASTA
- `1` — the maximal unitigs, their connectivities, and the input sequence
  tilings, in GFA 1.0
- `2` — the same, in GFA 2.0
- `3` — the maximal unitigs and the input sequence tilings, in
  [GFA-reduced](../output-formats/#cuttlefish-1-output)

## Cuttlefish 2 arguments

- `--read` and `--ref` are the *input type* arguments, for sequencing reads and
  reference sequences respectively.
- The (*k* + 1)-mer frequency threshold `c` defaults to **2** for read input and
  **1** for reference input.
- `--path-cover` constructs a maximal vertex-disjoint path cover of the de
  Bruijn graph instead of its compacted variant.

## Note on open files

The edge- and vertex-set generation step can produce a large number of
temporary files on disk — up to 2000. Failing to allow that many open files
produces errors like:

```text
Error: Cannot open temporary file ./kmc_00000.bin
```

Raise the limit for the user running the process:

```bash
ulimit -n 2048
```

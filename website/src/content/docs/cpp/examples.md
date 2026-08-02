---
title: Example usage
description: Worked Cuttlefish 1 and Cuttlefish 2 invocations on the repository's sample data.
---

These examples use *k* = 3 and 4 CPU threads, with a working directory named
`temp`, on the sample files in the repository's `data` directory.

Remember that the working directory must exist beforehand:

```bash
mkdir -p temp
```

## Using Cuttlefish 2

### From FASTQ

To construct the maximal unitigs of the example FASTQ file `reads.fq` with
frequency cutoff `c = 1`:

```bash
cuttlefish build -s reads.fq -k 3 -t 4 -o cdbg -w temp/ --read -c 1
```

### From FASTA

To construct the maximal unitigs of the example FASTA file `refs1.fa`:

```bash
cuttlefish build -s refs1.fa -k 3 -t 4 -o cdbg -w temp/ --ref
```

Each of these produces two output files: `cdbg.fa`, containing the maximal
unitigs, and `cdbg.json`, a metadata file with structural characteristics of
the graph.

Multiple sequence files, lists of them, or directories of them may also be
passed — see [Usage](../usage/#input).

## Using Cuttlefish 1

To output the compacted de Bruijn graph in GFA 2.0 for the example FASTA files
`refs1.fa` and `refs2.fa`:

```bash
cuttlefish build -s refs1.fa,refs2.fa -k 3 -t 4 -o cdbg.gfa2 -f 2 -w temp/
```

Note that no `--read` or `--ref` is passed — that is what selects Cuttlefish 1
rather than Cuttlefish 2.

Lists or directories of reference files may be provided as input in the same
way, as described in [Usage](../usage/#input).

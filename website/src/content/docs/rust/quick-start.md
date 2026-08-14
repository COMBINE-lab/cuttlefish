---
title: Quick start
description: Building your first uncolored, colored, and read graphs with Cuttlefish 3.
---

All three examples use the sample data in the repository's `data/` directory,
and assume `target/release/cuttlefish` is on your `PATH` as `cuttlefish`.

## An uncolored graph from references

```bash
mkdir -p work
cuttlefish build \
  --ref \
  --seq data/refs1.fa \
  --kmer-len 7 \
  --min-len 3 \
  --threads 16 \
  --work-dir work \
  --output graph
```

The maximal unitigs are written to `graph.fa`.

## A colored graph from a collection

```bash
cuttlefish build \
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
`colored-graph.fa` plus the color repository
`colored-graph.cf3rs.color-repository/` beside it. See [Output
formats](../output/) for how the colors are encoded, and read them back with
[`cuttlefish colors`](../colors/):

```bash
# unitigs carrying source 3 but not source 7
cuttlefish colors grep -r colored-graph.cf3rs.color-repository \
    -i colored-graph.fa --all-of 3 --none-of 7
```

## A graph from sequencing reads

Select `--read` instead of `--ref`; the default (*k* + 1)-mer frequency cutoff
then changes from 1 to 2, so singleton *k*-mers from sequencing error are
dropped:

```bash
cuttlefish build \
  --read \
  --seq reads.fastq.gz \
  --cutoff 2 \
  --work-dir work \
  --output reads-graph
```

## Checking two runs agree

Unitigs may be emitted in either orientation, and cyclic ones at any rotation,
so `diff` is not a valid comparison. Use the `compare` subcommand:

```bash
cuttlefish compare -a graph-a.fa -b graph-b.fa --kmer-len 31 --work-dir cmp
```

See [Comparing graphs](../comparing-graphs/).

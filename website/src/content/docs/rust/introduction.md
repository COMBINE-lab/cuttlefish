---
title: Introduction
description: What Cuttlefish 3 is, what it builds, and how it differs from the C++ Cuttlefish.
---

Cuttlefish 3 constructs uncolored and colored compacted de Bruijn graphs from
sequencing reads or reference sequences. It is a parallel, external-memory
implementation written in Rust, designed for collections that are too large to
keep entirely in RAM.

:::caution[Under active development]
Cuttlefish 3 has not had a stable release yet. The command-line interface and
the private intermediate formats may still change. The published `cuttlefish-rs`
and `cuttlefish-rs-cli` crates are at `0.0.x`, which reserve the names and prove
the release path — they are not a release.
:::

## What it builds

The **compacted de Bruijn graph** collapses every non-branching path of the de
Bruijn graph into a single vertex, a *maximal unitig*. Cuttlefish 3 emits those
unitigs as FASTA.

A **colored** build additionally records, for each position along each unitig,
which input sources cover it. Sources are the individual input files, numbered
in resolved input order.

## How it works

The build runs as four external-memory phases:

1. Parse input and emit weak super-*k*-mers into atlas-buffered partition
   buckets.
2. Build and contract independent local de Bruijn subgraphs in parallel.
3. Contract and expand the blocked discontinuity graph that connects the local
   paths.
4. Map local labels into maximal-unitig coordinate buckets and reduce them
   directly to FASTA and color runs.

Both uncolored and colored builds use this same pipeline. See
[Architecture](../architecture/) for more.

## Highlights

- FASTA, FASTQ, and gzip-compressed input
- Reference and sequencing-read graph construction
- Uncolored and positional colored compacted graphs
- Parallel partitioning, local contraction, discontinuity contraction, and
  final collation
- External-memory intermediates with adaptive file-descriptor fanout
- User-controlled worker count and soft memory budget
- Optional LZ4 compression for uncolored weak super-*k*-mer buckets
- Odd *k* values from 3 through 63

## How it relates to the C++ Cuttlefish

Cuttlefish 3 is a from-scratch Rust implementation. It is not a port of the C++
code, and it does not read or write the C++ intermediate formats.

If you are running Cuttlefish 1 or 2 today — from Bioconda, or built with
CMake — that is the [C++ documentation](../../cpp/overview/), and it remains the
right choice for published, reproducible work. See [which version should I
use?](../../about/which-version/) if you are deciding.

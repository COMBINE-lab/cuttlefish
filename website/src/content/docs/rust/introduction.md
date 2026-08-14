---
title: Introduction
description: What Cuttlefish 3 is, what it builds, and how it differs from the C++ Cuttlefish.
---

Cuttlefish 3 constructs uncolored and colored compacted de Bruijn graphs from
sequencing reads or reference sequences. It is a parallel, external-memory
implementation written in Rust, designed for collections that are too large to
keep entirely in RAM. This Rust implementation is the canonical,
forward-looking implementation of Cuttlefish 3, and lives on the repository's
default branch.

:::note[Release status]
Cuttlefish 3 is released as version **3.0.0**, feature-complete and validated
on reference and read inputs, uncolored and colored, for odd *k* from 3 to 63.
The major version tracks the product generation, so a backward-incompatible
change to what a user depends on — the output FASTA, the color repository
format, or the command line — bumps the *minor* version and is called out in
the changelog.
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

The Cuttlefish 3 algorithm was first carefully implemented in C++, preserved on
the [`cuttlefish3-cpp` branch](../../cpp/cuttlefish3/). This Rust
implementation succeeds it: it is a from-scratch rewrite — not a port, and it
does not read or write the C++ intermediate formats — and it is where all
Cuttlefish 3 development continues.

The earlier product generations, Cuttlefish 1 and 2, are separate C++ tools
with their own output formats; they live on the `cuttlefish-1-2` branch and are
covered by the [C++ documentation](../../cpp/overview/). See [which version
should I use?](../../about/which-version/) if you are deciding.

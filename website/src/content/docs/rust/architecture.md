---
title: Architecture
description: The four external-memory phases of a Cuttlefish 3 build.
---

A production build runs as four external-memory phases. Every phase writes its
output to the working directory and the next phase reads it back, which is what
keeps peak resident set a function of the worker set rather than of the graph.

## 1. Partition

Input is parsed and split into **weak super-*k*-mers** — maximal stretches of
sequence sharing a minimizer — which are emitted into atlas-buffered partition
buckets. Each bucket collects the super-*k*-mers belonging to one subgraph, so
the buckets partition the *k*-mer space into independent pieces.

Characters outside `A`, `C`, `G`, `T` split a record into independent
fragments here, so an ambiguity code ends a fragment rather than corrupting
one.

## 2. Local contraction

Each bucket is loaded and its **local de Bruijn subgraph** built and contracted
independently, in parallel. Local contraction produces maximal unitigs of the
subgraph, plus a record of every place where a path leaves the subgraph — a
*discontinuity*.

Because the buckets partition *k*-mer space, no two workers can touch the same
vertex, so this phase needs no coordination.

## 3. Discontinuity contraction

The discontinuities from every local subgraph form their own graph, connecting
the local paths into global ones. This **blocked discontinuity graph** is held
as an external edge matrix and contracted, then expanded to recover, for each
local path, the global unitig it belongs to and its position within it.

This is the phase that stitches subgraph-local results into a global answer,
and it is the reason the local phase can ignore everything outside its bucket.

## 4. Collation

Local labels are mapped into **maximal-unitig coordinate buckets** and reduced
directly to FASTA — and, for a colored build, to positional color runs written
alongside. Reducing straight to the output format avoids materializing the
whole graph anywhere.

## Colored builds

Colored and uncolored builds use the **same pipeline**. Color is carried
through as an additional payload on the super-*k*-mers, resolved during local
contraction, and written into the color repository during collation. There is
no separate colored code path to fall out of step with the uncolored one.

## Reading further

The repository's `docs/` directory holds the development notes behind these
decisions — the module decomposition, the external-memory container design, the
performance measurements, and the test harness.

API documentation for the library crate can be generated with:

```bash
cargo doc --workspace --no-deps --open
```

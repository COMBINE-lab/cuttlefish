---
title: Output formats
description: The unitig FASTA, the positional color encoding, and the color repository layout.
---

## Uncolored graphs

`<PREFIX>.fa` contains one maximal unitig per FASTA record. Headers currently
use `>0`; **record order and identifiers are not stable**, and should not be
used to refer to a unitig across runs.

A unitig may be emitted in either its forward or its reverse-complement
orientation. Any correct comparison of two graphs must canonicalize strand, and
cyclic unitigs must additionally account for rotation. The
[`compare`](../comparing-graphs/) subcommand does both.

## Colored graphs

A colored build writes the same unitig FASTA, plus **positional color runs** in
each FASTA header. A run is stored as a packed decimal integer holding its
offset along the unitig together with a coordinate into the color repository.

Positional means the color set may change *along* a unitig — Cuttlefish 3 does
not require a unitig to be monochromatic. This differs from the Cuttlefish 1
notion of color, where the color set is a property of the whole unitig; see
[the C++ colored output](../../cpp/output-formats/#colored-output-for-cuttlefish-1).

### The color repository

The repository is written to:

```text
<WORK_DIR>/<OUTPUT_NAME>.cf3rs.color-repository/
```

and contains:

| File | Contents |
| --- | --- |
| `metadata.tsv` | Graph parameters, source numbering, and source paths |
| `manifest.tsv` | Repository shard metadata |
| `NNN.colors` | Delta-coded source sets |

Source IDs are **one-based** and follow resolved input order — the order in
which `--seq`, `--list`, and `--dir` resolved to concrete files, which
`metadata.tsv` records explicitly so you never have to reconstruct it.

The repository format is versioned as `cf3rs-color-repository-v1`.

:::caution
The repository layout is a private intermediate format and may change before
the first stable release. Read it through `metadata.tsv` and `manifest.tsv`
rather than by assuming shard names.
:::

## The working directory

Beyond the color repository, the working directory holds partition buckets,
blocked edge matrices, path information, and collation buckets. These are
intermediates: they are removed as the build consumes them, and nothing outside
the color repository is meant to outlive the run.

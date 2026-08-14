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

The repository is written **beside the output FASTA**, not in the working
directory, because it is output: the FASTA cannot be interpreted without it.

```text
<OUTPUT_PREFIX>.cf3rs.color-repository/
```

It contains:

| File | Contents |
| --- | --- |
| `metadata.tsv` | Graph parameters, source numbering, and source paths |
| `manifest.tsv` | Repository shard metadata |
| `NNN.colors` | Delta-coded source sets |

Source IDs are **one-based** and follow resolved input order — the order in
which `--seq`, `--list`, and `--dir` resolved to concrete files, which
`metadata.tsv` records explicitly so you never have to reconstruct it.

The repository format is versioned as `cf3rs-color-repository-v2`: source sets
are stored in a hybrid Elias-delta encoding that picks per record between
sparse gap codes, a plain bitmap, and gap codes over the complement. A
backward-incompatible change to this format bumps Cuttlefish's minor version
and is called out in the changelog. Read the repository through
[`cuttlefish colors`](../colors/) — or through `metadata.tsv` and
`manifest.tsv` rather than by assuming shard names.

## Working directory

The working directory holds the build's external-memory intermediates —
partition buckets, blocked edge matrices, path information, and collation
buckets. They are removed as the build consumes them, so a successful run
leaves the directory empty, and nothing there is meant to outlive the run.

Everything a build creates there is named `<OUTPUT_NAME>.cf3rs.<SUFFIX>`:

| Name | What it holds | Phase |
|---|---|---|
| `.cf3rs.wsk` | weak super-*k*-mer buckets | partition |
| `.cf3rs.lmtig-labels` | concatenated local unitig labels | local contraction |
| `.cf3rs.lmtig-unitigs` | local unitig records into the labels | local contraction |
| `.cf3rs.lmtig-unitigs.edge-matrix` | blocked discontinuity edges | local contraction |
| `.cf3rs.local-unitig-buckets` | local unitigs, bucketed for collation | local contraction |
| `.cf3rs.colors`, `.cf3rs.color-runs` | positional color runs (colored builds) | local contraction |
| `.cf3rs.stitch-coords`, `.cf3rs.stitch-coords.cpp-expansion` | discontinuity path coordinates | expansion |
| `.cf3rs.final-unitigs` | final unitig buckets before the FASTA | collation |
| `.cf3rs.trivial.fa` | unitigs with no discontinuity exits, when not written straight to the output | local contraction |

A run that fails partway leaves these behind;
[`cuttlefish cleanup`](../cleanup/) removes them safely.

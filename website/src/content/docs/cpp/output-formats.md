---
title: Output formats
description: FASTA, GFA 1.0 and 2.0, the reduced GFA format, and the Cuttlefish 1 notion of color.
---

## Cuttlefish 2 output

The currently supported output format is:

- The set of maximal unitigs (non-branching paths) of the de Bruijn graph, in
  FASTA.

Other output formats are on the development roadmap.

## Cuttlefish 1 output

The supported output formats are:

- The set of maximal unitigs of the de Bruijn graph, in FASTA
- The compacted de Bruijn graph in [GFA 1.0](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md)
  and [GFA 2.0](https://github.com/GFA-spec/GFA-spec/blob/master/GFA2.md)
- The compacted de Bruijn graph in a *reduced* GFA format

### The reduced GFA format

The reduced format consists of two files:

#### `.cf_seg` — the segments

All the maximal unitig fragments of the graph — the `S`-tagged segment entries
of GFA — each with a unique id. The file is a list of pairs `<id segment>`.

#### `.cf_seq` — the sequence tilings

The tiling of each input sequence by maximal unitig fragments — the `P` entries
of GFA 1 or the ordered groups (`O`) of GFA 2. Each line is `<id tiling>`,
where `id` names the sequence and `tiling` is a space-separated list of unitig
ids that completely cover it. Each id carries a `+` or `-` depending on whether
that unitig appears in canonical or reverse-complemented form at that position.

For the example reference `refs1.fa` (in the repository's `data` directory)
containing:

```text
>ref1
CGACATGTCTTAG
```

the output files *may* look like:

```text title=".cf_seg"
1 CGA
3 ATGTC
6 CTAAGA
```

```text title=".cf_seq"
Reference:1_Sequence:ref1 1+ 3- 3+ 6-
```

#### What is missing, and why it does not matter

The only GFA information not stated explicitly is the links (GFA 1) or edges
and gaps (GFA 2) — the `L`, `E`, and `G` entries. These are recoverable from
the tilings: a tiling `u₀ u₁ … uₙ` corresponds to the edge-and-gap multiset

```text
{ (u₀, u₁), (u₁, u₂), …, (uₙ₋₁, uₙ) }
```

Whether a pair `(uᵢ, uᵢ₊₁)` is an edge or a gap follows from comparing the
length-(*k* − 1) suffix of `uᵢ` with the prefix of `uᵢ₊₁`, in the orientations
their `+`/`-` signs give. A gap is only possible when the sequence contains
characters outside `A`, `C`, `G`, `T`.

For moderate to large genomes this format is much preferable to full GFA, which
gets very verbose in this particular scenario while carrying effectively the
same information. For the 7-human-genome dataset at *k* = 31, the compacted
graph takes **112 GB in GFA2 but only 29.3 GB** in the reduced format.

## Orientation of the output

Cuttlefish works with the **canonical** representation of *k*-mers: a *k*-mer
and its reverse complement are the same vertex of the original graph.

The maximal unitig fragments — the *segments*, in GFA terminology — are always
output in canonical form, and the orientations are guaranteed identical across
identical executions.

## "Colored" output for Cuttlefish 1

In the GFA output formats, the graph is a list of vertices (the maximal
unitigs) and the adjacencies between them, plus a path tiling for each input
sequence. So the GFA output describes a **colored** de Bruijn graph, in the
sense that each vertex's color information is encoded in the `P` (GFA 1.0) or
`O` (GFA 2.0) entries — or in the tilings of the `.cf_seq` file for the reduced
output.

### What "color" means here

Throughout the [Cuttlefish 1
manuscript](https://academic.oup.com/bioinformatics/article/37/Supplement_1/i177/6319696),
"colored de Bruijn graph" refers to a **specific** definition of color. It is
intuitive and natural when building the compacted colored graph from a set of
reference genomes, but it is *not* the case that the Cuttlefish algorithm
permits arbitrary coloring of the *k*-mers.

Under the definition adopted there, the color set of a unitig is the subset of
input references in which the unitig appears. That information is implicitly
encoded in the path entries of the output.

A consequence: **all unitigs Cuttlefish produces are monochromatic** under this
coloring. A change of color set internal to a unitig would imply either a
branch — which would end the unitig — or the start or end of some reference and
a sentinel *k*-mer, which would also end it.

If you are building from raw sequencing reads or from highly fractured
assemblies, you may want a different notion of color, one where color sets vary
*along* a unitig. That is what [Cuttlefish 3's positional
colors](../../rust/output/#colored-graphs) provide.

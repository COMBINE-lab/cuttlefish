# A k-mer index for color queries: design notes

`cuttlefish colors grep` answers "which unitigs carry these sources" by
streaming the whole FASTA and decoding a source set per color run. That is the
right shape for the question it asks -- the predicate is over colors, and every
unitig has to be considered -- but it is the wrong shape for the question users
actually arrive with, which is over *sequence*: given this k-mer, or this read,
what colors does it have?

Answering that without an index means a linear scan per query. These notes
sketch what an index would look like, what it would cost, and what in the
current output it can reuse. Nothing here is built.

## What the query is

The primitive is: **k-mer to color set**. Everything else composes from it.

- *Membership.* Is this k-mer in the graph, and if so which sources contain it?
- *Read mapping.* For a read, intersect (or vote over) the color sets of its
  constituent k-mers. This is what a Fulgor-style pseudoalignment does, and it
  is the query that makes the index worth building.
- *Presence/absence matrix.* For a panel of sequences, one row per query.

The graph already stores the answer. A k-mer lies in exactly one unitig at
exactly one offset, and the unitig's color runs give the coordinate covering
that offset, which the repository decodes to a source set. The missing piece is
purely the locator: **k-mer to (unitig, offset)**.

## Why SSHash is the right class of structure

The locator has an unusual shape that a general hash table wastes. Consecutive
k-mers of a unitig overlap in k-1 bases, so storing them independently stores
each base k times over. SSHash exploits exactly this: it keeps the unitigs as a
concatenated 2-bit string and builds a minimal perfect hash over *super-k-mers*
(maximal runs sharing a minimizer) rather than over k-mers, so the string is
stored once and the index holds one entry per super-k-mer.

Three properties make it fit here rather than merely being available:

1. **The partitioner already computes the minimizers.** Cuttlefish's buckets
   *are* minimizer classes -- `canonical_minimizer` with the same `l` the build
   used. An index built at the same `l` inherits the build's minimizer
   distribution, including whatever skew that corpus has, and can be built
   bucket-parallel by the machinery that already exists.
2. **The unitig FASTA is already the concatenated string**, minus the headers.
   A 2-bit pack of the labels is a linear pass, and the colored header's runs
   are already offset-indexed into the same coordinate space the locator needs.
3. **Its failure mode is the one we can afford.** A minimal perfect hash over a
   set answers *wrongly* for keys outside it, so the index must verify by
   comparing the k-mer against the stored string. That verification is a random
   read into the 2-bit string, which is exactly one cache miss -- the same cost
   profile the local-contraction probe loop already lives with.

## Sketch

Four artifacts, all derivable from what a colored build already writes:

```
graph.fa                    unitig labels + packed color runs  (exists)
graph.cf3rs.color-repository  deduplicated source sets         (exists)
graph.cf3rs.kmer-index/
  strings.2bit              unitig labels, 2 bits per base, concatenated
  offsets                   unitig id -> start offset in strings.2bit
  minimizers.mphf           minimal perfect hash over distinct minimizers
  buckets                   minimizer -> super-k-mer list (offset, length)
```

A lookup is then: canonical-minimize the query k-mer; hash to its bucket; scan
the bucket's super-k-mers comparing against `strings.2bit`; on a hit, convert
the position to (unitig, offset) through `offsets`; read the color run covering
that offset; decode the set.

Sizing, for the 150k-genome Salmonella corpus (252 M unitigs, 16.4 G bases,
k = 31): `strings.2bit` is 16.4 G bases at 2 bits = **4.1 GB**, which is the
floor and dominates. `offsets` at 5 bytes per unitig is **1.3 GB**, compressible
to well under half that with Elias-Fano, since it is a monotone sequence --
the same encoding family the color repository already uses. The minimizer
structures are proportional to distinct minimizers, not k-mers: at l = 12 that
is at most 4^12 = 16.7 M, negligible. So an index is **roughly a third of the
FASTA it indexes**, and unlike the FASTA it must be resident to be useful.

## The decisions that need measuring, not guessing

- **Does it need to be an MPHF at all?** With l = 12 the minimizer alphabet is
  16.7 M -- small enough for a direct-addressed table of bucket pointers, no
  MPHF, no verification-by-construction subtlety. The MPHF earns its complexity
  only at larger `l`, and `l` is a build parameter here rather than an index
  parameter. Try direct addressing first.
- **Skew.** Minimizer buckets are famously uneven, and this corpus's skew is
  already measured for the partitioner (`max subgraph bucket` in the build
  log). A bucket holding a large fraction of the graph turns the scan into the
  linear search the index exists to avoid. SSHash's answer is to special-case
  the densest minimizers; whether that is needed here is an empirical question
  about *this* corpus's histogram, which the partitioner could emit for free.
- **Canonical form.** The graph stores canonical unitigs, so a query k-mer and
  its reverse complement must reach the same entry. Cuttlefish already
  canonicalizes minimizers strand-independently
  (`canonical_minimizer_is_strand_invariant_for_fixture_sequence` pins this),
  so the bucket is strand-free; only the in-bucket comparison must try both
  orientations.
- **Build cost against build cost.** The index is worth building during the
  graph build only if it costs meaningfully less than the build does. The 2-bit
  pack and the offset array are one linear pass over output that is already
  being written, so they are nearly free. The bucket construction is a sort by
  minimizer over super-k-mers -- the *same* sort the partitioner performs, on
  the same key. That suggests the honest design is not a separate index build
  at all, but an option on the existing one.

## What this would change about `colors grep`

With a locator, `grep` gains a second mode: instead of a color predicate
scanned over all unitigs, a *sequence* query resolved by lookup. The color
predicate mode stays as it is -- it is inherently a scan, and the index does
not help it -- but the sequence mode becomes proportional to the query rather
than to the graph, which is the difference between an interactive tool and a
batch job.

The order to build it in, if it is built: 2-bit strings and offsets first (they
are independently useful for any random access into unitigs, including making
`colors dump` seekable), then direct-addressed minimizer buckets, then measure
the bucket histogram before deciding whether SSHash's skew handling or an MPHF
is needed at all.

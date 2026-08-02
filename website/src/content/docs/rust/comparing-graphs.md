---
title: Comparing graphs
description: Why diff cannot compare two unitig FASTA files, and what cuttlefish compare does instead.
---

## Why `diff` does not work

Two *correct* builds of the same input can disagree record for record:

- A maximal unitig may be emitted in either its forward or its
  reverse-complement orientation.
- A **cyclic** unitig has no canonical starting point, so it may be emitted at
  any rotation.
- Record order is not stable, and headers are not stable identifiers.

So a byte-level comparison reports differences that are not differences.
Changing the thread count is enough to produce a file that differs everywhere
while describing exactly the same graph.

## `cuttlefish compare`

```bash
cuttlefish compare -a graph-a.fa -b graph-b.fa --kmer-len 31 --work-dir cmp
```

On agreement it prints the number of unitigs matched and exits 0:

```text
OK: 41028 strand-normalized unitigs match
```

On disagreement it reports the first differing record and exits 1.

### What it does

1. **Canonicalize** every record. A linear unitig is replaced by the
   lexicographically smaller of itself and its reverse complement. A cyclic
   unitig is additionally rotated to its lexicographically minimal rotation,
   over both strands.
2. **Sort** each side, in parallel, into chunks written to `--work-dir`.
3. **Merge** the two sorted streams and compare them as multisets.

Because the sort is external, memory is bounded by `--chunk-records` rather
than by the size of the graph — the comparison scales to graphs that do not fit
in RAM, which is exactly the case Cuttlefish is built for.

The `-k` you pass must be the *k* the graphs were built at. It is what
identifies a record as cyclic: a unitig is cyclic when its first and last
*k* − 1 bases coincide.

### Reporting every difference

By default the comparison stops at the first mismatch. `--full-diff` walks both
streams to the end and reports the totals on each side, with a few examples:

```bash
cuttlefish compare -a a.fa -b b.fa -k 31 -w cmp --full-diff
```

```text
cuttlefish compare: matching=13309865 A-only=1 B-only=1
cuttlefish compare: A-only len=16 hash=1349b48ccb96713f head=AAGTATAATAGACATG …
cuttlefish compare: B-only len=10 hash=21f582aad907d745 head=CATGTCTTAG …
```

### Re-running against the same input

Sorting dominates the cost, so the chunks can be kept and reused:

- `--reuse-chunks` reuses both sides' chunks under `--work-dir`.
- `--a-chunks DIR` / `--b-chunks DIR` reuse pre-sorted chunks for one side —
  useful when comparing many builds against one fixed reference graph.
- `--dump-mismatch DIR` writes the first differing pair of sequences to
  `a.seq` and `b.seq` for inspection.

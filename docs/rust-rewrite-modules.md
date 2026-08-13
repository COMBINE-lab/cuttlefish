# Rust rewrite module boundaries

This document defines safe decomposition boundaries for the Rust Cuttlefish 3
implementation. File splitting is a maintenance change, not an invitation to
change the algorithm. Every extraction must preserve topology, packed layouts,
phase ordering, worker scheduling, and external-memory behavior.

## Current pressure points

| file | size | responsibilities |
| --- | ---: | --- |
| `discontinuity.rs` | 14.2k lines | local streams, edge matrix, contraction, expansion, stitching, reduction |
| `buckets.rs` | 3.7k lines | bucket format, readers, serial/shared emitters, atlas buffering |
| `subgraph.rs` | 2.2k lines | local vertex tables (both widths), graph construction, unitig walks, colored contraction |
| `color.rs` | 1.9k lines | color runs, sidecars, concurrent repository, repository reader, metadata |
| `partition.rs` | 1.8k lines | input scheduling, minimizer scan, bucket packing, colored source windows |

`discontinuity.rs` came down from 18k with the removal of both stitch
subsystems and, later, the K > 31 arms that duplicated the streaming expansion
and blocked contraction. It is still four times the next-largest file and still
the first thing to split.

`discontinuity.rs` should be split first. The other files remain readable
enough that splitting them before the discontinuity graph would add churn
without addressing the main problem.

## Target discontinuity layout

The eventual module should have a small `discontinuity/mod.rs` facade and the
following private implementation modules:

1. `resource.rs`: file limits, live descriptor accounting, allocator trimming,
   and opt-in RSS reporting. This extraction is complete.
2. `local_stream.rs`: local-unitig ranges, label/color streams, and local bucket
   readers/writers.
3. `edge_matrix.rs`: matrix endpoints, packed edge codecs, blocked files, and
   concurrent append buffers.
4. `contract.rs`: diagonal compression, partition tables, high-to-low blocked
   contraction, and meta-vertex output.
5. `expand.rs`: path-info tables, reverse-order propagation, compressed
   diagonal expansion, and range-bucket emission.
6. `stitch/format.rs`: path-info, coordinate, label, and color record layouts.
7. `stitch/map.rs`: local-unitig/path-info merge and coordinate bucket mapping.
8. `stitch/reduce.rs`: coordinate sorting, maximal-unitig assembly,
   canonicalization, and direct FASTA writes.
9. `emit.rs`: parallel local-subgraph contraction and handoff into external
   discontinuity inputs.

Move modules in dependency order. `edge_matrix` must precede `contract`, and
`contract` plus `local_stream` must precede `expand`. The stitch modules should
move only after those interfaces are stable.

## Visibility policy

- Keep the existing public names re-exported from `discontinuity::mod` so
  downstream callers do not change.
- Internal cross-module APIs should be `pub(super)`, not globally public.
- Prefer concrete structs and generic functions. Do not introduce trait objects
  or boxed callbacks in graph hot paths to make a split easier.
- Keep packed record definitions beside their encoders, decoders, and size
  assertions.
- Keep worker-owned buffers concrete and reusable; a module boundary must not
  add cloning, `Arc`, or per-record allocation.

## Performance and correctness gates

Each extraction must satisfy all of the following before the next one starts:

1. `cargo test -p cuttlefish-rs --lib` and CLI tests pass.
2. `cargo doc --workspace --no-deps` reports no broken links.
3. Compile-time packed-size assertions remain unchanged.
4. The 1,000-genome colored topology comparator matches all canonical unitigs.
   Anything touching a width-dependent arm additionally needs the k = 33/63
   compat tests and a real-data comparison at k > 31: those arms have no
   coverage from the k = 31 workloads the other gates use, which is how three
   independently fatal defects survived in them (see the review record).
5. Release timings for the affected phase remain within ordinary run-to-run
   noise. Any repeatable regression requires reverting or profiling the split.
6. Peak RSS, intermediate bytes, bucket flushes, and descriptor peaks do not
   increase.

Moving Rust source between modules normally has no runtime cost. Changes to
visibility can, however, encourage accidental abstraction or inhibit expected
inlining when combined with other edits. For that reason, source extraction and
algorithmic optimization should be separate commits.

## Later splits

After discontinuity decomposition, `buckets.rs` can separate `format`,
`reader`, and `atlas_writer`; `subgraph.rs` can separate dense map primitives
from contraction walks; and `color.rs` can separate run manipulation from the
repository. These are lower priority and should use the same gate sequence.

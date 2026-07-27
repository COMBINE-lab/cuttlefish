# Rust rewrite architecture and performance

This document records the production architecture of the Rust rewrite and the
evidence required for performance-parity claims. Benchmark results are valid
only when C++ and Rust use the same inputs, `k`, minimizer length, input mode,
thread limit, storage, and host allocation.

## Production architecture

The implementation follows the Cuttlefish 3 external-memory algorithm instead
of constructing a global in-memory de Bruijn graph.

1. Input records are parsed without copying sequence fragments. Rolling packed
   k-mers and the C++ wyhash minimizer specialization produce weak super-kmers.
2. Partitioning size-sorts the complete source set and schedules it dynamically
   across the requested workers, matching the C++ partitioner. Each colored
   worker drains its atlas tails after a source; uncolored workers defer atlas
   tails until the source set is exhausted.
3. A 128-subgraph atlas is source-partitioned in place. The active worker pool
   drains worker-atlas chunks directly in worker order, matching the reference
   handoff and avoiding a second full counting-sort allocation.
4. Compressed buckets are decoded through borrowed packed records. Local
   subgraphs use a dense state arena indexed by a fast, non-adversarial hash
   table, then emit compact local-unitig, label, color-run, and blocked-edge
   streams.
5. The discontinuity graph uses a true blocked external edge matrix. Partition
   contraction follows C++'s high-to-low dependency order and emits prepared
   blocked edges directly from the concurrent table scan.
6. Expansion uses the C++ path-info propagation order and compact generation-
   tagged slots. Final mapping scatters path information into max-unitig
   coordinate buckets, followed by parallel bucket reduction and FASTA output.
7. Color sets are deduplicated in a concurrent repository using a fast hash,
   compact coordinates, varint sidecars, and positional writes.

The binary layouts are Rust-specific where a smaller or cheaper representation
is available; phase boundaries and dependency ordering remain equivalent to
C++.

## Adaptive resource policy

The Rust CLI derives its work directory from `std::env::temp_dir()` and, like
C++, defaults to one quarter of the available hardware threads (with an
eight-thread fallback when the host count is unavailable). `--threads` is an
upper bound used by every production phase; there are no dataset-size
crossovers or fixed 64-worker caps.

`-m, --max-memory <GiB>` is a Rust-only soft memory budget; the current C++
Cuttlefish 3 CLI has no corresponding option. It reduces partition,
local-contraction, and post-local concurrency according to their replicated
state estimates. Memory-induced worker reductions are rounded down to leave
headroom for workload-sized shared tables.
Largest-bucket graph tables, output buffers, and the colored repository impose
a workload-dependent minimum, so this option is not an RSS hard limit.
An explicit budget that admits every requested worker preserves the exact
request, including non-power-of-two counts. Only a phase whose estimate exceeds
the budget is rounded down. A regression test covers 96 requested workers at
both 24 GiB (96/96/96 workers) and 12 GiB (96/32/96 workers).

On the 1000-genome HumGut2 control, every point below produced 40,370,131
unitigs and 2,726,597,540 bases. Wall time and peak RSS were measured on the
same host and local filesystem:

| mode | threads | C++ wall | Rust wall | Rust/C++ | C++ RSS KiB | Rust RSS KiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| uncolored | 32 | 19.81 s | 20.67 s | 1.043 | 9,260,764 | 9,258,312 |
| uncolored | 64 | 17.70 s | 18.62 s | 1.052 | 12,754,100 | 9,689,004 |
| uncolored | 128 | 21.50 s | 22.05 s | 1.026 | 17,995,892 | 10,495,064 |
| uncolored | 256 | 31.33 s | 22.35 s | 0.713 | 25,856,860 | 15,981,980 |
| colored | 32 | 23.03 s | 23.34 s | 1.013 | 9,964,708 | 10,100,644 |
| colored | 64 | 21.81 s | 19.93 s | 0.914 | 17,229,760 | 12,262,212 |
| colored | 128 | 24.31 s | 22.38 s | 0.921 | 31,056,452 | 20,184,668 |
| colored | 256 | 35.85 s | 25.50 s | 0.711 | 37,863,660 | 23,918,524 |

Uncolored is within 5.2% or faster throughout this matrix. Colored is within
1.3% at 32 threads and faster at 64, 128, and 256 threads. The colored result
uses C++-style direct worker-tail handoff: final worker-atlas tails retain
worker order instead of undergoing an extra Rust-only source counting sort.

The current soft-memory checks used 256 requested threads and the 1000-genome
HumGut2 control. Correct output was preserved at every point tested:

| mode | budget | partition/local/post workers | wall | peak RSS KiB |
| --- | ---: | ---: | ---: | ---: |
| uncolored | 8 GiB | 128/16/64 | 21.11 s | 7,958,716 |
| colored | 8 GiB | 128/16/64 | 25.63 s | 8,651,520 |
| colored | 16 GiB | 256/64/128 | 19.72 s | 14,657,080 |

Additional non-power-of-two controls used 96 requested threads. They catch a
previous policy bug which rounded 96 to 64 merely because `--max-memory` was
present, even when the estimate admitted all requested workers:

| mode | budget | partition/local/post workers | implementation | wall | peak RSS KiB |
| --- | ---: | ---: | --- | ---: | ---: |
| colored | none | 96/96/96 | C++ | 21.19 s | 23,838,576 |
| colored | 24 GiB | 96/96/96 | Rust | 18.63 s | 8,592,332 |
| colored | 12 GiB | 96/32/96 | Rust | 20.29 s | 7,540,308 |
| uncolored | none | 96/96/96 | C++ | 18.04 s | 16,139,028 |
| uncolored | 12 GiB | 96/32/96 | Rust | 17.43 s | 7,350,200 |

All four Rust controls emitted 40,370,131 unitigs and 2,726,597,540 bases.
The 12 GiB controls remain faster than C++ despite deliberately reducing only
local contraction to 32 workers. The external strand-normalizing comparator
also matched all 40,370,131 unitigs from the colored 12 GiB run to C++.

The colored 8 GiB run exceeds the requested value by 5.6%; its shared color
table and largest graph tables are not reducible by lowering worker count.
Uncolored remains below the requested limit. These checks validate budget
adaptation, not C++ memory-limit parity because Cuttlefish 3 exposes no such
limit.

## Engineering decisions

- Use phase-local bounded worker pools. The compact `k <= 31` path does not keep
  an idle global Rayon pool competing with explicit phase workers.
- Use fast deterministic hash functions for packed integer keys. Rust's default
  collision-resistant hasher is intentionally not used in graph hot paths.
- Grow worker-local graph maps from observed insertions and reuse them across
  subgraphs, matching the reference. There is no HumGut-derived vertices-per-
  weak-super-kmer preallocation factor.
- Keep labels, colors, edges, and path information external. Memory is bounded
  by the largest active bucket/table/window rather than total graph size.
- Preserve the C++ atlas order. Overflowed worker fragments and final worker
  tails are appended directly; the production path performs no extra global
  source counting sort.
- Use the active worker pool for atlas draining and bucket finalization, as the
  reference does; no Rust-only 16- or 4-worker cap remains.
- Compare unitigs after canonicalizing strand. A maximal unitig may legally be
  emitted in forward or reverse-complement orientation.

## HumGut2 acceptance result

The first authoritative colored parity run used 1000 genomes from
`/scratch3/tmp/humgut2/profile1000.list`, `k=31`, minimizer length 12, reference
mode, and a requested 256 threads. Rust bounds phase-local pools internally.

| implementation | process wall time | peak RSS | unitigs | bases |
| --- | ---: | ---: | ---: | ---: |
| C++ baseline | 20.60 s | 17,234,540 KiB | 40,370,131 | 2,726,597,540 |
| Rust | 19.92 s | 13,840,356 KiB | 40,370,131 | 2,726,597,540 |

The comparator confirmed that all 40,370,131 strand-normalized unitigs match.
Rust was 3.3% faster and used 19.7% less peak RSS on this run.

## HumGut2 scale results

A matched 5,000-genome run on the first 5,000 entries of `genomes.list` used
64 threads, `k=31`, minimizer length 12, cutoff 1, reference mode, and local
`/scratch3/tmp` storage. All implementations produced 179,767,076 unitigs and
11,337,993,544 bases. Internal local-unitig, edge-matrix, meta-vertex, and
path-info counts also matched exactly.

| path | implementation | wall time | peak RSS |
| --- | --- | ---: | ---: |
| colored | C++ | 73.28 s | 29,495,880 KiB |
| colored | Rust, direct worker tails | 75.89 s | 14,228,856 KiB |
| uncolored | C++ | 61.70 s | 26,122,016 KiB |
| uncolored | Rust, deferred atlas emission | 64.17 s | 14,231,404 KiB |

Colored Rust is 3.6% slower and uses 51.8% less memory at this rung. Uncolored
Rust is 4.0% slower and uses 45.5% less memory. Both satisfy the current 5%
scale-parity gate while preserving exact unitig and base counts.

At 1,000 genomes the same deferred atlas path completed in 18.81 s and used
13,893,160 KiB, versus 16.32 s and 23,481,744 KiB for C++. Partitioning took
1.62 s and performed exactly 16,384 flushes. This confirms the partition
optimization independently of the noisier downstream 5,000-genome phases.

### Current high-thread acceptance

The latest 5,000-genome runs used the same local filesystem and a process file
limit of 1,024 for Rust. C++ required a limit of 65,536. Rust preserves C++'s
1,024 maximal-unitig coordinate buckets when descriptors permit it and reduces
the fanout from the live descriptor count otherwise. Local contraction applies
the same accounting to its unitig buckets, weak-super-kmer readers, color
writers, and transient blocked-edge appends.

| mode | threads | C++ wall | Rust wall | Rust/C++ | C++ RSS KiB | Rust RSS KiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| uncolored | 128 | 63.10 s | 60.94 s | 0.966 | 30,692,100 | 25,799,636 |
| uncolored | 256 | 67.66 s | 62.21 s | 0.919 | 41,323,272 | 45,648,484 |
| colored | 128 | 73.14 s | 73.21 s | 1.001 | 42,448,572 | 47,550,812 |
| colored | 256 | 79.56 s | 71.20 s | 0.895 | 59,645,080 | 51,103,548 |

Rust is within 0.1% or faster at every high-thread point. The 256-thread
colored Rust run emitted the exact 179,767,076 unitigs and 11,337,993,544
bases. The 1,000-genome external comparator also matched all 40,370,131
strand-normalized unitigs after descriptor adaptation.

The 128-thread colored point was repeated immediately after its paired C++
run because an earlier Rust sample was 7.6% slower while the host load was
falling. The repeat was 73.21 seconds against C++'s 73.14 seconds. The
256-thread Rust run was 10.5% faster than C++ under the same scale workload.

On the 1,000-genome colored control at 256 threads, a matched descriptor A/B
showed that adaptation is not a hidden performance tax. With `RLIMIT_NOFILE`
1024, Rust selected 128 final coordinate buckets and completed in 20.50 seconds
at 14,332,044 KiB RSS. With a 65,536 limit, it retained C++'s 1,024 buckets and
completed in 21.95 seconds at 23,046,624 KiB RSS. Both emitted exact output;
the default-limit path was 6.6% faster and used 37.8% less peak memory. C++
still requires the raised descriptor limit on this workload.

The final no-preallocation uncolored binary also produced a 74.05-second
outlier while host load averaged 73; its unchanged post-local expansion alone
rose by about ten seconds. An immediate repeat completed in 62.21 seconds, in
line with the preceding 61.98-second run, while local contraction improved to
12.12 seconds. The outlier is not used as the acceptance sample.

C++ output at 128 and 256 threads is nondeterministic on this build. For
example, the 256-thread colored run emitted 179,764,872 records, and repeated
128-thread uncolored runs emitted different counts below the expected
179,767,076. These C++ timings remain useful performance baselines, but their
outputs cannot serve as high-thread correctness references. Rust's uncolored
5,000-genome output was checked against the known-good 64-thread C++ topology.

Subsequent profiling found that reverse-front local-unitig assembly compiled
the checked ASCII base conversion into range checks and indirect jump tables.
A valid-ASCII complement lookup reduced the 5,000-genome local contraction
from 22.77 s to 19.98 s, slightly faster than the C++ 20.46 s phase. Folding
the blocked contractor's lock into the unused high bit of its 62-bit k-mer key
and overlapping diagonal compression with non-diagonal block reads reduced a
matched 1,000-genome graph build from 18.13 s to 17.60 s (19.72 s process
wall, 13,376,800 KiB peak RSS). Counts remained exact at 40,370,131 unitigs
and 2,726,597,540 bases. A 5,000-genome validation reached local contraction
in 19.35 s and 48 of 128 discontinuity partitions in 2.8 s, but the execution
wrapper stopped the run before final output; do not treat that partial timing
as a completed parity result.

A later blocked-contraction migration replaced the Rust read-all-files barrier
with dynamically scheduled, aligned 1 MiB reads that immediately decode,
contract, and emit edges through reusable worker buffers. On the matched 5,000
genome input, contraction fell from 13.76 s to 9.02 s while preserving the
179,767,076-unitig / 11,337,993,544-base result. C++ takes 7.15 s, so this
closes 4.74 s of the 6.61 s contraction deficit but does not yet reach parity.
The same run completed in 77.53 s process wall time, down from 82.88 s; C++
remains at 61.70 s. A fused path-info `pread` experiment won on 1,000 genomes
but regressed at 5,000 and was rejected. Final collation now uses C++-style
32-bit bucket-local label offsets and packs flags into the 30-bit color-count
word, reducing its hot in-memory coordinate from 40 to 32 bytes. Scale timing
for that isolated layout change was collected under heavy host contention and
is not suitable as a new total-time acceptance result.

The uncolored local path now follows C++ for trivial maximal unitigs: unitigs
with no discontinuity exits are canonicalized and emitted during local
contraction instead of being written to local-coordinate buckets and reread by
the global collator. On 5,000 HumGut2 genomes this bypassed 60,463,337 unitigs,
reduced intermediate writes by about 5.8 GB, and lowered path-info mapping from
4.00 s to 3.62 s. The matched run completed in 77.16 s with exact counts, a
small improvement over 77.53 s; local contraction absorbed the canonicalization
cost. The 1,000-genome strand-normalized comparator still matched all
40,370,131 C++ unitigs exactly.

The materialized coordinate buckets now use the same native 32-byte layout as
the reducer's in-memory coordinate record. This removes the field-by-field
decode and narrowing pass that C++ avoids by loading `Unitig_Coord` records
directly. On the 1,000-genome control, reducer wall time fell from 1.51 to 1.33
seconds and the largest group's load time fell from 225 to 145 milliseconds.
The full build completed in 20.54 seconds and the topology comparator matched
all 40,370,131 C++ unitigs.

Expansion vertex path-info uses a native aligned 24-byte temporary record for
`k <= 31`, with rank and flags already packed for the compact lookup table.
The compact table uses byte-sized per-slot epochs instead of borrowing the two
unused high bits of a 62-bit `k=31` key. The old tags wrapped every four graph
partitions and repeatedly cleared the full table; epochs make all 128 normal
partition clears constant-time. A direct 5,000-genome A/B measured 34.12
seconds for epoch expansion and 73.32 seconds for the graph build, versus 37.84
and 77.44 seconds with the tagged-key table. Peak RSS was 37.79 GB versus 38.13
GB, and both runs emitted 179,767,076 unitigs and 11,337,993,544 bases.

Earlier profiling of the 5,000-genome colored partition identified a different
barrier.
Rust spends about 3.6-4.2 s scanning and packing, followed by 5.8-7.3 s of
window-level atlas collation and compression. C++ drains 64 KiB worker-local
atlas chunks during scanning. Whole-window pipelining was rejected because it
left partition time near 10 s and raised peak RSS to 36.2 GB. An initial attempt
to remove source collation was rejected because it retained other Rust-only
layout differences and moved work into local contraction.

The subsequent implementation drains 64 KiB worker-atlas chunks after source
boundaries and writes subgraph buffers once they reach 64 KiB. Exact source
sets are recovered with the reference algorithm's reusable bit-vector
deduplication. At 5,000 genomes this reduced colored partitioning to 6.16 s and
the measured full-run peak RSS to 14,563,704 KiB. The full run took 97.45 s,
but unchanged contraction, expansion, and reduction phases were substantially
slower during that shared-system run, so it is not a valid total-time A/B
against the earlier 78.70 s result. Keeping all 16,384 bucket descriptors open
was tested and rejected: partitioning regressed to 6.53 s.

The uncolored partitioner now uses the same bounded hierarchy as the C++
atlas: each worker retains one interleaved 64 KiB chunk per atlas, including
the destination graph ID, and drains records into 64 KiB subgraph buffers.
This replaces the graph-separated source-window representation, which created
over one million independently growing vectors at 64 workers. With the memory
bound moved to worker-atlas chunks, uncolored inputs use one globally balanced
source range rather than 2,048-source windows. On 5,000 genomes partitioning
took 3.943 s versus 3.838 s for C++, and partition peak RSS was 3.50 GB. The
complete Rust run took 69.16 s and 15,211,244 KiB peak RSS, produced the exact
179,767,076-unitig / 11,337,993,544-base totals, and remained 12.1% slower than
the 61.70 s C++ baseline while using 41.8% less peak memory. On the 1,000-genome
control, it completed in 16.49 s and 10,355,644 KiB; all 40,370,131
strand-normalized unitigs matched C++ exactly.

The CLI uses jemalloc by default. A feature-controlled mimalloc A/B on the
same 1,000-genome binary and input took 19.74 s and 13,380,700 KiB, versus
16.49 s and 10,355,644 KiB with jemalloc. Mimalloc is retained only as an
experimental build feature; it is 19.7% slower and uses 29.2% more memory on
this workload.

For `k <= 31`, blocked discontinuity edges now use the reference widths:
two 64-bit endpoints, a 16-bit weight, a 32-bit local-unitig coordinate, one
flags byte, and a 16-bit coordinate-bucket ID. This reduces each external edge
record from 31 to 25 bytes. The contraction table likewise stores its weight
as `u16`, reducing an atomic slot from 32 to 24 bytes. On 1,000 genomes the
change reduced graph construction from 15.176 to 14.640 seconds and
contraction from 2.848 to 2.707 seconds while preserving all 40,370,131
strand-normalized unitigs. The 5,000-genome run remained exact and reduced
peak RSS from 15,211,244 to 14,388,616 KiB; its 69.65-second wall time was
effectively flat on a noisy shared host, so the layout is retained for its
lower external I/O and memory footprint rather than claiming a scale runtime
win from that single sample.

The final path-info map now also matches the reference buffering condition:
each worker/bucket flush waits for independent 8 KiB coordinate and label
buffers (and an independent 8 KiB color buffer when colored), rather than
flushing when their combined size reaches 8 KiB. Materialized coordinates use
the reference `uint16_t` rank, local-unitig length, and per-unitig color count,
reducing the sortable record from 32 to 24 bytes. On the 5,000-genome
uncolored workload this reduced path-info reduction from 5.57 to 4.96 seconds,
graph construction from 64.25 to 62.92 seconds, process wall time from 68.36
to 67.00 seconds, and external writes by about 11.4 GB. The run emitted the
exact 179,767,076 unitigs and 11,337,993,544 bases. A 1,000-genome colored run
also completed without a width overflow and emitted the exact 40,370,131
unitigs and 2,726,597,540 bases in 17.42 seconds of graph time with 12,478,300
KiB peak RSS.

Diagonal contraction now preserves its endpoint hash-table allocation across
all 128 partitions, matching the lifetime of C++'s
`ankerl::unordered_dense` workspace. Reconstructing and reserving a new Rust
table per partition was insignificant at 1,000 genomes but costly at scale.
On 5,000 genomes, reuse reduced measured diagonal compression from about 3.9
to 2.93 seconds and total discontinuity contraction to 8.65 seconds. The run
retained the exact 119,303,739 meta-vertex and 713,577,963 expanded path-info
counts and produced the exact final unitig/base totals.

Path-info reduction now retains one materialized coordinate, label, and color
buffer per worker and grows each buffer only to that worker's largest group.
This mirrors the reference reducer's max-bucket workspaces and removes repeated
allocation and initialization for every group. On the 5,000-genome uncolored
workload, reduction fell from 5.11 to 2.88 seconds and graph construction took
61.93 seconds, compared with 61.70 seconds for C++. An external canonical
comparison confirmed that all 179,767,076 unitigs match C++ topology. The
1,000-genome colored workload completed in 16.44 seconds and emitted the exact
40,370,131 unitigs and 2,726,597,540 bases. Multi-shard tests cover label and
color-offset rebasing, and the generated two-source fixture passes
source-derived color validation.

Colored materialized buckets now match C++'s coordinate/color split: the
24-byte coordinate stores both the color offset and count, while the color
sidecar is one contiguous native `UnitigColor` array. The earlier Rust format
wrote an extra count word per colored coordinate and reconstructed spans while
loading. On 5,000 genomes, raw color arrays reduced mapping from 4.74 to 4.29
seconds, reduction from 3.67 to 3.34 seconds, graph construction from 69.00 to
67.41 seconds, and filesystem output by about 2.8 GB. Process wall time was
73.75 seconds and peak RSS was 14,387,348 KiB, versus C++ at 73.28 seconds and
29,495,880 KiB. The run emitted the exact 179,767,076 unitigs and
11,337,993,544 bases; the 1,000-genome topology and source-derived colored
fixture also pass.

The 5,000-genome topology can be compared exactly with bounded memory:

```bash
python3 scripts/compare_colored_fasta.py \
  --cpp-fasta /scratch3/tmp/humgut2/scale5000/cpp-uncolored3.fa \
  --rust-fasta /scratch3/tmp/humgut2/scale5000/rust-uncolored.fa \
  --topology-only --external-work-dir /scratch3/tmp/cf3-compare -k 31
```

Colored local contraction now resolves source sets directly into the final
8-byte `UnitigColor` coordinate runs. The production path no longer builds a
source-bearing `LocalUnitigColorRun` vector for every local unitig and then
allocates a second vector to repack it in the output sink; that richer form is
retained only by the test/helper API. On 5,000 HumGut2 genomes, local
contraction fell from 33.04 to 31.07 seconds and peak RSS from 14,229,424 to
14,067,268 KiB. On 10,000 genomes, local contraction fell from 63.57 to 55.41
seconds, process wall time from 127.73 to 126.20 seconds, and peak RSS from
16,102,060 to 15,674,340 KiB. The run preserved the exact 277,146,851-unitig /
16,550,721,590-base result. The 1,000-genome external comparator matched all
40,370,131 strand-normalized C++ unitigs, and all core color/source tests pass.
Unchanged blocked contraction and expansion were slower in the accepted 10,000
genome run, absorbing most of the 8.17-second local-stage gain; those phases
are the next profiling target.

The accepted colored partitioner removes source windows entirely. Inspection
of the active C++ small-source path showed that it does not run a final source
counting sort: overflowed worker-atlas chunks and final worker tails are
appended in worker order. Migrating that behavior directly reduced the
1,000-genome 64-thread run from 29.06 to 19.93 seconds and the 5,000-genome run
from 80.73 to 75.89 seconds. The latter is within 3.6% of C++ and uses 51.8%
less peak RSS.

The colored local scheduler no longer applies its historical 64-worker cap.
That cap addressed allocator and sink contention in the older source-bearing
color-run path, but direct packed color resolution and bucketed local output
removed those bottlenecks. In a paired 10,000-genome run with 256 requested
threads, Rust local contraction fell from 53.57 to 34.61 seconds and total
process wall time fell from 121.69 to 99.59 seconds. Peak RSS increased from
15,699,600 to 25,846,328 KiB, still 62% below the paired C++ peak of
68,087,012 KiB. C++ completed in 109.51 seconds, so Rust was 9.1% faster on
the paired run. Rust and C++ produced the same 1,100,356,370 local unitigs,
182,024,288 meta-vertices, 1,005,346,617 expanded path-info records, and final
277,146,851-unitig / 16,550,721,590-base result. The 1,000-genome external
comparator also matched all 40,370,131 strand-normalized unitigs under the
256-worker scheduler.

Local contraction now emits each completed unitig directly into the local
output representation instead of first retaining every `LocalUnitig` and then
making a second pass to classify unitigs and prepare blocked edges. Colored
contraction retains only the compact pending color runs needed for resolution.
Subgraph groups are scheduled by compressed bytes, matching the C++
longest-bucket-first policy more closely than the previous logical-record
ordering. On the 5,000-genome uncolored workload, direct emission reduced local
contraction from 14.64 to 12.54 seconds and preserved all 179,767,076 canonical
unitigs. A repeated colored run completed in 68.25 seconds and 22,119,004 KiB
RSS, versus C++ at 73.28 seconds and 29,495,880 KiB.

Earlier builds used a two-billion-exit crossover to cap post-local work at 64
threads. That dataset-specific policy has been removed. Every unrestricted
phase now receives the user-requested thread count, while an explicit soft
memory budget can reduce phase concurrency. The 32/64/128/256 control matrix
above verifies parity across thread counts; the 5,000-genome 64-thread
uncolored run completed in 64.17 seconds versus C++ at 61.70 seconds and
emitted the exact 179,767,076-unitig topology.

The Rust command was:

```bash
CF3_RS_PROFILE_RSS=1 /usr/bin/time -v \
  target/release/cuttlefish3-rs build \
  --list /scratch3/tmp/humgut2/profile1000.list \
  --kmer-len 31 --min-len 12 --ref --color --threads 256 \
  --output /scratch3/tmp/humgut2/profile1000-inplace16.fa \
  --work-dir /scratch3/tmp/humgut2/profile1000-inplace16-work
```

Topology was checked with:

```bash
python3 scripts/compare_colored_fasta.py \
  --cpp-fasta /scratch3/tmp/humgut2/cpp-profile1000-current.fa.fa \
  --rust-fasta /scratch3/tmp/humgut2/profile1000-inplace16.fa.fa \
  --topology-only -k 31
```

## Scale-validation protocol

Use at least three rungs: the 1000-genome control, a larger subset, and the full
31,225-genome HumGut2 list when resources permit. Run colored and uncolored
paths separately. For every rung:

1. Build the same ordered input list for both programs.
2. Set `PARLAY_NUM_THREADS` for C++ and pass `--threads` to Rust.
3. Raise `RLIMIT_NOFILE` equally and use separate work directories on the same
   local filesystem.
4. Capture complete stdout, stderr, `/usr/bin/time -v`, binary revision, command
   line, and intermediate/output byte counts.
5. Require equal unitig count and base count.
6. Require exact strand-normalized FASTA multiset equality. For colored runs,
   validate source-derived colors on a tractable subset as well.
7. Report process wall time and peak RSS, not selective phase timers.

The authoritative full-corpus comparisons used all 31,225 genomes, `k=31`,
minimizer length 12, reference mode, 256 requested threads, and paired runs on
`/scratch3/tmp`. The direct-emission colored path completed in 3:18.25 with
37,574,892 KiB peak RSS; C++ completed in 3:47.25 with 88,006,804 KiB peak
RSS. Rust was 12.8% faster and used 57.3% less peak memory. The direct-emission
uncolored path completed in 2:49.11 with 42,467,428 KiB peak RSS; C++ completed
in 3:12.36 with 62,787,908 KiB peak RSS. Rust was 12.1% faster and used 32.4%
less peak memory.

Both Rust modes produced 1,987,425,531 local unitigs, 1,780,367,364 blocked
edges, 360,512,617 meta-vertices, 1,780,691,705 expanded path-info records, and
the final 567,570,784 unitigs / 31,235,773,275 bases. The native external-memory
comparator streamed both full FASTA files and matched all 567,570,784
strand-normalized unitigs with zero one-sided records. The 1,000-genome
source-derived fixture supplies the corresponding color-content check.

Repeated full C++ runs were not topology-deterministic at high thread counts;
an experimental slot-lock change reduced but did not eliminate that mismatch
and was not retained. Full-scale Rust colored/uncolored equality is therefore
the deterministic topology gate, while tractable matched C++ runs remain the
cross-implementation topology gate.

The best pre-parity full-corpus Rust implementation took 8:40.59 and
22,354,288 KiB RSS. The current direct packed-color path is 2.63 times faster;
its higher peak memory is the deliberate consequence of honoring all 256
requested local workers rather than applying the obsolete 64-worker cap.

## Read-mode performance record

Reference workloads drove every measurement above. Read workloads were never
benchmarked, and they exposed a structural limit that reference input hides.

`BuildParams::partition_workers` clamps workers to the input *file* count, and
the production partitioner gives each worker one whole file. Sequencing-read
input is a handful of very large files, so partitioning ran on two cores of 256
for a two-file input, and its cost was identical at `--threads 4` and
`--threads 256`. On `SRR105788_{1,2}.fastq.gz` (4.4 GB gzip, `k=31`, `l=15`),
partitioning took 48.3 s of a 63.9 s build at either thread count.

Two changes address this. The gzip backend moved from `flate2`'s default
`miniz_oxide` to `zlib-rs`, the pure-Rust port of zlib-ng, which is enabled with
`flate2`'s `zlib-rs` feature and introduces no C code. The partitioner then
gained a streamed path: when the input presents fewer files than requested
workers, one reader thread per file decompresses and assembles records into a
bounded pool of recyclable fragment batches, and every requested worker drains
that pool to run the minimizer scan and bucket packing. This mirrors the C++
partitioner's chunk queue (`Graph_Partitioner::read_chunks`). Inputs with at
least as many files as workers keep the direct one-file-per-worker path
unchanged, so reference-mode behavior and tuning are untouched.
`CF3_RS_DIRECT_PARTITION=1` forces the direct path for A/B measurement.

| stage | wall | partition | peak RSS |
| --- | ---: | ---: | ---: |
| baseline | 63.91 s | 48.32 s | 13.0 GB |
| `zlib-rs` | 56.88 s | 40.94 s | 12.6 GB |
| streamed partitioner | 32.07 s | 15.82 s | 12.2 GB |
| C++ | 68.54 s | — | 19.8 GB |

Partitioning is 3.05 times faster and the build 1.99 times faster; against C++
the build is 2.14 times faster using 38% less memory. Every run emitted the
exact 13,309,867 unitigs and 1,223,963,745 bases, and `cf3-compare-fasta`
matched all 13,309,867 strand-normalized unitigs against both the previous Rust
output and a C++ reference.

Read rungs at 256 threads with the accepted build, all counts exact:

| dataset | files | gzip bytes | before | after | C++ |
| --- | ---: | ---: | ---: | ---: | ---: |
| `SRR105788` pair | 2 | 4.4 GB | 63.91 s | 30.87 s | 68.54 s |
| `ggallus` | 12 | 23 GB | 92.62 s | 55.19 s | 91.33 s |
| `SRR622461_1` | 1 | 2.4 GB | — | 52.60 s | — |
| `mbal` | 4 | 48 GB | — | 167.68 s | — |

`mbal` is the remaining limit rather than a win. Two of its four files hold 96%
of the data, so only four reader threads can be used and 145.9 s of its 167.7 s
is partitioning. Single-member gzip cannot be split without a speculative or
indexed decoder, so reader-side parallelism is bounded by the file count. Inputs
that are one enormous gzip member will stay decompression-bound.

## Descriptor and fanout policy

Maximal-unitig coordinate buckets were previously fixed at C++'s 1024 and
reduced only when `RLIMIT_NOFILE` forced it. Two experiments revised that.

Making `OpenMaterializedWriterCache` a real descriptor budget — its `limit` was
set to the bucket count, so its eviction path never ran — was rejected. A
mapping worker scatters every batch across all buckets, so a writer cache
smaller than the bucket count thrashes. Capping open writers to 768, 256, 128,
and 64 cost 8.7%, 24%, 29%, and 31% of wall time on the colored 10,000-genome
workload at 256 threads. Under `ulimit -n 1024` the evicting plan took 52.59 s
against 39.70 s for the existing fanout-reducing plan.

Reducing the fanout itself was the actual lever. On the same workload:

| buckets | wall | peak RSS | descriptor peak |
| ---: | ---: | ---: | ---: |
| 1024 | 42.17 s | 35.8 GB | 3562 |
| 512 | 40.53 s | 25.4 GB | 2057 |
| 256 | 39.56 s | 28.3 GB | 1285 |
| 128 | 38.95 s | 24.4 GB | 978 |
| 64 | 39.04 s | 22.8 GB | 974 |

Two opposing costs decide the fanout. Each mapping worker keeps staging and a
writer per bucket, which favours fewer buckets as workers increase. But the
reducer sizes its per-worker workspaces from the *largest* bucket, so halving
the fanout doubles that workspace in every worker, which favours more buckets as
the graph grows.

A thread-count rule captured only the first effect and was wrong at scale. On
149,998 Salmonella assemblies at 256 threads, narrowing to 256 buckets cost 33%
more peak memory than 1024 — 74.9 GB against 50.3 GB — for no time difference,
exactly inverting the 10,000-genome result. A 50,000-genome guard had not been
A/B'd on fanout and did not catch it.

The fanout is therefore chosen from bucket *size*, using the local-unitig base
count already carried in the inputs: at 128 threads or more it narrows to 256
buckets only while each would hold at most 64 MiB of local-unitig bases, and
otherwise keeps Cuttlefish's 1024. Measured colored runs at 256 threads:

| genomes | local unitig bases | fanout | wall | peak RSS |
| ---: | ---: | ---: | ---: | ---: |
| 10,000 | 1.21e10 | 256 | 39.94 s | 26.9 GB |
| 50,000 | — | 1024 | 1:54.63 | 39.8 GB |
| 149,998 | 4.38e10 | 1024 | 4:43.26 | 49.9 GB |

Below 128 threads the full fanout is always used, so that path is unchanged.
`CF3_RS_MCOORD_BUCKETS` overrides the fanout for measurement.

A related observation: before this change, running with a *lower* descriptor
limit was already both faster and lighter than running with a generous one
(39.70 s and 22.5 GB at `ulimit -n 1024` versus 42.58 s and 33.5 GB unrestricted).
The tool was leaving that on the table whenever descriptors were plentiful.

## Rejected experiments

- **Link-time optimization.** `lto = "fat"` with `codegen-units = 1` was 1.4% to
  2.7% slower on reference workloads across interleaved repeats, for about 1.8%
  on read workloads, at 4.6 times the build time. Reverted.
- **Materialized writer eviction.** Described above.

## Colored repository at scale

The fixed color table is capped at 2^26 slots and diverts to an `scc` overflow
map once three quarters full, after which every lookup on a saturated table
spins in `wait_until_quiescent`. On 10,000 and 50,000 Salmonella assemblies the
table reported 67,108,864 slots and **zero** overflow entries, so this cliff is
not reached on that corpus and the redundancy of a single-species collection
keeps the distinct-color-set count far below the cap.

At 149,998 assemblies the table does saturate, reporting 67,108,864 slots and
36,810,127 overflow entries. Admitting all of them into the primary table — by
raising both the estimate ceiling and the slot cap, which produced a
134,217,728-slot table with zero overflow — changed nothing measurable: the
colored full-corpus build took 4:42.34 and 75.0 GB against 4:41.22 and 75.7 GB
for the saturated configuration. Overflow is a slope, not a cliff, and the
`scc` map absorbs tens of millions of colours at the cost of the flat table.
The default ceiling is therefore unchanged; `CF3_RS_EXPECTED_COLORS` and
`CF3_RS_COLOR_TABLE_SLOTS` expose both bounds for further measurement.

`AtomicColorTable::entry` now resolves an already-published key with a plain
probe before touching `active_insertions`, which previously performed a
`fetch_add` and `fetch_sub` on one shared cache line for every call including
pure hits. `wait_until_quiescent` latches once the counter drains, since
saturation stops further increments. Both measured neutral on the workloads
tested here — colored 10,000 genomes at 256 threads was 40.52 s median either
way — and are retained for the contention they remove rather than a measured
win. The 2,097,151-source ceiling is now a `TooManySources` error raised during
partitioning instead of an assertion inside a contraction worker.

## Full-corpus Salmonella acceptance

149,998 of the 149,999 assemblies in `/scratch4/rob/blackwell_salmonella`,
`k=31`, minimizer length 12, reference mode, 256 requested threads. One input,
`SAMEA2674605.contigs.fa.gz`, is zero bytes; `gzip -t` and the pre-change binary
reject it identically, and it is excluded. `--skip-unreadable` now allows such a
corpus to build without curating the list by hand.

| mode | wall | peak RSS | descriptor peak | unitigs | bases |
| --- | ---: | ---: | ---: | ---: | ---: |
| uncolored | 3:42.82 | 65,858,656 KiB | 1033 | 252,487,658 | 16,417,233,428 |
| colored | 4:43.26 | 52,325,000 KiB | 3593 | 252,487,658 | 16,417,233,428 |

Colored and uncolored agree exactly, which is the deterministic topology gate at
this scale.

Matched against C++ on the same corpus, colored, same host and filesystem:

| threads | impl | wall | peak RSS | descriptor peak | unitigs |
| ---: | --- | ---: | ---: | ---: | ---: |
| 16 | Rust | 21:44.64 | 10.2 GB | 3113 | 252,487,658 |
| 16 | C++ | 24:33.02 | 10.9 GB | 27,915 | 252,487,658 |
| 32 | Rust | 11:48.29 | 11.9 GB | 3145 | 252,487,658 |
| 32 | C++ | 13:12.56 | 17.1 GB | 27,980 | 252,487,658 |
| 64 | Rust | 7:21.32 | 15.7 GB | 3209 | 252,487,658 |
| 64 | C++ | 7:33.78 | 28.3 GB | 28,109 | 252,487,658 |
| 256 | Rust | 4:43.08 | 44.6 GB | 3588 | 252,487,658 |
| 256 | C++ | 8:34.00 | 69.4 GB | 28,878 | 252,481,862 |

Rust is ahead on wall time and peak memory at every point. Its margin is not
monotonic — 11.4%, 10.6%, 2.7%, then 81.6% — because 64 threads is C++'s best
scaling step (1.75 times from 32) after which it degrades, running 1.13 times
*slower* at 256 than at 64. Rust improves throughout but sub-linearly: 1.84,
1.60, and 1.56 times for the 16-to-32, 32-to-64, and 64-to-256 steps.

Memory parity is closest at 16 threads, where the two are 6.4% apart. Rust's
advantage comes largely from per-worker replicated state that barely exists at
low worker counts, so the gap widens with thread count rather than being a
constant factor.

The 16-thread point is worth noting on its own: a 150,000-genome colored graph
builds in 10.2 GB, which fits a modest machine.

Below 128 threads C++ is deterministic and both implementations emit identical
counts at full corpus scale, so only the 256-thread C++ output is unusable.

Smaller reference rungs at low thread counts, all counts matching:

| genomes | mode | threads | Rust | C++ | Rust RSS | C++ RSS |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1,000 | uncolored | 16 | 23.78 s | 26.03 s | 4.2 GB | 5.8 GB |
| 1,000 | uncolored | 32 | 15.99 s | 17.57 s | 5.3 GB | 8.6 GB |
| 1,000 | colored | 16 | 28.43 s | 31.39 s | 6.2 GB | 6.9 GB |
| 1,000 | colored | 32 | 18.97 s | 22.51 s | 7.9 GB | 9.3 GB |
| 10,000 | uncolored | 16 | 1:28.12 | 1:30.74 | 4.8 GB | 6.8 GB |
| 10,000 | uncolored | 32 | 52.32 s | 54.79 s | 5.6 GB | 10.8 GB |
| 10,000 | colored | 16 | 1:56.01 | 2:05.89 | 7.0 GB | 8.6 GB |
| 10,000 | colored | 32 | 1:08.59 | 1:15.46 | 8.8 GB | 11.7 GB |

`cf3-compare-fasta` matched all 252,487,658 strand-normalized unitigs between
the two 64-thread outputs, so the implementations are verified identical at full
corpus scale, not merely equal in count. C++'s 5,796-unitig deficit at 256
threads is therefore a function of thread count rather than workload.

## Dense path-info arrays dominated peak memory

Tracing RSS across phases showed the peak was not where the earlier fanout work
assumed. At 64 threads it arrived during the path-info map, which added 25.5 GB
on top of 13.7 GB; at 256 threads the reduce phase added a further 31.5 GB.

The map phase builds a dense array with one 16-byte `DenseLocalPathInfo` per
local unitig *in the bucket being processed*, and keeps one per worker. Local
unitig buckets were capped at the worker count, so the total was

    workers * (local_unitigs / buckets) * 16  =  local_unitigs * 16

independent of thread count — 18.2 GB for the 1,136,040,479 local unitigs of
149,998 assemblies, which is exactly the amount the phase trace could not
otherwise account for.

Each worker now owns a private contiguous span of several buckets and rotates
through them, so bucket ownership stays worker-exclusive and the writer mutexes
stay uncontended, while each bucket — and therefore each dense array — shrinks
by the oversubscription factor. Total buckets are capped at 1024, because the
extra files are nearly free at 64 workers but cost 3.3% of wall time at 256.

| workload | before | after |
| --- | ---: | ---: |
| 149,998 colored, 64 threads | 7:13.75, 39.2 GB | 7:21.32, 15.7 GB |
| 149,998 colored, 256 threads | 4:43.26, 49.9 GB | 4:43.08, 44.6 GB |
| 10,000 colored, 64 threads | 49.12 s, 16.0 GB | 49.00 s, 12.0 GB |
| 10,000 colored, 256 threads | 39.94 s, 26.9 GB | 40.18 s, 25.9 GB |

Peak memory falls 60% at 64 threads and 25% at 10,000 genomes, for 1.7% of wall
time at the full-corpus 64-thread point and nothing measurable elsewhere.

## Allocation and inner-loop cleanups

Four changes, measured together on interleaved repeats:

- `complement_ascii` is now a 256-entry lookup table built at compile time from
  `Base::from_ascii(..).complement().to_ascii()`, so it is equivalent by
  construction (a test asserts this for all 256 bytes). Reverse-strand label
  assembly complements one base at a time, and the enum path compiled to range
  checks and a jump table.
- The direct-FASTA writer reused a header buffer instead of allocating a `Vec`
  per unitig, and stopped cloning a `PathBuf` per record purely to build an
  error that is almost never constructed.
- Expansion takes `compressed_diagonal_edges` from the contraction instead of
  deep-copying it. The caller drops the contraction immediately afterwards, so
  the copy was pure waste.
- `add_prepared_edges` counts phi and diagonal edges inside its existing loop
  rather than making two further passes over the batch.

| workload | before | after |
| --- | ---: | ---: |
| Salmonella 10,000, colored, 256 threads | 40.64 s | 39.73 s |
| Salmonella 10,000, uncolored, 256 threads | 34.83 s | 34.44 s |
| `SRR105788` read, 256 threads | 29.56 s | 28.98 s |

Read-mode peak RSS also fell from about 9.7 GB to 9.0 GB. All counts were
preserved and `cf3-compare-fasta` matched all 51,644,203 strand-normalized
unitigs against the deterministic 64-thread C++ reference.

## The single-large-gzip limit

`mbal` is the workload this does not help. Its four files are 21.4 GB, 2.2 GB,
21.4 GB, and 2.2 GB, so only four reader threads apply and the two large files
dominate: 145.9 s of a 167.7 s build is partitioning. That works out to about
147 MB/s of compressed input per reader, which is the single-thread `zlib-rs`
inflate rate. Reader parallelism is bounded by file count because a gzip member
is a single LZ77 stream: block boundaries are not byte-aligned, are not indexed,
and back-references reach up to 32 KiB into already-decoded output.

Surveyed options, none currently adoptable under a pure-Rust dependency policy:

| approach | status |
| --- | --- |
| `rapidgzip` crate | Rust bindings over the vendored C++ decoder; violates the no-C-code policy |
| `gzp` | parallel *compression*, and decompression only for block formats (BGZF, mgzip); no arbitrary-gzip parallel decode, and its backends are not pure Rust |
| `flate2` + `miniz_oxide`/`zlib-rs` | pure Rust, single-threaded decode only |
| pugz | C++, and restricted to text whose bytes lie in 9-126 |

There is no pure-Rust parallel gzip decoder. The `flate2` maintainers have
discussed adopting the rapidgzip approach but nothing exists.

Three routes remain, in increasing order of effort:

1. **Exploit block-structured input when present.** Implemented. BGZF
   concatenates independent gzip members and records each member's compressed
   size in a `BC` extra-field subfield, so members are located by reading headers
   alone and inflated concurrently. `crate::bgzf::ParallelBgzfReader` deals whole
   blocks round-robin to a worker pool and returns them in the same rotation, so
   the decompressed byte stream is identical to a serial decode and parsing is
   untouched. Detection reads the first 64 bytes; plain gzip falls back to the
   serial decoder at no cost.

   On a 1.66 GB BGZF transcode of `SRR197985_1` (6.3 GB of payload) at 256
   threads, against the same data as plain gzip:

   | input | wall | partition | graph build |
   | --- | ---: | ---: | ---: |
   | plain gzip | 20.80 s | 11.65 s | 9.04 s |
   | BGZF | 16.81 s | 7.34 s | 9.35 s |

   Partitioning is 37% faster and the build 19% faster, with the graph phase
   unchanged, confirming the gain is confined to decompression. Both emitted
   9,714,930 unitigs and 682,165,789 bases. The remaining serial work in a
   reader is record assembly, which bounds the speedup well below the worker
   count. Sequencer and SRA FASTQ is usually single-member, so this helps
   callers who control their own compression — including anyone willing to
   transcode once with `bgzip` or `rebgzf`.
2. **Persist a gzip index.** One serial pass records the bit offset and 32 KiB
   history window at chosen block boundaries; later runs over the same file seek
   and inflate in parallel. This does not help a first run, which is the common
   case for a one-shot assembly, and adds a sidecar file to manage.
3. **Implement rapidgzip-style speculative parallel inflate in Rust.** Workers
   start at arbitrary offsets, find plausible deflate block headers by trial and
   error, decode with an unknown initial window, and reconcile once the true
   window arrives from the preceding chunk. This is the only route that speeds
   up a first pass over an arbitrary single-member gzip, and the published
   speedups are large, but it is a substantial standalone project and a
   correctness-critical one.


## Closing the uncolored partition gap

Uncolored full-corpus builds trailed C++ at 16, 32, and 64 threads while
leading everywhere else. Bisecting both implementations at the same phase
boundaries located the cause precisely, after several plausible explanations
failed to survive measurement.

Diagnostic switches make the phases separable: `CF3_RS_EXIT_AFTER_PARTITION`
and `CF3_EXIT_AFTER_SPLIT` stop each implementation after partitioning,
`CF3_RS_SCAN_ONLY` and `CF3_SCAN_ONLY` suppress record emission, and
`CF3_RS_PARSE_ONLY` additionally suppresses the minimizer scan.

On 10,000 genomes as plain FASTA, the original split was:

| stage | Rust | C++ |
| --- | ---: | ---: |
| parse | 0.708e12 | — |
| minimizer scan | 4.350e12 | — |
| scan total (parse + minimizer) | 5.257e12 | 3.274e12 |
| emit | 1.790e12 | 3.024e12 |
| partition total | 7.634e12 | 6.298e12 |

Rust's emit path was already 41% cheaper than C++'s; the entire deficit, and
more, was in scanning. Four changes addressed it:

| change | partition instructions |
| --- | ---: |
| baseline | 7.634e12 |
| `PEXT` label packing | 7.141e12 |
| single bulk record write | 7.047e12 |
| vectorized line append | 6.848e12 |
| fused scan loop | 5.763e12 |

The last is the substantive one. The scan ran as a per-base window walker
invoking a visitor closure, plus a closure tracking super-k-mer boundaries.
Disassembly showed the closure was never inlined, so its state lived in a
closure environment and was reloaded and stored through memory on every base,
on top of a call. `#[inline(always)]` did not persuade the compiler; merging
the loops did, and the visitor now runs once per weak super-k-mer.

Rust's partition is now cheaper than C++'s: 5.763e12 instructions against
6.298e12, and 106.6 s against 113.5 s on the full corpus at 64 threads.

| threads | Rust before | Rust after | C++ |
| ---: | ---: | ---: | ---: |
| 16 | 15:15.10 | 14:05.75 | 14:42.67 |
| 32 | 8:21.99 | 7:48.26 | 8:02.13 |
| 64 | 5:12.76 | 4:57.06 | 4:48.37 |

Rust leads at 16 and 32 threads and trails by 3.0% at 64. Colored improved to
7:00.71 at 64 threads against C++'s 7:33.78, and read mode to 29.51 s. All
counts exact.

### Rejected along the way

Each was measured and reverted: atlas-lock scope (partition 4% worse), bucket
compression as an uncolored default (7-8% slower at 32-64 threads), retained
bucket-file handles (2% at 32 threads for eight times the descriptors),
jemalloc decay tuning (11.6% slower), a smaller per-file read buffer (no
effect), and a seed-0 wyhash specialization (already constant-folded). Two
hypotheses were disproved outright: C++'s ISA-L inflate accounts for only 8%
of the gap, and page faults concentrate downstream of partitioning rather than
in file reading.

## Closing the collation gap

Collation trailed C++ at every thread count, unlike local contraction which
only differed at 64. Splitting the phase against matched C++ runs on the full
Salmonella corpus showed the deficit was not where it looked:

| sub-phase | C++ t64 | Rust t64 | Δ |
| --- | ---: | ---: | ---: |
| discontinuity contraction | 9.15 | 11.38 | +2.24 |
| expansion | 31.97 | 29.12 | −2.85 |
| path-info map | 4.06 | 5.43 | +1.37 |
| path-info reduce | 4.68 | 4.28 | −0.40 |
| unattributed | — | 4.08 | +4.08 |

Rust already led on expansion and reduce. The largest single item was time no
timer accounted for, and it was flat across thread counts — 4.09 s at 16, 4.56 s
at 32, 4.26 s at 64 — which rules out anything that scales with workers.

Instrumenting the untimed statements between phases found it: `remove_dir_all`
on the edge-matrix block directory cost **3.835 s**, single-threaded, on the
critical path. A second `remove_dir_all` over the expansion buckets cost about
0.9 s and sat *inside* the map timer, which is most of why map appeared slower.
The neighbouring suspects were nothing: dropping the contraction table took
0.004 s and `malloc_trim` 0.000 s.

Nothing downstream reads either directory. `spawn_background_dir_removal`
unlinks them on eight threads — the work is syscall-bound, so a handful of
threads saturates the filesystem without competing for the cores map and reduce
are using — and the caller joins before returning so the work directory is
still clean when the build finishes. The join measures 0.000 s at every thread
count: the removal completes entirely within the phases it now overlaps.

| threads | collation before | collation after | residual before | residual after |
| ---: | ---: | ---: | ---: | ---: |
| 16 | 90.10 | 86.77 | 4.09 | 0.02 |
| 32 | 65.44 | 57.61 | 4.56 | 0.03 |
| 64 | 54.30 | 50.00 | 4.26 | 0.06 |

Post-local total against C++ at 32 threads went from +9.29 s to +1.46 s. Wall
time is 7:38.51 at 32 threads against C++'s 8:02.13, and 4:55.00 at 64 against
4:48.37. Topology is unchanged: 252,487,658 strand-normalized unitigs match the
predecessor binary exactly, on a byte-identical FASTA.

Wall-time deltas are noisier than the phase timers here. The 16-thread run
recorded 14:06.81 against a 14:05.75 baseline despite collation improving 3.3 s,
because partition and local contraction — phases this change does not touch —
happened to run 4.4 s slower on that pass. Run-to-run spread on the earlier
phases is roughly ±4 s at this scale, so the collation timer is the signal.

Trivial (exit-free) local unitigs now go straight into the output FASTA instead
of a temporary file that collation copied in. This was pursued as the suspected
residual and was not it: the copy measured 0.185 s for 3.4 GB because
`std::io::copy` reaches `copy_file_range`. It is kept only because it removes
3.4 GB of redundant temporary write and read, and it is neutral on wall time.

The remaining collation gap at 64 threads is +1.03 s, now dominated by
discontinuity contraction at 12.43 s against C++'s 9.15 s. Its scan and
diagonal halves run concurrently under a `join`, so the two reported timers
overlap and cannot be added.

## Anatomy of the remaining contraction gap

Discontinuity contraction is the last sub-phase where Rust trails: 12.43 s
against C++'s 9.15 s at 64 threads, 13.13 s against 10.61 s at 32. Aligning the
two implementations' sub-timers shows the deficit is narrower than it looks.

| | t16 | t32 | t64 |
| --- | ---: | ---: | ---: |
| Rust table scan | 14.41 | 6.43 | 4.85 |
| C++ non-diagonal contraction | 14.36 | 6.90 | 4.50 |
| Rust diagonal | 4.54 | 4.25 | 4.74 |
| C++ diagonal (computation + contraction) | 4.45 | 2.38 | 3.11 |

Rust's scan matches C++ and scales 2.97x from 16 to 64 threads. The diagonal
does not scale at all. `compress_diagonal_block_with_ends` is a sequential
chain walk: each edge's result depends on `ends` insertions made by earlier
edges, and it runs once per partition on one thread. It cannot be hoisted or
run across partitions either, because `block(p, p)` only becomes final after
every partition above `p` has reinserted into it.

The raw timer overstates the cost. The diagonal runs as one half of a
`rayon::join` against the column scan, so what matters is how much the scan
hides. Timing the join's wall clock alongside both halves:

| threads | scan | diagonal | join wall | exposed |
| ---: | ---: | ---: | ---: | ---: |
| 32 | 7.37 | 4.42 | 8.54 | 1.18 |
| 64 | 4.69 | 4.79 | 6.27 | 1.58 |

Only 1.2-1.6 s of the diagonal is on the critical path; the scan absorbs the
rest. Removing it entirely would leave contraction around 10.8 s at 64 threads,
still above C++'s 9.15 s.

Instrumenting the statements the phase timers skipped found the other half of
the gap, the same way the collation residual was found. Per-partition setup was
negligible at 0.04 s, but concatenating each scan task's meta-vertex vector into
one cost **1.58 s**, single-threaded. Pre-sizing that vector — a partition
contributes on the order of a million meta-vertices, so growth by doubling
copies gigabytes across the run — takes it to 1.23 s, and leaves 0.37 s
unattributed in the phase.

C++ avoids both costs by construction: its workers append diagonal chain edges
to per-worker buckets *during* the parallel non-diagonal scan, then flatten them
with a prefix sum and a parallel memcpy and contract them with a
`parlay::parallel_for`. There is no serial diagonal pass and no serial gather.

Matching that is the remaining work, and it is not cheap. Eliminating the gather
copy means threading nested per-task vectors through five call sites and the
meta-vertex writer, for about 1.2 s. Eliminating the exposed diagonal means
restructuring the contraction core so chains are collected under the atomic
table's slot locks during the scan, for about 1.6 s. Together they are worth
roughly 2.8 s of a 295 s run — under 1% — against real risk to a
correctness-critical path, so they are recorded here rather than attempted.

## The 4-thread read-mode memory gap

An earlier reading of the record had Rust at 1.55x C++'s peak RSS on
`SRR105788` at 4 threads. That number came from a stale binary. Re-measured on
the current build:

| | wall | peak RSS |
| --- | ---: | ---: |
| Rust, as recorded | 1:21.86 | 5.24 GB |
| Rust, current | 1:20.62 | 3.71 GB |
| C++ | 1:55.87 | 3.38 GB |

The real figure is 1.10x, against a 30% lead on wall time. The phase profile
puts the whole peak in local contraction — partitioning peaks at 1.46 GB and
collation never exceeds 0.88 GB — climbing with bucket progress and collapsing
to 93 MB the moment the phase ends.

Three explanations were tested and all three are wrong:

- **Per-block edge write buffers.** `ConcurrentBlockedEdgeWriters` keeps a
  buffer per block, `clear()` retains capacity, and the block count is
  `DEFAULT_VERTEX_PARTITIONS` squared regardless of threads — 16,384 blocks at
  256 KiB is 4 GB of latent staging that a 4-thread run would show and a
  256-thread run would hide. Sweeping the total budget over a 16x range moved
  peak RSS from 3.66 to 3.84 GB, in the wrong direction. Not it.
- **Vertex-map high-water retention.** Each worker carries its map to the next
  subgraph, so it grows to the largest subgraph seen. Disabling reuse entirely
  changed peak RSS from 4.10 to 4.11 GB. Not it.
- **jemalloc page retention.** Rust links jemalloc and C++ does not, which fits
  a diffuse excess with no single structure behind it. Forcing immediate return
  with `dirty_decay_ms:0,muzzy_decay_ms:0` *raised* peak RSS to 3.86 GB from
  3.76 GB and cost 1.7 s. Not it.

The remaining 0.3 GB is unattributed. It is 10% on the one configuration where
Rust is 30% faster, so it is recorded rather than chased further.

### Vertex-map reuse is workload-dependent

The investigation turned up a larger effect than the memory gap. Removing the
reuse entirely made 4-thread read mode *faster*: local contraction 44.67 s to
33.45 s, wall 1:20.91 to 1:09.14, on identical unitig counts. `clear_and_reserve`
calls `HashTable::clear`, which walks the table's capacity rather than its
length, so a map still sized for the largest bucket makes every later, smaller
subgraph pay for that size.

It does not generalize. On the full Salmonella corpus reuse is worth 8.2 s of
local contraction at 64 threads and 7.7 s at 16, and at 256 threads in read mode
the whole phase is 4.4 s and the policy is irrelevant. The effect only appears
where one worker processes thousands of subgraphs in sequence.

A carried-map policy keyed on `LocalBucketGroup::stored_bytes` — reuse only when
the held map is not sized more than N times beyond the incoming subgraph — was
implemented and swept:

| slack | read t4 local | ref t64 local |
| ---: | ---: | ---: |
| always reuse | 44.67 | 136.72 |
| 64 | 44.08 | 135.52 |
| 16 | 34.87 | 140.92 |
| 4 | 35.22 | 139.76 |
| never reuse | 33.45 | 144.93 |

No threshold serves both: the read win needs 16 or below, and those same values
cost 3-4 s on the reference corpus. That `stored_bytes` fails to separate them
suggests it is a poor proxy — bytes per vertex differs sharply between long
low-coverage reference contigs and short high-coverage reads. The change was
reverted rather than landed with a regression.

### A self-calibrating criterion works

Keying on the map's actual capacity against a running mean of the vertex counts
each worker has recently seen removes the need for a per-workload constant.
`LocalVertexMap::capacity` reports what a clear would walk; the carried map
records an exponentially-weighted mean of `len()` across subgraphs, and is kept
only while its capacity stays within eight times that mean. The mean is held
separately from the map so discarding an oversized table does not also discard
what the worker has learned. Small maps below 4096 slots are always carried,
since clearing them costs less than the branch declining to.

Uniform reference buckets keep capacity near the mean and never trip the test.
Read data trips it on the subgraph after an outlier, which is exactly where the
old policy was paying to clear a table sized for a bucket it would not see
again.

| workload | metric | before | after |
| --- | --- | ---: | ---: |
| read, 4 threads | local contraction | 43.87 / 44.15 | 34.59 / 34.90 |
| read, 4 threads | wall | 1:19.05 / 1:21.17 | 1:10.81 / 1:10.34 |
| Salmonella, 64 threads | local contraction | 138.75 / 137.45 / 137.41 | 138.60 / 138.57 / 136.98 |
| Salmonella, 16 threads | local contraction | 408.74 | 408.83 |
| Salmonella colored, 64 threads | wall | 6:57.60 | 6:54.27 |

Read mode at 4 threads improves 21% in the phase and 11.9% end to end. The
reference corpus is flat: 0.13% at 64 threads over three interleaved pairs,
against a 1.3 s spread within the baselines themselves, and 0.02% at 16.
Colored is slightly faster. Topology was checked against the C++ output
directly — 13,309,867 strand-normalized unitigs match — and unitig and
discontinuity-exit counts are identical across every run above.

Two caveats on the record. The 16-thread and colored rows are single pairs
rather than three. And wall time on the reference corpus sits 0.4-0.6% above
baseline while the phase the change actually touches is flat, which is
consistent with variance in the surrounding phases rather than an effect, but
is not separable at this sample count.

## Colour-set encoding

The colour repository is the largest artifact a colored build produces, and it
was four times the size of C++'s: 176.0 GB against 43.6 GB on 149,998
Salmonella assemblies.

Two things about it were initially misread and are worth recording, because
both look like defects and neither is. It is not residue -- it carries a
`metadata.tsv` declaring format, k, encoding and source count, the CLI prints
its path, and it is the colour index that accompanies the FASTA. And the
primary colour table's overflow, 36,810,130 entries past a 2^26 ceiling, has
nothing to do with it: that table is 1 GiB of memory, and admitting every
overflowed colour was already measured to leave wall time and peak RSS
unchanged.

The size came from the encoding. Each set was a varint count followed by varint
gaps, which spends a byte per member. In a pangenome of near-identical
assemblies most k-mers are core, so a typical set names almost every source and
costs over a hundred kilobytes.

C++ avoids this with `fulgor::color_set_builder`, vendored at
`include/color_sets/hybrid.hpp`, which picks one of three regimes by how much of
the source set a colour covers:

| coverage | encoding |
| --- | --- |
| below a quarter | Elias-delta gaps between members |
| a quarter to three quarters | a plain bitmap, one bit per source |
| above three quarters | Elias-delta gaps between *absences* |

The last regime is what matters at pangenome scale, and Elias delta rather than
varints is what makes the first and last cheap: a run of adjacent sources costs
a few bits instead of a byte each.

Implementing the same scheme, with Fulgor's LSB-first bit order:

| | before | after | C++ |
| --- | ---: | ---: | ---: |
| 1,000 assemblies | 135 MB | 58 MB | — |
| 149,998 assemblies | 176.0 GB | 43.4 GB | 43.6 GB |

Source IDs are one-based here, so the thresholds and the bitmap width are sized
by the largest ID plus one, matching C++'s `max_source_id + 1`. Sizing them by
the source *count* silently drops the top source from the bitmap; the colored
compat test caught it.

The cost was 3.4% of wall time, measured over two interleaved pairs at 64
threads -- 6:31.13 against 6:44.29, with the two baselines within 1.1 s of each
other, so an effect rather than noise. Most of it was in the bit writer rather
than the encoding: it drained its accumulator a byte at a time, so a full
64-bit append cost eight `Vec::push` calls. Buffering a whole word and
appending it in one go, plus dropping a per-absence `Option` check from the
complement loop, took two further pairs from 6:44.67 to 6:35.50, a 2.3%
recovery. The net cost of the new encoding is 1.1%.

Fusing the three appends of an Elias delta code into a single register-built
word was tried and reverted: it is 2% *slower* at 1,000 assemblies across two
pairs, despite issuing a third as many appends. The complement walk itself was
also examined and left alone -- it already iterates only the breaks in the
sorted source list, and scanning that list is a lower bound on finding them.

Against C++'s 7:33.78 the build is 12.8% faster while matching its output size.
Peak RSS moves by about 1%.

The repository also now defaults to `<output_prefix>.cf3rs.color-repository`,
beside the FASTA. It was written into `--work-dir`, which callers treat as
scratch; C++ likewise writes its own colour files next to the FASTA as
`<output>.fa.col.N`.

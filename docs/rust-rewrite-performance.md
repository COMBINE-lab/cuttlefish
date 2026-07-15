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
2. Colored partitioning processes a bounded window of 4096 sources. Files are
   scheduled largest-first across at most 64 scan workers. Each worker retains
   atlas-local records until the source barrier.
3. A 128-subgraph atlas is source-partitioned in place. Sixteen long-lived
   consumers merge the sorted worker streams and reuse compression buffers.
   This avoids a second full counting-sort allocation and avoids leaving one
   compression allocation in every worker arena.
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

## Engineering decisions

- Use phase-local bounded worker pools. The compact `k <= 31` path does not keep
  an idle global Rayon pool competing with explicit phase workers.
- Use fast deterministic hash functions for packed integer keys. Rust's default
  collision-resistant hasher is intentionally not used in graph hot paths.
- Keep labels, colors, edges, and path information external. Memory is bounded
  by the largest active bucket/table/window rather than total graph size.
- Preserve source order at the atlas barrier. Spilling worker fragments before
  the barrier both duplicates memory and breaks colored source collation.
- Reuse buffers on long-lived workers. Limiting atlas consumers to 16 lowers
  allocator-arena high-water usage without reducing scan parallelism.
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
| colored | Rust, 2048-source windows | 78.70 s | 25,006,660 KiB |
| uncolored | C++ | 61.70 s | 26,122,016 KiB |
| uncolored | Rust, deferred atlas emission | 84.84 s | 32,879,540 KiB |

Colored Rust is 7.4% slower and uses 15.2% less memory at this rung. This is
near parity, not strict time parity. Uncolored scale parity is not achieved:
the deferred worker-local atlas path reduced partitioning from 10.64 s to
6.59 s and bucket flushes to 49,152, but the measured process remained 37.5%
slower and used 25.9% more memory. Its largest deficits are the local graph
construction, blocked contraction/expansion, and final path-info map/reduce
phases. Enabling
the optional uncolored LZ4 bucket format reduced bucket bytes but did not fix
the deficit (82.59 s, 36,949,396 KiB), because it emits many small compressed
blocks. Do not cite the 1,000-genome result as evidence of uncolored scale
parity.

At 1,000 genomes the same deferred atlas path completed in 18.81 s and used
13,893,160 KiB, versus 16.32 s and 23,481,744 KiB for C++. Partitioning took
1.62 s and performed exactly 16,384 flushes. This confirms the partition
optimization independently of the noisier downstream 5,000-genome phases.

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

Profiling the 5,000-genome colored partition identified a different barrier.
Rust spends about 3.6-4.2 s scanning and packing, followed by 5.8-7.3 s of
window-level atlas collation and compression. C++ drains 64 KiB worker-local
atlas chunks during scanning. Whole-window pipelining was rejected because it
left partition time near 10 s and raised peak RSS to 36.2 GB. Removing source
collation was also rejected: sorting extracted source relationships moved work
into local contraction and increased total time.

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

Colored atlas windows were subsequently increased from 2,048 to 4,096 sources.
Worker-atlas tails already drain at 64 KiB after each source, so the larger
window does not change that memory bound; it removes two of five all-atlas tail
barriers on the 10,000-genome input. An unbounded-window experiment improved
5,000-genome partitioning but was rejected at 10,000 genomes because broader
source mixing enlarged compressed buckets by 2.6% and slowed local contraction.
The 4,096-source compromise kept bucket growth to 0.27%, reduced partitioning
from 10.35 to 9.84 seconds, local contraction from 55.41 to 53.57 seconds, and
process wall time from 126.20 to 121.69 seconds. Peak RSS remained effectively
flat at 15,699,600 KiB and the exact 277,146,851-unitig /
16,550,721,590-base result was preserved. Against the matched 114.63-second C++
run, the remaining colored gap is 7.06 seconds (6.2%).

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

The post-local graph has a measured parallelism crossover. With 1,188,492,662
uncolored discontinuity exits, using 64 workers for the blocked
contraction/expansion/collation path reduced the 5,000-genome process time from
69.06 to 56.41 seconds and RSS from 22,302,284 to 18,897,360 KiB. C++ required
61.70 seconds and 26,122,016 KiB. At the full-corpus 2,840,034,051-exit scale,
64 workers regressed Rust from 169.11 to 182.77 seconds. The uncolored driver
therefore uses at most 64 post-local workers below two billion exits and honors
the requested count above that threshold. Both sides of the policy emit the
same topology: the 5,000-genome 64-worker run matched all 179,767,076 unitigs,
and the full 256-worker run matched all 567,570,784 unitigs against the colored
path.

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

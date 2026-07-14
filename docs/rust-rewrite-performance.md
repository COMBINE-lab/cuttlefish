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
2. Colored partitioning processes a bounded window of 2048 sources. Files are
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

The 5,000-genome topology can be compared exactly with bounded memory:

```bash
python3 scripts/compare_colored_fasta.py \
  --cpp-fasta /scratch3/tmp/humgut2/scale5000/cpp-uncolored3.fa \
  --rust-fasta /scratch3/tmp/humgut2/scale5000/rust-uncolored.fa \
  --topology-only --external-work-dir /scratch3/tmp/cf3-compare -k 31
```

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

Historical full-corpus Rust runs on all 31,225 genomes produced 567,570,784
unitigs and 31,235,773,275 bases. The best pre-parity implementation took
8:40.59 and 22,354,288 KiB RSS. These runs establish scale correctness counts,
but not current performance parity because they predate the final atlas design
and no matched C++ timing log was retained.

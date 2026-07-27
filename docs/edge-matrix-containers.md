# Edge-matrix container files

## The problem

A 150k uncolored build creates roughly 36,600 physical files, and two classes
are 90% of them: 16,384 weak-super-k-mer buckets and 16,641 edge-matrix blocks
(`vertex_partitions + 1` squared, 129x129). Peak descriptor use is only 2,185 --
already twelve times better than C++'s 26,892 -- so this is a *file count*
problem, not a descriptor one. What it costs is inode and directory-metadata
pressure, which dominates on network filesystems, plus open/close churn and
unlink time. Removing the edge-matrix directory alone measured 3.835 s.

Every block access today is an open, a read, and a close. `BucketFile` is worse
still: `open_existing` re-reads and revalidates a header on every 64 KiB flush.

## The access pattern, which is not what it first appears

The edge matrix is a *single* object. Local contraction writes it; discontinuity
contraction reads it and reinserts into it; expansion then reads the
accumulated result. The two readers traverse it orthogonally:

- contraction reads a **column**: blocks `(0..=c, c)`, partitions in reverse order
- expansion reads a **row**: blocks `(r, r+1..N)`

So every block is read exactly once per axis, and no single container axis makes
both sequential. The matrix is upper-triangular, so row `r` holds `N - r` blocks
and column `c` holds `c + 1`.

## Why the axis choice is sharper than it looks

Extents inside a container are interleaved by write order, because blocks
sharing a container are written concurrently. That does not blunt the axis
choice -- it sharpens it.

The favoured phase wants *every* block in the container it is reading, so it can
stream the container front to back in one linear pass and demultiplex by extent;
interleaving costs it nothing. The other phase wants one block from each of ~129
containers, so it does one `pread` per block.

The choice is therefore between *one streaming read* and *one pread per block* --
not between two shades of partial contiguity.

## Options considered

| option | contraction | expansion | cost |
| --- | --- | --- | --- |
| column containers | streams | pread per block | — |
| **row containers** | pread per block | **streams** | — |
| tiled, sqrt(N) x sqrt(N) | partial | partial | gives up streaming entirely |
| write each edge twice | streams | streams | ~+50 GB writes, double write CPU |
| transpose between phases | streams | streams | extra ~50 GB read and write pass |

Tiling is the weakest: once streaming is the property that matters, partial
contiguity buys much less than full contiguity, and it forfeits both.

## Chosen design

Row containers, with the axis switchable through `CF3_RS_EDGE_CONTAINER_AXIS`
so both can be measured rather than argued.

Expansion is the heavier phase -- 28.8 s against contraction's 12.8 s at 64
threads -- so it should be the one that streams. And the phase that loses is
still better off than today: its scattered access becomes a `pread` on a
permanently open descriptor, against today's open, header read, revalidate,
seek and close per block. Neither phase regresses; the choice is only about
which one gains more.

Duplicate writes and the transpose both buy the second phase's streaming at a
real cost, and only pay off if scattered preads turn out to be expensive. That
is measurable once row containers land, and duplicate-writing is additive
rather than a redesign, so it is the fallback rather than the starting point.

## Shape

```rust
struct BlockExtent { offset: u64, len: u32 }

struct EdgeContainers {
    files: Vec<EdgeContainerFile>,   // one per row (or column)
    axis: EdgeContainerAxis,
    partition_count: usize,
}
struct EdgeContainerFile { path: PathBuf, file: File, cursor: AtomicU64 }
```

`BlockedEdgeBlock` drops its `path` and gains `extents: Vec<BlockExtent>`.
Appends reserve with `cursor.fetch_add` and then `write_all_at`, so no
container-wide lock is held and nothing is opened or closed per flush; the
existing per-block buffer mutex still serialises a block's own extent order.

## Call sites

The core is `BlockedEdgeMatrix::create`, `flush_blocked_edge_block`,
`append_blocked_bytes` (deleted), and `read_flushed_block`. Beyond those, 28
sites need converting, in three clusters: the `ConcurrentBlockedAppend` writers
and their five append paths, contraction's chunked task-building, which turns
`BlockedReadTask::File { path, offset, len }` into container-relative extents,
and the two expansion read paths that call `fs::read(&block.path)`.

Cleanup keeps working unchanged: the background directory removal added for
collation now unlinks 129 files instead of 16,641, which should take the
3.835 s measured for that step to near zero.

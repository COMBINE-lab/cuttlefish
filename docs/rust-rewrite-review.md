# Rust rewrite: code review and opportunity inventory

*2026-08-11. A systematic review of `crates/cuttlefish-rs` (~30.5k lines) and
its harness, covering correctness/robustness observations, ranked performance
opportunities, and simplification candidates. Every claim below was verified
against source at review time; file:line references are to the tree at the
review commit. Items are numbered for selection — nothing in section B or C is
implemented until chosen, except where a decision is already recorded.*

**Outcome status.** Section C landed in full (seven commits, ~5,300 lines
removed, no measured regression across uncolored/colored/mbal interleaved
pairs — see "The dead-code removal campaign" in `rust-rewrite-performance.md`).
B1's rapidgzip half landed and measured a 52% partition / 45% wall win on
`mbal` ("Parallel plain-gzip inflation", same document). The BGZF-reader
replacement measured parity and landed: `bgzf.rs` is deleted and all gzip
flows through rapidgzip. The thread-broker half was integrated, measured, and
turned **off by default**: on this host the static split's oversubscription
beats every broker split, and one run in three misconverged ("A closed-loop
split evaluated", same document); it stays behind
`CF3_RS_DYNAMIC_INFLATE_SPLIT` as the evaluation harness. B2 landed: the
per-unitig label allocation is gone (borrowed-label emit into a reused
scratch), measuring -3.0% uncolored / -0.9% colored local contraction with
non-overlapping arms; B5's trivial-path allocation went with it, and the
restructure absorbed C6's remaining dispatch fold. A4 was removed with C3.
B2–B5, A1–A3, A5–A6, and the remaining C5 walker-family consolidation are
open for selection.

The review deliberately does not re-propose anything in the rejected-approaches
record of `rust-rewrite-performance.md`: LTO, mimalloc, container axis changes,
stage-3 containers, fused path-info preads, Elias-append fusion, retained
bucket descriptors, thread-count-derived coordinate fanout, and the rest listed
there remain rejected.

---

## A. Correctness and robustness observations

None of these is a known-wrong-output bug; they are latent-hazard and
hygiene findings.

**A1. Poison-swallowing locks (14 sites in `discontinuity.rs`, 1 in
`color.rs`).** Every `Mutex` acquisition uses
`unwrap_or_else(|poison| poison.into_inner())`. For the edge-block buffer
mutexes this means a worker that panics mid-append leaves a half-written
buffer that the next writer extends — the record-alignment invariant
(`sum(extents) + buffer.len() == edges * record_len`) would be silently
violated and surface three phases later as `MalformedCoordBucket`. In practice
a worker panic already aborts the build (`panic = "abort"` in the release
profile), so the release binary cannot reach the poisoned state; the hazard is
real only for debug/test builds and for any future move away from
`panic = "abort"`. **Recommendation:** leave the sites alone, but add one
comment at the first site explaining that `panic = "abort"` is what makes the
pattern sound, so the invariant is recorded where it can be violated.

**A2. `unsafe impl Sync` justified only by runtime convention.**
`ConcurrentStitchedCoordWriters` (`discontinuity.rs:13317`) hands each rayon
worker an `UnsafeCell` slot indexed by `rayon::current_thread_index()`. The
safety argument is correct today but lives entirely in a comment; nothing
stops a future call site outside the pool. Low priority: convert to
`thread_local`-style per-worker ownership or add a debug assertion on the
index. Same class: `JoinNeighborWriter` (`:10322`), `PathInfoSlot` (`:15662`).

**A3. Unbounded spin-waits in the color table** (`color.rs:479, 508, 547`).
Readers finding a `PENDING_VALUE` slot spin with `std::hint::spin_loop()` and
no escalation; `wait_until_quiescent` likewise. Under oversubscription (or if
an inserting thread is descheduled) this burns cores the inserter needs.
**Remediation candidate:** the libradicl policy
(`libradicl/src/readers.rs:70` @ 9e74b82) — `crossbeam_utils::Backoff`:
`snooze()` escalating spin→yield, then 50 µs sleeps once `is_completed()`.
`crossbeam-utils` is already in the tree transitively via rayon; using it
directly adds no new code, only a manifest line. The windows are normally
nanoseconds, so this is robustness, not measured time — verify no regression
with one colored interleaved pair.

**A4. `partition::source_hash` footgun** (`partition.rs:1238`). A `pub fn`
with the same name as `state::source_hash` but a different algorithm
(`hash_u64` vs `xxh3_64` over 3 bytes) and **zero callers**. If anything ever
calls the wrong one, color signatures silently diverge. Zero-codegen removal;
should go regardless of what else is chosen.

**A5. Uncommented hard cap** — `threads.min(256)` on color-repository workers
(`discontinuity.rs:17045`). The cap is load-bearing (worker index is packed
into 8 bits of `ColorCoordinate`, `state.rs:203-209`), but the connection is
stated nowhere near the cap. One comment.

**A6. `Vec::set_len`-before-read bucket loads** (`discontinuity.rs:14392`
region). Sound as written (`read_exact` fills before any read), but the
pattern is fragile under refactoring. Candidate: `spare_capacity_mut` +
`set_len` after the read, or `read_buf` when stabilized. Zero expected
performance change; do it only if touching the area anyway.

---

## B. Performance opportunities, ranked

Ranking is by expected end-to-end value on the three gate workloads
(150k-uncolored t64, 150k-colored t64, `mbal` read mode). Every item carries
its measurement plan; the interleaved-pair discipline (≥2 pairs, 3 if close,
`--mem-limit 64G`, phase timers not wall) applies throughout.

**B1. Decompression parallelism (the big one — separate design, Stage 2).**
Plain gzip is single-threaded per file (`MultiGzDecoder`, `input.rs:150`);
`mbal` spends 145.9 s of 167.7 s in partitioning at ~147 MB/s compressed per
reader. Colored partitioning is the worst case: `workers = min(threads,
files)` and `inflate_workers = 1` hardwired. The plan: `rapidgzip-core`
`stream_reader()` for plain gzip behind the existing `Box<dyn Read>` seam in
`parse_fragments_borrowed_with` (`input.rs:118-166`); `thread-broker` to
replace the static `inflate_workers = workers/readers` split
(`partition.rs:322`) with closed-loop reallocation between decompression and
sharding; parallel inflate for colored readers (invariant-safe — see the
source-grouping note below). Adoption is decided **conditional on winning**
on `mbal`/`SRR105788`/`ggallus` without regressing either 150k gate.

*Constraint recorded for any colored-partitioning change:* the color-signature
hash requires each source's records to be **contiguous per bucket**
(`add_source_hashed` dedups only against the last source; XOR-combining makes
group order free but interleaving fatal). One worker per source + per-source
drains is what enforces it today. Parallel decompression *within* a file
preserves it; mixed-source batch pools do not.

**B2. Per-unitig label allocation in local contraction.**
`extract_maximal_unitig_compact` (`subgraph.rs:1224`) allocates a fresh
`Vec<u8>` per local unitig — 1.14 G allocations on 150k Salmonella — whose
contents are immediately copied into the shared buffer by
`DiscontinuityInputs::push_unitig` (`discontinuity.rs:331-352`) and freed.
Cyclic unitigs pay twice (`back_walk.label.clone()` + rotation buffer,
`subgraph.rs:1214`). Local contraction is the largest phase (124–132 s at
t64), and this is its hottest allocation. **Candidate:** thread a per-worker
scratch `Vec` through the emit path (`contract_compact_with` already takes a
closure; change `LocalUnitig.label` to borrow or pass `&[u8]` to the sink).
Touches the `LocalUnitig` type used by tests. **Measure:** `local_s` on both
150k gates; allocation count via `CF3_RS_PROFILE_RSS` phases won't show it —
use an interleaved pair and accept only a clear `local_s` win.

**B3. Per-fragment `Instant::now()` on the direct and colored partition
paths** (`partition.rs:524/527` and `:636/639`). Two `clock_gettime` per
sequence fragment; the streamed path already times per batch. On 150k
Salmonella (many small fragments, colored path) this is millions of syscalls
inside the parse callback. **Candidate:** time per file (the loop already has
per-file boundaries), matching the streamed path's granularity. Trivial,
low-risk; measure with one colored pair (`partition_s`).

**B4. Lock-scope items — evidence before action.** Two sites hold a lock
across real work: LZ4 compression runs inside the per-atlas mutex (shared
`CompressionScratch`, `buckets.rs:462-465` — deliberate, commented, and the
lock *is* the serialization design), and container `write_all_at` runs under
the per-block buffer mutex (`spill`, `discontinuity.rs:1385-1396`, ~14.4 M
spills per 150k build). Contention was not measured by this review, and both
designs have recorded justifications. **Do not touch without an off-CPU /
lock-contention profile first** (`perf sched` or eBPF off-CPU on a gated run);
if contention shows, the fixes are per-worker scratch (atlas) and
reserve-then-write-outside-lock (spill, already how the cursor works).

**B5. Per-unitig allocations in the path-info map.**
`map_external_cpp_path_info_range_bucket_owned` allocates
`reverse_complement_label(&label)` or `label.clone()` per unitig plus a
`Vec::new()` for `direct_colors` (`discontinuity.rs:12223-12238`). Path-info
map is 5.4 s at t64 — bounded upside; fold into B2's scratch-reuse pattern if
B2 is done, otherwise skip.

**B6. K>31 pool churn** (note only). `flush_all_with_threads`,
`load_column_into`, `read_flushed_row` spawn a fresh thread set per partition
— 128× per phase — on the K>31 branch only. k≤31 is the overwhelmingly common
case; record and move on.

**B7. Colored bucket source-ID run-length encoding** (note only). ~2 B/record
of staging memory and pre-compression volume; LZ4 already collapses the
on-disk cost because sources are grouped. Format change not justified by the
evidence; revisit only if colored staging memory becomes a limit.

---

## C. Design and simplification inventory

Hard criterion: no noticeable runtime/memory regression. Items marked
**[zero-codegen]** cannot regress anything by construction (they compile to
nothing today).

**C1. Stitch subsystems (decision recorded: remove both).**
(i) The in-memory stitch behind three hardcoded `false` locals
(`discontinuity.rs:7447, 9935, 9936` — verified still present) plus the four
zero-caller `pub` collation entries: ~15-function closure. **[zero-codegen]**
(ii) The endpoint-stitch pipeline rooted at
`stitch_external_discontinuity_paths_to_final_buckets` (`:7547`), its
~20-function closure, and 6 env switches (`CF3_RS_ENDPOINT_STITCH` +
`DIRECT_STITCH_LABEL_READS`, `ORDERED_STITCH_LABELS`,
`NEIGHBOR_STITCH_PATH_INFO`, `DISABLE_INMEM_STITCH_PATH_INFO`,
`INMEM_STITCH_PATH_INFO_BYTES`). **[zero-codegen]** while the switch is unset.
(iii) The `compact_unitigs` plumbing (23 sites through the live
local-contraction sink) — **not** zero-codegen; gets its own interleaved pair.
Total ~2,500–3,000 lines.

**C2. `#[allow(dead_code)]` clusters** — 39 documented items in three groups
(retired partition-table strategies incl. unannotated transitive members
`scan_partition_column_flat`/`_owned`; superseded materialization variants;
superseded write/encode forms) plus 4 orphans. All **[zero-codegen]**. The
two documented as "tests reach this" (`ConcurrentBlockedEdgeWriters::add`,
`add_prepared_edges_absolute`) in fact have zero test references — verify
once more at removal time.

**C3. API tightening. [zero-codegen]** 14 `pub` items with zero references
anywhere (incl. `partition::source_hash` = A4, `VertexState::is_discontinuity`
— re-verified zero refs, `LocalSubgraph::from_bucket_path`,
`coalesce_bucket_manifest`, four `SerialUncoloredCollator::collate_stitched*`
entries which fall with C1); ~19 `pub` → `pub(crate)` candidates; unread
`BuildParams` fields (`vertex_partitions`, `lmtig_buckets`, `gmtig_buckets` —
initialized, validated, never read; `discontinuity.rs:17016` hardcodes
`DEFAULT_VERTEX_PARTITIONS` instead).

**C4. `CF3_RS_LEGACY_UNCOLORED_PARTITION` path** (`partition.rs:145-254`).
Spawns **one OS thread per input file** — 150k threads on the full corpus if
anyone ever sets it. It predates the streamed path and survives only as a
comparison arm. Candidate: remove, or at minimum clamp its producer count.

**C5. Duplicated helpers.** Byte-identical `minimal_rotation` /
`least_rotation_start` (`subgraph.rs:558/567` = `discontinuity.rs:15595/15604`);
`reverse_complement_label` ×3; the triplicated `DirectedKmer` walker family
(`uncolored.rs:201`, `subgraph.rs:439`, `StitchCycleKmer` — the third copy
falls with C1); test-only `BucketEmitter` (production uses `SharedBucketSink`);
unreferenced `contract_compact` (`subgraph.rs:725`). Consolidation is safe for
the rotation/RC helpers (leaf functions); the walker family sits in
contraction inner loops — consolidate only with `#[inline]` verified in
disassembly or skip. `minimizer.rs` is confirmed test-only (one compat-test
consumer); keep, but mark it as reference/test code at module level.

**C6. Diagnostic-switch hygiene.** 24 live `CF3_RS_*` switches remain after
C1 removes 6. The phase-truncation set (`PARSE_ONLY`, `SCAN_ONLY`,
`DECODE_ONLY`, `EXIT_AFTER_PARTITION`) earns its keep for measurement (B1
depends on the first two). `SORT_LOCAL_VERTICES` duplicates an identical
dispatch block in two functions (`subgraph.rs:867-905` vs `:1053-1088`) —
fold into one helper when next touching either.

---

## Harness note

`bench_variant.sh` scraped only the uncolored local-contraction message, so
every colored TSV row had an empty `local_s`; fixed as part of this review
(alternate pattern for "colored local contraction completed in"). The
`CF3_RS_PARSE_ONLY` / `CF3_RS_SCAN_ONLY` diagnostics were re-verified live —
they provide the decompress/parse vs shard attribution B1's measurements need
with no new instrumentation.

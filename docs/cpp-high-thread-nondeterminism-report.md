# Cuttlefish 3: nondeterministic missing unitigs at high thread counts

## Summary

Cuttlefish 3 sometimes emits fewer unitigs from identical input when run with
128 or 256 worker threads. Repeated runs with the same binary, input, and
arguments can produce different record counts, including both incorrect and
exact results. This is silent output corruption: the process exits normally,
but some unitigs are absent.

This is not explained by unitig orientation or output order. Unitigs were
compared after canonicalizing forward and reverse-complement orientations.

## Test configuration

- Input: first 5,000 genomes from HumGut2
- Input list: `/scratch3/tmp/humgut2/scale5000/genomes.list`
- Graph type: reference compacted de Bruijn graph
- `k`: 31
- Minimizer length: 12
- Cutoff: 1
- File descriptor limit: 65,536
- Thread counts showing the problem: 128 and 256

The known-correct result, established by lower-thread C++ runs and independent
Rust output comparison, is:

- Unitigs: `179,767,076`
- Bases: `11,337,993,544`

## Observed results

Identical high-thread C++ runs have produced the following unitig counts:

| Mode | Threads | Unitigs | Missing from expected |
| --- | ---: | ---: | ---: |
| uncolored | 128 | 179,766,461 | 615 |
| uncolored | 128 | 179,767,049 | 27 |
| uncolored | 128 | 179,767,076 | 0 |
| uncolored | 256 | 179,766,543 | 533 |
| colored | 256 | 179,764,872 | 2,204 |
| colored | 128 | both incorrect and exact runs observed | varies |

The variation between repetitions indicates a concurrency bug rather than a
deterministic algorithmic disagreement.

## Reproduction

```bash
ulimit -n 65536
PARLAY_NUM_THREADS=128 build-cpp-compare/src/cuttlefish build \
  --list /scratch3/tmp/humgut2/scale5000/genomes.list \
  --output out \
  --work-dir work \
  --kmer-len 31 \
  --min-len 12 \
  --cutoff 1 \
  --ref
```

Count the emitted records:

```bash
rg -c '^>' out.fa
```

Repeat the build several times with a fresh output and work directory. Add
`--color` to reproduce the colored case. Using 256 workers increases the
probability or magnitude of the mismatch in the runs observed so far.

## Investigation update (2026-07-25)

The `Concurrent_Hash_Table` hypothesis below was tested and **ruled out as the
cause**, though it did identify real undefined behavior that has now been fixed.

### A one-minute reproducer

The 5,000-genome HumGut2 case is not needed. 10,000 Salmonella assemblies from
`/scratch4/rob/blackwell_salmonella` reproduce the defect in about a minute:

```bash
scripts/bench_variant.sh --variant x --impl cpp \
  --bin build-cpp-compare/src/cuttlefish \
  --list /scratch4/rob/cf3-bench/lists/salmonella-10000.list \
  --mode ref --threads 256 --color
```

The correct result is 51,644,203 unitigs / 4,130,164,227 bases. Three identical
colored runs at 256 threads emitted 51,644,119, 51,644,068, and 51,644,038.
Uncolored at 256 threads emitted 51,644,192; read-mode `ggallus` at 256 threads
emitted 33,694,839 against the correct 33,694,954. The defect is not specific to
colored builds or to reference input.

### 64 threads is deterministic and correct

Three repeated uncolored runs at 64 threads each emitted exactly 51,644,203 /
4,130,164,227, and `cuttlefish compare` matched all 51,644,203 strand-normalized
unitigs against Rust. **Low-thread C++ is therefore a valid correctness
reference**; only 128+ thread runs are unusable for topology.

### The discontinuity graph is not where output is lost

Phase counters are bit-identical across repeated 256-thread runs that produced
different final outputs:

| counter | value in every run |
| --- | ---: |
| trivial maximal unitigs | 16,290,319 |
| expected further non-DCC maximal unitigs | 35,353,474 |
| meta-vertices | 35,353,884 |
| ICCs | 410 |
| phantoms | 24,316 |

`16,290,319 + 35,353,474 + 410 = 51,644,203`, the correct total. Contraction and
expansion agree on the graph every time; the discrepancy appears afterwards.

### The defect corrupts collation, not just the final write

Comparing a correct 64-thread run against a 256-thread run is not a pure
record-deletion difference. `cuttlefish compare` reports differing unitig
*content*: a 1,685-base unitig in the good run appears as 1,695 bases in the bad
one, with a 1,463-base common prefix and 223-base common suffix.

What does vary between runs is bucket assignment — "Maximum maximal unitig
bucket size" was 314,777 / 314,678 / 314,824 in three runs, and the edge-bucket
and label-length maxima varied likewise. The collator's own workspace sizing
(`Unitig_Collator::reduce`, `src/Unitig_Collator.cpp:208-224`) is a serial
`std::for_each` and is not itself racy. The remaining suspect is path-ID
assignment during expansion splitting one maximal unitig's pieces across
buckets, which would explain both a changed record count and changed record
content. This has not yet been confirmed.

### What was fixed

`include/Concurrent_Hash_Table.hpp` and `src/Discontinuity_Graph_Contractor.cpp`
were corrected for the undefined behavior described below:

- Values are now written before their key, and the key is published with a
  release store (`store_key`) and read with an acquire load (`load_key`),
  replacing plain concurrent reads and writes of slot keys.
- `insert(..., slot)` / `find(key, slot)` plus `lock_slot` / `unlock_slot` let a
  caller hold the slot lock across a read-modify-write. The contractor's
  non-diagonal and diagonal-chain sites now update `Other_End` under that lock
  and snapshot the previous value, so `form_meta_vertex` and `G.add_edge` still
  run outside the critical section. The diagonal-chain site takes its two slot
  locks one at a time, so no lock-ordering deadlock is possible.

The header's claim that this is safe because "a key is accessed only twice" is
false: the `z.in_same_part()` branch reassigns and re-processes a slot, re-arming
it for a third access.

This changed neither the failure nor the runtime (fixed 65.7-68.3 s versus
original 65.2-68.5 s on the 10,000-genome colored workload), so it is retained on
correctness grounds alone.

## Likely race in `Concurrent_Hash_Table`

The strongest code-level suspect is `Concurrent_Hash_Table`, used by
`Discontinuity_Graph_Contractor`.

### Non-atomic key publication

In `include/Concurrent_Hash_Table.hpp`, the `atomic_read` path reads a plain
`T_key_` outside the slot lock while another worker can write that key under
the lock:

```cpp
if (T[i].key == empty_key_) {
    lock[i].lock();
    if (T[i].key == empty_key_)
        T[i].key = key, T[i].val = val;
    lock[i].unlock();
}

if (T[i].key == key) {
    // ...
}
```

An aligned 64-bit load may be indivisible on the target machine, but a
concurrent plain read and write is still a data race under the C++ memory model
and therefore undefined behavior. Hardware atomicity alone does not establish
the required synchronization or publication ordering.

### Mutable values used after releasing their lock

Both `insert(..., T_val_*& val_add)` and `find()` acquire the slot lock, return
a pointer to the value, and release the lock before the caller accesses it:

```cpp
lock[i].lock();
val_add = &T[i].val;
lock[i].unlock();
return false;
```

`Discontinuity_Graph_Contractor` then mutates the returned value without
holding the slot lock:

```cpp
auto& z = *p_z;
z = Other_End(...);
z.process();
```

The diagonal-chain path similarly obtains mutable pointers with `find()` and
then modifies them in parallel:

```cpp
auto& m_x = *M.find(e.x());
auto& m_y = *M.find(e.y());
// ...
m_x.process();
m_y.process();
```

The hash-table implementation comments that this is safe because a key is
accessed only twice. That invariant should be verified explicitly for every
non-diagonal, compressed-diagonal, false-phantom, and isolated-cordless-cycle
case. If multiple workers can receive the same value pointer, updates can be
lost or partially observed. Even if the twice-accessed invariant holds, the
plain concurrent key reads and writes remain undefined behavior.

## Suggested investigation

1. Make slot keys proper atomics, with release publication after initializing
   the value and acquire loads before reading it.
2. Keep the slot lock held while a caller reads or mutates a value, or replace
   the pointer-returning API with a callback/update operation executed under
   the lock.
3. Assert or instrument the claimed two-access invariant and record the number
   of insertions/updates per key at 128 and 256 threads.
4. Run a smaller reproducer under ThreadSanitizer if the Parlay and dependency
   build supports it.
5. Compare phase-level counts across repeated runs, especially meta-vertices,
   phantom edges, expanded path-info records, and final mapped edge counts, to
   locate the first divergent phase.

## Separate resolved issue

This report is distinct from the earlier C++ segfault caused by reading
external buckets before their output streams had been flushed. That issue was
fixed in commit `f9ed7e7` (`Flush external buckets before concurrent reads`).
The nondeterministic record-count problem described here was observed after
that fix.


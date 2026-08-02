---
title: Resource control
description: Threads, the soft memory budget, working-directory sizing, and the file-descriptor limit.
---

## Threads

`--threads` is an **upper bound** for every parallel phase, not a target. By
default Cuttlefish uses one quarter of the machine's available hardware
threads.

A phase may use fewer workers than you allow. Input file count, the amount of
work actually available, and `--max-memory` can each reduce phase-local
concurrency — a build with two input files has no way to spend sixty-four
parsing threads.

## Memory

`--max-memory` is a **soft** budget in GiB. It bounds replicated per-worker
state, which is the part of the footprint that scales with the worker count.

It is not an operating-system-enforced limit, and it is not a guarantee.
Workload-sized shared tables can impose a higher minimum resident set than the
budget you ask for: if the graph needs a table of a given size, that table is
allocated regardless. Treat `--max-memory` as a lever on worker replication,
not as a cap you can rely on for scheduling.

## Working directory

`--work-dir` holds the external-memory intermediates: partition buckets,
blocked edge matrices, path information, and collation buckets.

Required space depends on input redundancy and graph structure rather than on
input size alone. Put it on **fast local storage** where possible — the
intermediates are written and re-read at every phase boundary, so working-
directory bandwidth is often what bounds a large run.

The default is the system temporary directory, which on many clusters is small
or on a shared filesystem. Set it explicitly for anything but a toy run.

## File descriptors

The implementation adapts its bucket fanout to the process file-descriptor
limit. Raising `ulimit -n` can expose more I/O parallelism on large runs.

You should not normally have to: at startup Cuttlefish raises its own soft
limit to the hard limit, which needs no privileges, and reports it when it
does:

```text
cuttlefish: raised open-file limit from 1024 to 262144
```

If the hard limit itself is low, the build narrows its container plan rather
than failing, so a constrained machine produces the same graph more slowly
rather than producing no graph.

## Bucket compression

`--compress-buckets` (the default) LZ4-compresses uncolored partition buckets,
reducing working-directory footprint. Whether it improves wall time depends on
storage bandwidth against spare CPU: on a fast NVMe with a busy CPU,
`--no-compress-buckets` can be quicker; on a shared network filesystem,
compression usually wins.

---
title: The C++ Cuttlefish 3
description: Where the C++ implementation of Cuttlefish 3 lives, and what its status is.
---

There was a C++ implementation of Cuttlefish 3 before the Rust one. It is no
longer developed, but it has not been discarded.

## Where it is

The `cuttlefish3-cpp` branch:

```bash
git clone https://github.com/COMBINE-lab/cuttlefish.git
cd cuttlefish
git checkout cuttlefish3-cpp
```

It is built the same way as [Cuttlefish 1 and 2](../installation/#from-source),
with CMake.

## Why it is a branch

During the Rust rewrite, the two implementations were developed side by side so
their outputs could be compared directly. A number of **partial correctness
fixes** were made to the C++ Cuttlefish 3 over that period.

Those fixes were incomplete — they were made to keep the comparison honest, not
to finish the C++ implementation — but they were worth keeping. The
`cuttlefish3-cpp` branch preserves that work exactly as it stood when the Rust
implementation took over.

The `rust-rewrite` branch, where Cuttlefish 3 development continues, has the C++
sources removed entirely.

## Should you use it?

Almost certainly not.

- For **published, reproducible work**, use [Cuttlefish 1 or
  2](../overview/) — they are released, packaged on Bioconda, and described by
  peer-reviewed papers.
- For **Cuttlefish 3 features** — colored graphs, positional colors,
  collection-scale external-memory construction — use the [Rust
  implementation](../../rust/introduction/).

The C++ Cuttlefish 3 is documented here so the branch is not a mystery, and so
that anyone reproducing the rewrite's comparisons can find what those
comparisons ran against.

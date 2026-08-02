---
title: The C++ Cuttlefish 3
description: Where the C++ implementation of Cuttlefish 3 lives, and how it relates to the Rust one.
---

There is a C++ implementation of Cuttlefish 3 alongside the Rust one. The two
are developed in parallel, with the eventual plan being for the Rust
implementation to become the standard one.

## Where it is

The `cuttlefish3-cpp` branch:

```bash
git clone https://github.com/COMBINE-lab/cuttlefish.git
cd cuttlefish
git checkout cuttlefish3-cpp
```

It is built the same way as [Cuttlefish 1 and 2](../installation/#from-source),
with CMake.

## Why it is a separate branch

Developing the two side by side means their outputs can be compared directly,
which is how each keeps the other honest — a disagreement between them is a bug
in one of them, and that is a far sharper signal than either implementation
produces alone. A number of **correctness fixes** were made to the C++
Cuttlefish 3 over this period, and the `cuttlefish3-cpp` branch carries them.

The two live on separate branches simply because they are separate codebases
with separate build systems. The `rust-rewrite` branch, where the Rust
Cuttlefish 3 is developed, has the C++ sources removed entirely so that a
checkout of either branch is one implementation and one build.

## Which should you use?

- For **published, reproducible work**, use [Cuttlefish 1 or
  2](../overview/) — they are released, packaged on Bioconda, and described by
  peer-reviewed papers. Neither Cuttlefish 3 has a release yet.
- For **Cuttlefish 3 features** going forward — colored graphs, positional
  colors, collection-scale external-memory construction — the [Rust
  implementation](../../rust/introduction/) is the one heading for standard,
  and is where new work lands.
- The **C++ Cuttlefish 3** is the right thing to reach for if you are
  reproducing or extending the cross-implementation comparisons, or if you need
  to check a result the Rust implementation produced against an independent
  one.

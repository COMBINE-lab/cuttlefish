---
title: Licenses
description: Cuttlefish's own license and those of its bundled dependencies.
---

Cuttlefish itself is **Revised BSD** (BSD 3-Clause) licensed, in both the C++
and the Rust implementations.

## C++ dependencies

The C++ Cuttlefish bundles or links the following:

| Library | License |
| --- | --- |
| [BBHash](https://github.com/rizkg/BBHash) | MIT |
| [Boost C++ Metaprogramming](https://www.boost.org/doc/libs/1_31_0/libs/mpl/doc/index.html) | Boost Software License |
| [compact_vector](https://github.com/gmarcais/compact_vector) | MIT |
| [cxxopts](https://github.com/jarro2783/cxxopts) | MIT |
| [fmt](https://github.com/fmtlib/fmt) | MIT |
| [KMC](https://github.com/refresh-bio/KMC) | GNU GPL 3 |
| [kseq](http://lh3lh3.users.sourceforge.net/kseq.shtml) | MIT |
| [spdlog](https://github.com/gabime/spdlog) | MIT |
| [xxHash](https://github.com/Cyan4973/xxHash) | BSD |

:::caution[KMC and the GPL]
KMC is GPL 3 licensed. That is the most restrictive license in the C++
dependency set, and it governs what you may do with a distributed binary of the
C++ Cuttlefish. The Rust implementation does not use KMC.
:::

## Rust dependencies

The Rust Cuttlefish 3 depends only on crates.io packages, all permissively
licensed (MIT and/or Apache-2.0):

| Crate | Purpose |
| --- | --- |
| `flate2` (with `zlib-rs`) | gzip decompression of input |
| `hashbrown` | hash tables |
| `lz4_flex` | bucket compression |
| `libc` | `fallocate` hole punching, descriptor limits |
| `rayon` | data parallelism |
| `scc` | concurrent containers |
| `xxhash-rust` | hashing |
| `tikv-jemallocator` / `mimalloc` | optional allocators |

For the exact versions and their full license texts, see `Cargo.lock` and run:

```bash
cargo tree
```

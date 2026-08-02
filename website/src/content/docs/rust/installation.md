---
title: Installation
description: Building and installing the Rust Cuttlefish 3 from source.
---

## Requirements

- A 64-bit Linux or macOS system
- Rust 1.85 or newer
- A standard linker and the native build tools for your Rust target
- Enough temporary disk space for the external-memory intermediates

No C++ compiler, CMake, or system bioinformatics libraries are needed. The Rust
implementation has no C++ dependency.

## From source

```bash
git clone https://github.com/COMBINE-lab/cuttlefish.git
cd cuttlefish
git checkout rust-rewrite
cargo build --release
```

The executable is written to `target/release/cuttlefish`. To install it under a
Cargo-style prefix:

```bash
cargo install --path crates/cuttlefish-rs-cli
```

## Allocators

The default build links **jemalloc**, which is what the performance figures were
measured with. The system allocator and mimalloc are also available, though
neither generally performs as well:

```bash
# system allocator
cargo build --release -p cuttlefish-rs-cli --no-default-features

# mimalloc
cargo build --release -p cuttlefish-rs-cli --no-default-features --features mimalloc
```

## As a library

The graph construction pipeline is a library crate, so it can be driven from
Rust directly rather than through the CLI:

```toml
[dependencies]
cuttlefish-rs = "0.0.1"
```

API documentation can be generated locally:

```bash
cargo doc --workspace --no-deps --open
```

:::note
`cuttlefish-rs 0.0.1` on crates.io is a placeholder that holds the name and
validates the publish path. Depend on the git repository until the first real
release.
:::

## Verifying the build

```bash
cargo test --workspace
```

The compatibility test crate asserts exact unitig and base counts against
fixtures checked into the repository, canonicalizing unitig strand and cyclic
rotation before comparing.

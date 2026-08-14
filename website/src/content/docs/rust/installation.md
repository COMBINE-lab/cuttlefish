---
title: Installation
description: Installing Cuttlefish 3 from a prebuilt binary, bioconda, crates.io, or source.
---

## Prebuilt binaries

Every release attaches prebuilt binaries for Linux and macOS, on x86-64 and
arm64, to the [GitHub releases
page](https://github.com/COMBINE-lab/cuttlefish/releases), along with a shell
installer script that picks the right one:

```bash
curl --proto '=https' --tlsv1.2 -LsSf \
    https://github.com/COMBINE-lab/cuttlefish/releases/latest/download/cuttlefish-rs-cli-installer.sh | sh
```

Cuttlefish is also on
[bioconda](https://bioconda.github.io/recipes/cuttlefish/README.html):

```bash
conda install -c bioconda cuttlefish
```

## From crates.io

With a Rust toolchain installed:

```bash
cargo install cuttlefish-rs-cli
```

This installs the single binary `cuttlefish`.

## From source

Requirements:

- A 64-bit Linux or macOS system
- Rust 1.91 or newer
- A standard linker and the native build tools for your Rust target
- Enough temporary disk space for the external-memory intermediates

No C++ compiler, CMake, or system bioinformatics libraries are needed. The Rust
implementation has no C++ dependency.

```bash
git clone https://github.com/COMBINE-lab/cuttlefish.git
cd cuttlefish
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
cuttlefish-rs = "=3.0.0"
```

API documentation can be generated locally:

```bash
cargo doc --workspace --no-deps --open
```

:::note[Pin an exact version]
The library API is not covered by semver: Cuttlefish's major version tracks the
product generation, and a breaking change bumps the minor version — which
cargo's default caret ranges would accept silently. Pin exactly (`=3.0.0`) and
review the changelog before moving the pin.
:::

## Verifying the build

```bash
cargo test --workspace
```

The compatibility test crate asserts exact unitig and base counts against
fixtures checked into the repository, canonicalizing unitig strand and cyclic
rotation before comparing.

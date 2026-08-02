---
title: Overview
description: What the C++ Cuttlefish is, and how Cuttlefish 1 and Cuttlefish 2 relate.
---

Cuttlefish is a fast, parallel, and very low-memory tool for constructing the
compacted de Bruijn graph from sequencing reads or reference sequences. It is
highly scalable in the size of the input data.

[![Anaconda-Server Badge](https://anaconda.org/bioconda/cuttlefish/badges/version.svg)](https://anaconda.org/bioconda/cuttlefish)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/cuttlefish/README.html)

This section documents the **C++ implementation** — Cuttlefish 1 and
Cuttlefish 2, which share one `cuttlefish build` executable. For the Rust
implementation, see [Cuttlefish 3](../../rust/introduction/).

## Cuttlefish 1 and Cuttlefish 2

Both live in the same binary, and which one runs is selected by the arguments:

- Passing **`--read` or `--ref`** invokes **Cuttlefish 2**.
- Passing **neither** invokes **Cuttlefish 1**.

They differ in scope:

| | Cuttlefish 1 | Cuttlefish 2 |
| --- | --- | --- |
| Reference sequences | yes | yes |
| Sequencing reads | no | yes |
| FASTA output | yes | yes |
| GFA 1.0 / GFA 2.0 output | yes | not yet |
| Reduced GFA (tilings) output | yes | not yet |

See [Differences between 1 & 2](../differences/) for the full comparison.

## The papers

The work is described in two papers:

- [Cuttlefish (original)](https://academic.oup.com/bioinformatics/article/37/Supplement_1/i177/6319696)
- [Cuttlefish 2](https://doi.org/10.1186/s13059-022-02743-6)

Please [cite the appropriate one](../../about/citations/) if Cuttlefish
contributed to your work.

## Dependencies

If you install from [Bioconda](../installation/), nothing else is required.
Building from source needs:

- [GCC](https://gcc.gnu.org/) **or** [Clang](https://clang.llvm.org), for C++17 and C11
- [CMake](https://cmake.org/) version 3.14 or newer
- [zlib](https://zlib.net/)
- [bzip2](https://www.sourceware.org/bzip2/)

These are usually already present, and otherwise available from your package
manager:

```bash
# Linux
sudo apt-get install build-essential cmake zlib1g-dev libbz2-dev

# macOS
brew install --with-toolchain llvm
brew install cmake zlib bzip2
```

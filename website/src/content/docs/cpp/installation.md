---
title: Installation
description: Installing the C++ Cuttlefish from Bioconda or from source.
---

## From Bioconda

```bash
conda install -c bioconda cuttlefish
```

The Conda package supports *k* values up to **127**. For larger *k*, install
from source and see [Larger *k*-mer sizes](../large-k/).

## From source

```bash
git clone https://github.com/COMBINE-lab/cuttlefish.git
cd cuttlefish/
mkdir build && cd build/
cmake -DCMAKE_INSTALL_PREFIX=../ ..
make -j 8 install
cd ..
```

Replace `8` in `make -j 8` with the number of threads you want to build with.

This installs Cuttlefish into a `bin` subdirectory inside the project root. To
choose a different location, pass its path as `-DCMAKE_INSTALL_PREFIX`:

```bash
cmake -DCMAKE_INSTALL_PREFIX=<custom_path>/ ..
```

The executable then lands in `<custom_path>/bin/`. Omitting
`-DCMAKE_INSTALL_PREFIX` entirely installs to `/usr/local/bin/`, which may
require `sudo make -j 8 install`.

:::note
A source installation supports *k* up to **63** by default — lower than the
Conda package's 127. To raise it, see [Larger *k*-mer sizes](../large-k/).
:::

## Which branch

Cuttlefish 1 and 2 live on the `cuttlefish-1-2` branch (their integration
history is on `develop-legacy`):

```bash
git clone --branch cuttlefish-1-2 https://github.com/COMBINE-lab/cuttlefish.git
```

The initial C++ implementation of Cuttlefish 3 is preserved on the
`cuttlefish3-cpp` branch — see [The C++ Cuttlefish 3](../cuttlefish3/).

The repository's default branch, `main`, contains only the Rust Cuttlefish 3;
the C++ sources were removed from it, so a checkout of any of these branches
gives you one implementation and one build system.

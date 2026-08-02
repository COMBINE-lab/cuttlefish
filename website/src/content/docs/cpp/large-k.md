---
title: Larger k-mer sizes
description: Raising the maximum k-mer size supported by a source build of the C++ Cuttlefish.
---

The default maximum *k*-mer size for a source installation is **63**. (The
Bioconda package is built for 127.)

## Raising the limit

To set the maximum *k*-mer size capacity to some `MAX_K`, pass
`-DINSTANCE_COUNT=<instance_count>` to `cmake`, where `<instance_count>` is the
number of *k* values Cuttlefish should support:

```text
instance_count = (MAX_K + 1) / 2
```

For example, to support a `MAX_K` of 127:

```bash
cmake -DINSTANCE_COUNT=64 ..
```

The division by two is because Cuttlefish supports only **odd** *k* values,
for theoretical reasons — an even *k* admits a *k*-mer that is its own reverse
complement, which breaks the canonical-vertex correspondence.

## Limits

`MAX_K` is currently supported up to **255**. Please contact the authors if you
need support for a larger one.

## Cost

There is none for smaller *k*. Cuttlefish uses only as many bytes as a *k*-mer
actually requires, rounded up to a multiple of 8. Raising the maximum *k*-mer
size capacity therefore does **not** affect performance at smaller *k* — it
only adds template instantiations to the build.

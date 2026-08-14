---
title: Which version?
description: Choosing between the Rust Cuttlefish 3 and the C++ Cuttlefish 1 and 2.
---

Three implementations carry the Cuttlefish name. This page is the short answer
to which one you want.

## Use Cuttlefish 3 (Rust) if…

- You are starting fresh: Cuttlefish 3 is the **current generation and the
  canonical implementation**, released, on crates.io, and installable from
  Bioconda.
- You need **colored** compacted graphs, with colors that can vary *along* a
  unitig rather than being a property of the whole unitig.
- Your input is a **collection too large to hold in RAM** — Cuttlefish 3 keeps
  its intermediates in external memory throughout, and typically runs in a
  fraction of the memory of Cuttlefish 2.

→ [Rust documentation](../../rust/introduction/)

## Use the C++ Cuttlefish 1 or 2 if…

- You need **GFA output** or **sequence tilings**. Only Cuttlefish 1 emits
  these.
- You need ***k* above 63** — Cuttlefish 1 and 2 support *k* up to 127 from
  Bioconda, or 255 built from source.
- You are **reproducing published results** that used them, your own or
  someone else's.

→ [C++ documentation](../../cpp/overview/)

## At a glance

| | Cuttlefish 1 | Cuttlefish 2 | Cuttlefish 3 (Rust) |
| --- | --- | --- | --- |
| Language | C++ | C++ | Rust |
| Released | yes | yes | yes (3.0.0) |
| Bioconda | yes | yes | yes |
| References (FASTA) | yes | yes | yes |
| Reads (FASTQ) | no | yes | yes |
| FASTA output | yes | yes | yes |
| GFA output | yes | no | no |
| Colored graphs | per-unitig, via GFA paths | no | positional |
| Max *k* | 127 (255 from source) | 127 (255 from source) | 63 |

## Versioning

Cuttlefish's **major version tracks the product generation** — 3.x is
Cuttlefish 3, the successor to Cuttlefish 2, not an ordinary semver major. A
backward-incompatible change to what a user depends on (the output FASTA, the
color repository format, or the command line) bumps the *minor* version and is
called out in the changelog. Rust library dependents should
[pin an exact version](../../rust/installation/#as-a-library).

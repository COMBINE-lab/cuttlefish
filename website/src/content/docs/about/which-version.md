---
title: Which version?
description: Choosing between the Rust Cuttlefish 3 and the C++ Cuttlefish 1 and 2.
---

Three implementations carry the Cuttlefish name. This page is the short answer
to which one you want.

## Use the C++ Cuttlefish 1 or 2 if…

- You need a **released, citable, packaged** tool. Cuttlefish 1 and 2 are
  described by peer-reviewed papers and installable in one line from Bioconda.
- You need **GFA output** or **sequence tilings**. Only Cuttlefish 1 emits
  these.
- You are **reproducing published results**, your own or someone else's.

→ [C++ documentation](../../cpp/overview/)

## Use the Rust Cuttlefish 3 if…

- You need **colored** compacted graphs, with colors that can vary *along* a
  unitig rather than being a property of the whole unitig.
- Your input is a **collection too large to hold in RAM** — Cuttlefish 3 keeps
  its intermediates in external memory throughout.
- You are willing to build from source and to track a moving target.

→ [Rust documentation](../../rust/introduction/)

## At a glance

| | Cuttlefish 1 | Cuttlefish 2 | Cuttlefish 3 (Rust) |
| --- | --- | --- | --- |
| Language | C++ | C++ | Rust |
| Released | yes | yes | not yet |
| Bioconda | yes | yes | no |
| References (FASTA) | yes | yes | yes |
| Reads (FASTQ) | no | yes | yes |
| FASTA output | yes | yes | yes |
| GFA output | yes | no | no |
| Colored graphs | per-unitig, via GFA paths | no | positional |
| Max *k* | 127 (255 from source) | 127 (255 from source) | 63 |

## A caution about Cuttlefish 3

Cuttlefish 3 has no stable release. Its command-line interface and its private
intermediate formats may still change, and the `0.0.x` crates on crates.io hold
the names rather than marking a release.

If your work needs to be reproducible by someone else next year, that argues
for Cuttlefish 2 today.

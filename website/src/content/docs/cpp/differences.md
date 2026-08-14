---
title: Differences between 1 & 2
description: What Cuttlefish 2 added over Cuttlefish 1, and what Cuttlefish 1 still does that 2 does not.
---

Cuttlefish 1 and Cuttlefish 2 ship in the same `cuttlefish build` executable.
Which one runs is decided by the arguments you pass.

## Selecting one

- Passing **`--read` or `--ref`** invokes **Cuttlefish 2**.
- Passing **neither** invokes **Cuttlefish 1**.

## What differs

- **Input types.** Cuttlefish 1 applies only to assembled reference sequences.
  Cuttlefish 2 applies to both sequencing reads and reference sequences.
- **Output formats.** For reference sequences, Cuttlefish 1 supports the GFA
  output formats; Cuttlefish 2 does not support them *yet*. Both emit FASTA.

## At a glance

| | Cuttlefish 1 | Cuttlefish 2 |
| --- | --- | --- |
| Selected by | passing neither `--read` nor `--ref` | passing `--read` or `--ref` |
| Reference sequences (FASTA) | yes | yes |
| Sequencing reads (FASTQ) | no | yes |
| FASTA output | yes | yes |
| GFA 1.0 / 2.0 output | yes | not yet |
| Reduced GFA output | yes | not yet |
| Sequence tilings | yes | not yet |
| Frequency cutoff (`-c`) | n/a | yes |
| Path cover (`--path-cover`) | n/a | yes |

## And Cuttlefish 3?

The Rust [Cuttlefish 3](../../rust/introduction/) — the current generation and
the canonical implementation — is a separate implementation again, with colored
graph construction and a positional notion of color that neither C++ version
has. It does not currently emit GFA. See [which version should I
use?](../../about/which-version/).

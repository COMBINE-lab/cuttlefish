---
title: Cleaning up
description: Removing the intermediates an interrupted build left behind with cuttlefish cleanup.
---

A build stages its external-memory intermediates under `--work-dir` and
unlinks them as it consumes them, so a **successful run leaves the working
directory empty**. A run that is killed, fills the disk, or fails partway
leaves them behind, and they are large: hundreds of gigabytes for a
collection-scale graph. `cuttlefish cleanup` reclaims that space safely.

```text
Usage:
  cuttlefish cleanup -w DIR [OPTION...]

 options:
  -w, --work-dir <arg>   the --work-dir the build was given
  -p, --prefix <arg>     only this build's artifacts, by output name
  -n, --dry-run          report what would be removed, remove nothing
      --include-repository
                         also remove <name>.cf3rs.color-repository, which a
                         completed colored build needs to interpret its FASTA
  -h, --help             print this help
```

```bash
cuttlefish cleanup -w /path/to/work-dir --dry-run   # report, remove nothing
cuttlefish cleanup -w /path/to/work-dir             # remove them
```

## What it will and will not touch

A working directory is usually shared scratch, so "delete everything" is not
an option. Cleanup removes **only** names matching the exact
`<name>.cf3rs.<suffix>` shape that cuttlefish itself produces, against a fixed
suffix table; anything else in the directory is reported and left alone.
Removals are listed largest-first with human-readable sizes. Use
`-p`/`--prefix` to restrict it to one build's artifacts when several builds
share a directory.

Two things are never removed by default:

- **The output FASTA.** After a bailed run it is partial, and whether to keep
  it is yours to decide.
- **The color repository** (`<name>.cf3rs.color-repository/`). It is durable
  output, not an intermediate: a completed colored build's FASTA cannot be
  interpreted without it. Cleanup prints a `keeping dir …` line instead; pass
  `--include-repository` if you really mean to remove one. (The repository is
  written beside the output prefix, so it only appears here at all when the
  output prefix points into the working directory.)

## Inspecting a run instead

Setting the `CF3_RS_KEEP_INTERMEDIATES` environment variable makes a build
retain all of its intermediates — the thing to use when inspecting a run
rather than cleaning up after one. See the [working-directory
table](../output/#working-directory) for what each `.cf3rs.*` file holds.

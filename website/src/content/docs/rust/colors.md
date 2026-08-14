---
title: Querying colors
description: Reading a colored build back with cuttlefish colors — dump, sets, and grep.
---

A colored build leaves two artifacts behind: the unitig FASTA, whose headers
carry packed positional color runs, and the
`<PREFIX>.cf3rs.color-repository/` directory holding the deduplicated,
delta-coded source sets those runs point into. Neither is meant to be read by
hand — `cuttlefish colors` is the supported way in.

```text
Usage:
  cuttlefish colors dump -r DIR -i FASTA [OPTION...]
  cuttlefish colors sets -r DIR [OPTION...]
  cuttlefish colors grep -r DIR -i FASTA (--any-of|--all-of|--none-of LIST) [OPTION...]

 modes:
  dump                   every color run of every unitig, as
                         unitig, vertex range, sources
  sets                   every distinct source set in the repository
  grep                   unitigs with a run matching a source predicate

 common options:
  -r, --repository <arg> the .cf3rs.color-repository directory
  -i, --input <arg>      the colored unitig FASTA (dump, grep)
  -o, --output <arg>     write here instead of stdout
  -z, --gzip             gzip the output (implied by an .gz output name)
      --no-gzip          never gzip, whatever the output is named
      --names            print source paths instead of 1-based ids

 grep options:
      --any-of <arg>     comma-separated source ids; keep runs containing any
      --all-of <arg>     keep runs containing all of these
      --none-of <arg>    drop runs containing any of these
  -c, --count            print only how many unitigs matched
```

All three modes stream: output is TSV with a `#` header line, optionally
gzipped, and closing the pipe early (`| head`) is fine. **Sources are numbered
from 1** in resolved input order; `--names` prints the source paths recorded in
the repository metadata instead. The source-to-path mapping is also in the
repository's `metadata.tsv`.

## `colors sets` — what does the repository hold?

One line per **distinct** source set, in repository order — the natural first
look at a colored build. It needs no FASTA:

```bash
cuttlefish colors sets -r graph.cf3rs.color-repository --names
```

Columns: `#coordinate  worker  index  size  sources`. The coordinate is what
the FASTA's packed runs refer to.

## `colors dump` — every run of every unitig

One line per color run per unitig: `#unitig  vertex_start  vertex_end
sources`. This is the largest output the tool can produce — for a real corpus
a dump is **bigger than the graph itself** — so write it gzipped:

```bash
cuttlefish colors dump -r graph.cf3rs.color-repository -i graph.fa \
    -o dump.tsv.gz
```

## `colors grep` — unitigs by source predicate

Finds unitigs with a color run satisfying a predicate over sources:

```bash
# unitigs carrying source 3 but not source 7
cuttlefish colors grep -r graph.cf3rs.color-repository -i graph.fa \
    --all-of 3 --none-of 7

# how many unitigs does source 12 touch?
cuttlefish colors grep -r graph.cf3rs.color-repository -i graph.fa \
    --any-of 12 --count
```

At least one of `--any-of`, `--all-of`, `--none-of` is required; they combine
conjunctively. For each match it prints the run's own substring rather than
the whole unitig (which can be megabases): `#unitig  vertex_start  vertex_end
label  sources`. Without `--count`, a `N of M unitig(s) matched` summary goes
to stderr.

The predicate is evaluated **once per distinct source set**, not once per run:
grep streams the repository in one sequential pass, records which sets match
in a small per-coordinate bitset, and only then scans the FASTA. On a
150,000-genome build that turns 41 GB of random repository reads into one
sequential pass plus an 11 MB bitset.

## Queries by sequence?

`colors grep` selects by *color*, not by *k*-mer. Querying which sources
contain a given sequence would need a k-mer-to-unitig locator, which is
sketched in the engineering record but not built.

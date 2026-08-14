# cuttlefish-rs-cli

The Cuttlefish 3 command-line program: parallel, external-memory construction
of uncolored and colored compacted de Bruijn graphs from reference sequences or
sequencing reads.

```bash
cargo install cuttlefish-rs-cli
```

This installs one binary, `cuttlefish`. It carries `build`, `compare`,
`colors` (`dump`/`sets`/`grep`), and `cleanup`, plus `help` and `version`:

```bash
# construct a compacted graph
cuttlefish build --ref --seq refs.fa -k 31 --work-dir work --output graph

# construct a colored graph: each listed path is one source color
cuttlefish build --ref --list genomes.list -k 31 --work-dir work \
    --output graph --color

# decide whether two unitig FASTA files describe the same graph, up to
# strand and cyclic rotation
cuttlefish compare -a graph.fa -b other.fa -k 31 --work-dir cmp

# query a colored build's repository: dump runs, list source sets, or
# select unitigs by source predicate
cuttlefish colors grep -r graph.cf3rs.color-repository -i graph.fa --all-of 3

# remove the intermediates an interrupted build left behind
cuttlefish cleanup -w work --dry-run
```

The default build links jemalloc. `--no-default-features` selects the system
allocator, and `--features mimalloc` selects mimalloc.

For the library, see [`cuttlefish-rs`][lib]. For full usage, resource controls,
output formats, and citation, see the [documentation][docs] and the
[repository README][repo].

> Cuttlefish's major version tracks the product generation: a
> backward-incompatible change to the outputs or the command line bumps the
> *minor* version and is called out in the changelog.

Distributed under the BSD 3-Clause License.

[lib]: https://crates.io/crates/cuttlefish-rs
[docs]: https://combine-lab.github.io/cuttlefish/
[repo]: https://github.com/COMBINE-lab/cuttlefish

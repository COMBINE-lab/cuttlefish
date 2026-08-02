# cuttlefish-rs-cli

The Cuttlefish 3 command-line program: parallel, external-memory construction
of uncolored and colored compacted de Bruijn graphs from reference sequences or
sequencing reads.

```bash
cargo install cuttlefish-rs-cli
```

This installs two binaries: `cuttlefish3-rs`, the graph constructor, and
`cf3-compare-fasta`, which compares two unitig FASTA files up to strand and
cyclic rotation.

```bash
cuttlefish3-rs build --ref --seq refs.fa -k 31 --work-dir work --output graph
```

The default build links jemalloc. `--no-default-features` selects the system
allocator, and `--features mimalloc` selects mimalloc.

For the library, see [`cuttlefish-rs`][lib]. For full usage, resource controls,
output formats, and citation, see the [repository README][repo].

> Cuttlefish 3 is under active development. The command-line interface and the
> private intermediate formats may still change before a stable release.

Distributed under the BSD 3-Clause License.

[lib]: https://crates.io/crates/cuttlefish-rs
[repo]: https://github.com/COMBINE-lab/cuttlefish

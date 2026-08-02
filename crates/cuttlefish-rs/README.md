# cuttlefish-rs

The Cuttlefish 3 library: parallel, external-memory construction of uncolored
and colored compacted de Bruijn graphs from reference sequences or sequencing
reads, including collections too large to keep entirely in RAM.

The pipeline partitions input into weak super-k-mer buckets, contracts
independent local subgraphs in parallel, resolves their discontinuities through
a blocked external graph, and reduces the result to maximal unitigs.

For the command-line program, see [`cuttlefish-rs-cli`][cli]. For build
instructions, usage, output formats, and citation, see the [repository
README][repo].

> Cuttlefish 3 is under active development. The API and the private
> intermediate formats may still change before a stable release.

Distributed under the BSD 3-Clause License.

[cli]: https://crates.io/crates/cuttlefish-rs-cli
[repo]: https://github.com/COMBINE-lab/cuttlefish

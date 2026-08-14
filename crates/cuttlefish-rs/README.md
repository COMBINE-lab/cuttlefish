# cuttlefish-rs

The Cuttlefish 3 library: parallel, external-memory construction of uncolored
and colored compacted de Bruijn graphs from reference sequences or sequencing
reads, including collections too large to keep entirely in RAM.

The pipeline partitions input into weak super-k-mer buckets, contracts
independent local subgraphs in parallel, resolves their discontinuities through
a blocked external graph, and reduces the result to maximal unitigs.

For the command-line program, see [`cuttlefish-rs-cli`][cli]. For build
instructions, usage, output formats, and citation, see the
[documentation][docs] and the [repository README][repo].

> **Pin an exact version.** The library API is not covered by semver:
> Cuttlefish's major version tracks the product generation, and a breaking
> change bumps the minor version — which cargo's default caret ranges would
> accept silently. Depend on `cuttlefish-rs = "=X.Y.Z"` and review the
> changelog before moving the pin.

Distributed under the BSD 3-Clause License.

[cli]: https://crates.io/crates/cuttlefish-rs-cli
[docs]: https://combine-lab.github.io/cuttlefish/
[repo]: https://github.com/COMBINE-lab/cuttlefish

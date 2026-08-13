# Changelog

## 0.1.0

First release of the Rust implementation of Cuttlefish 3. It constructs
uncolored and colored compacted de Bruijn graphs from reference sequences or
sequencing reads, in external memory, for odd `k` from 3 through 63.

### Highlights

- **Reference and read graphs**, uncolored and colored, from FASTA, FASTQ, and
  gzip-compressed input.
- **Exact positional colors.** A vertex's color set is the set of inputs that
  contain it; the implementation does not introduce approximation beyond a
  genuine hash collision between distinct observed color sets.
- **Parallel decompression within the thread budget.** The whole gzip family is
  handled by one decoder: BGZF members decode block-parallel, and a plain
  single-member `.gz` decodes through speculative mid-stream splits, so a large
  compressed input is not bounded by one decompressing thread. Decode workers
  come out of `--threads`, never in addition to it.
- **`cuttlefish colors`** reads a colored build's repository back: `dump` for
  every color run of every unitig, `sets` for the distinct source sets, and
  `grep` for unitigs matching a source predicate. All three stream and can gzip
  their output, since a dump is larger than the graph it describes.
- **`cuttlefish cleanup`** removes the intermediates an interrupted build left
  behind, matching only names cuttlefish itself produces and reporting anything
  else it finds.
- **`cuttlefish compare`** decides whether two unitig FASTA files describe the
  same graph, canonicalizing strand and cyclic rotation, sorting to disk so the
  comparison is bounded by chunk size rather than by the graph.

### Performance

Measured against the C++ implementation on the same host and inputs, 64
threads, 149,998 Salmonella assemblies at k = 31:

| | Rust | C++ |
| --- | ---: | ---: |
| uncolored, wall | 4:27.6 | 4:57.6 |
| uncolored, peak RSS | 10.0 GB | 23.3 GB |
| colored, partition phase | 112.6 s | 119.5 s |

Roughly half the peak memory and half the intermediate disk of the C++
implementation, with identical unitig and base counts. At k = 55 on 10,000
assemblies the uncolored build is 40.8 s against 42.2 s, and colored is at
parity.

### Notes

- Odd `k` only, up to 63. The k-mer representation switches from one word to
  two above k = 31; both widths are covered by the test suite.
- Uncolored partition buckets are LZ4-compressed by default;
  `--no-compress-buckets` turns that off. Colored buckets are always compressed.
- A successful build leaves `--work-dir` empty. An interrupted one does not; see
  `cuttlefish cleanup`.
- The color repository format is `cf3rs-color-repository-v2` and sits beside the
  output FASTA, which cannot be interpreted without it.

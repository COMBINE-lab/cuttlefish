# Cuttlefish

[![Anaconda-Server Badge](https://anaconda.org/bioconda/cuttlefish/badges/version.svg)](https://anaconda.org/bioconda/cuttlefish)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/cuttlefish/badges/platforms.svg)](https://anaconda.org/bioconda/cuttlefish)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/cuttlefish/badges/license.svg)](https://anaconda.org/bioconda/cuttlefish)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/cuttlefish/README.html)

Cuttlefish is a fast, parallel, and very lightweight memory tool to construct the compacted de Bruijn graph from sequencing reads or reference sequences. It is highly scalable in terms of the size of the input data.

## Table of contents

- [Overview](#overview)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Output formats](#output-formats)
  - [''Colored'' output for Cuttlefish 1](#colored-output-for-cuttlefish-1)
- [Example usage](#example-usage)
- [Larger _k_-mer sizes](#larger-k-mer-sizes)
- [Differences between Cuttlefish 1 & 2](#differences-between-cuttlefish-1--2)
- [Citations & Acknowledgement](#citations--acknowledgement)
- [Licenses](#licenses)

## Overview

Cuttlefish is a program to produce the compacted de Bruijn graph from sequencing reads or reference sequences.

The papers describing the work are: [Cuttlefish (original)](https://academic.oup.com/bioinformatics/article/37/Supplement_1/i177/6319696) and [Cuttlefish 2](https://doi.org/10.1186/s13059-022-02743-6).

## Dependencies

Cuttlefish can be installed using Bioconda (check [Installation](#installation)).
If installing _from source_, the following are required:

- [GCC](https://gcc.gnu.org/) **or** [Clang](https://clang.llvm.org) compilers for C++17 and C11
- [CMake](https://cmake.org/) (version >= 3.14)
- [zlib](https://zlib.net/)
- [bzip2](https://www.sourceware.org/bzip2/)

These should already be available in your platform; and if not, then these can be easily installed from their sources.
Besides, these should also be available via some package manager for your operating system:

- **Linux**
  
  ```bash
  sudo apt-get install build-essential cmake zlib1g-dev libbz2-dev
  ```

- **MacOS**
  
  ```bash
  brew install --with-toolchain llvm
  brew install cmake zlib bzip2
  ```

## Installation

- From [Bioconda](https://bioconda.github.io/user/install.html):

  ```bash
  conda install -c bioconda cuttlefish
  ```

  The Conda package supports _k_ values up-to 127.
  To use larger _k_ values, please install Cuttlefish from the source.

- From source:

  ```bash
  git clone https://github.com/COMBINE-lab/cuttlefish.git
  cd cuttlefish/
  mkdir build && cd build/
  cmake -DCMAKE_INSTALL_PREFIX=../ ..
  make -j 8 install
  cd ..
  ```

  You may replace `8` in `make -j 8` with the preferred count of threads to use in the installation process.

  This installs Cuttlefish in a sub-directory named `bin`, inside the project root directory.
  To specify a different installation directory, its path may be passed as the value of `-DCMAKE_INSTALL_PREFIX` with the `cmake` command, i.e. you may use `cmake -DCMAKE_INSTALL_PREFIX=<custom_path>/ ..` .
  Then the installed Cuttlefish executable will be found in `<custom_path>/bin/`.
  Skipping `-DCMAKE_INSTALL_PREFIX` entirely will install Cuttlefish in `/usr/local/bin/`, for which `sudo` access might be required (i.e. `sudo make -j 8 install`).

  This installation supports _k_ values up-to `63`.
  To ensure support for larger values, please compile the source with the slight modification described in [Larger _k_-mer sizes](#larger-k-mer-sizes).

## Usage

`cuttlefish build --help` displays the following message (the default `threads` argument is machine-configuration specific):

```txt
Efficiently construct the compacted de Bruijn graph from sequencing reads or reference sequences
Usage:
  cuttlefish build [OPTION...]

 common options:
  -s, --seq arg            input files
  -l, --list arg           input file lists
  -d, --dir arg            input file directories
  -k, --kmer-len arg       k-mer length (default: 27)
  -t, --threads arg        number of threads to use (default: 22)
  -o, --output arg         output file
  -w, --work-dir arg       working directory (default: .)
  -m, --max-memory arg     soft maximum memory limit in GB (default: 3)
      --unrestrict-memory  do not impose memory usage restriction
  -h, --help               print usage

 cuttlefish_1 options:
  -f, --format arg  output format (0: FASTA, 1: GFA 1.0, 2: GFA 2.0, 3:
                    GFA-reduced)

 cuttlefish_2 options:
      --read        construct a compacted read de Bruijn graph (for FASTQ
                    input)
      --ref         construct a compacted reference de Bruijn graph (for
                    FASTA input)
  -c, --cutoff arg  frequency cutoff for (k + 1)-mers (default: refs: 1,
                    reads: 2)
      --path-cover  extract a maximal path cover of the de Bruijn graph

 debug options:
      --vertex-set arg  set of vertices, i.e. k-mers (KMC database) prefix
                        (default: "")
      --edge-set arg    set of edges, i.e. (k + 1)-mers (KMC database) prefix
                        (default: "")

 specialized options:
      --save-mph       save the minimal perfect hash (BBHash) over the vertex
                       set
      --save-buckets   save the DFA-states collection of the vertices
      --save-vertices  save the vertex set of the graph

```

It supports GNU style arguments, `--` for long options, and `-` for short options.
Long options `opt` taking a parameter can be written as `--opt=parameter` or as `--opt parameter`.
Short options `o` taking a parameter is written as `-o parameter`.

The common arguments (for Cuttlefish 1 and 2) are set as following.

- The input files can be passed in any of the following ways (and the options may be mixed together).
  - `-s <data files>`
  - `-l <newline-separated list files of data files>`
  - `-d <directories containing only the data files>`

  Multiple values for each option can be passed as `--seq=s1,s2,...`, `--seq s1 --seq s2 ...`, `-s s1,s2 ...`, or `-s s1 -s s2` (similarly for `list` and `dir`).

  In case of using sequencing reads as input, the files should be in the FASTQ format.
  For reference sequences, those should be in the FASTA format.
  The input files can also be gzipped.
- The _k_-mer length `k` must be odd and within `127` (and `63` if installed from source; see [Larger _k_-mer sizes](#larger-k-mer-sizes) to increase the _k_-mer size capacity beyond these).
The default value is `27`.
- The number of threads `t` is set to a quarter of the number of concurrent threads supported, by default.
The use of high-enough values is recommended.
- Cuttlefish generates two output files:
  - A FASTA / GFA1 / GFA2 file containing the maximal unitigs of the de Bruijn graph (with the extension `.fa` / `.gfa1` / `.gfa2`).
  The GFA output formats are exclusive for Cuttlefish 1.
  - A metadata file containing some structural characteristics of the de Bruijn graph and its compacted form (with the extension `.json`).
- The working directory `w` is used for temporary files created by the process—it is not created by Cuttlefish, and must exist beforehand.
The current directory is set as the default working directory.
- A soft maximum memory-limit `m` (in GB) can be provided to trade-off the RAM usage for faster execution time;
this will only be adhered to if the provided limit is at least the minimum required memory for Cuttlefish, determined internally.
- Memory-usage restrictions can be lifted by using `unrestrict-memory`, trading off extra RAM usage for faster execution time.

Cuttlefish 1 specific arguments are set as following.

- The output formats (`f`) are —
  - `0`: only the maximal unitig (non-branching path) fragments, in FASTA;
  - `1`: the maximal unitigs, their connectivities, and the input sequence tilings, in GFA 1.0;
  - `2`: the maximal unitigs, their connectivities, and the input sequence tilings, in GFA 2.0; and
  - `3`: the maximal unitigs and the input sequence tilings, in GFA-reduced (see [I/O formats](#io-formats)).

Cuttlefish 2 specific arguments are set as following.

- `read` and `ref` are ''input type'' arguments, based on whether you are providing sequencing reads or reference sequences as input, respectively.
- The frequency threshold `c` (of (k + 1)-mers) is set to `2` for read inputs, and `1` for reference inputs, by default.
- `path-cover` is used to construct a maximal vertex-disjoint path cover of the de Bruijn graph, instead of its compacted variant.

### Note

The edge- and / or the vertex-set generation step could produce a high number of temporary files in disk, up-to 2000.
Failure to ensure the capability of opening this many files could produce error messages of the following form:
> Error: Cannot open temporary file ./kmc_00000.bin

The concurrently open file-handle limit for the user running the process can be raised with the following command:

  ```bash
  ulimit -n 2048
  ```

## Output formats

### Cuttlefish 2 output

The currently supported output format is

- The set of the maximal unitigs (non-branching paths) of the de Bruijn graph, in FASTA

Other output formats are currently in the development roadmap.

### Cuttlefish 1 output

The currently supported output formats are —

- The set of the maximal unitigs (non-branching paths) of the de Bruijn graph, in FASTA
- The compacted de Bruijn graph in the [GFA 1.0](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) and the [GFA 2.0](https://github.com/GFA-spec/GFA-spec/blob/master/GFA2.md) formats
- The compacted de Bruijn graph in a ''reduced'' GFA format. It consists of two files, with the extensions: `.cf_seg` and `.cf_seq`:
  - The `.cf_seg` file contains all the maximal unitig fragments of the graph (the segment outputs from GFA, i.e. the `S`-tagged entries), each one with a unique id.
  This file is a list of pairs `<id segment>`.
  - The `.cf_seq` file contains the ''tiling'' of each input sequence, made by the maximal unitig fragments (the paths in GFA 1 / ordered groups in GFA 2, i.e. the `P`- / `O`-tagged entries).
  Each line of the file is of the format `<id tiling>`, where `id` is a unique identifier (name) of this sequence, and `tiling` is a space-separated list of the unitig ids, completely covering the sequence.
  Each unitig id also has a `+` or `-` sign following it, depending on whether the corresponding unitig is present in the canonical or the reverse-complemented form in this tiling order.

  For the example reference file `refs1.fa` (provided in the `data` directory), the output files _may_ look like the following.

  <table>
  <tr><th>Reference file</th></tr>
  <td>
    >ref1

    CGACATGTCTTAG
  </td>
  </table>

  <table>
  <tr><th>Segments file</th><th>Sequence tilings file</th></tr>
  <tr>
  <td>

    1 CGA  
    3 ATGTC  
    6 CTAAGA
  
  </td>
  <td>

    Reference:1_Sequence:ref1 1+ 3- 3+ 6-

  </td></tr>
  </table>

  The only GFA information missing _explictly_ in this format is the links (GFA 1) / edges and gaps (GFA 2), i.e. the `L`- or the `E`- and the `G`-tagged entries.
  These can be readily inferred from the sequence-tilings.
  For example, a tiling \\(u_0 u_1 \dots u_n\\) corresponds to the edge and gap multi-set \\( \{ (u_0, u_1), (u_1, u_2), \dots , (u_{n-1}, u_n)\} \\).
  Whether a pair \\( (u_i, u_{i+1}) \\) is an edge or a gap can be inferred by checking the suffix and the prefix (of length `k - 1`) of the unitigs 
  \\( u_i \\) and \\( u_{i+1} \\), respectively (in their correct orientations, based on their following `+`/`-` signs).
  Note that, a gap is possible in a sequence-tiling only if the sequence contains characters outside of `A`, `C`, `G`, and `T`.
  
  For moderate to large sized genomes, this output format is preferrable to the GFA ones as the GFA formats can be quite verbose for this particular scenario, while the reduced representation provides effitively the same information, while taking much less space.
  For example, for the 7-human genome dataset (experimented with in the manuscripts) and using `k = 31`, the compacted graph takes 112 GB in GFA2, but only 29.3 GB in this reduced format.

### Orientation of the output

Cuttlefish works with the canonical representations of the _k_-mers, i.e. each _k_-mer and its reverse complement are treated as the same vertex in the original graph.
The maximal unitig fragments (the ''segments'' in the GFA-terminology) are always output in their canonical forms—the orientations are guaranteed to be the same across identical executions.

### ''Colored'' output for Cuttlefish 1

In the [GFA](https://github.com/GFA-spec/GFA-spec) output formats for the compacted de Bruijn graph, the graph is represented as a list of the vertices (i.e. the maximal unitigs) and the adjacencies between them.
The output also includes a path-tiling for each individual sequence in the input references, i.e. an ordered list of the maximal unitig ids that completely tile that sequence.
Put differently, the GFA outputs describe a colored de Bruijn graph in the sense that the color information for each vertex (maximal unitig) is encoded in the `P` (GFA 1.0) or the `O` (GFA 2.0) entries (or the tilings in the `.cf_seq` file, in the reduced output).

Throughout the [manuscript (Cuttlefish 1)](https://academic.oup.com/bioinformatics/article/37/Supplement_1/i177/6319696), when we mention the colored de Bruijn graph, we refer to a specific definition of colors.
While this definition is intuitive and natural when constructing the compacted colored de Bruijn graph from a set of reference genomes, it is not the case that the Cuttlefish algorithm allows arbitrary coloring of the _k_-mers in the de Bruijn graph.
Specifically, in the definition adopted herein, the color set of a unitig is the subset of input references <code>s<sub>i<sub>1</sub></sub>, s<sub>i<sub>2</sub></sub>, ..., s<sub>i<sub>l</sub></sub></code> in which the unitig appears.
This color information is implicitly encoded in the path entries of the output GFA files (the `P` entries in GFA 1.0 and the `O` entries in GFA 2.0).
As a result, all unitigs produced by Cuttlefish are ''monochromatic'' under this coloring definition, as a change to the color set internally to a unitig would imply either a branch (which would terminate the unitig) or the start or end of some reference string and a sentinel _k_-mer (which would also terminate the unitig).
If one were constructing the compacted colored de Bruijn graph from raw sequencing reads or from highly-fractured assemblies, then one may wish to adopt a different notion of color, wherein color sets may vary across an individual unitig.

## Example usage

We use _k_ = 3, and 4 CPU threads, with a working directory named `temp` in the following examples.

### Using Cuttlefish 2

- **From FASTQ files**

To construct the maximal unitigs of the example FASTQ file `reads.fq` (provided in the `data` directory) with frequency cutoff `c = 1`, the following may be used.

```bash
cuttlefish build -s reads.fq -k 3 -t 4 -o cdbg -w temp/ --read -c 1
```

- **From FASTA files**

To construct the maximal unitigs of the example FASTA file `refs1.fa` (provided in the `data` directory), the following may be used.

```bash
cuttlefish build -s refs1.fa -k 3 -t 4 -o cdbg -w temp/ --ref
```

These executions will produce two output files each: `cdbg.fa`, containing the maximal unitigs of the graph; and `cdbg.json`, a metadata file with some structural characteristics of the graph.

Multiple seq-files, lists of seq-files, or directories of seq-files may also be passed, as described in [Usage](#usage).

### Using Cuttlefish 1

To output the compacted de Bruijn graph (in GFA 2.0) for the example FASTA files `refs1.fa` and `refs2.fa` (provided in the `data` directory), the following may be used:

```bash
cuttlefish build -s refs1.fa,refs2.fa -k 3 -t 4 -o cdbg.gfa2 -f 2 -w temp/
```

You may also provide lists or directories of reference files as input, as described in [Usage](#usage).

## Larger _k_-mer sizes

The default maximum _k_-mer size supported with the installation from source is `63`.
To set the maximum _k_-mer size capacity to some `MAX_K`, add `-DINSTANCE_COUNT=<instance_count>` with the `cmake` command—where `<instance_count>` is the number of `k`-values that are to be supported by Cuttlefish, and should be set to `(MAX_K + 1) / 2`.
For example, to support a `MAX_K` of 127, use the following:

```bash
cmake -DINSTANCE_COUNT=64 ..
```

Cuttlefish supports only the odd `k` values within `MAX_K` due to theoretical reasons.
Currently, `MAX_K` is supported upto 255.
Please contact the authors if support for a larger `MAX_K` is required.

Note that, Cuttlefish uses only as many bytes as required (rounded up to multiples of 8) for a _k_-mer. Thus, increasing the maximum _k_-mer size capacity through setting large values for `MAX_K` does not affect the performance for smaller _k_-mer sizes.

## Differences between Cuttlefish 1 & 2

- Cuttlefish 1 is applicable only for assembled reference sequences.
Whereas Cuttlefish 2 is applicable for both sequencing reads and reference sequences.
- For reference sequences, Cuttlefish 1 supports outputting the compacted graph in the GFA formats, whereas Cuttlefish 2 does not support this _yet_.
- Cuttlefish 2 can be used by passing either one of the following arguments to the `cuttlefish build` command: `--read` or `--ref`.
Passing neither of these invokes Cuttlefish 1.

## Citations & Acknowledgement

If you use Cuttlefish or Cuttlefish 2 in your work, please include the following citations, as appropriate:

### [Cuttlefish (original)](https://doi.org/10.1093/bioinformatics/btab309)

> Jamshed Khan, Rob Patro, Cuttlefish: fast, parallel and low-memory compaction of de Bruijn graphs from large-scale genome collections, Bioinformatics, Volume 37, Issue Supplement_1, July 2021, Pages i177–i186, <https://doi.org/10.1093/bioinformatics/btab309>

### [Cuttlefish 2](https://doi.org/10.1186/s13059-022-02743-6)

> Khan, J., Kokot, M., Deorowicz, S. et al. Scalable, ultra-fast, and low-memory construction of compacted de Bruijn graphs with Cuttlefish 2. Genome Biol 23, 190 (2022). <https://doi.org/10.1186/s13059-022-02743-6>

This work is supported by _NIH R01 HG009937_, and by _NSF CCF-1750472_, and _CNS-1763680_.

## Licenses

- The [BBHash](https://github.com/rizkg/BBHash) library is MIT licensed.
- The [Boost C++ Metaprogramming](https://www.boost.org/doc/libs/1_31_0/libs/mpl/doc/index.html) library is Boost Software licensed.
- The [compact_vector](https://github.com/gmarcais/compact_vector) library is MIT licensed.
- The [cxxopts](https://github.com/jarro2783/cxxopts) library is MIT licensed.
- The [fmt](https://github.com/fmtlib/fmt) library is MIT licensed.
- The [KMC](https://github.com/refresh-bio/KMC) software is GNU GPL 3 licensed.
- The [kseq](http://lh3lh3.users.sourceforge.net/kseq.shtml) library is MIT licensed.
- The [spdlog](https://github.com/gabime/spdlog) library is MIT licensed.
- The [xxHash](https://github.com/Cyan4973/xxHash) library is BSD licensed.
- Cuttlefish itself is Revised BSD licensed.

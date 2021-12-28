# Cuttlefish

[![Anaconda-Server Badge](https://anaconda.org/bioconda/cuttlefish/badges/version.svg)](https://anaconda.org/bioconda/cuttlefish)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/cuttlefish/badges/platforms.svg)](https://anaconda.org/bioconda/cuttlefish)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/cuttlefish/badges/license.svg)](https://anaconda.org/bioconda/cuttlefish)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/cuttlefish/README.html)

Cuttlefish is a fast, parallel, and very lightweight memory tool to construct the compacted de Bruijn graph from sequencing reads or reference sequences, which is highly scalable in terms of the size of the input data.

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

The papers describing the work are: [Cuttlefish 1](https://academic.oup.com/bioinformatics/article/37/Supplement_1/i177/6319696) and [Cuttlefish 2](https://doi.org/10.1101/2021.12.14.472718) (pre-print).

## Dependencies

<!-- The earlier version, Cuttlefish 1, can be installed using Bioconda (check [Installation](#installation)). -->
To install Cuttlefish from source, the following are required:

- [GCC](https://gcc.gnu.org/) compilers for C++14 and C11
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

Cuttlefish also makes use of the [KMC 3](https://github.com/refresh-bio/KMC) tool.
If you are installing Cuttlefish from the source, then it will be automatically installed.
To use with Cuttlefish 1 while it is installed using `conda`, you may use the following to install KMC 3:

- From [Bioconda](https://bioconda.github.io/user/install.html):

  ```bash
  conda install -c bioconda kmc
  ```

## Installation

- From [Bioconda](https://bioconda.github.io/user/install.html) (only for Cuttlefish 1, _for now_):

  ```bash
  conda install -c bioconda cuttlefish
  ```

  The Conda package supports _k_ values up-to 127.
  To use larger _k_ values, please install Cuttlefish from the source.

- From source (works for both Cuttlefish 1 and 2):

  ```bash
  git clone https://github.com/COMBINE-lab/cuttlefish.git
  cd cuttlefish/ && mkdir build && cd build/
  cmake -DCMAKE_INSTALL_PREFIX=../ ..
  make -j 8 install
  cd ..
  ```

  You may replace `8` in `make -j 8` with the preferred count for threads to use in the installation process.

  This installs Cuttlefish in a sub-directory named `bin`, inside the project root directory.
  To specify a different installation directory, its path may be passed as the value of `-DCMAKE_INSTALL_PREFIX` with the `cmake` command, i.e. you may use `cmake -DCMAKE_INSTALL_PREFIX=<custom_path>/ ..` .
  Then the installed Cuttlefish executable will be found in `<custom_path>/bin/`.
  Skipping `-DCMAKE_INSTALL_PREFIX` entirely will install Cuttlefish in `/usr/local/bin/`, for which `sudo` access might be required (i.e. `sudo make -j 8 install`).

  This installation supports _k_ values up-to `63`.
  To ensure support for larger values, please compile the source with the slight modification described in [Larger _k_-mer sizes](#larger-k-mer-sizes).

## Usage

`cuttlefish build --help` displays the following message:

```bash
Efficiently construct the compacted de Bruijn graph from sequencing reads or reference sequences
Usage:
  cuttlefish build [OPTION...]

 common options:
  -r, --refs arg      input files (default: "")
  -l, --lists arg     input file lists (default: "")
  -d, --dirs arg      input file directories (default: "")
  -k, --kmer-len arg  k-mer length (default: 25)
  -t, --threads arg   number of threads to use (default: 1)
  -o, --output arg    output file (default: "")
  -w, --work-dir arg  working directory (default: .)
  -h, --help          print usage

 cuttlefish 1.0 options:
  -s, --kmc-db arg  set of vertices, i.e. k-mers (KMC database) prefix
                    (default: .)
  -f, --format arg  output format (0: txt, 1: GFA 1.0, 2: GFA 2.0, 3:
                    GFA-reduced) (default: 0)
      --rm          remove the KMC database

 cuttlefish 2.0 options:
      --read               construct a compacted read de Bruijn graph
      --ref                construct a compacted reference de Bruijn graph
  -c, --cutoff arg         frequency cutoff for (k + 1)-mers (default: 2)
  -m, --max-memory arg     soft maximum memory limit (in GB) (default: 3)
      --unrestrict-memory  do not impose memory usage restriction
      --path-cover         extract a maximal path cover of the de Bruijn
                           graph

 debug options:
  -e, --edge-db arg  set of edges, i.e. (k + 1)-mers (KMC database) prefix
                     (default: "")

 specialized options:
      --mph arg        minimal perfect hash (BBHash) file (optional)
                       (default: "")
      --buckets arg    hash table buckets (cuttlefish) file (optional)
                       (default: "")
      --save-vertices  save the vertex set of the graph

```

### Cuttlefish 2

To construct a compacted de Bruijn graph, use Cuttlefish as following:

```bash
cuttlefish build <input_type> <input_files> -k <k-mer_length> -c <cutoff_frequency> -o <output_prefix> -t <thread_count> -w <working_directory>
```

The arguments are set as following:

- The `<input_type>` argument should be either `--read` or `--ref`, based on whether you are providing sequencing reads or reference sequences as input, respectively.
- The input files `<input_files>` can be passed in any of the following ways (and the options may be mixed together).
  - `-r <comma-separated list of data files>`
  - `-l <newline-separated list file of data files>`
  - `-d <directory containing only the data files>`

  In case of using sequencing reads as input, the files should be in the FASTQ format.
  For reference sequences, those should be in the FASTA format.
  The input files can also be possibly gzipped.
- The _k_-mer length `k` must be odd and within `63` (see [Larger _k_-mer sizes](#larger-k-mer-sizes) to increase the _k_-mer size capacity beyond this).
The default value is `25`.
- The frequency threshold `c` is set to `2` by default, and should be set to `1` when passing reference sequences as input.
- Cuttlefish 2 generates two output files with the prefix `<output_prefix>`:
  - A FASTA file containing the maximal unitigs of the de Bruijn graph (with the extension `.fa`).
  - A metadata file containing some structural characteristics of the de Bruijn graph and its compacted form (with the extension `.json`).
- The number of threads `t` is set to `1` by default, and the use of higher values is recommended.
- The working directory `w` is used for temporary files created by the process—it is not created by Cuttlefish, and must exist beforehand.
The current directory is set as the default working directory.

Some other useful arguments:

- `--path-cover` to construct a maximal vertex-disjoint path cover of the de Bruijn graph, instead of its compacted variant
- `-m <max_memory>` to pass a soft maximum memory-limit (in GB) to trade-off RAM usage for faster execution time; this will only be adhered to if the provided limit is larger than the minimum required memory for Cuttlefish, determined internally
- `--unrestrict-memory` to not impose any memory-usage restriction, trading off RAM usage for faster execution time

### Cuttlefish 1

Unlike Cuttlefish 2, Cuttlefish 1 does not execute KMC 3 by itself (_for now_).
To produce the _k_-mer set from an individual input reference sequence using KMC 3, the following may be used:

```bash
kmc -k<k-mer_length> -fm -ci1 -t<thread_count> <input_reference> <output_set> <working_directory>
```

If working with multiple references, you may use:

```bash
kmc -k<k-mer_length> -fm -ci1 -t<thread_count> @<input_reference_list> <output_set> <working_directory>
```

The input file `<input_reference>` or the files listed in `<input_reference_list>` should be in the FASTA format, possibly gzipped.
The `k` value should be odd (required by Cuttlefish), and is `25` by default.
Having executed, KMC 3 will produce two files with the same prefix `<output_database>`, and extensions `.kmc_pre` and `.kmc_suf`.
When working within strict memory limits, you should add the arguments `-m<max_memory> -sm` with these invocations, where `<max_memory>` is your memory budget in gigabytes.

<!-- It supports GNU style arguments, -- for long options, and - for short options. Long options that take a parameter can be written with an equals sign as --long=parameter or without as --long parameter, and short options can be written together, with the last option able to take a parameter, such as in tar -xvf file.tar.
Follow https://github.com/jarro2783/cxxopts/wiki -->

Then to build the compacted de Bruijn graph, use Cuttlefish as following:

```bash
cuttlefish build <input_references> -k <k-mer_length> -s <k-mer_set> -o <output_file> -f <output_format> -t <thread_count> -w <working_directory>
```

The arguments are set as following:

- The input references can be passed in any of the following ways (and the options may be mixed together).
  - `-r <comma-separated list of refs>`
  - `-l <newline-separated list file of refs>`
  - `-d <directory containing only the refs>`
  
  Each input reference should be in the FASTA format, possibly gzipped.
- The _k_-mer length `k` must be odd and within `63` (or, `127` if you install Cuttlefish using `conda`; see [Larger _k_-mer sizes](#larger-k-mer-sizes) to increase the _k_-mer size capacity beyond these).
The default value is `25`.
- The _k_-mer set prefix `s` must match exactly the output path used in the `kmc` invocation, i.e. it should be the `<output_set>` argument from the `kmc` invocation.
- The output formats (`f`) are —
  - `0`: only the maximal unitig (non-branching path) fragments;
  - `1`: GFA 1.0;
  - `2`: GFA 2.0; and
  - `3`: GFA-reduced (see [I/O formats](#io-formats)).
- The number of threads `t` is set to `1` by default, and the use of higher values is recommended.
- The working directory `-w` is used for temporary files created by the process—it is not created by Cuttlefish, and must exist beforehand.
The current directory is set as the default working directory.

## Output formats

### Cuttlefish 2

The currently supported output format is

- The set of the maximal unitigs (non-branching paths) from the original de Bruijn graph, in FASTA

Other output formats are currently in the development roadmap.

### Cuttlefish 1

The currently supported output formats are —

- The set of the maximal unitigs (non-branching paths) from the original de Bruijn graph, in plain text
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
  For example, a tiling <code><seq_id u<sub>0</sub> u<sub>1</sub> ... u<sub>n</sub>></code> corresponds to the edge and gap multi-set <code>{(u<sub>0</sub>, u<sub>1</sub>), (u<sub>1</sub> u<sub>2</sub>), ... , (u<sub>n-1</sub>, u<sub>n</sub>)}</code>.
  Whether a pair <code>(u<sub>i</sub>, u<sub>i+1</sub>)</code> is an edge or a gap can be inferred by checking the suffix and the prefix (of length `k - 1`) of the unitigs <code>u<sub>i</sub></code> and <code>u<sub>i+1</sub></code>, respectively (in their correct orientations, based on their following `+`/`-` signs).
  Note that, a gap is possible in a sequence-tiling only if the sequence contains characters outside of `A`, `C`, `G`, and `T`.
  
  For moderate to large sized genomes, this output format is preferrable to the GFA ones—the GFA formats can be quite verbose for this particular scenario, while the reduced representation provides effitively the same information, while taking much lesser space.
  For example, for the 7-human genomes (experimented with in the manuscripts) and using `k = 31`, the compacted graph takes 112 GB in GFA2, while 29.3 GB in this reduced format.

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

### Cuttlefish 2

_To be completed_

### Cuttlefish 1

Please use the `kmc` and the `cuttlefish` executables from their respective paths in the following examples.
We use _k_ = 3, and 4 CPU threads, with a working directory named `temp` in the following examples.

- **For individual input genome reference**
  
  To output the compacted de Bruijn graph (in GFA 1.0) for the example FASTA file `refs1.fa` (provided in the `data` directory), the following may be used:

  - Generate the _k_-mer set:

    ```bash
    kmc -k3 -fa -ci1 -t4 refs1.fa kmers temp/
    ```

  - Output the compacted graph (in GFA 1.0):

    ```bash
    cuttlefish build -r refs1.fa -k 3 -s kmers -t 4 -o cdbg.gfa1 -f 1 -w temp/
    ```

    To get only the maximal unitig fragments (which is `-f 0` by default):

    ```bash
    cuttlefish build -r refs1.fa -k 3 -s kmers -t 4 -o cdbg.txt -w temp/
    ```
  
- **For multiple input genome references**

  To output the compacted de Bruijn graph (in GFA 2.0) for the example FASTA files `refs1.fa` and `refs2.fa` (provided in the `data` directory), the following may be used:

  - Produce a newline-separated list of the paths of the input references. For example,

    ```bash
    readlink -f refs1.fa > refs.lst
    readlink -f refs2.fa >> refs.lst
    ```

  - Generate the _k_-mer set:

    ```bash
    kmc -k3 -fa -ci1 -t4 @refs.lst kmers temp/
    ```

  - Output the compacted graph (in GFA 2.0):

    ```bash
    cuttlefish build -l refs.lst -k 3 -s kmers -t 4 -o cdbg.gfa2 -f 2 -w temp/
    ```

    Or,

    ```bash
    cuttlefish build -r refs1.fa,refs2.fa -k 3 -s kmers -t 4 -o cdbg.gfa2 -f 2 -w temp/
    ```

## Larger _k_-mer sizes

The default maximum _k_-mer size supported with the installation from source is `63`.
To set the maximum _k_-mer size capacity to some `MAX_K`, add `-DINSTANCE_COUNT=<instance_count>` with the `cmake` command—where `<instance_count>` is the number of `k`-values that are to be supported by Cuttlefish, and should be set to `(MAX_K + 1) / 2`.
For example, to support a `MAX_K` of 127, use the following:

```bash
cmake -DINSTANCE_COUNT=64 ..
```

Cuttlefish supports only the odd `k` values within `MAX_K` due to theoretical reasons.
Currently, KMC3 supports a `MAX_K` of `255`.
<!-- Also, to increase `MAX_K` beyond 255, the maximum supported _k_ should also be updated with KMC3. -->
<!-- Please note that, the second step of the pipeline, i.e. the construction of a minimal perfect hash function (using [BBHash](https://github.com/rizkg/BBHash)) gets less efficient (time-wise) with increasing _k_, due to disk-read throughput bottlenecks associated with reading the _k_-mers from the KMC database. -->

Note that, Cuttlefish uses only as many bytes as required (rounded up to multiples of 8) for a _k_-mer as necessary—thus increasing the maximum _k_-mer size capacity through setting large values for `MAX_K` does not affect the performance for smaller _k_-mer sizes.

<!-- ## Intermediate disk usage

The Cuttlefish pipeline uses a non-trivial amount of intermediate disk space, in the forms of — the _k_-mer set produced by KMC3, and temporary files produced during the minimal perfect hash construction and the GFA output constructions. The produced KMC3 database (the `.kmc_pre` and the `.kmc_suf` extension files) is not removed automatically (unless the flag `--rm` is passed to the graph build), and can be safely removed by the user after the Cuttlefish execution. -->

## Differences between Cuttlefish 1 & 2

- Cuttlefish 1 is applicable only for (whole-genome or transcriptome) reference sequences.
Whereas Cuttlefish 2 is applicable for both sequencing reads and reference sequences.
- For reference sequences, Cuttlefish 1 supports outputting the compacted graph in the GFA formats, whereas Cuttlefish 2 does not support this _yet_.
- Cuttlefish 2 can be used by passing either one of the following arguments to the `cuttlefish build` command: `--read` or `--ref`.
Passing neither of these invokes Cuttlefish 1.

## Citations & Acknowledgement

```bibtex
@article{10.1093/bioinformatics/btab309,
  author = {Khan, Jamshed and Patro, Rob},
  title = "{Cuttlefish: fast, parallel and low-memory compaction of de Bruijn graphs from large-scale genome collections}",
  journal = {Bioinformatics},
  volume = {37},
  number = {Supplement\_1},
  pages = {i177-i186},
  year = {2021},
  month = {07},
  issn = {1367-4803},
  doi = {10.1093/bioinformatics/btab309},
  url = {https://doi.org/10.1093/bioinformatics/btab309},
}
```

```bibtex
@article{Khan2021.12.14.472718,
  author = {Khan, Jamshed and Kokot, Marek and Deorowicz, Sebastian and Patro, Rob},
  title = "{Scalable, ultra-fast, and low-memory construction of compacted de Bruijn graphs with Cuttlefish 2}",
  elocation-id = {2021.12.14.472718},
  year = {2021},
  doi = {10.1101/2021.12.14.472718},
  publisher = {Cold Spring Harbor Laboratory},
  URL = {https://www.biorxiv.org/content/early/2021/12/16/2021.12.14.472718},
  journal = {bioRxiv}
}
```

This work is supported by _NIH R01 HG009937_, and by _NSF CCF-1750472_, and _CNS-1763680_.

## Licenses

- The [BBHash](https://github.com/rizkg/BBHash) library is MIT licensed.
- The [Boost C++ Metaprogramming](https://www.boost.org/doc/libs/1_31_0/libs/mpl/doc/index.html) library is Boost Software licensed.
- The [compact_vector](https://github.com/gmarcais/compact_vector) library is MIT licensed.
- The [cxxopts](https://github.com/jarro2783/cxxopts) library is MIT licensed.
- The [fmt](https://github.com/fmtlib/fmt) library is MIT licensed.
- The [(ghc) filesystem](https://github.com/gulrak/filesystem) library is MIT licensed.
- The [KMC](https://github.com/refresh-bio/KMC) software is GNU GPL 3 licensed.
- The [kseq](http://lh3lh3.users.sourceforge.net/kseq.shtml) library is MIT licensed.
- The [spdlog](https://github.com/gabime/spdlog) library is MIT licensed.
- The [xxHash](https://github.com/Cyan4973/xxHash) library is BSD licensed.
- Cuttlefish is Revised BSD licensed.

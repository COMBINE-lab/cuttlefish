# Cuttlefish

Building the compacted de Bruijn graph efficiently from set of references.

## Dependencies

Cuttlefish uses the [KMC3](https://github.com/refresh-bio/KMC) tool, which is a disk-based _k_-mer counter.

## Compile

(out-of-source build)

```C++
mkdir build
cd build
cmake ..
make install
```

## Usage

Using KMC3 to produce the _k_-mers set:

```Bash
kmc -k<k-mer_len> -m<max_mem> -fm -ci0 -t<thread_count> <input_file_name> <output_file_name> <working_directory>
```

If working with colored de Bruijn graphs (i.e. multi-references), use:

```Bash
kmc -k<k-mer_len> -m<max_mem> -fm -ci0 -t<thread_count> <@input_file_list> <output_file_name> <working_directory>
```

Then to build the (colored) compacted de Bruijn graph, use Cuttlefish as follows:

(from project directory)

```Bash
./bin/cuttlefish build <options>
```

**Options**:

```Bash
Usage:
  cuttlefish build [OPTION...]

  -r, --refs arg      reference files (optional)
  -l, --lists arg     reference file lists (optional)
  -d, --dirs arg      reference file directories (optional)
  -k, --kmer_len arg  k-mer length
  -s, --kmc_db arg    set of k-mers (KMC database) prefix
  -t, --threads arg   number of threads to use (default: 1)
  -o, --output arg    output file
  -f, --format arg    output format (0: txt, 1: GFAv1, 2: GFAv2) (default: 0)
  -w, --work_dir arg  working directory (default: .)
      --mph arg       minimal perfect hash (BBHash) file (optional)
      --buckets arg   hash table buckets (cuttlefish) file (optional)
  -h, --help          print usage
```

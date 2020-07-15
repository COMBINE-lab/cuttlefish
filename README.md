# Cuttlefish

Building the compacted de Bruijn graph efficiently from set of reference.

## Compile
(out-of-source build)
```
$ mkdir build
$ cd build
$ cmake ..
$ make install
```

## Usage
(from project directory)
```
./bin/cuttlefish -r <reference-file> -k <k> -d <KMC-database-prefix> -t <thread-count> -o <output-file> -b <bbhash-file>
```

Options:
```
efficiently constructed the compacted de Bruijn graph.
Usage:
  cuttlefish [OPTION...]

  -r, --ref arg       the refrence sequence (FASTA)
  -k, --klen arg      the k-mer length
  -d, --database arg  the KMC database
  -t, --threads arg   the number of threads to use
  -o, --output arg    the output file
  -b, --bbhash arg    the BBHash file
  -h, --help          print usage
```

# kmc-super-kmers-iterator
This repo is for iterating super-k-mers from bins in kmc.

# Generating bins
This repo contains a precompiled modified version of KMC (in release page, click [here](https://github.com/refresh-bio-sandbox/kmc-super-kmers-iterator/releases/download/v0.0.1/kmc.tar.gz) to download directly) that prepares its intermediate files in a form that may be read using the API described here.
KMC version used here does not stop after generating these intermediate files, but it continues till the end of computations. 
Because this is based on a version that we used to build partial unitigs it stores counter k-mers of each bin in a separate binary file. This is not important at this moment.
Because the number of files may be quite huge I would recommend storing them in separate directories.
Let's say we have a fastq file (input.fq) that we want to use to generate bins. 
To generate bins one may use:
```
mkdir -p bins
mkdir -p out_dir
./kmc -k27 -t8 input.fq out_dir/out bins
```
The content of `out_dir` may be ignored.
The content of `bins` dir is important. There are `2 * n_bins + 2` files inside. `n_bins` is a number of bins that may be controlled with `-n` parameter.
All these files are in binary format and understanding them is not very important to use the API (the API just wants a path to this directory).
Below I describe them just for documentation purposes:
 - `bins.global` - some global configuration, contains length of k-mers, number of bins and flag if bins are zstd compressed
 - `stbm.bin` - a binary file containing a mapping of signatures values to bin_id, not used now, but may be needed in the future to determine the id of bin on a linking character end
 - `kmc_<bin_id>.bin` - this is the binary data containing k-mers - this is exactly the same file as the unmodified version of KMC produces
 - `kmc_<bin_id>.bin.meta` - this is the binary that contains additional metadata required to decode bin (in unmodified KMC this metadata is kept in RAM

# overlaps of super-k-mers like in ggcat
To generate super-k-mers with overlaps like in ggact one should add a flag `--super-kmers-with-overlaps` to kmc run.

# overlaps inside the same bin
If it happens that two consecutive super-k-mers are from the same bin, the overlap will be of (k-1) symbols. This way we avoid overrepresentation of k-mers, but as a consequence not all (k+1)-mers from the original input are represented in super-k-mers.
To have all (k+1)-mers one may use an additional flag `--inter-bin-k-overlap`. This will of course cause the overrepresentation of some k-mers.

# The API
The whole API is enclosed in a single class: `IterateSuperKmers` (file `iterate_super_kmers.h`) with the following methods:
 - `IterateSuperKmers(std::string bins_path, size_t bin_id, size_t queue_size)` - constructor
   -  `bins_path` - the path where kmc bins are stored (`bins` in the example command above)
   -  `bin_id` - which bin do we want to iterate
   -  `queue_size` - this is used to determine the size of one of the internal queues, I would recommend setting it to the number of threads one wants to use but probably may be lower - I have not experimented with this a lot
 - `void AddConsumer(SUPER_KMER_CALLBACK_T&& super_kmer_callback)` - the most important method of a class. **Each call of this method creates a new thread** that iterates super-k-mers from some data pack and for each super-k-mer calls a callback. The callback may be anything callable that takes two parameters: `const unsigned long long**` and `size_t` and returns `void` (I would recommend using lambda).
Here is a simple example of how to call it.
```
IterateSuperKmers iterate("bins", 0, 4);
iterate.AddConsumer([](const unsigned long long* super_kmer, size_t super_kmer_len_symbols, bool left_flag, bool right_flag) {
 //this part of the code will be called for super-k-mers
 //super_kmer is binary encoded and the format will be described below
 //super_kmer_len_symbols is length of the super-k-mer in symbols
 // the size of const unsigned long long* super_kmer is always iterate.GetSuperKmerDataLen()
 //left_flag is `true` if the left-most symbol of super-k-mer is a linking character (i.e. left-most k-mer exists also in a different bin)
 //right_flag is `true` if the right-most symbol of super-k-mer is a linking character (i.e. right-most k-mer exists also in a different bin)
 //if it happens that super-k-mer was broken into two parts but both are stored in the same bin these flags are not set although the character was appended to one of the two super-k-mers (it must be appended because in the opposite case we would lose a k-mer, and it must be appended only to a one of the two because in the opposite case, it would be overrepresented in a bin)
});
```
 - `void WaitForAll()` - this is a barrier, it must be called before the destruction of the object it joins threads
 - `size_t GetSuperKmerDataLen()` - the callback of `AddConsumer` takes two parameters, `GetSuperKmerDataLen()` returns the size of an array that is the first parameter

# Parallelization
Currently, the single internal pack of K-mers may be quite big so for small datasets it may be impossible to actually use more threads. For bigger datasets, it should be fine.
To lower the size of these packs one may run KMC with a smaller amount of memory (`-m` switch). This may increase kmc computation time, but for testing may be OK.
If needed we may redesign the internals of the API a little to exploit parallelization also for smaller datasets.

# Just generate super-k-mers in text format
This repo contains also a prebuild version of KMC (`kmc_text_dump`) that does all the above but also dumps super-k-mers in FASTA format for each bin. It results in additional files `kmers_<bin_id>.fa`.
It makes the computation much longer but may be used for some verifications (I have been using it to verify the API - keep in mind that the order of super-k-mers may be different in these files and in the API, which is a result of multithreading - if one thread is used in kmc and in API the order should be the same).
# ZSTD bins
It is possible to apply zstd compression to bins which may be very profitable for some datasets (like Salmonella from ggcat paper). To do this use `--bin-storage-modezstd` as a KMC flag.
The decompression may be done in parallel (and this is how its done in KMC), but the API in this repo currently does it single-threaded.
There are also some other improvements related to zstd reading that may be done here. If you notice any slowdowns because of this let us know.


# The example
This repository contains also an example application (`example.cpp`) that just converts the binary representation of a bin to fasta format. The conversion of super-k-mer binary representation to string is not optimized (and is in fact very slow), but this is to keep the code simple and demonstrate how symbols are stored in the raw memory.
Class `super_kmer_decoder` shows how to access symbols of super-k-mers. It may be compiled using `Makefile` which will result in a binary `bin/example`. This binary may be used as follows:
```
bin/example bins 0 4 > 0.fa # convert bin 0 (where all bins of kmc run are in bins directory) using 4 threads into fasta with super-k-mers
```

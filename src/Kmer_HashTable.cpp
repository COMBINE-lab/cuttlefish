
#include "Kmer_HashTable.hpp"

#include <fstream>
#include <vector>
#include <chrono>


void Kmer_HashTable::construct(const std::string& kmers_file, const uint64_t kmer_count)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Open the k-mers collection file.
    std::ifstream input(kmers_file.c_str(), std::ifstream::in);
    if(!input)
    {
        std::cerr << "Error opening k-mer collection file " << kmers_file << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    
    // Get the k-mers collection into memory.
    std::vector<uint64_t> kmers;
    kmers.reserve(kmer_count);
    
    uint64_t kmer;
    while(input >> kmer)
        kmers.push_back(kmer);
    

    // Build the MPHF.
    mph = new boomphf::mphf<u_int64_t, hasher_t>(kmer_count, kmers, 16, gamma_factor);

    std::cout << "MPH table size in bits / elem: " << (float)(mph->totalBitSize()) / kmer_count;


    // Remove the k-mers from memory forcefully.
    kmers.clear();
    kmers.shrink_to_fit();


    // Allocate the hash table.
    hash_table.resize(kmer_count);


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "\nDone constructing the MPH table. Time taken = " << elapsed_seconds << " seconds.\n";
}


void Kmer_HashTable::clear()
{
    // hash.clear();
    delete mph;
    
    hash_table.clear();
    hash_table.shrink_to_fit();
}

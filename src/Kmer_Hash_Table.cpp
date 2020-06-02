
#include "Kmer_Hash_Table.hpp"

#include <fstream>
#include <vector>
#include <chrono>


void Kmer_Hash_Table::construct(const std::string& kmers_file, const uint64_t kmer_count)
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
    // TODO: Replace this mechanism totally.
    std::vector<uint64_t> kmers;
    kmers.reserve(kmer_count);
    
    uint64_t kmer;
    while(input >> kmer)
        kmers.push_back(kmer);


    // Number of threads to use in the MPH build.
    // TODO: fix this.
    const uint32_t thread_count = 16;

    // Build the MPHF.
    mph = new boomphf::mphf<u_int64_t, hasher_t>(kmer_count, kmers, thread_count, gamma_factor);


    const uint64_t total_bits = mph->totalBitSize();
    std::cout << "Total MPH size (in MB): " << total_bits / (8 * 1024 * 1024) << "\n";
    std::cout << "MPH table size in bits / elem: " << (float)(total_bits) / kmer_count << "\n";


    // Remove the k-mers collection from memory forcefully.
    kmers.clear();
    kmers.shrink_to_fit();


    // Allocate the hash table.
    hash_table.resize(kmer_count);


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done constructing the MPH table. Time taken = " << elapsed_seconds << " seconds.\n";
}


void Kmer_Hash_Table::clear()
{
    // hash.clear();
    delete mph;
    
    hash_table.clear();
    hash_table.shrink_to_fit();
}

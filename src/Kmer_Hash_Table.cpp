
#include "Kmer_Hash_Table.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_Iterator.hpp"

#include <fstream>
#include <vector>
#include <chrono>


void Kmer_Hash_Table::construct(const std::string& kmc_file_name, const uint16_t thread_count)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Open a container over the k-mer database.
    Kmer_Container kmer_container(kmc_file_name);

    // Build the MPHF.
    uint64_t kmer_count = kmer_container.size();
    auto data_iterator = boomphf::range(kmer_container.begin(), kmer_container.end());
    mph = new boomphf::mphf<cuttlefish::kmer_t, Kmer_Hasher> (kmer_count, data_iterator, thread_count, gamma_factor);

    const uint64_t total_bits = mph->totalBitSize();
    std::cout << "Total MPH size (in MB): " << total_bits / (8 * 1024 * 1024) << "\n";
    std::cout << "MPH table size in bits / elem: " << (float)(total_bits) / kmer_count << "\n";

    // Allocate the hash table.
    hash_table.resize(kmer_count);


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done constructing the MPH table. Time taken = " << elapsed_seconds << " seconds.\n";
}


void Kmer_Hash_Table::clear()
{
    delete mph;
    
    hash_table.clear();
    hash_table.shrink_to_fit();
}

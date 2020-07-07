
#include "Kmer_Hash_Table.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_Iterator.hpp"

#include <fstream>
#include <vector>
#include <chrono>
#include <sys/stat.h>


void Kmer_Hash_Table::construct(const std::string& kmc_file_name, const std::string& bbhash_file_name, const uint16_t thread_count)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Open a container over the k-mer database.
    Kmer_Container kmer_container(kmc_file_name);

    uint64_t kmer_count = kmer_container.size();


    // The serialized BBHash file (saved from some earlier execution) exists.
    struct stat buffer;
    if(stat(bbhash_file_name.c_str(), &buffer) == 0)
    {
        std::cout << "Loading the MPH function from file " << bbhash_file_name << "\n";
        
        std::ifstream input(bbhash_file_name.c_str(), std::ifstream::in);
        mph = new boomphf::mphf<cuttlefish::kmer_t, Kmer_Hasher>();
        mph->load(input);
        input.close();
        
        std::cout << "Loaded the MPH function into memory.\n";
    }
    else    // No BBHash file exists. Build and save one now.
    {
        // Build the MPHF.
        std::cout << "Building the MPH function from the k-mer database " << kmc_file_name << "\n";

        auto data_iterator = boomphf::range(kmer_container.begin(), kmer_container.end());
        mph = new boomphf::mphf<cuttlefish::kmer_t, Kmer_Hasher> (kmer_count, data_iterator, thread_count, gamma_factor);

        std::cout << "Built the MPH function in memory.\n";
        

        // Save the MPHF.
        std::cout << "Saving the MPH function in file " << bbhash_file_name << "\n";

        std::ofstream output(bbhash_file_name.c_str(), std::ofstream::out);
        mph->save(output);
        output.close();

        std::cout << "Saved the MPH function in disk.\n";
    }

    const uint64_t total_bits = mph->totalBitSize();
    std::cout << "Total MPH size (in MB): " << total_bits / (8 * 1024 * 1024) << "\n";
    std::cout << "MPH table size in bits / elem: " << (float)(total_bits) / kmer_count << "\n";

    // Allocate the hash table.
    hash_table.resize(kmer_count);
    hash_table.clear_mem();
    std::cout << "Allocated hash table for " << kmer_count << " k-mers. Total memory taken (in MB): "
                << hash_table.bytes() / (1024 * 1024) << "\n";


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done constructing the MPH table. Time taken = " << elapsed_seconds << " seconds.\n";
}


void Kmer_Hash_Table::clear()
{
    delete mph;
    
    // hash_table.clear();
    hash_table.resize(0);
}

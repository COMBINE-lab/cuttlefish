
#include "Kmer_Hash_Table.hpp"
#include "Kmer_Iterator.hpp"

#include <fstream>
#include <vector>
#include <chrono>
#include <sys/stat.h>


void Kmer_Hash_Table::build_mph_function(const Kmer_Container& kmer_container, const std::string& bbhash_file_name, const uint16_t thread_count)
{
    // The serialized BBHash file (saved from some earlier execution) exists.
    struct stat buffer;
    if(!bbhash_file_name.empty() && stat(bbhash_file_name.c_str(), &buffer) == 0)
    {
        std::cout << "Loading the MPH function from file " << bbhash_file_name << "\n";

        load_mph_function(bbhash_file_name);

        std::cout << "Loaded the MPH function into memory.\n";
    }
    else    // No BBHash file name provided, or does not exist. Build and save (if specified) one now.
    {
        // Build the MPHF.
        std::cout << "Building the MPH function from the k-mer database " << kmer_container.container_location() << "\n";

        auto data_iterator = boomphf::range(kmer_container.begin(), kmer_container.end());
        mph = new boomphf::mphf<cuttlefish::kmer_t, Kmer_Hasher> (kmer_container.size(), data_iterator, thread_count, gamma_factor);

        std::cout << "Built the MPH function in memory.\n";
        

        // Save the MPHF if specified.
        if(!bbhash_file_name.empty())
        {
            std::cout << "Saving the MPH function in file " << bbhash_file_name << "\n";

            save_mph_function(bbhash_file_name);

            std::cout << "Saved the MPH function in disk.\n";
        }
    }
}


void Kmer_Hash_Table::load_mph_function(const std::string& file_name)
{
    std::ifstream input(file_name.c_str(), std::ifstream::in);
    if(input.fail())
    {
        std::cerr << "Error opening file " << file_name << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    mph = new cuttlefish::mphf_t();
    mph->load(input);

    input.close();
}


void Kmer_Hash_Table::save_mph_function(const std::string& file_name) const
{
    std::ofstream output(file_name.c_str(), std::ofstream::out);
    if(output.fail())
    {
        std::cerr << "Error writing to file " << file_name << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    mph->save(output);
    
    output.close();
}


void Kmer_Hash_Table::construct(const std::string& kmc_file_name, const std::string& bbhash_file_name, const uint16_t thread_count)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Open a container over the k-mer database.
    Kmer_Container kmer_container(kmc_file_name);
    uint64_t kmer_count = kmer_container.size();
    
    lock_range_size = uint64_t(std::ceil(double(kmer_count) / lock_count));


    // Build the minimal perfect hash function.
    build_mph_function(kmer_container, bbhash_file_name, thread_count);
    const uint64_t total_bits = mph->totalBitSize();
    std::cout << "Total MPH size (in MB): " << total_bits / (8 * 1024 * 1024) << "\n";
    std::cout << "MPH table size in bits / elem: " << (float)(total_bits) / kmer_count << "\n";

    // Allocate the hash table buckets.
    hash_table.resize(kmer_count);
    hash_table.clear_mem();
    std::cout << "Allocated hash table buckets for " << kmer_count << " k-mers. ";
    std::cout << "Total memory taken (in MB): " << hash_table.bytes() / (1024 * 1024) << "\n";

    uint64_t total_mem = (total_bits / 8) + hash_table.bytes();   // in bytes
    std::cout << "Total memory usage by the hash table: " << total_mem / (1024 * 1024)  << " MB.\n";
    std::cout << "Bits per k-mer: " << (total_mem * 8.0) / kmer_count << "\n";


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done allocating the hash table. Time taken = " << elapsed_seconds << " seconds.\n";
}


void Kmer_Hash_Table::clear()
{
    delete mph;
    
    // hash_table.clear();
    hash_table.resize(0);
}

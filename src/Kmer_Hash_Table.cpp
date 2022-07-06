
#include "Kmer_Hash_Table.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "Build_Params.hpp"
#include "utility.hpp"

#include <fstream>
#include <algorithm>
#include <chrono>
#include <cstdio>



template <uint16_t k, uint8_t BITS_PER_KEY> constexpr double Kmer_Hash_Table<k, BITS_PER_KEY>::gamma_min;
template <uint16_t k, uint8_t BITS_PER_KEY> constexpr double Kmer_Hash_Table<k, BITS_PER_KEY>::gamma_max;
template <uint16_t k, uint8_t BITS_PER_KEY> constexpr double Kmer_Hash_Table<k, BITS_PER_KEY>::min_bits_per_hash_key;
template <uint16_t k, uint8_t BITS_PER_KEY> constexpr double Kmer_Hash_Table<k, BITS_PER_KEY>::bits_per_gamma[];
template <uint16_t k, uint8_t BITS_PER_KEY> constexpr double Kmer_Hash_Table<k, BITS_PER_KEY>::gamma_resolution;


template <uint16_t k, uint8_t BITS_PER_KEY>
Kmer_Hash_Table<k, BITS_PER_KEY>::Kmer_Hash_Table(const std::string& kmc_db_path): Kmer_Hash_Table(kmc_db_path, Kmer_Container<k>::size(kmc_db_path))
{}


template <uint16_t k, uint8_t BITS_PER_KEY>
Kmer_Hash_Table<k, BITS_PER_KEY>::Kmer_Hash_Table(const std::string& kmc_db_path, const uint64_t kmer_count):
    gamma(gamma_min),
    kmc_db_path(kmc_db_path),
    kmer_count(kmer_count),
    sparse_lock(kmer_count, lock_count)
{}


template <uint16_t k, uint8_t BITS_PER_KEY>
Kmer_Hash_Table<k, BITS_PER_KEY>::Kmer_Hash_Table(const std::string& kmc_db_path, const uint64_t kmer_count, const std::size_t max_memory): Kmer_Hash_Table(kmc_db_path, kmer_count)
{
    set_gamma(max_memory);
}


template <uint16_t k, uint8_t BITS_PER_KEY>
Kmer_Hash_Table<k, BITS_PER_KEY>::Kmer_Hash_Table(const std::string& kmc_db_path, const uint64_t kmer_count, const std::size_t max_memory, const double gamma):
    Kmer_Hash_Table(kmc_db_path, kmer_count)
{
    if(gamma > 0)
        this->gamma = std::min(std::max(gamma, gamma_min), gamma_max);
    else
        set_gamma(max_memory);
}


template <uint16_t k, uint8_t BITS_PER_KEY>
void Kmer_Hash_Table<k, BITS_PER_KEY>::set_gamma(const std::size_t max_memory)
{
    const std::size_t max_memory_bits = max_memory * 8U;
    const std::size_t min_memory_bits = static_cast<std::size_t>(kmer_count * (min_bits_per_hash_key + BITS_PER_KEY));
    if(max_memory_bits > min_memory_bits)
    {
        const double max_bits_per_hash_key = (static_cast<double>(max_memory_bits) / kmer_count) - BITS_PER_KEY;
        const std::size_t gamma_idx = (std::upper_bound(bits_per_gamma, bits_per_gamma + (sizeof(bits_per_gamma) / sizeof(*bits_per_gamma)), max_bits_per_hash_key) - 1) - bits_per_gamma;
        gamma = gamma_idx * gamma_resolution;
    }
}


template <uint16_t k, uint8_t BITS_PER_KEY>
void Kmer_Hash_Table<k, BITS_PER_KEY>::build_mph_function(const uint16_t thread_count, const std::string& working_dir_path, const std::string& mph_file_path)
{
    // The serialized BBHash file (saved from some earlier execution) exists.
    if(!mph_file_path.empty() && file_exists(mph_file_path))
    {
        std::cout << "Found the MPHF at file " << mph_file_path << ".\n";
        std::cout << "Loading the MPHF.\n";

        load_mph_function(mph_file_path);

        std::cout << "Loaded the MPHF into memory.\n";
    }
    else    // No BBHash file name provided, or does not exist. Build one now.
    {
        // Open a container over the k-mer database.
        const Kmer_Container<k> kmer_container(kmc_db_path);

        // Build the MPHF.
        std::cout << "Building the MPHF from the k-mer database " << kmer_container.container_location() << ".\n";

        // auto data_iterator = boomphf::range(kmer_container.buf_begin(), kmer_container.buf_end());
        const auto data_iterator = boomphf::range(kmer_container.spmc_begin(thread_count), kmer_container.spmc_end(thread_count));
        std::cout << "Using gamma = " << gamma << ".\n";
        mph = new mphf_t(kmer_count, data_iterator, working_dir_path, thread_count, gamma);

        std::cout << "Built the MPHF in memory.\n";

        // std::cout << "Total data copy time to BBHash buffers " << mph->data_copy_time << "\n\n";
    }
}


template <uint16_t k, uint8_t BITS_PER_KEY>
void Kmer_Hash_Table<k, BITS_PER_KEY>::load_mph_function(const std::string& file_path)
{
    std::ifstream input(file_path.c_str(), std::ifstream::in);
    if(input.fail())
    {
        std::cerr << "Error opening file " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    mph = new mphf_t();
    mph->load(input);

    input.close();
}


template <uint16_t k, uint8_t BITS_PER_KEY>
void Kmer_Hash_Table<k, BITS_PER_KEY>::save_mph_function(const std::string& file_path) const
{
    std::ofstream output(file_path.c_str(), std::ofstream::out);
    if(output.fail())
    {
        std::cerr << "Error writing to file " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    mph->save(output);
    
    output.close();
}


template <uint16_t k, uint8_t BITS_PER_KEY>
void Kmer_Hash_Table<k, BITS_PER_KEY>::save_hash_buckets(const std::string& file_path) const
{
    std::ofstream output(file_path.c_str(), std::ofstream::out);
    if(output.fail())
    {
        std::cerr << "Error writing to file " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    hash_table.serialize(output);
    
    output.close();
}


template <uint16_t k, uint8_t BITS_PER_KEY>
void Kmer_Hash_Table<k, BITS_PER_KEY>::load_hash_buckets(const std::string& file_path)
{
    std::ifstream input(file_path.c_str(), std::ifstream::in);
    if(input.fail())
    {
        std::cerr << "Error opening file " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    input.close();


    hash_table.deserialize(file_path);
}


template <uint16_t k, uint8_t BITS_PER_KEY>
void Kmer_Hash_Table<k, BITS_PER_KEY>::save(const Build_Params& params) const
{
    save_mph_function(params.mph_file_path());
    save_hash_buckets(params.buckets_file_path());
}


template <uint16_t k, uint8_t BITS_PER_KEY>
void Kmer_Hash_Table<k, BITS_PER_KEY>::load(const Build_Params& params)
{
    load_mph_function(params.mph_file_path());
    load_hash_buckets(params.buckets_file_path());
}


template <uint16_t k, uint8_t BITS_PER_KEY>
void Kmer_Hash_Table<k, BITS_PER_KEY>::remove(const Build_Params& params) const
{
    const std::string mph_file_path = params.mph_file_path();
    const std::string buckets_file_path = params.buckets_file_path();

    if( (file_exists(mph_file_path) && std::remove(mph_file_path.c_str()) != 0) ||
        (file_exists(buckets_file_path) && std::remove(buckets_file_path.c_str()) != 0))
    {
        std::cerr << "Error removing the hash table files from disk. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k, uint8_t BITS_PER_KEY>
void Kmer_Hash_Table<k, BITS_PER_KEY>::construct(const uint16_t thread_count, const std::string& working_dir_path, const std::string& mph_file_path, const bool save_mph)
{
    // std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    std::cout << "Total number of k-mers in the set (KMC database): " << kmer_count << ".\n";


    // Build the minimal perfect hash function.
    build_mph_function(thread_count, working_dir_path, mph_file_path);

    if(save_mph)
    {
        save_mph_function(mph_file_path);
        std::cout << "Saved the hash function at " << mph_file_path << "\n";
    }

    const uint64_t total_bits = mph->totalBitSize();
    std::cout <<    "\nTotal MPHF size: " << total_bits / (8 * 1024 * 1024) << " MB."
                    " Bits per k-mer: " << static_cast<double>(total_bits) / kmer_count << ".\n";

    // Allocate the hash table buckets.
    hash_table.resize(kmer_count, State().code);
    std::cout << "Allocated hash table buckets for the k-mers. Total size: " <<
                hash_table.bytes() / (1024 * 1024) << " MB.\n";

    const uint64_t total_mem = (total_bits / 8) + hash_table.bytes();   // in bytes
    std::cout <<    "Total memory usage by the hash table: " << total_mem / (1024 * 1024)  << " MB."
                    " Bits per k-mer: " << (total_mem * 8.0) / kmer_count << ".\n";


    // std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    // const double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    // std::cout << "Done allocating the hash table. Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k, uint8_t BITS_PER_KEY>
void Kmer_Hash_Table<k, BITS_PER_KEY>::clear()
{
    if(mph != NULL)
        delete mph;

    mph = NULL;

    
    // hash_table.clear();
    hash_table.resize(0);
}


template <uint16_t k, uint8_t BITS_PER_KEY>
Kmer_Hash_Table<k, BITS_PER_KEY>::~Kmer_Hash_Table()
{
    clear();
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BIT, Kmer_Hash_Table)

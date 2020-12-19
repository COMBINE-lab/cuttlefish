
#include "Kmer_Hash_Table.hpp"
//#include "Kmer_Iterator.hpp"
#include "Kmer_Buffered_Iterator.hpp"

#include <fstream>
#include <vector>
#include <chrono>
#include <sys/stat.h>


template <uint16_t k>
void Kmer_Hash_Table<k>::build_mph_function(const Kmer_Container<k>& kmer_container, const uint16_t thread_count, const std::string& working_dir_path, const std::string& mph_file_path)
{
    // The serialized BBHash file (saved from some earlier execution) exists.
    struct stat buffer;
    if(!mph_file_path.empty() && stat(mph_file_path.c_str(), &buffer) == 0)
    {
        std::cout << "Found the MPHF at file " << mph_file_path << ".\n";
        std::cout << "Loading the MPHF.\n";

        load_mph_function(mph_file_path);

        std::cout << "Loaded the MPHF into memory.\n";
    }
    else    // No BBHash file name provided, or does not exist. Build and save (if specified) one now.
    {
        // Build the MPHF.
        std::cout << "Building the MPHF from the k-mer database " << kmer_container.container_location() << ".\n";

        auto data_iterator = boomphf::range(kmer_container.buf_begin(), kmer_container.buf_end());
        mph = new mphf_t(kmer_container.size(), data_iterator, working_dir_path, thread_count, GAMMA_FACTOR);

        std::cout << "Built the MPHF in memory.\n";

        // std::cout << "Total data copy time to BBHash buffers " << mph->data_copy_time << "\n\n";
        

        // Save the MPHF if specified.
        if(!mph_file_path.empty())
        {
            std::cout << "Saving the MPHF in file " << mph_file_path << ".\n";

            save_mph_function(mph_file_path);

            std::cout << "Saved the MPHF in disk.\n";
        }
    }
}


template <uint16_t k>
void Kmer_Hash_Table<k>::load_mph_function(const std::string& file_path)
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


template <uint16_t k>
void Kmer_Hash_Table<k>::save_mph_function(const std::string& file_path) const
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


template <uint16_t k>
void Kmer_Hash_Table<k>::save_hash_buckets(const std::string& file_path) const
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


template <uint16_t k>
void Kmer_Hash_Table<k>::load_hash_buckets(const std::string& file_path)
{
    std::ifstream input(file_path.c_str(), std::ifstream::in);
    if(input.fail())
    {
        std::cerr << "Error opening file " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    input.close();


    hash_table.deserialize(file_path, false);
}


template <uint16_t k>
void Kmer_Hash_Table<k>::construct(const std::string& kmc_db_path, const uint16_t thread_count, const std::string& working_dir_path, const std::string& mph_file_path)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Open a container over the k-mer database.
    Kmer_Container<k> kmer_container(kmc_db_path);
    const uint64_t kmer_count = kmer_container.size();
    std::cout << "Total number of k-mers in the set (KMC database): " << kmer_count << ".\n";
    
    lock_range_size = uint64_t(std::ceil(double(kmer_count) / lock_count));


    // Build the minimal perfect hash function.
    build_mph_function(kmer_container, thread_count, working_dir_path, mph_file_path);

    const uint64_t total_bits = mph->totalBitSize();
    std::cout <<    "\nTotal MPHF size: " << total_bits / (8 * 1024 * 1024) << " MB."
                    " Bits per k-mer: " << (float)(total_bits) / kmer_count << ".\n";

    // Allocate the hash table buckets.
    hash_table.resize(kmer_count);
    hash_table.clear_mem();
    std::cout << "Allocated hash table buckets for the k-mers. Total size: " <<
                hash_table.bytes() / (1024 * 1024) << " MB.\n";

    uint64_t total_mem = (total_bits / 8) + hash_table.bytes();   // in bytes
    std::cout <<    "Total memory usage by the hash table: " << total_mem / (1024 * 1024)  << " MB."
                    " Bits per k-mer: " << (total_mem * 8.0) / kmer_count << ".\n";


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done allocating the hash table. Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
void Kmer_Hash_Table<k>::clear()
{
    if(mph != NULL)
        delete mph;

    mph = NULL;
    
    // hash_table.clear();
    hash_table.resize(0);
}



// Template instantiations for the required specializations.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Kmer_Hash_Table)

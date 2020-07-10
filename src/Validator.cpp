
#include "Validator.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_Iterator.hpp"
#include "Kmer_Hasher.hpp"
#include "BBHash/BooPHF.h"
#include "Directed_Kmer.hpp"

#include <fstream>
#include <sys/stat.h>


Validator::Validator(const std::string& ref_file_name, const uint16_t k, const std::string& kmc_db_name, const std::string& cdbg_file_name):
    ref_file_name(ref_file_name), k(k), kmc_db_name(kmc_db_name), cdbg_file_name(cdbg_file_name), mph(NULL)
{
    Kmer::set_k(k);
}


void Validator::build_mph_function(const std::string& bbhash_file_name, const uint16_t thread_count)
{
    const Kmer_Container kmer_container(kmc_db_name);

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
        std::cout << "Building the MPH function from the k-mer database " << kmer_container.container_location() << "\n";

        auto data_iterator = boomphf::range(kmer_container.begin(), kmer_container.end());
        mph = new boomphf::mphf<cuttlefish::kmer_t, Kmer_Hasher> (kmer_container.size(), data_iterator, thread_count, 2);

        std::cout << "Built the MPH function in memory.\n";
        

        // Save the MPHF.
        std::cout << "Saving the MPH function in file " << bbhash_file_name << "\n";

        std::ofstream output(bbhash_file_name.c_str(), std::ofstream::out);
        mph->save(output);
        output.close();

        std::cout << "Saved the MPH function in disk.\n";
    }
}


bool Validator::validate_kmer_set() const
{
    Kmer_Container kmer_container(kmc_db_name);
    uint64_t kmer_count = kmer_container.size();

    std::cout << "Number of k-mers in the database: " << kmer_count << "\n";

    std::vector<bool> is_present(kmer_count);
    uint64_t kmers_seen = 0;
    uint64_t kmers_repeated = 0;
    uint64_t unitigs_processed = 0;
    uint64_t kmers_invalid = 0;

    std::string unitig;
    std::ifstream input(cdbg_file_name.c_str(), std::ifstream::in);
    while(input >> unitig)
    {
        cuttlefish::kmer_t first_kmer(unitig, 0);
        Directed_Kmer kmer(first_kmer);

        for(size_t kmer_idx = 0; kmer_idx <= unitig.length() - k; ++kmer_idx)
        {
            uint64_t hash_val = mph->lookup(kmer.canonical);

            if(hash_val >= kmer_count)
            {
                std::cout << "Invalid k-mer encountered.\n";
                kmers_invalid++;
            }
            else if(!is_present[hash_val])
                is_present[hash_val] = true;
            else
            {
                std::cout << "Repeated k-mer encountered.\n";
                kmers_repeated++;
            }

            if(kmer_idx < unitig.length() - k)
                kmer.roll_to_next_kmer(unitig[kmer_idx + k]);
        }


        kmers_seen += unitig.length() - k + 1;

        if(++unitigs_processed % 1000000 == 0)
            std::cout << "Processed " << unitigs_processed << " unitigs.\n";
    }


    std::cout << "Total number of repeated k-mers: " << kmers_repeated << "\n";
    std::cout << "Total number of invalid k-mers: " << kmers_invalid << "\n";
    std::cout << "Total number of k-mers seen: " << kmers_seen << "\n";
    std::cout << "Total number of k-mers expected: " << kmer_count << "\n";

    input.close();

    return !kmers_repeated && !kmers_invalid && kmers_seen == kmer_count;
}


bool Validator::validate(const std::string& bbhash_file_name, const uint16_t thread_count)
{
    build_mph_function(bbhash_file_name, thread_count);

    return validate_kmer_set();
}


#include "Validator.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_Iterator.hpp"
#include "BBHash/BooPHF.h"
#include "spdlog/sinks/stdout_color_sinks.h"

#include <fstream>
#include <sys/stat.h>


template <uint16_t k>
void Validator<k>::build_mph_function()
{
    const std::string& kmc_db_path = params.kmc_db_path();
    const uint16_t thread_count = params.thread_count();
    const std::string& mph_file_path = params.mph_file_path();

    const Kmer_Container<k> kmer_container(kmc_db_path);

    // The serialized BBHash file (saved from some earlier execution) exists.
    struct stat buffer;
    if(stat(mph_file_path.c_str(), &buffer) == 0)
    {
        console->info("Loading the MPH function from file {}\n", mph_file_path);
        
        std::ifstream input(mph_file_path.c_str(), std::ifstream::in);
        mph = new mphf_t();
        mph->load(input);
        input.close();
        
        console->info("Loaded the MPH function into memory.\n");
    }
    else    // No BBHash file exists. Build and save one now.
    {
        // Build the MPHF.
        console->info("Building the MPH function from the k-mer database {}\n", kmer_container.container_location());

        auto data_iterator = boomphf::range(kmer_container.begin(), kmer_container.end());
        mph = new mphf_t(kmer_container.size(), data_iterator, thread_count, GAMMA_FACTOR);

        console->info("Built the MPH function in memory.\n");
        

        // Save the MPHF.
        console->info("Saving the MPH function in file {}\n", mph_file_path);

        std::ofstream output(mph_file_path.c_str(), std::ofstream::out);
        mph->save(output);
        output.close();

        console->info("Saved the MPH function in disk.\n");
    }
}


template<uint16_t k>
void Validator<k>::clear()
{
    if(mph != NULL)
        delete mph;

    mph = NULL;
}



// Template instantiations for the required specializations.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Validator)

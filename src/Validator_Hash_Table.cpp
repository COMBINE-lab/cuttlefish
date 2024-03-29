
#include "Validator.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "utility.hpp"
#include "spdlog/sinks/stdout_color_sinks.h"

#include <fstream>


template <uint16_t k>
void Validator<k>::build_mph_function()
{
    const std::string& kmc_db_path = params.kmc_db_path();
    const uint16_t thread_count = params.thread_count();
    const std::string& working_dir_path = params.working_dir_path();
    const std::string& mph_file_path = params.mph_file_path();

    const Kmer_Container<k> kmer_container(kmc_db_path);

    // The serialized BBHash file (saved from some earlier execution) exists.
    if(file_exists(mph_file_path))
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

        // auto data_iterator = boomphf::range(kmer_container.begin(), kmer_container.end());
        auto data_iterator = boomphf::range(kmer_container.spmc_begin(thread_count), kmer_container.spmc_end(thread_count));
        mph = new mphf_t(kmer_container.size(), data_iterator, working_dir_path, thread_count, GAMMA_FACTOR);

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



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Validator)

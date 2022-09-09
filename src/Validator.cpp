
#include "Validator.hpp"
#include "DNA_Utility.hpp"

#include <thread>


template <uint16_t k>
Validator<k>::Validator(const Validation_Params& params, logger_t console):
    params(params), console(console)
{}


template <uint16_t k>
Validator<k>::Validator(const Validation_Params& params):
    params(params), console(spdlog::stdout_color_mt("Validator"))
{}


template <uint16_t k>
bool Validator<k>::validate()
{
    build_mph_function();


    bool is_valid_sequence;
    std::thread sequence_validator(&Validator::validate_sequence_completion, this, std::ref(is_valid_sequence));

    bool is_valid_kmer_set;
    std::thread kmer_set_validator(&Validator::validate_kmer_set, this, std::ref(is_valid_kmer_set));


    if(!sequence_validator.joinable())
    {
        console->error("Early termination faced for a worker thread. Aborting.\n");
        std::exit(EXIT_FAILURE);
    }

    sequence_validator.join();
    console->info("{} validation of complete coverage of the sequence by the produced unitigs.\n",
                    is_valid_sequence ? "Passed" : "Failed");

    
    if(!kmer_set_validator.joinable())
    {
        console->error("Early termination faced for a worker thread. Aborting.\n");
        std::exit(EXIT_FAILURE);
    }

    kmer_set_validator.join();
    console->info("{} validation of the k-mer set.\n", is_valid_kmer_set ? "Passed" : "Failed");


    clear();

    return is_valid_kmer_set && is_valid_sequence;
}


template <uint16_t k>
size_t Validator<k>::search_valid_kmer(const char* const seq, const size_t seq_len, const size_t start_idx) const
{
    size_t valid_start_idx;
    uint16_t base_count;

    size_t idx = start_idx;
    while(idx <= seq_len - k)
    {
        // Go over the contiguous subsequence of 'N's.
        for(; idx <= seq_len - k && DNA_Utility::is_placeholder(seq[idx]); idx++);

        // Go over the contiguous subsequence of non-'N's.
        if(idx <= seq_len - k)
        {
            valid_start_idx = idx;
            base_count = 0;

            for(; idx < seq_len && !DNA_Utility::is_placeholder(seq[idx]); ++idx)
                if(++base_count == k)
                    return valid_start_idx;
        }
    }


    return seq_len;
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Validator)

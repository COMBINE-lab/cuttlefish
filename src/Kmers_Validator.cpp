
#include "Validator.hpp"
#include "Directed_Kmer.hpp"
#include "Kmer_Container.hpp"
#include "spdlog/sinks/stdout_color_sinks.h"

#include <fstream>


template <uint16_t k>
void Validator<k>::validate_kmer_set(bool& result) const
{
    console->info("Testing validation of the uniqueness of the k-mers and completeness of the k-mer set in the produced unitigs.\n");

    const std::string& kmc_db_path = params.kmc_db_path();
    const std::string& cdbg_file_path = params.cdbg_file_path();

    const Kmer_Container<k> kmer_container(kmc_db_path);
    const uint64_t kmer_count = kmer_container.size();

    console->info("Number of k-mers in the database: {}\n", kmer_count);

    std::vector<bool> is_present(kmer_count);
    uint64_t kmers_seen = 0;
    uint64_t kmers_repeated = 0;
    uint64_t unitigs_processed = 0;
    uint64_t kmers_invalid = 0;

    // Scan through the unitigs one-by-one.
    std::string unitig;
    std::ifstream input(cdbg_file_path.c_str(), std::ifstream::in);
    while(input >> unitig)
    {
        const Kmer<k> first_kmer(unitig, 0);
        Directed_Kmer<k> kmer(first_kmer);

        // Scan through the k-mers one-by-one.
        for(size_t kmer_idx = 0; kmer_idx <= unitig.length() - k; ++kmer_idx)
        {
            uint64_t hash_val = mph->lookup(kmer.canonical());

            // Encountered a k-mer that is absent at the k-mer database and hashes outside of the valid range.
            if(hash_val >= kmer_count)
            {
                console->error("Invalid k-mer encountered.\n");
                kmers_invalid++;
            }
            // Encountered a k-mer for the first time that is either a k-mer present at the database, or is
            // absent there but hashes to the value of a present one (and that present one hasn't been seen yet).
            else if(!is_present[hash_val])
                is_present[hash_val] = true;
            // A repeated k-mer is seen.
            else
            {
                console->info("Repeated k-mer encountered.\n");
                kmers_repeated++;
            }

            if(kmer_idx < unitig.length() - k)
                kmer.roll_to_next_kmer(unitig[kmer_idx + k]);
        }


        kmers_seen += unitig.length() - k + 1;
        unitigs_processed++;

        if(unitigs_processed % PROGRESS_GRAIN_SIZE == 0)
            console->info("Validated {}M unitigs.\n", unitigs_processed / 1000000);
    }


    console->info("Total number of repeated k-mers: {}\n", kmers_repeated);
    console->info("Total number of invalid k-mers: {}\n", kmers_invalid);
    console->info("Total number of k-mers seen: {}\n", kmers_seen);
    console->info("Total number of k-mers expected: {}\n", kmer_count);

    input.close();

    result = (!kmers_repeated && !kmers_invalid && kmers_seen == kmer_count);
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Validator)

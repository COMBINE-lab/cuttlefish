
#ifndef KMER_ENUMERATOR_HPP
#define KMER_ENUMERATOR_HPP



#include "Build_Params.hpp"
#include "kmer_Enumeration_Stats.hpp"
#include "kmc_runner.h"


template <uint16_t k> class kmer_Enumeration_Stats;


// Class to enumerate all the k-mers for some provided input collection.
template <uint16_t k>
class kmer_Enumerator
{
private:

    static constexpr std::size_t min_memory = 3;    // In GB; set as per the KMC3 library requirement.
    static constexpr uint16_t bin_count = 2000;
    static constexpr uint16_t signature_len = 11;
    static constexpr uint64_t counter_max = 1;  // The `-cs` argument for KMC3; we're not interested in the counts and `cs = 1` will trigger skipping the counts.

    KMC::Stage1Params stage1_params;    // Parameters collection for the k-mer statistics approximation step of KMC3.
    KMC::Stage1Results stage1_results;  // Results of the k-mer statistics approximation.
    KMC::Stage2Params stage2_params;    // Parameters collection for the actual KMC3 execution (some execution parameters are absent and is present in `stage1_params`).
    KMC::Stage2Results stage2_results;  // Results of the actual k-mer set enumeration.

    KMC::Runner kmc;    // The KMC3 executor.


    // Returns the count of solid k-mers, i.e. k-mers occuring at least `cutoff` number of times,
    // estimated through KMC3's approximation step.
    uint64_t solid_kmer_count_approx(uint16_t cutoff) const;

    // Returns the strict memory limit (in GB) for the actual KMC3 execution, based on the number of
    // unique k-mers `unique_kmer_count` (typically approximated earlier) and the expected bits/kmer
    // `bits_per_kmer` requested.
    std::size_t memory_limit(uint64_t unique_kmer_count, double bits_per_kmer) const;


public:

    static constexpr uint16_t small_k_threshold = 13;   // KMC's internal threshold to switch modes for small-k optimizations.

    // Enumerates the k-mers from the sequences (of type `input_file_type`) present is `seqs`, that
    // are present at least `cutoff` times. Employs `thread_count` number of processor threads and
    // uses a soft memory-cap of `max_memory` (in GB). If `strict_memory` is `true`, then the memory
    // usage is attempted to be kept within a limit—the max of `max_memory` and the estimated memory
    // to be used by the downstream stages of Cuttlefish. This memory estimation is made only if
    // `estimate_mem_usage` is `true`, otherwise `max_memory` is the limit. Temporary files are
    // written to `working_dir_path`. The output database is stored at path prefix `output_db_path`.
    // Returns summary statistics of the enumeration.
    kmer_Enumeration_Stats<k> enumerate(
        KMC::InputFileType input_file_type, const std::vector<std::string>& seqs, uint32_t cutoff, uint16_t thread_count,
        std::size_t max_memory, bool strict_memory, bool estimate_mem_usage, double bits_per_kmer,
        const std::string& working_dir_path, const std::string& output_db_path);
};



// A class to display progress of the k-mer enumeration execution.
class FunnyProgress: public KMC::IPercentProgressObserver
{
    std::string funnChars = "/-\\|";
    int current = 0;

    void SetLabel(const std::string& label) override
    {
        //ignore
        (void)label;
    }
    void ProgressChanged(int newValue) override
    {
        if(newValue == 100)
            std::cerr << "\rDone.\n";
        else
            std::cerr << "\r" << funnChars[current++ % funnChars.size()];
        
    }
};



#endif

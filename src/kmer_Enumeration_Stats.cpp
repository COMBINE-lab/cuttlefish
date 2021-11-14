
#include <iostream>

#include "kmer_Enumeration_Stats.hpp"
#include "globals.hpp"


template <uint16_t k>
kmer_Enumeration_Stats<k>::kmer_Enumeration_Stats(const KMC::Stage1Results& stage1_results, const KMC::Stage2Results& stage2_results, const std::size_t max_memory, const std::size_t db_size):
    stage1_results(stage1_results),
    stage2_results(stage2_results),
    max_memory_(max_memory),
    db_size_(db_size)
{}


template <uint16_t k>
uint64_t kmer_Enumeration_Stats<k>::seq_count() const
{
    return stage1_results.nSeqences;
}


template <uint16_t k>
uint64_t kmer_Enumeration_Stats<k>::seq_len() const
{
    return total_kmer_count() + (seq_count() * (k - 1));
}


template <uint16_t k>
uint64_t kmer_Enumeration_Stats<k>::total_kmer_count() const
{
    return stage2_results.nTotalKmers;
}


template <uint16_t k>
uint64_t kmer_Enumeration_Stats<k>::unique_kmer_count() const
{
    return stage2_results.nUniqueKmers;
}


template <uint16_t k>
uint64_t kmer_Enumeration_Stats<k>::below_min_cutoff_kmer_count() const
{
    return stage2_results.nBelowCutoffMin;
}


template <uint16_t k>
uint64_t kmer_Enumeration_Stats<k>::above_max_cutoff_kmer_count() const
{
    return stage2_results.nAboveCutoffMax;
}


template <uint16_t k>
uint64_t kmer_Enumeration_Stats<k>::counted_kmer_count() const
{
    return unique_kmer_count() - (below_min_cutoff_kmer_count() + above_max_cutoff_kmer_count());
}


template <uint16_t k>
std::size_t kmer_Enumeration_Stats<k>::max_memory() const
{
    return max_memory_;
}


template <uint16_t k>
std::size_t kmer_Enumeration_Stats<k>::temp_disk_usage() const
{
    return stage2_results.maxDiskUsage;
}


template <uint16_t k>
std::size_t kmer_Enumeration_Stats<k>::db_size() const
{
    return db_size_;
}


template <uint16_t k>
void kmer_Enumeration_Stats<k>::log_stats() const
{
    std::cout << k << "-mer enumeration statistics:\n";

    std::cout << "\tNumber of sequences:\t" << seq_count() << ".\n";
    std::cout << "\tTotal sequence length:\t" << seq_len() << ".\n";
    std::cout << "\tTotal number of " << k << "-mers:\t" << total_kmer_count() << ".\n";
    std::cout << "\tNumber of unique " << k << "-mers:\t" << unique_kmer_count() << ".\n";
    std::cout << "\tNumber of solid " << k << "-mers:\t" << counted_kmer_count() << ".\n";
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_ALL, kmer_Enumeration_Stats)


#ifndef KMER_ENUMERATION_STATS_HPP
#define KMER_ENUMERATION_STATS_HPP



#include "kmc_runner.h"

#include <cstdint>
#include <cstddef>


// A class to wrap summary statistics of k-mer enumeration by `kmer_Enumerator`.
template <uint16_t k>
class kmer_Enumeration_Stats
{
private:

    const KMC::Stage1Results stage1_results;    // Results stats of KMC stage 1 execution.
    const KMC::Stage2Results stage2_results;    // Results stats of KMC stage 2 execution.
    const std::size_t max_memory_;  // Maximum memory usage allowed for the KMC executions.
    const std::size_t db_size_; // Size of the output KMC database size in bytes.


public:

    // Constructs a a k-mer enumeration stats wrapper object for a KMC execution with
    // first stage results in `stage1_results`, second stage results in `stage2_results`,
    // maximum allowed memory usage to be `max_memory` (in GB), and output database size
    // of `db_size`.
    kmer_Enumeration_Stats(const KMC::Stage1Results& stage1_results, const KMC::Stage2Results& stage2_results, std::size_t max_memory, std::size_t db_size);

    // Returns the number of sequences in the execution input.
    uint64_t seq_count() const;

    // Returns the total length of the sequences in the execution input.
    uint64_t seq_len() const;

    // Returns the total number of k-mers in the execution input.
    uint64_t total_kmer_count() const;

    // Returns the number of unique k-mers (irrespective of the cutoff frequency used) in the
    // execution input.
    uint64_t unique_kmer_count() const;
    
    // Returns the number of unique k-mers in the execution input that have frequency below
    // the minimum cutoff frequency used.
    uint64_t below_min_cutoff_kmer_count() const;

    // Returns the number of unique k-mers in the execution input that have frequency above
    // the maximum cutoff frequency used.
    uint64_t above_max_cutoff_kmer_count() const;

    // Returns the number of unique k-mers in the execution input that have frequencies within
    // the min and max cutoff frequencies used.
    uint64_t counted_kmer_count() const;

    // Returns the maximum memory (in GB) allowed for the execution.
    std::size_t max_memory() const;

    // Returns the temporary disk usage (in bytes) used by the execution.
    std::size_t temp_disk_usage() const;

    // Returns the size of the output KMC database size in bytes.
    std::size_t db_size() const;

    // Logs a summary statistics of the execution.
    void log_stats() const;
};



#endif

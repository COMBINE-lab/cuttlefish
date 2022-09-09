
#ifndef UNIPATHS_META_INFO_HPP
#define UNIPATHS_META_INFO_HPP



#include "Maximal_Unitig_Scratch.hpp"

#include <cstddef>
#include <cstdint>


// A class to track meta-information over maximal unipaths extracted by some worker thread.
template <uint16_t k>
class Unipaths_Meta_info
{
private:

    uint64_t unipath_count_;    // Total number of maximal unitigs.
    uint64_t kmer_count_;   // Total number of k-mers in the maximal unitigs.
    std::size_t max_len_;   // Length of the longest maximal unitig.
    std::size_t min_len_;   // Length of the shortest maximal unitig.
    uint64_t sum_len_;  // Sum length of the maximal unitigs.

    uint64_t dcc_count_;    // Total number of DCCs (Detached Chordless Cycles).
    uint64_t dcc_kmer_count_;   // Total number of k-mers in the DCCs.
    uint64_t dcc_sum_len_;  // Sum length of the DCCs.


public:

    // Constructs a meta-information tracker for maximal unitigs.
    Unipaths_Meta_info();

    // Adds information of the maximal unitig at the scratch space `maximal_unitig`
    // to the tracker.
    void add_maximal_unitig(const Maximal_Unitig_Scratch<k>& maximal_unitig);

    // Adds information of a maximal unitig with vertex count `size` to the tracker.
    void add_maximal_unitig(std::size_t size);

    // Aggregates the information of the tracker `other` to this tracker.
    void aggregate(const Unipaths_Meta_info<k>& other);

    // Returns the total number of maximal unitigs.
    uint64_t unipath_count() const;

    // Returns the total number of k-mers in the extracted maximal unitigs.
    uint64_t kmer_count() const;

    // Returns the length of the longest maximal unitig.
    std::size_t max_len() const;

    // Returns the length of the shortest maximal unitig.
    std::size_t min_len() const;

    // Returns the sum length of the maximal unitigs.
    uint64_t sum_len() const;

    // Returns the average length of the maximal unitigs.
    uint64_t avg_len() const;

    // Returns the total number of DCCs (Detached Chordless Cycles).
    uint64_t dcc_count() const;

    // Returns the total number of k-mers in the DCCs.
    uint64_t dcc_kmer_count() const;
    
    // Returns the sum length of the DCCs.
    uint64_t dcc_sum_len() const;

    // Prints the tracked information to the standard output.
    void print() const;
};


template <uint16_t k>
inline void Unipaths_Meta_info<k>::add_maximal_unitig(const std::size_t size)
{
    unipath_count_++;

    const std::size_t vertex_count = size;
    const std::size_t unipath_size = vertex_count + (k - 1);
    kmer_count_ += vertex_count;

    if(max_len_ < unipath_size)
        max_len_ = unipath_size;

    if(min_len_ > unipath_size)
        min_len_ = unipath_size;

    sum_len_ += unipath_size;
}


template <uint16_t k>
inline void Unipaths_Meta_info<k>::add_maximal_unitig(const Maximal_Unitig_Scratch<k>& maximal_unitig)
{
    add_maximal_unitig(maximal_unitig.size());

    if(maximal_unitig.is_cycle())
    {
        dcc_count_++;
        dcc_kmer_count_ += maximal_unitig.size();
        dcc_sum_len_ += maximal_unitig.size() + (k - 1);
    }
}



#endif

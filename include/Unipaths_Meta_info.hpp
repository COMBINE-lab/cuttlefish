
#ifndef UNIPATHS_META_INFO_HPP
#define UNIPATHS_META_INFO_HPP



#include <cstdint>
#include <limits>


// A class to track meta-information over maximal unipaths extracted by some worker thread.
template <uint16_t k>
class Unipaths_Meta_info
{
private:

    uint64_t unipath_count; // Total number of maximal unitigs.
    uint64_t kmer_count;    // Total number of k-mers in the maximal unitigs.
    std::size_t max_len;    // Length of the longest maximal unitig.
    std::size_t min_len;    // Length of the shortest maximal unitig.
    uint64_t sum_len;   // Sum length of the maximal unitigs.


public:

    // Constructs a meta-information tracker for maximal unitigs.
    Unipaths_Meta_info();

    // Adds information of the maximal unitig `unipath` to the tracker.
    template <typename T_container_>
    void add_maximal_unitig(const T_container_& unipath);

    // Aggregates the information of the tracker `other` to this tracker.
    void aggregate(const Unipaths_Meta_info<k>& other);

    // Prints the tracked information to the standard output.
    void print() const;
};


template <uint16_t k>
template <typename T_container_>
inline void Unipaths_Meta_info<k>::add_maximal_unitig(const T_container_& unipath)
{
    unipath_count++;
    
    kmer_count += unipath.size() - (k - 1);

    if(max_len < unipath.size())
        max_len = unipath.size();

    if(min_len > unipath.size())
        min_len = unipath.size();

    sum_len += unipath.size();
}



#endif

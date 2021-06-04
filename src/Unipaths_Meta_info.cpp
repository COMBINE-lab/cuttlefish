
#include "Unipaths_Meta_info.hpp"
#include "globals.hpp"

#include <algorithm>
#include <iostream>


template <uint16_t k>
Unipaths_Meta_info<k>::Unipaths_Meta_info():
    unipath_count(0),
    kmer_count(0),
    max_len(0),
    min_len(std::numeric_limits<std::size_t>::max()),
    sum_len(0)
{}


template <uint16_t k>
void Unipaths_Meta_info<k>::aggregate(const Unipaths_Meta_info& other)
{
    unipath_count += other.unipath_count;
    kmer_count += other.kmer_count;
    
    max_len = std::max(max_len, other.max_len);
    min_len = std::min(min_len, other.min_len);
    sum_len += other.sum_len;
}


template <uint16_t k>
void Unipaths_Meta_info<k>::print() const
{
    std::cout << "Number of maximal unitigs: " << unipath_count << ".\n";
    std::cout << "Number of k-mers in the maximal unitigs: " << kmer_count << ".\n";
    std::cout << "Length of the longest maximal unitig (in bases):  " << max_len << ".\n";
    std::cout << "Length of the shortest maximal unitig (in bases): " << min_len << ".\n";
    std::cout << "Sum length of the maximal unitigs (in bases): " << sum_len << ".\n";
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Unipaths_Meta_info)

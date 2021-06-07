
#include "Unipaths_Meta_info.hpp"
#include "globals.hpp"

#include <algorithm>
#include <iostream>


template <uint16_t k>
Unipaths_Meta_info<k>::Unipaths_Meta_info():
    unipath_count_(0),
    kmer_count_(0),
    max_len_(0),
    min_len_(std::numeric_limits<std::size_t>::max()),
    sum_len_(0)
{}


template <uint16_t k>
void Unipaths_Meta_info<k>::aggregate(const Unipaths_Meta_info& other)
{
    unipath_count_ += other.unipath_count_;
    kmer_count_ += other.kmer_count_;
    
    max_len_ = std::max(max_len_, other.max_len_);
    min_len_ = std::min(min_len_, other.min_len_);
    sum_len_ += other.sum_len_;
}


template <uint16_t k>
uint64_t Unipaths_Meta_info<k>::kmer_count() const
{
    return kmer_count_;
}


template <uint16_t k>
void Unipaths_Meta_info<k>::print() const
{
    std::cout << "Number of maximal unitigs: " << unipath_count_ << ".\n";
    std::cout << "Number of k-mers in the maximal unitigs: " << kmer_count_ << ".\n";
    std::cout << "Length of the longest maximal unitig (in bases):  " << max_len_ << ".\n";
    std::cout << "Length of the shortest maximal unitig (in bases): " << min_len_ << ".\n";
    std::cout << "Sum length of the maximal unitigs (in bases): " << sum_len_ << ".\n";
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Unipaths_Meta_info)

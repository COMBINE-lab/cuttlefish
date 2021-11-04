
#include "Unipaths_Meta_info.hpp"

#include "nlohmann/json.hpp"

#include <algorithm>
#include <iostream>
#include <cmath>


template <uint16_t k>
Unipaths_Meta_info<k>::Unipaths_Meta_info():
    unipath_count_(0),
    kmer_count_(0),
    max_len_(0),
    min_len_(std::numeric_limits<std::size_t>::max()),
    sum_len_(0),
    dcc_count_(0),
    dcc_kmer_count_(0),
    dcc_sum_len_(0)
{}


template <uint16_t k>
void Unipaths_Meta_info<k>::aggregate(const Unipaths_Meta_info& other)
{
    unipath_count_ += other.unipath_count_;
    kmer_count_ += other.kmer_count_;
    
    max_len_ = std::max(max_len_, other.max_len_);
    min_len_ = std::min(min_len_, other.min_len_);
    sum_len_ += other.sum_len_;

    dcc_count_ += other.dcc_count_;
    dcc_kmer_count_ += other.dcc_kmer_count_;
    dcc_sum_len_ += other.dcc_sum_len_;
}


template <uint16_t k>
uint64_t Unipaths_Meta_info<k>::unipath_count() const
{
    return unipath_count_;
}


template <uint16_t k>
uint64_t Unipaths_Meta_info<k>::kmer_count() const
{
    return kmer_count_;
}


template <uint16_t k>
std::size_t Unipaths_Meta_info<k>::max_len() const
{
    return max_len_;
}


template <uint16_t k>
std::size_t Unipaths_Meta_info<k>::min_len() const
{
    return min_len_;
}


template <uint16_t k>
uint64_t Unipaths_Meta_info<k>::sum_len() const
{
    return sum_len_;
}


template <uint16_t k>
uint64_t Unipaths_Meta_info<k>::avg_len() const
{
    return static_cast<uint64_t>(std::round(static_cast<double>(sum_len_) / unipath_count_));
}


template <uint16_t k>
uint64_t Unipaths_Meta_info<k>::dcc_count() const
{
    return dcc_count_;
}


template <uint16_t k>
uint64_t Unipaths_Meta_info<k>::dcc_kmer_count() const
{
    return dcc_kmer_count_;
}


template <uint16_t k>
uint64_t Unipaths_Meta_info<k>::dcc_sum_len() const
{
    return dcc_sum_len_;
}


template <uint16_t k>
void Unipaths_Meta_info<k>::print() const
{
    std::cout << "Number of maximal unitigs: " << unipath_count_ << ".\n";
    std::cout << "Number of k-mers in the maximal unitigs: " << kmer_count_ << ".\n";
    std::cout << "Length of the longest maximal unitig (in bases):  " << max_len_ << ".\n";
    std::cout << "Length of the shortest maximal unitig (in bases): " << min_len_ << ".\n";
    std::cout << "Sum length of the maximal unitigs (in bases): " << sum_len_ << ".\n";

    if(dcc_count_ > 0)
    {
        std::cout << "\nThere are Detached Chordless Cycles (DCC) present in the graph:\n";

        std::cout << "DCC count: " << dcc_count_ << ".\n";
        std::cout << "Number of vertices in the DCCs: " << dcc_kmer_count_ << ".\n";
        std::cout << "Sum length of the DCCs (in bases): " << dcc_sum_len_ << ".\n";
    }
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Unipaths_Meta_info)

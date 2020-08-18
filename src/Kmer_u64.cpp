
#include "Kmer_u64.hpp"

#include <algorithm>
#include <cassert>


// Static fields initialization.
uint16_t Kmer_u64::k = 0;
uint64_t Kmer_u64::bitmask_MSN = 0;


Kmer_u64::Kmer_u64(const std::string& label):
    Kmer_u64(label.c_str(), 0)
{}


Kmer_u64::Kmer_u64(const std::string& label, const size_t kmer_idx):
    Kmer_u64(label.c_str(), kmer_idx)
{}


void Kmer_u64::set_k(uint16_t k)
{
    Kmer_u64::k = k;
    Kmer_u64::bitmask_MSN = ~(uint64_t(0b11) << (2 * (k - 1)));
}


Kmer_u64 Kmer_u64::canonical() const
{
    return canonical(reverse_complement());
}


std::string Kmer_u64::string_label() const
{
    uint64_t kmer = this->kmer;
    char* label = new char[k + 1];
    
    for(uint16_t idx = 0; idx < k; ++idx)
    {
        switch(kmer & 0b11)
        {
        case DNA::A:
            label[idx] = 'A';
            break;

        case DNA::C:
            label[idx] = 'C';
            break;

        case DNA::G:
            label[idx] = 'G';
            break;

        case DNA::T:
            label[idx] = 'T';
            break;

        default:
            label[idx] = 'N';
        }

        
        kmer >>= 2;
    }


    label[k] = '\0';

    std::reverse(label, label + k);
    std::string str_label(label);
    
    delete label;

    return str_label;
}


uint16_t Kmer_u64::get_k()
{
    return k;
}


std::ostream& operator<<(std::ostream& out, const Kmer_u64& kmer)
{
    out << kmer.string_label();
    
    return out;
}


#include "Kmer_u64.hpp"

#include <cassert>


// Static fields initialization.
uint16_t Kmer_u64::k = 0;
uint64_t Kmer_u64::bitmask_MSN = 0;


Kmer_u64::Kmer_u64(const std::string& label)
{
    assert(label.length() == k);

    kmer = 0;

    for(const char p: label)
    {
        uint8_t nucleotide = map_nucleotide(p);
        // assert(nucleotide < DNA_Base::N);

        kmer = (kmer << 2) | nucleotide;
    }
}


Kmer_u64::Kmer_u64(const std::string& label, const size_t kmer_idx)
{
    kmer = 0;

    for(size_t idx = kmer_idx; idx < kmer_idx + k; ++idx)
    {
        uint8_t nucleotide = map_nucleotide(label[idx]);
        kmer = (kmer << 2) | nucleotide;
    }
}


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
    uint64_t kmer = this -> kmer;
    char* label = new char[k + 1];
    
    for(uint16_t idx = 0; idx < k; ++idx)
    {
        switch(kmer & 0b11)
        {
        case DNA_Base::A:
            label[idx] = 'A';
            break;

        case DNA_Base::C:
            label[idx] = 'C';
            break;

        case DNA_Base::G:
            label[idx] = 'G';
            break;

        case DNA_Base::T:
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

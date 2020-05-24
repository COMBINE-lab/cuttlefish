
#include "Kmer.hpp"

#include <cassert>
#include <cstdlib>
#include <algorithm>


uint16_t Kmer::k = 0;


inline uint8_t Kmer::map_nucleotide(const char nucleotide)
{
    switch(nucleotide)
    {
    case 'A':
        return DNA::A;
    
    case 'C':
        return DNA::C;

    case 'G':
        return DNA::G;

    case 'T':
        return DNA::T;

    default:
        return DNA::N;
    }
}


inline uint8_t Kmer::complement_nucleotide(const uint8_t nucleotide)
{
    switch(nucleotide)
    {
    case DNA::A:
        return DNA::T;

    case DNA::C:
        return DNA::G;

    case DNA::G:
        return DNA::C;

    case DNA::T:
        return DNA::A;

    default:
        return DNA::N;
    }
}


Kmer::Kmer(const std::string& label)
{
    assert(label.length() == k);

    kmer = 0;

    for(const char p: label)
    {
        uint8_t nucleotide = map_nucleotide(p);
        assert(nucleotide < DNA::N);

        kmer = (kmer << 2) | nucleotide;
    }
}


Kmer::Kmer(const char* label, const uint32_t kmer_idx)
{
    kmer = 0;

    for(uint32_t idx = kmer_idx; idx < kmer_idx + k; ++idx)
    {
        uint8_t nucleotide = map_nucleotide(label[idx]);

        // Placeholder rule to handle `N` nucleotides.
        // TODO: Need to make an informed rule for this.
        if(nucleotide == DNA::N)
            nucleotide = DNA::A;

        kmer = (kmer << 2) | nucleotide;
    }
}


Kmer::Kmer(const uint64_t kmer): kmer(kmer)
{}


void Kmer::set_k(uint16_t k)
{
    Kmer::k = k;
}


Kmer Kmer::reverse_complement() const
{
    uint64_t kmer = this -> kmer;
    uint64_t rev_comp = 0;

    for(uint16_t idx = 0; idx < k; ++idx)
    {
        rev_comp = ((rev_comp << 2) | complement_nucleotide(kmer & 0b11));

        kmer >>= 2;
    }

    return Kmer(rev_comp);
}


bool Kmer::operator <(const Kmer& rhs) const
{
    return kmer < rhs.kmer;
}


Kmer Kmer::canonical() const
{
    return std::min(*this, reverse_complement());
}


bool Kmer::operator ==(const Kmer& rhs) const
{
    return kmer == rhs.kmer;
}


cuttlefish::kmer_dir_t Kmer::direction(const Kmer& kmer_hat) const
{
    return *this == kmer_hat;
}


bool Kmer::is_same_kmer(const Kmer& rhs) const
{
    return canonical() == rhs.canonical();
}


std::string Kmer::string_label() const
{
    uint64_t kmer = this -> kmer;
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


uint64_t Kmer::int_label() const
{
    return kmer;
}


std::ostream& operator <<(std::ostream& out, const Kmer& kmer)
{
    out << kmer.string_label();
    
    return out;
}

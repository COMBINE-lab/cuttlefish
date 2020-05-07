
#include "Kmer.hpp"
#include "globals.hpp"

#include <algorithm>


Kmer Kmer::reverse_complement() const
{
    Kmer rev_comp(label);

    std::reverse(rev_comp.label.begin(), rev_comp.label.end());

    for(std::string::iterator p = rev_comp.label.begin(); p != rev_comp.label.end(); ++p)
        *p = complement(*p);

    
    return rev_comp;
}


bool Kmer::operator <(const Kmer& rhs) const
{
    return label < rhs.label;
}


bool Kmer::operator ==(const Kmer& rhs) const
{
    return label == rhs.label;
}

Kmer Kmer::min(const Kmer& lhs, const Kmer& rhs)
{
    return lhs < rhs ? lhs : rhs;
}


Kmer Kmer::canonical() const
{
    return Kmer::min(*this, this -> reverse_complement());
}


cuttlefish::kmer_dir_t Kmer::direction(const Kmer& kmer_hat) const
{
    return *this == kmer_hat;
}


bool Kmer::is_same_kmer(const Kmer& rhs) const
{
    return canonical() == rhs.canonical();
}


std::ostream& operator <<(std::ostream& out, const Kmer& kmer)
{
    out << kmer.label;

    return out;
}

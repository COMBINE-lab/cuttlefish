
#include "Kmer_Str.hpp"
#include "globals.hpp"

#include <algorithm>


Kmer_Str Kmer_Str::reverse_complement() const
{
    Kmer_Str rev_comp(label);

    std::reverse(rev_comp.label.begin(), rev_comp.label.end());

    for(std::string::iterator p = rev_comp.label.begin(); p != rev_comp.label.end(); ++p)
        *p = complement(*p);

    
    return rev_comp;
}


bool Kmer_Str::operator <(const Kmer_Str& rhs) const
{
    return label < rhs.label;
}


bool Kmer_Str::operator ==(const Kmer_Str& rhs) const
{
    return label == rhs.label;
}

Kmer_Str Kmer_Str::min(const Kmer_Str& lhs, const Kmer_Str& rhs)
{
    return lhs < rhs ? lhs : rhs;
}


Kmer_Str Kmer_Str::canonical() const
{
    return Kmer_Str::min(*this, this -> reverse_complement());
}


cuttlefish::kmer_dir_t Kmer_Str::direction(const Kmer_Str& kmer_hat) const
{
    return *this == kmer_hat;
}


bool Kmer_Str::is_same_kmer(const Kmer_Str& rhs) const
{
    return canonical() == rhs.canonical();
}


std::ostream& operator <<(std::ostream& out, const Kmer_Str& kmer)
{
    out << kmer.label;

    return out;
}

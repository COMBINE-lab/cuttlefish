
#ifndef KMER_STR_HPP
#define KMER_STR_HPP


#include <string>
#include <iostream>

#include "globals.hpp"


class Kmer_Str
{
private:
    std::string label;

    bool operator ==(const Kmer_Str& rhs) const;

    static Kmer_Str min(const Kmer_Str& lhs, const Kmer_Str& rhs);

public:
    Kmer_Str()
    {}

    Kmer_Str(const std::string& label):
        label(label)
    {}

    Kmer_Str reverse_complement() const;

    bool operator <(const Kmer_Str& rhs) const;

    Kmer_Str canonical() const;

    cuttlefish::kmer_dir_t direction(const Kmer_Str& kmer_hat) const;

    bool is_same_kmer(const Kmer_Str& rhs) const;

    friend std::ostream& operator <<(std::ostream& out, const Kmer_Str& kmer);
};


#endif

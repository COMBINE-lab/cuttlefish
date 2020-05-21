
#ifndef KMER_HPP
#define KMER_HPP


#include <string>
#include <iostream>

#include "globals.hpp"


class Kmer
{
private:
    std::string label;

    bool operator ==(const Kmer& rhs) const;

    static Kmer min(const Kmer& lhs, const Kmer& rhs);

public:
    Kmer()
    {}

    Kmer(const std::string& label):
        label(label)
    {}

    Kmer reverse_complement() const;

    bool operator <(const Kmer& rhs) const;

    Kmer canonical() const;

    cuttlefish::kmer_dir_t direction(const Kmer& kmer_hat) const;

    bool is_same_kmer(const Kmer& rhs) const;

    friend std::ostream& operator <<(std::ostream& out, const Kmer& kmer);
};


#endif

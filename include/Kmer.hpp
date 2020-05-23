
#ifndef KMER_HPP
#define KMER_HPP


#include "globals.hpp"

#include <cstdint>
#include <string>
#include <iostream>


class Kmer
{
private:
    static uint16_t k;

    // A = 0, C = 1, G = 2, T = 3
    enum DNA
    {
        A = 0b00, C = 0b01, G = 0b10, T = 0b11, N = 0b100
    };

    uint64_t kmer = 0;

    // Returns the mapping integer value of the given character `nucleotide`.
    static uint8_t map_nucleotide(const char nucleotide);

    // Returns the complement of `nucleotide`.
    static uint8_t complement_nucleotide(const uint8_t nucleotide);

    bool operator ==(const Kmer& rhs) const;

    static Kmer min(const Kmer& lhs, const Kmer& rhs);


public:
    Kmer()
    {}

    Kmer(const std::string& label);

    Kmer(const uint64_t kmer);

    static void set_k(const uint16_t k);

    Kmer reverse_complement() const;

    bool operator <(const Kmer& rhs) const;

    Kmer canonical() const;

    cuttlefish::kmer_dir_t direction(const Kmer& kmer_hat) const;

    bool is_same_kmer(const Kmer& rhs) const;

    std::string string_label() const;

    uint64_t int_label() const;

    friend std::ostream& operator <<(std::ostream& out, const Kmer& kmer);
};


#endif 

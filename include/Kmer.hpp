
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


    // Constructs a k-mer with the provided encoded value `kmer`.
    Kmer(const uint64_t kmer);

    // Returns the mapping integer value of the given character `nucleotide`.
    static uint8_t map_nucleotide(const char nucleotide);

    // Returns the complement of `nucleotide`.
    static uint8_t complement_nucleotide(const uint8_t nucleotide);

    bool operator ==(const Kmer& rhs) const;

    static Kmer min(const Kmer& lhs, const Kmer& rhs);


public:
    Kmer()
    {}

    // Constructs a k-mer from the provided string `label`.
    Kmer(const std::string& label);

    // Constructs a k-mer from the provided characters at
    // `label[kmer_idx,...,kmer_idx + k - 1]`.
    Kmer(const char* label, const uint32_t kmer_idx);

    // Set the value of the `k` parameter across the `Kmer` class.
    static void set_k(const uint16_t k);

    // Returns the reverese complement of the k-mer.
    Kmer reverse_complement() const;

    bool operator <(const Kmer& rhs) const;

    // Returns the canonical version of the k-mer.
    Kmer canonical() const;

    // Returns the direction of the k-mer relative to its canonical version.
    cuttlefish::kmer_dir_t direction(const Kmer& kmer_hat) const;

    // Returns true iff this k-mer and the provided k-mer `rhs` are actually
    // the same k-mer irrespective of the strands they originate from, i.e.
    // their canonical versions are the same.
    bool is_same_kmer(const Kmer& rhs) const;

    // Returns the string label of the k-mer.
    std::string string_label() const;

    // Returns the 64-bit encoding of the k-mer.
    uint64_t int_label() const;

    // For debugging purposes.
    friend std::ostream& operator <<(std::ostream& out, const Kmer& kmer);
};


#endif 

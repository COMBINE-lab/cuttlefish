
#ifndef KMER_HPP
#define KMER_HPP


#include "globals.hpp"

#include <cstdint>
#include <string>
#include <iostream>


class Kmer
{
private:
    static uint16_t k;  // The k-parameter.
    uint64_t kmer = 0;  // The 64-bit encoding of the underlying k-mer.
    static uint64_t bitmask_MSN;    // Bitmask used to clear the most significant nucleotide character, i.e. the first nucleotide of the k-mer which is at the bits `2k-1 : 2k-2`.

    // A = 0, C = 1, G = 2, T = 3
    enum DNA_Base
    {
        A = 0b00,   // 0b00
        C = 0b01,   // 0b01
        G = 0b10,   // 0b11
        T = 0b11,   // 0b11
        N = 0b100   // 0b100
    };


    // Constructs a k-mer with the provided encoded value `kmer`.
    Kmer(const uint64_t kmer);

    // Returns the mapping integer value of the given character `nucleotide`.
    static DNA_Base map_nucleotide(const char nucleotide);

    // Returns the mapping integer value of the complement of `nucleotide`.
    static DNA_Base complement_nucleotide(const DNA_Base nucleotide);

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

    bool operator ==(const Kmer& rhs) const;

    // Returns the canonical version of the k-mer.
    Kmer canonical() const;

    // Returns the direction of the k-mer relative to its canonical version.
    cuttlefish::kmer_dir_t direction(const Kmer& kmer_hat) const;

    // Returns true iff this k-mer and the provided k-mer `rhs` are actually
    // the same k-mer irrespective of the strands they originate from, i.e.
    // their canonical versions are the same.
    bool is_same_kmer(const Kmer& rhs) const;

    // Transforms this k-mer by chopping off the first nucleotide and
    // appending the next nucleotide `next_nucl` to the end, i.e.
    // rolls the k-mer by one nucleotide. Also sets the passed reverse
    // complement `rev_compl` of the k-mer accordingly.
    void roll_to_next_kmer(const cuttlefish::nucleotide_t next_nucl, cuttlefish::kmer_t& rev_compl);

    // Returns the canonical version of the k-mer, comparing it to its
    // reverse complement `rev_compl`.
    Kmer canonical(const Kmer& rev_compl) const;

    // Returns the string label of the k-mer.
    std::string string_label() const;

    // Returns the 64-bit encoding of the k-mer.
    uint64_t int_label() const;

    // For debugging purposes.
    friend std::ostream& operator <<(std::ostream& out, const Kmer& kmer);
};



inline Kmer::DNA_Base Kmer::map_nucleotide(const char nucleotide)
{
    switch(nucleotide)
    {
    case 'A':
        return DNA_Base::A;
    
    case 'C':
        return DNA_Base::C;

    case 'G':
        return DNA_Base::G;

    case 'T':
        return DNA_Base::T;

    default:
        // Placeholder rule to handle `N` nucleotides.
        // TODO: Need to make an informed rule for this.
        // Current: As per the rule used by the KMC tool.
        
        std::cerr << "Encountered invalid nucleotide " << nucleotide << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


inline Kmer::DNA_Base Kmer::complement_nucleotide(const DNA_Base nucleotide)
{
    switch(nucleotide)
    {
    case DNA_Base::A:
        return DNA_Base::T;

    case DNA_Base::C:
        return DNA_Base::G;

    case DNA_Base::G:
        return DNA_Base::C;

    case DNA_Base::T:
        return DNA_Base::A;

    default:
        // Placeholder rule to handle `N` nucleotides.
        // TODO: Need to make an informed rule for this.
        // Current: As per the rule used by the KMC tool.
        
        std::cerr << "Encountered invalid DNA_Base. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


inline uint64_t Kmer::int_label() const
{
    return kmer;
}


inline void Kmer::roll_to_next_kmer(const cuttlefish::nucleotide_t next_nucl, cuttlefish::kmer_t& rev_compl)
{
    const DNA_Base mapped_nucl = map_nucleotide(next_nucl);

    kmer = ((kmer & bitmask_MSN) << 2) | mapped_nucl;
    rev_compl.kmer = (rev_compl.kmer >> 2) | (uint64_t(complement_nucleotide(mapped_nucl)) << (2 * (k - 1)));
}


inline Kmer Kmer::canonical(const Kmer& rev_compl) const
{
    return std::min(*this, rev_compl);
}


#endif

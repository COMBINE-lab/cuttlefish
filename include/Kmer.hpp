
#ifndef KMER_HPP
#define KMER_HPP



#include "globals.hpp"
#include "kmc_api/kmc_file.h"

#include <string>
#include <iostream>


class Kmer_Hasher;


class Kmer
{
    friend class Kmer_Hasher;

private:

    static uint16_t k;  // The k-parameter.
    uint64_t kmer = 0;  // The 64-bit encoding of the underlying k-mer.
    static uint64_t bitmask_MSN;    // Bitmask used to clear the most significant nucleotide character, i.e. the first nucleotide of the k-mer which is at the bits `2k-1 : 2k-2`.

    // A = 0, C = 1, G = 2, T = 3.
    // Note that, this is not possible to change this mapping w/o modifications to the
    // interfacing of our code with the KMC api. This mapping is essential for some
    // crucial performance hacks in the interfacing.
    enum DNA_Base: uint8_t
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

    // Returns the 64-bit encoding of the k-mer.
    uint64_t to_u64() const;
    

public:

    Kmer(): kmer(0)
    {}

    // Constructs a k-mer from the provided string `label`.
    Kmer(const std::string& label);

    // Constructs a k-mer from the provided characters at
    // `label[kmer_idx,...,kmer_idx + k - 1]`.
    Kmer(const char* label, const size_t kmer_idx);

    // Constructs a k-mer from the provided characters at
    // `label[kmer_idx,...,kmer_idx + k - 1]`.
    Kmer(const std::string& label, const size_t kmer_idx);

    // Constructs a k-mer from `kmer_api` which is a k-mer object built from KMC.
    Kmer(const CKmerAPI& kmer_api);

    // Copy constructs the k-mer from another k-mer `rhs`.
    Kmer(const Kmer& rhs);

    // Copy assignment operator
    Kmer& operator=(const Kmer& rhs) = default;

    // Sets the value of the `k` parameter across the `Kmer` class.
    static void set_k(const uint16_t k);

    // Returns the DNA-complement character of the nucleotide character `nucl`.
    static cuttlefish::nucleotide_t complement(const cuttlefish::nucleotide_t nucl);

    // Returns the reverese complement of the k-mer.
    Kmer reverse_complement() const;

    // Returns true iff the bitwise encoding of this k-mer is lesser to the
    // encoding of the other k-mer `rhs`.
    bool operator<(const Kmer& rhs) const;

    // Returns true iff the bitwise encoding of this k-mer is equal to the
    // encoding of the other k-mer `rhs`.
    bool operator==(const Kmer& rhs) const;

    // Returns the direction of the k-mer relative to its canonical version.
    cuttlefish::dir_t direction(const Kmer& kmer_hat) const;

    // Transforms this k-mer by chopping off the first nucleotide and
    // appending the next nucleotide `next_nucl` to the end, i.e.
    // rolls the k-mer by one nucleotide. Also sets the passed reverse
    // complement `rev_compl` of the k-mer accordingly.
    void roll_to_next_kmer(const cuttlefish::nucleotide_t next_nucl, cuttlefish::kmer_t& rev_compl);

    // Returns the canonical version of the k-mer, comparing it to its
    // reverse complement `rev_compl`.
    Kmer canonical(const Kmer& rev_compl) const;

    // Returns the canonical version of the k-mer.
    Kmer canonical() const;

    // Returns the string label of the k-mer.
    std::string string_label() const;

    // Gets the k-mer from the KMC api object `kmer_api`.
    void from_CKmerAPI(const CKmerAPI& kmer_api);

    // Returns the k-parameter.
    static uint16_t get_k();
    
    // For debugging purposes.
    friend std::ostream& operator<<(std::ostream& out, const Kmer& kmer);
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
        // Placeholder rule to handle `N` nucleotides. Currently, as per the rule used by the KMC tool.
        
        std::cerr << "Encountered invalid nucleotide " << nucleotide << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


inline Kmer::Kmer(const char* label, const size_t kmer_idx)
{
    kmer = 0;

    for(size_t idx = kmer_idx; idx < kmer_idx + k; ++idx)
    {
        uint8_t nucleotide = map_nucleotide(label[idx]);
        kmer = (kmer << 2) | nucleotide;
    }
}


inline Kmer::Kmer(const CKmerAPI& kmer_api)
{
    kmer = 0;

    for(uint16_t idx = 0; idx < k; ++idx)
    {
        // uint8_t nucleotide = map_nucleotide(kmer_api.get_asci_symbol(idx));
        uint8_t nucleotide = kmer_api.get_num_symbol(idx);  // Works as long as our nucleotide-to-integer mapping is the same as KMC.

        kmer = (kmer << 2) | nucleotide;
    }
}


inline Kmer::Kmer(const Kmer& rhs): kmer(rhs.kmer)
{}
/*
inline Kmer& Kmer::operator=(const Kmer& rhs) {
    kmer = rhs.kmer;
    return *this; 
}
*/

inline Kmer Kmer::reverse_complement() const
{
    uint64_t kmer_val = kmer;
    Kmer rev_comp;

    for(uint16_t idx = 0; idx < k; ++idx)
    {
        rev_comp.kmer = ((rev_comp.kmer << 2) | complement_nucleotide(DNA_Base(kmer_val & 0b11)));

        kmer_val >>= 2;
    }

    return rev_comp;
}


inline void Kmer::from_CKmerAPI(const CKmerAPI& kmer_api)
{
    kmer = 0;
    kmer_api.to_u64(kmer);
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
        // Placeholder rule to handle `N` nucleotides. Currently, as per the rule used by the KMC tool.
        
        std::cerr << "Encountered invalid DNA_Base. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


inline cuttlefish::nucleotide_t Kmer::complement(const cuttlefish::nucleotide_t nucl)
{
    switch (nucl)
    {
    case 'A':
        return 'T';

    case 'C':
        return 'G';

    case 'G':
        return 'C';

    case 'T':
        return 'A';
    
    default:
        // Placeholder rule to handle `N` nucleotides.
        // TODO: Need to make an informed rule for this.
        
        std::cerr << "Invalid nucleotide " << nucl << " encountered. Aborting.";
        std::exit(EXIT_FAILURE);
    }
}


inline uint64_t Kmer::to_u64() const
{
    return kmer;
}


inline bool Kmer::operator<(const Kmer& rhs) const
{
    return kmer < rhs.kmer;
}


inline bool Kmer::operator==(const Kmer& rhs) const
{
    return kmer == rhs.kmer;
}


inline cuttlefish::dir_t Kmer::direction(const Kmer& kmer_hat) const
{
    return this->operator==(kmer_hat);
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

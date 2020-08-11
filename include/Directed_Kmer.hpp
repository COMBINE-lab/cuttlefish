
#ifndef DIRECTED_KMER_HPP
#define DIRECTED_KMER_HPP



#include "globals.hpp"


// K-mer and its reverse complement, canonical form, and direction.
template <uint16_t k>
class Directed_Kmer
{
protected:

    Kmer<k> kmer_;
    Kmer<k> rev_compl_;
    Kmer<k> canonical_;
    cuttlefish::dir_t dir_;


public:

    Directed_Kmer()
    {}

    // Constructs a k-mer with its reverse complement, canonical form, and direction.
    Directed_Kmer(const Kmer<k>& kmer);

    // Copy constructs the directed k-mer from `rhs`.
    Directed_Kmer(const Directed_Kmer<k>& rhs) = default;

    // Transforms this k-mer by chopping off the first nucleotide and
    // appending the next nucleotide `next_nucl` to the end, i.e.
    // rolls the k-mer by one nucleotide and sets all the relevant
    // information accordingly.
    void roll_to_next_kmer(cuttlefish::nucleotide_t next_nucl);

    void operator=(const Directed_Kmer<k>& rhs);

    // Returns the k-mer.
    Kmer<k> kmer() const;

    // Returns the reverse complement.
    Kmer<k> rev_compl() const;

    // Returns the canonical form of the k-mer.
    Kmer<k> canonical() const;

    // Returns the direction of the k-mer.
    cuttlefish::dir_t dir() const;
};


template <uint16_t k>
inline Directed_Kmer<k>::Directed_Kmer(const Kmer<k>& kmer):
    kmer_(kmer)
{
    rev_compl_ = kmer.reverse_complement();
    canonical_ = kmer.canonical(rev_compl_);
    dir_ = kmer.in_forward(canonical_);
}


template <uint16_t k>
inline void Directed_Kmer<k>::roll_to_next_kmer(const cuttlefish::nucleotide_t next_nucl)
{
    kmer_.roll_to_next_kmer(next_nucl, rev_compl_);
    
    canonical_ = kmer_.canonical(rev_compl_);
    dir_ = kmer_.in_forward(canonical_);
}


template <uint16_t k>
inline void Directed_Kmer<k>::operator=(const Directed_Kmer<k>& rhs)
{
    kmer_ = rhs.kmer_;
    rev_compl_ = rhs.rev_compl_;
    canonical_ = rhs.canonical_;
    dir_ = rhs.dir_;
}


template <uint16_t k>
inline Kmer<k> Directed_Kmer<k>::kmer() const
{
    return kmer_;
}


template <uint16_t k>
inline Kmer<k> Directed_Kmer<k>::rev_compl() const
{
    return rev_compl_;
}


template <uint16_t k>
inline Kmer<k> Directed_Kmer<k>::canonical() const
{
    return canonical_;
}


template <uint16_t k>
inline cuttlefish::dir_t Directed_Kmer<k>::dir() const
{
    return dir_;
}



#endif

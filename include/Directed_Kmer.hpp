
#ifndef DIRECTED_KMER_HPP
#define DIRECTED_KMER_HPP


#include "globals.hpp"
#include "Kmer.hpp"


class Annotated_Kmer;


// K-mer and its reverse complement, canonical form, and direction.
class Directed_Kmer
{
    friend class Annotated_Kmer;

private:

    cuttlefish::kmer_t kmer_;
    cuttlefish::kmer_t rev_compl_;
    cuttlefish::kmer_t canonical_;
    cuttlefish::kmer_dir_t dir_;


public:

    Directed_Kmer()
    {}

    // Constructs a k-mer with its reverse complement, canonical form, and direction.
    Directed_Kmer(const cuttlefish::kmer_t& kmer);

    // Copy constructs the directed k-mer from `rhs`.
    Directed_Kmer(const Directed_Kmer& rhs) = default;

    // Transforms this k-mer by chopping off the first nucleotide and
    // appending the next nucleotide `next_nucl` to the end, i.e.
    // rolls the k-mer by one nucleotide and sets all the relevant
    // information accordingly.
    void roll_to_next_kmer(const cuttlefish::nucleotide_t next_nucl);

    void operator=(const Directed_Kmer& rhs);

    // Returns the k-mer.
    cuttlefish::kmer_t kmer() const;

    // Returns the reverse complement.
    cuttlefish::kmer_t rev_compl() const;

    // Returns the canonical form of the k-mer.
    cuttlefish::kmer_t canonical() const;

    // Returns the direction of the k-mer.
    cuttlefish::kmer_dir_t dir() const;
};



inline Directed_Kmer::Directed_Kmer(const cuttlefish::kmer_t& kmer):
    kmer_(kmer)
{
    rev_compl_ = kmer.reverse_complement();
    canonical_ = kmer.canonical(rev_compl_);
    dir_ = kmer.direction(canonical_);
}


inline void Directed_Kmer::roll_to_next_kmer(const cuttlefish::nucleotide_t next_nucl)
{
    kmer_.roll_to_next_kmer(next_nucl, rev_compl_);
    
    canonical_ = kmer_.canonical(rev_compl_);
    dir_ = kmer_.direction(canonical_);
}


inline void Directed_Kmer::operator=(const Directed_Kmer& rhs)
{
    kmer_ = rhs.kmer_;
    rev_compl_ = rhs.rev_compl_;
    canonical_ = rhs.canonical_;
    dir_ = rhs.dir_;
}


inline cuttlefish::kmer_t Directed_Kmer::kmer() const
{
    return kmer_;
}


inline cuttlefish::kmer_t Directed_Kmer::rev_compl() const
{
    return rev_compl_;
}


inline cuttlefish::kmer_t Directed_Kmer::canonical() const
{
    return canonical_;
}


inline cuttlefish::kmer_dir_t Directed_Kmer::dir() const
{
    return dir_;
}


#endif

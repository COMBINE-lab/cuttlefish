
#ifndef DIRECTED_KMER_HPP
#define DIRECTED_KMER_HPP


#include "globals.hpp"
#include "Kmer.hpp"


class Directed_Kmer
{
public:
    cuttlefish::kmer_t kmer;
    cuttlefish::kmer_t rev_compl;
    cuttlefish::kmer_t canonical;
    cuttlefish::kmer_dir_t dir;


    Directed_Kmer()
    {}

    Directed_Kmer(const cuttlefish::kmer_t& kmer):
        kmer(kmer)
    {
        rev_compl = kmer.reverse_complement();
        canonical = kmer.canonical(rev_compl);
        dir = kmer.direction(canonical);
    }


    void roll_to_next_kmer(const cuttlefish::nucleotide_t next_nucl)
    {
        kmer.roll_to_next_kmer(next_nucl, rev_compl);
        canonical = kmer.canonical(rev_compl);
        dir = kmer.direction(canonical);
    }


    void operator =(const Directed_Kmer& rhs)
    {
        kmer = rhs.kmer;
        rev_compl = rhs.rev_compl;
        canonical = rhs.canonical;
        dir = rhs.dir;
    }
};


#endif

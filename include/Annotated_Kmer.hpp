
#ifndef ANNOTATED_KMER_HPP
#define ANNOTATED_KMER_HPP



#include "Directed_Kmer.hpp"
#include "Kmer_Hash_Table.hpp"


class CdBG;


// Complete k-mer information: k-mer itself and its reverse complement,
// canonical form, direction, index in corresponding sequence, and its class.
class Annotated_Kmer: public Directed_Kmer
{
    friend class CdBG;

private:

    size_t idx_;
    cuttlefish::Vertex_Class vertex_class_;


    Annotated_Kmer()
    {}

    // Constructs an annotated k-mer with its complete information.
    Annotated_Kmer(const cuttlefish::kmer_t& kmer, const size_t kmer_idx, const Kmer_Hash_Table& hash);

    // Copy constructs the annotated k-mer from `rhs`.
    Annotated_Kmer(const Annotated_Kmer& rhs) = default;

    // Transforms this k-mer by chopping off the first nucleotide and
    // appending the next nucleotide `next_nucl` to the end, i.e.
    // rolls the k-mer by one nucleotide, sets all the relevant k-mer
    // information accordingly (k-mer state is set using the `hash`).
    void roll_to_next_kmer(const cuttlefish::nucleotide_t next_nucl, const Kmer_Hash_Table& hash);

    void operator=(const Annotated_Kmer& rhs);

    // Returns the index of the k-mer.
    size_t idx() const;

    // Returns the vertex class of the k-mer.
    cuttlefish::Vertex_Class vertex_class() const;
};



inline Annotated_Kmer::Annotated_Kmer(const cuttlefish::kmer_t& kmer, const size_t kmer_idx, const Kmer_Hash_Table& hash):
    Directed_Kmer(kmer), idx_(kmer_idx), vertex_class_(hash[canonical_].vertex_class())
{}


inline void Annotated_Kmer::roll_to_next_kmer(const cuttlefish::nucleotide_t next_nucl, const Kmer_Hash_Table& hash)
{
    Directed_Kmer::roll_to_next_kmer(next_nucl);

    idx_++;
    vertex_class_ = hash[canonical_].vertex_class();
}


inline void Annotated_Kmer::operator=(const Annotated_Kmer& rhs)
{
    kmer_ = rhs.kmer_;
    rev_compl_ = rhs.rev_compl_;
    canonical_ = rhs.canonical_;
    dir_ = rhs.dir_;
    
    idx_ = rhs.idx_;
    vertex_class_ = rhs.vertex_class_;
}


inline size_t Annotated_Kmer::idx() const
{
    return idx_;
}


inline cuttlefish::Vertex_Class Annotated_Kmer::vertex_class() const
{
    return vertex_class_;
}



#endif
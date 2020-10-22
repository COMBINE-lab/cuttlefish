
#ifndef ANNOTATED_KMER_HPP
#define ANNOTATED_KMER_HPP



#include "Directed_Kmer.hpp"
#include "Kmer_Hash_Table.hpp"


// Complete k-mer information: k-mer itself and its reverse complement,
// canonical form, direction, index in corresponding sequence, and its class.
template <uint16_t k>
class Annotated_Kmer: public Directed_Kmer<k>
{
private:

    size_t idx_;
    cuttlefish::Vertex_Class vertex_class_;


public:

    Annotated_Kmer()
    {}

    // Constructs an annotated k-mer with its complete information.
    Annotated_Kmer(const Kmer<k>& kmer, size_t kmer_idx, const Kmer_Hash_Table<k>& hash);

    // Copy constructs the annotated k-mer from `rhs`.
    Annotated_Kmer(const Annotated_Kmer& rhs) = default;

    // Transforms this k-mer by chopping off the first base and
    // appending the next base `next_base` to the end, i.e. rolls
    // the k-mer by one base, sets all the relevant k-mer
    // information accordingly (k-mer state is set using the `hash`).
    void roll_to_next_kmer(char next_base, const Kmer_Hash_Table<k>& hash);

    void operator=(const Annotated_Kmer<k>& rhs);

    // Returns the index of the k-mer.
    size_t idx() const;

    // Returns the vertex class of the k-mer.
    cuttlefish::Vertex_Class vertex_class() const;
};


template <uint16_t k>
inline Annotated_Kmer<k>::Annotated_Kmer(const Kmer<k>& kmer, const size_t kmer_idx, const Kmer_Hash_Table<k>& hash):
    Directed_Kmer<k>(kmer), idx_(kmer_idx), vertex_class_(hash[this->canonical_].vertex_class())
{}


template <uint16_t k>
inline void Annotated_Kmer<k>::roll_to_next_kmer(const char next_base, const Kmer_Hash_Table<k>& hash)
{
    Directed_Kmer<k>::roll_to_next_kmer(next_base);

    idx_++;
    vertex_class_ = hash[this->canonical_].vertex_class();
}


template <uint16_t k>
inline void Annotated_Kmer<k>::operator=(const Annotated_Kmer<k>& rhs)
{
    this->kmer_ = rhs.kmer_;
    this->rev_compl_ = rhs.rev_compl_;
    this->canonical_ = rhs.canonical_;
    this->dir_ = rhs.dir_;
    
    idx_ = rhs.idx_;
    vertex_class_ = rhs.vertex_class_;
}


template <uint16_t k>
inline size_t Annotated_Kmer<k>::idx() const
{
    return idx_;
}


template <uint16_t k>
inline cuttlefish::Vertex_Class Annotated_Kmer<k>::vertex_class() const
{
    return vertex_class_;
}



#endif

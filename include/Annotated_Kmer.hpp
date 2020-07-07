
#ifndef ANNOTATED_KMER_HPP
#define ANNOTATED_KMER_HPP


#include "Directed_Kmer.hpp"
#include "Kmer_Hash_Table.hpp"


// Complete k-mer information: k-mer itself and its reverse complement,
// canonical form, direction, index in corresponding sequence, and type / state.
class Annotated_Kmer: public Directed_Kmer
{
public:
    uint32_t idx;
    cuttlefish::Vertex_Class vertex_class;



    Annotated_Kmer()
    {}


    // Constructs an annotated k-mer with its complete information.
    Annotated_Kmer(const cuttlefish::kmer_t& kmer, const uint32_t kmer_idx, const Kmer_Hash_Table& hash):
        Directed_Kmer(kmer), idx(kmer_idx)
    {
        vertex_class = hash[canonical].vertex_class();
    }


    // Transforms this k-mer by chopping off the first nucleotide and
    // appending the next nucleotide `next_nucl` to the end, i.e.
    // rolls the k-mer by one nucleotide, sets all the relevant k-mer
    // information accordingly (k-mer state is set using the `hash`).
    void roll_to_next_kmer(const cuttlefish::nucleotide_t next_nucl, const Kmer_Hash_Table& hash)
    {
        Directed_Kmer::roll_to_next_kmer(next_nucl);

        idx++;
        vertex_class = hash[canonical].vertex_class();
    }

    Annotated_Kmer(const Annotated_Kmer& rhs) : Directed_Kmer(rhs), idx(rhs.idx), vertex_class(rhs.vertex_class)
    {}

    void operator =(const Annotated_Kmer& rhs)
    {
        kmer = rhs.kmer;
        rev_compl = rhs.rev_compl;
        canonical = rhs.canonical;
        dir = rhs.dir;
        
        idx = rhs.idx;
        vertex_class = rhs.vertex_class;
    }
};


#endif
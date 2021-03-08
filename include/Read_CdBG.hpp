
#ifndef READ_CDBG_HPP
#define READ_CDBG_HPP



#include "Build_Params.hpp"
#include "Kmer_Hash_Table.hpp"


// Read de Bruijn graph class to support the compaction algorithm.
template <uint16_t k>
class Read_CdBG
{
private:

    const Build_Params params;  // Required parameters (wrapped inside).
    Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER> hash_table;  // Hash table for the vertices (canonical k-mers) of the graph.


public:

    // Constructs a `Read_CdBG` object with the parameters required for the construction of
    // the compacted representation of the underlying read de Bruijn graph wrapped in `params`.
    Read_CdBG(const Build_Params& params);

    // Constructs the compacted read de Bruijn graph, employing the parameters received
    // with the object-constructor.
    void construct();
};



#endif

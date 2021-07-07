
#ifndef READ_CDBG_HPP
#define READ_CDBG_HPP



#include "globals.hpp"
#include "Build_Params.hpp"
#include "Kmer_Hash_Table.hpp"

#include "nlohmann/json.hpp"


// Read de Bruijn graph class to support the compaction algorithm.
template <uint16_t k>
class Read_CdBG
{
private:

    const Build_Params params;  // Required parameters (wrapped inside).
    Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER> hash_table;  // Hash table for the vertices (canonical k-mers) of the graph.

    cuttlefish::json_t dBg_info;    // JSON object to store structural information over the de Bruijn graph.


    // Writes the structural information about the de Bruijn graph — obtained from the algorithm
    // execution — to disk.
    void dump_dBg_info() const;

public:

    // Constructs a `Read_CdBG` object with the parameters required for the construction of
    // the compacted representation of the underlying read de Bruijn graph wrapped in `params`.
    Read_CdBG(const Build_Params& params);

    // Constructs the compacted read de Bruijn graph, employing the parameters received
    // with the object-constructor.
    void construct();
};



#endif

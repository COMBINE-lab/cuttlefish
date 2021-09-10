
#ifndef READ_CDBG_HPP
#define READ_CDBG_HPP



#include "globals.hpp"
#include "Build_Params.hpp"
#include "Kmer_Hash_Table.hpp"
#include "dBG_Info.hpp"

#include <memory>


// Read de Bruijn graph class to support the compaction algorithm.
template <uint16_t k>
class Read_CdBG
{
private:

    const Build_Params params;  // Required parameters (wrapped inside).
    std::unique_ptr<Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>> hash_table; // Hash table for the vertices (canonical k-mers) of the graph.

    dBG_Info<k> dbg_info;   // Wrapper object for structural information of the graph.


    // Computes the states of the automata, i.e. the vertices in the graph.
    void compute_DFA_states();

    // Extracts the maximal unitigs from the graph.
    void extract_maximal_unitigs();


public:

    // Constructs a `Read_CdBG` object with the parameters required for the construction of
    // the compacted representation of the underlying read de Bruijn graph wrapped in `params`.
    Read_CdBG(const Build_Params& params);

    // Constructs the compacted read de Bruijn graph, employing the parameters received
    // with the object-constructor.
    void construct();

    // Returns `true` iff the compacted de Bruijn graph to be built from the parameters
    // collection `params` had been constructed in an earlier execution.
    // NB: only the existence of the output meta-info file is checked for this purpose.
    static bool is_constructed(const Build_Params& params);
};



#endif


#ifndef READ_CDBG_EXTRACTPR_HPP
#define READ_CDBG_EXTRACTOR_HPP



#include "globals.hpp"
#include "Kmer_Hash_Table.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "Build_Params.hpp"
#include "Spin_Lock.hpp"
#include "Thread_Pool.hpp"


// A class to extract the vertices from a compacted de Bruin graph — which are the maximal unitigs of some ordinary de Bruijn graph.
template <uint16_t k>
class Read_CdBG_Extractor
{
    friend class Thread_Pool<k>;

private:

    const Build_Params params;  // Required parameters (wrapped inside).
    Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table; // Hash table for the vertices (i.e. canonical k-mers) of the original (uncompacted) de Bruijn graph.

    // Members required to keep track of the total number of vertices processed across different worker (i.e. extractor) threads.
    mutable Spin_Lock lock;
    mutable uint64_t vertices_processed = 0;


    // Distributes the maximal unitigs extraction task — disperses the graph vertices (i.e. k-mers)
    // parsed by the parser `vertex_parser` to the worker threads in the thread pool `thread_pool`,
    // for the unitpath-flanking vertices to be identified and the corresponding unipaths to be extracted.
    void distribute_unipaths_extraction(Kmer_SPMC_Iterator<k>* vertex_parser, Thread_Pool<k>& thread_pool);

    // Processes the vertices provided to the thread with id `thread_id` from the parser `vertex_parser`,
    // i.e. for each vertex `v` provided to that thread, identifies whether it is a unipath-flanking
    // vertex, and if it is, then piece-wise constructs the corresponding unipath.
    void process_vertices(Kmer_SPMC_Iterator<k>* vertex_parser, uint16_t thread_id);

    // Returns `true` iff some vertex `v` with the provided state `state` is a flanking vertex for
    // the maximal unitig containing it. If it is, then stores the side of `v` that is opposite to
    // the flanking side, i.e. the side extending the unipath `p` containing `v`, to `unipath_side`.
    // If `p` is trivial, then the side stored is implementation-specific.
    // NB: unless for flanking vertices that are branching, this function cannot possibly be defined
    // — for non-branching vertices, it's not possible to compute whether they are flanking solely
    // from their states. Nevertheless, given that the information-discarding heuristic for branching
    // k-mers has been implemented in the DFA states computation phase, this method correctly computes
    // the flanking-status for some vertex from just its state.
    static bool is_flanking_state(State_Read_Space state, cuttlefish::side_t& unipath_side);

    // Returns `true` iff the vertex-side `side` for a vertex with state `state` flanks the maximal
    // unitig containing the vertex.
    // NB: this method is only applicable when information-discarding is propagated to the neighbors
    // from branching vertices. In absence of the implementation for the heuristic, this function
    // cannot possibly be defined from solely the parameters `state` and `side` — for non-branching
    // vertices, it's not possible to compute whether some side of them they are flanking solely
    // from their states — a vertex `v` is also required.
    static bool is_flanking_side(State_Read_Space state, cuttlefish::side_t side);


public:

    // Constructs a vertex-extractor object for some compacted read de Bruijn graph, with the required
    // parameters wrapped inside `params`, and uses the Cuttlefish hash table `hash_table`.
    Read_CdBG_Extractor(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table);

    // Extracts the maximal unitigs of the de Bruijn graph.
    void extract_maximal_unitigs();
};


template <uint16_t k>
inline bool Read_CdBG_Extractor<k>::is_flanking_state(const State_Read_Space state, cuttlefish::side_t& unipath_side)
{
    if(is_flanking_side(state, cuttlefish::side_t::front))
    {
        unipath_side = cuttlefish::side_t::back;
        return true;
    }

    if(is_flanking_side(state, cuttlefish::side_t::back))
    {
        unipath_side = cuttlefish::side_t::front;
        return true;
    }


    return false;
}


template <uint16_t k>
inline bool Read_CdBG_Extractor<k>::is_flanking_side(const State_Read_Space state, const cuttlefish::side_t side)
{
    const cuttlefish::edge_encoding_t edge = state.edge_at(side);

    return edge == cuttlefish::edge_encoding_t::N || edge == cuttlefish::edge_encoding_t::E;
}



#endif


#ifndef READ_CDBG_EXTRACTPR_HPP
#define READ_CDBG_EXTRACTOR_HPP



#include "globals.hpp"
#include "Kmer_Hash_Table.hpp"
#include "Directed_Vertex.hpp"
#include "Build_Params.hpp"
#include "Spin_Lock.hpp"
#include "Async_Logger_Wrapper.hpp"
#include "Output_Sink.hpp"
#include "Unipaths_Meta_info.hpp"

#include <limits>
#include <fstream>


// Forward declarations.
template <uint16_t k> class Kmer_SPMC_Iterator;
template <uint16_t k> class Thread_Pool;


// A class to extract the vertices from a compacted de Bruin graph — which are the maximal unitigs of some ordinary de Bruijn graph.
template <uint16_t k>
class Read_CdBG_Extractor
{
    friend class Thread_Pool<k>;

private:

    const Build_Params params;  // Required parameters (wrapped inside).
    Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table; // Hash table for the vertices (i.e. canonical k-mers) of the original (uncompacted) de Bruijn graph.

    // typedef std::ofstream sink_t;
    typedef Async_Logger_Wrapper sink_t;
    Output_Sink<sink_t> output_sink;    // Sink for the output maximal unitigs.

    // TODO: give these limits more thoughts, especially their exact impact on the memory usage.
    static constexpr std::size_t BUFF_SZ = 100 * 1024ULL;   // 100 KB (soft limit) worth of maximal unitigs can be retained in memory, at most, before flushing.
    static constexpr std::size_t SEQ_SZ = 5 * 1024ULL * 1024ULL;    // 5 MB (soft limit) sized maximal unitig, at most, is constructed at a time.

    mutable uint64_t vertices_processed = 0;    // Total number of vertices scanned from the database.
    mutable Spin_Lock lock; // Mutual exclusion lock to access various unique resources by threads spawned off this class' methods.
    
    Unipaths_Meta_info<k> unipaths_meta_info;   // Meta-information over the extracted maximal unitigs.


    // Distributes the maximal unitigs extraction task — disperses the graph vertices (i.e. k-mers)
    // parsed by the parser `vertex_parser` to the worker threads in the thread pool `thread_pool`,
    // for the unitpath-flanking vertices to be identified and the corresponding unipaths to be extracted.
    void distribute_unipaths_extraction(Kmer_SPMC_Iterator<k>* vertex_parser, Thread_Pool<k>& thread_pool);

    // Processes the vertices provided to the thread with id `thread_id` from the parser `vertex_parser`,
    // i.e. for each vertex `v` provided to that thread, identifies whether it is a unipath-flanking
    // vertex, and if it is, then piece-wise constructs the corresponding unipath.
    void process_vertices(Kmer_SPMC_Iterator<k>* vertex_parser, uint16_t thread_id);

    // Extracts the maximal unitig `p` that is flanked by the vertex `v_hat` and connects to `v_hat`
    // through its side `s_v_hat`. Returns `true` iff the extraction is successful, which happens when
    // the maximal unitig is encountered and attempted for output-marking _first_, by some thread. If
    // the attempt is successful, then the maximal unitig is extracted in its canonical form into
    // `unipath` (it is overwritten); also, a unique ID for it is put in `id`. If not, `unipath` may
    // contain partial form of the unitig, and `id` is unaltered.
    bool extract_maximal_unitig(const Kmer<k>& v_hat, cuttlefish::side_t s_v_hat, uint64_t& id, std::vector<char>& unipath);

    // Marks the vertex `v` as outputted. Returns `true` iff `v` has not been marked yet and the hash
    // table update is successful.
    bool mark_vertex(const Directed_Vertex<k>& v);

    // Marks the two endpoint vertices of a maximal unitig `p` as outputted: the first vertex in the
    // canonical form of `p`, `sign_vertex`, and the last vertex in the form, `cosign_vertex`. Returns
    // `true` iff the vertices have not been marked yet and the corresponding hash table updates are
    // successful.
    bool mark_flanking_vertices(const Directed_Vertex<k>& sign_vertex, const Directed_Vertex<k>& cosign_vertex);

    // Initializes the output sink.
    void init_output_sink();

    // Closes the output sink.
    void close_output_sink();

    // Replaces the character sequence `seq` in-place with its reverse complement.
    template <typename T_container_>
    static void reverse_complement(T_container_& seq);

    // Note: The following methods are only applicable when the heuristic of information-discarding
    // from branching vertices to their neighbors has been implemented in the DFA states computation
    // phase. In the general case, these functions with their specified input parameters and their
    // intended output values can not possibly be defined. Given just the state of a vertex, unless
    // for flanking vertices that are branching, it's not possible to determine —
    //      i.   whether it is a maximal unitig flanking vertex;
    //      ii.  whether some specific side of it is a flanking side;
    //      iii. and which side of it may connect to the containing maximal unitig.
    // The vertex itself, along with the hash table containing the DFA states, are required in general
    // — so that the neighboring vertices can also be probed to answer these queries. Nevertheless,
    // given that the heuristic of propagation of information-discarding from branching vertices has
    // been implemented in the DFA states computation phase, the following method definitions can
    // correctly respond to the queries given just the state. This is because the heuristic transforms
    // the flanking non-branching vertices into branching ones. Thus, although their states are not
    // technically "correct" as per the theoretical model — we are throwing away more information from
    // the model than it already is doing — this does not affect the output for the purposes of maximal
    // unitigs extraction.

    // Returns `true` iff some vertex `v` with the provided state `state` is a flanking vertex for the
    // maximal unitig `p` containing it. If yes, then stores the side of `v` to `unipath_side` through
    // which `v` is connected to `p`. If `p` is trivial, i.e. `p = v`, then the returned side is `back`,
    // —  necessitated by an optimization for the extraction of unipaths in their canonical forms.
    // NB: this method is only applicable if the heuristic of information-propagation from branching vertices
    // has been implemented in the DFA states computation phase. See the detailed comment in the class body.
    static bool is_flanking_state(State_Read_Space state, cuttlefish::side_t& unipath_side);

    // Returns `true` iff the vertex-side `side` for a vertex with state `state` flanks the maximal unitig
    // containing the vertex.
    // NB: this method is only applicable if the heuristic of information-propagation from branching vertices
    // has been implemented in the DFA states computation phase. See the detailed comment in the class body.
    static bool is_flanking_side(State_Read_Space state, cuttlefish::side_t side);


public:

    // Constructs a vertex-extractor object for some compacted read de Bruijn graph, with the required
    // parameters wrapped inside `params`, and uses the Cuttlefish hash table `hash_table`.
    Read_CdBG_Extractor(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table);

    // Extracts the maximal unitigs of the de Bruijn graph.
    void extract_maximal_unitigs();
};


template <uint16_t k>
inline bool Read_CdBG_Extractor<k>::mark_vertex(const Directed_Vertex<k>& v)
{
    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket = hash_table[v.hash()];
    State_Read_Space& state = bucket.get_state();

    if(state.is_outputted())
        return false;

    state.mark_outputted();
    return hash_table.update(bucket);
}


template <uint16_t k>
inline bool Read_CdBG_Extractor<k>::mark_flanking_vertices(const Directed_Vertex<k>& sign_vertex, const Directed_Vertex<k>& cosign_vertex)
{
    return mark_vertex(sign_vertex) && (sign_vertex.hash() == cosign_vertex.hash() || mark_vertex(cosign_vertex));
}


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


template <uint16_t k>
template <typename T_container_>
inline void Read_CdBG_Extractor<k>::reverse_complement(T_container_& seq)
{
    assert(!seq.empty());

    auto fwd = seq.begin();
    auto bwd = seq.end() - 1;

    for(; fwd < bwd; ++fwd, --bwd)
    {
        std::swap(*fwd, *bwd);
        
        *fwd = DNA_Utility::complement(*fwd),
        *bwd = DNA_Utility::complement(*bwd);
    }

    if(fwd == bwd)
        *fwd = DNA_Utility::complement(*fwd);
}



#endif

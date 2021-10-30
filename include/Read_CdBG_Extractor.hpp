
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
#include "Progress_Tracker.hpp"

#include <limits>
#include <fstream>


// Forward declarations.
template <uint16_t k> class Kmer_SPMC_Iterator;
template <uint16_t k> class Thread_Pool;
template <uint16_t k> class dBG_Info;


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
    static constexpr std::size_t SEQ_SZ = 1 * 1024ULL * 1024ULL;    // 1 MB (soft limit) sized maximal unitig, at most, is constructed at a time.

    mutable uint64_t vertices_scanned = 0;    // Total number of vertices scanned from the database.
    mutable Spin_Lock lock; // Mutual exclusion lock to access various unique resources by threads spawned off this class' methods.

    mutable uint64_t vertices_marked = 0;   // Total number of vertices marked as present in maximal unitigs; used for the extraction of detached chordless cycle(s), if any.
    
    Unipaths_Meta_info<k> unipaths_meta_info_;  // Meta-information over the extracted maximal unitigs.

    Progress_Tracker progress_tracker;  // Progress tracker for the maximal unitigs extraction task.


    // Distributes the maximal unitigs extraction task — disperses the graph vertices (i.e. k-mers)
    // parsed by the parser `vertex_parser` to the worker threads in the thread pool `thread_pool`,
    // for the unitpath-flanking vertices to be identified and the corresponding unipaths to be extracted.
    void distribute_unipaths_extraction(Kmer_SPMC_Iterator<k>* vertex_parser, Thread_Pool<k>& thread_pool);

    // Scans the vertices provided to the thread with id `thread_id` from the parser `vertex_parser`
    // for potential unipath-flanking vertices, i.e. for each vertex `v` provided to that thread,
    // identifies whether it is a unipath-flanking vertex, and if it is, then piece-wise constructs
    // the corresponding unipath.
    void scan_vertices(Kmer_SPMC_Iterator<k>* vertex_parser, uint16_t thread_id);

    // Extracts the maximal unitig `p` that is flanked by the vertex `v_hat` and connects to `v_hat`
    // through its side `s_v_hat`. Returns `true` iff the extraction is successful, which happens when
    // the maximal unitig is encountered and attempted for output-marking _first_, by some thread. If
    // the attempt is successful, then the maximal unitig is extracted in its canonical form into
    // `unipath` (it is overwritten), a unique ID for it is put in `id`, and the hashes of the vertices
    // constituting the path overwrites `path_hashes` (when the user-option is specified). If not,
    // `unipath` and `path_hashes` may contain partial form of the path, and `id` is unaltered.
    bool extract_maximal_unitig(const Kmer<k>& v_hat, cuttlefish::side_t s_v_hat, uint64_t& id, std::vector<char>& unipath, std::vector<uint64_t>& path_hashes);

    // Marks all the vertices which have their hashes present in `path_hashes` as outputted.
    void mark_path(const std::vector<uint64_t>& path_hashes);

    // Marks all the vertices that are present in the maximal unitigs of the graph with its vertex
    // set being present at the path prefix `vertex_db_path`.
    void mark_maximal_unitig_vertices(const std::string& vertex_db_path);

    // Scans the vertices provided to the thread with id `thread_id` from the parser `vertex_parser`
    // for potential unipath-flanking vertices. If a vertex `v` is found to be a flanking one, then
    // piece-wise constructs the corresponding (partial) maximal unitig starting the traversal from
    // `v`, and marks the vertices along the way. Premature halts before traversing the entire unitig
    // `p` is possible, in cases when some other thread is concurrently constructing `p`, but from
    // the opposite flank — the halt happens at the threads' meeting-point.
    void mark_maximal_unitig_vertices(Kmer_SPMC_Iterator<k>* vertex_parser, uint16_t thread_id);

    // Marks (partially) the vertices of the maximal unitig `p` that is flanked by the vertex `v_hat`
    // from one side and connects to `v_hat` through its side `s_v_hat`. `p` might not be marked
    // completely by this thread if some other thread is concurrently traversing `p`, but from the
    // opposite flank. However, together these two threads mark `p` completely. Also returns the
    // number of vertices marked in this execution.
    std::size_t mark_maximal_unitig(const Kmer<k>& v_hat, cuttlefish::side_t s_v_hat);

    // Extracts all the detached chordless cycles present in the graph with its vertex set being
    // present at the path prefix `vertex_db_path`, into the output file at `output_file_path`.
    void extract_detached_chordless_cycles(const std::string& vertex_db_path, const std::string& output_file_path);

    // Scans the vertices provided to the thread with id `thread_id` from the parser `vertex_parser`
    // for potential detached chordless cycles. If a vertex `v` is found to be not marked as present
    // in the earlier extracted maximal unitigs, then it implies — by definition from the Cuttlefish
    // algorithm — that `v` belongs to a chordless cycle that is detached completely from the rest of
    // the graph. The method piece-wise constructs the cycle, starting the traversal from `v`.
    void extract_detached_chordless_cycles(Kmer_SPMC_Iterator<k>* vertex_parser, uint16_t thread_id);

    // Extracts the detached chordless cycle `p` that contains the vertex `v_hat`. Returns `true` iff
    // the extraction is successful, which happens when the cycle is encountered and attempted for
    // output-marking _first_ by this thread. If the attempt is successful, then the cycle is extracted
    // in its literal form that starts with `v_hat`'s canonical representation, into `cycle` (it is
    // overwritten); also, a unique ID for it is put in `id`. If not, `cycle` may contain partial form
    // of the unitig, and `id` is unaltered. The index of the lexicographically lowest (canonical) k-mer
    // in `cycle` is recorded at `pivot`.
    bool extract_cycle(const Kmer<k>& v_hat, uint64_t& id, std::vector<char>& cycle, std::size_t& pivot);

    // Marks the vertex `v` as outputted. Returns `true` iff `v` has not been marked yet and the hash
    // table update is successful.
    bool mark_vertex(const Directed_Vertex<k>& v);

    // Marks the two endpoint vertices of a maximal unitig `p` as outputted: the first vertex in the
    // canonical form of `p`, `sign_vertex`, and the last vertex in the form, `cosign_vertex`. Returns
    // `true` iff the vertices have not been marked yet and the corresponding hash table updates are
    // successful.
    bool mark_flanking_vertices(const Directed_Vertex<k>& sign_vertex, const Directed_Vertex<k>& cosign_vertex);

    // Initializes the output sink, corresponding to the file `output_file_path`.
    void init_output_sink(const std::string& output_file_path);

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

    // Extracts the maximal unitigs of the de Bruijn graph with the vertex set at path prefix `vertex_db_path`,
    // into the output file at `output_file_path`.
    void extract_maximal_unitigs(const std::string& vertex_db_path, const std::string& output_file_path);

    // Extracts the chordless cycles from the de Bruijn graph that are completely disconnected from the
    // rest of the graph. The graph is to contain its vertex set at the path prefix `vertex_db_path`,
    // and the cycles are appeneded to the output file at `output_file_path`. `dbg_info` is used to
    // determine whether the compacted graph had been constructed earlier—in which case some data
    // structures are re-used from the earlier construction.
    void extract_detached_cycles(const std::string& vertex_db_path, const std::string& output_file_path, const dBG_Info<k>& dbg_info);

    // Returns the parameters collection for the compacted graph construction.
    const Build_Params& get_params() const;

    // Returns a wrapper over the meta-information of the extracted unitigs.
    const Unipaths_Meta_info<k>& unipaths_meta_info() const;

    // Returns the number of vertices in the underlying graph.
    uint64_t vertex_count() const;

    // Returns `true` iff the de Bruijn graph has DCCs (Detached Chordless Cycles).
    bool has_dcc() const;

    // Returns the number of vertices present in maximal unitigs (excluding the DCCs).
    uint64_t unipaths_vertex_count() const;
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


template <uint16_t k>
inline void Read_CdBG_Extractor<k>::mark_path(const std::vector<uint64_t>& path_hashes)
{
    for(const uint64_t hash: path_hashes)
        hash_table.update(hash, State_Read_Space::get_outputted_state());
}



#endif


#ifndef READ_CDBG_CONSTRUCTOR_HPP
#define READ_CDBG_CONSTRUCTOR_HPP



#include "globals.hpp"
#include "Kmer_Hash_Table.hpp"
#include "State_Read_Space.hpp"
#include "Edge.hpp"
#include "Endpoint.hpp"
#include "Build_Params.hpp"
#include "Progress_Tracker.hpp"

#include <cstdint>
#include <string>


template <uint16_t k> class Kmer_SPMC_Iterator;
template <uint16_t k> class Thread_Pool;


// A class to construct compacted read de Bruijn graphs.
template <uint16_t k>
class Read_CdBG_Constructor
{
    friend class Thread_Pool<k>;

private:

    const Build_Params params;  // Required parameters (wrapped inside).
    Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table; // Hash table for the vertices (canonical k-mers) of the graph.

    uint64_t edge_count_;    // Number of edges in the underlying graph.
    
    // Members required to keep track of the total number of edges processed across different threads.
    mutable Spin_Lock lock;
    mutable uint64_t edges_processed = 0;
    
    Progress_Tracker progress_tracker;  // Progress tracker for the DFA states computation task.


    // Distributes the DFA-states computation task — disperses the graph edges (i.e. (k + 1)-mers)
    // parsed by the parser `edge_parser` to the worker threads in the thread pool `thread_pool`,
    // for the edges to be processed by making appropriate state transitions for their endpoints.
    void distribute_states_computation(Kmer_SPMC_Iterator<k + 1>* edge_parser, Thread_Pool<k>& thread_pool);

    // Processes the edges provided to the thread with id `thread_id` from the parser `edge_parser`,
    // based on the end-purpose of extracting either the maximal unitigs or a maximal path cover.
    void process_edges(Kmer_SPMC_Iterator<k + 1>* edge_parser, uint16_t thread_id);

    // Processes the edges provided to the thread with id `thread_id` from the parser `edge_parser`,
    // i.e. makes state-transitions for the DFA of the vertices `u` and `v` for each bidirected edge
    // `(u, v)` provided to that thread, in order to construct a CdBG.
    void process_cdbg_edges(Kmer_SPMC_Iterator<k + 1>* edge_parser, uint16_t thread_id);

    // Processes the edges provided to the thread with id `thread_id` from the parser `edge_parser`,
    // i.e. makes state-transitions for the DFA of the vertices `u` and `v` for each bidirected edge
    // `(u, v)` provided to that thread, to construct a maximal path cover of the dBG.
    void process_path_cover_edges(Kmer_SPMC_Iterator<k + 1>* edge_parser, uint16_t thread_id);

    // Adds the information of an incident edge `e` to the side `s` of some vertex `v`, all wrapped
    // inside the edge-endpoint object `endpoint` — making the appropriate state transitions for the
    // DFA of `v`. Returns `false` iff an attempted state transition failed.
    bool add_incident_edge(const Endpoint<k>& endpoint);
    
    // Adds the information of an incident loop that connects the two different endpoints of some
    // vertex `v`, wrapped inside the edge-endpoint object `endpoint` — making the appropriate state
    // transition for the DFA of `v`. Returns `false` iff an attempted state transition failed.
    bool add_crossing_loop(const Endpoint<k>& endpoint);

    // Adds the information of an incident loop for some vertex `v` that connects its side `s` to
    // the side itself, all wrapped inside the edge-endpoint object `endpoint` — making the
    // appropriate state transition for the DFA of `v`. Returns `false` iff an attempted state
    // transition failed.
    bool add_one_sided_loop(const Endpoint<k>& endpoint);

    // Adds the information of the edge `e = {u, v}` to its endpoint vertices `u` and `v` iff this
    // edge connects sides of `u` and `v` that do not have any edges added yet, which ensures that
    // neither of the vertices belong to two different paths in a path cover of the graph; and makes
    // the appropriate state transitions for the DFAs of `u` and `v`. Returns `false` iff the edge
    // could not be added as such.
    bool add_path_cover_edge(const Edge<k>& e);


public:

    // Consructs a read-CdBG builder object, with the required parameters wrapped in `params`, and uses
    // the Cuttlefish hash table `hash_table`.
    Read_CdBG_Constructor(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table);

    // Computes the states of the DFA in the de Bruijn graph with the edge set at path prefix `edge_db_path`.
    void compute_DFA_states(const std::string& edge_db_path);

    // Returns the number of distinct vertices in the underlying graph.
    uint64_t vertex_count() const;

    // Returns the number of distinct edges in the underlying graph.
    uint64_t edge_count() const;


// The following methods are not used anymore with the current algorithm.
/*
private:

    // Adds the information of an incident edge `e` to the side `s` of some vertex `v`, all wrapped
    // inside the edge-endpoint object `endpoint` — making the appropriate state transitions for the
    // DFA of `v`. Also stores the edge encodings of the incidence information of the side `s` before
    // and after to this addition, in `e_old` and `e_new` respectively. Returns `false` iff an
    // attempted state transition failed.
    bool add_incident_edge(const Endpoint<k>& endpoint, cuttlefish::edge_encoding_t& e_old, cuttlefish::edge_encoding_t& e_new);

    // Adds the information of an incident loop that connects the two different endpoints of some
    // vertex `v`, wrapped inside the edge-endpoint object `endpoint` — making the appropriate state
    // transition for the DFA of `v`. Also stores the edge encodings of the incidence information of
    // the front and the back sides before this addition, in `e_front` and `e_back` respectively.
    // Returns `false` iff an attempted state transition failed.
    bool add_crossing_loop(const Endpoint<k>& endpoint, cuttlefish::edge_encoding_t& e_front, cuttlefish::edge_encoding_t& e_back);

    // Adds the information of an incident loop for some vertex `v` that connects its side `s` to
    // the side itself, all wrapped inside the edge-endpoint object `endpoint` — making the
    // appropriate state transition for the DFA of `v`. Also stores the edge encoding of the incidence
    // information of the side `s` before this addition, in `e_old`. Returns `false` iff an attempted
    // state transition failed.
    bool add_one_sided_loop(const Endpoint<k>& endpoint, cuttlefish::edge_encoding_t& e_old);

    // If the endpoint object `v_end` connects to some neighboring endpoint `w_end` through a unique
    // edge encoded with `e`, then discards the incidence information of `w_end` — making the
    // appropriate state transition for the corresponding neighboring vertex `w`.
    void propagate_discard(const Endpoint<k>& v_end, cuttlefish::edge_encoding_t e);

    // For two neighboring endpoints `u_end` and `v_end`, discards the incidence information from
    // `v_end`. Also discards information from any other neighboring side of `u_end` that may have
    // connected to it through a unique edge encoded with `e`. Makes the appropriate state transitions
    // for these neighbors of `u_end`.
    void propagate_discard(const Endpoint<k>& u_end, const Endpoint<k>& v_end, cuttlefish::edge_encoding_t e);

    // Discards the incidence information of the endpoint `v_end`. Returns `false` iff an attempted
    // state transition failed.
    bool discard_side(const Endpoint<k>& v_end);

    // Discards the incidence information of some endpoint `w_end` that connects to the endpoint
    // `v_end` through the unique edge encoded with `e` — making the appropriate state transition.
    void discard_neighbor_side(const Endpoint<k>& v, cuttlefish::edge_encoding_t e);
*/
};


template <uint16_t k>
inline bool Read_CdBG_Constructor<k>::add_incident_edge(const Endpoint<k>& endpoint)
{
    // Fetch the hash table entry for the vertex associated to the endpoint.

    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket = hash_table[endpoint.hash()];
    State_Read_Space& state = bucket.get_state();
    const cuttlefish::edge_encoding_t e_curr = state.edge_at(endpoint.side());

    // If we've already discarded the incidence information for this side, then a self-transition happens.
    if(e_curr == cuttlefish::edge_encoding_t::N)
        return true;    // The side has already been determined to be branching—nothing to update here anymore.

    cuttlefish::edge_encoding_t e_new = endpoint.edge();
    if(e_curr != cuttlefish::edge_encoding_t::E)    // The side is not empty.
    {
        // We can get away without updating the same value again, because — (1) even if this DFA's state changes
        // in the hash table by the time this method completes, making no updates at this point is theoretically
        // equivalent to returning instantaneously as soon as the hash table value had been read; and also (2) the
        // ordering of the edges processed does not matter in the algorithm.
        if(e_new == e_curr)
            return true;

        // This side has been visited earlier, but with a different edge—discard the incidence information.
        e_new = cuttlefish::edge_encoding_t::N;
    }

    state.update_edge_at(endpoint.side(), e_new);
    return hash_table.update(bucket);
}


template <uint16_t k>
inline bool Read_CdBG_Constructor<k>::add_crossing_loop(const Endpoint<k>& endpoint)
{
    // Fetch the hash table entry for the DFA of vertex associated to the endpoint.
    
    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket = hash_table[endpoint.hash()];
    State_Read_Space& state = bucket.get_state();

    const State_Read_Space state_curr = state;

    if(state.edge_at(cuttlefish::side_t::front) != cuttlefish::edge_encoding_t::N)   // Discard the front-incidence information, if not done already.
        state.update_edge_at(cuttlefish::side_t::front, cuttlefish::edge_encoding_t::N);
    
    if(state.edge_at(cuttlefish::side_t::back) != cuttlefish::edge_encoding_t::N)    // Discard the back-incidence information, if not done already.
        state.update_edge_at(cuttlefish::side_t::back, cuttlefish::edge_encoding_t::N);

    // We can get away without updating the same value again: see detailed comment in `add_incident_edge`.
    return state == state_curr ? true : hash_table.update(bucket);
}


template <uint16_t k>
inline bool Read_CdBG_Constructor<k>::add_one_sided_loop(const Endpoint<k>& endpoint)
{
    // Fetch the hash table entry for the vertex associated to the endpoint.

    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket = hash_table[endpoint.hash()];
    State_Read_Space& state = bucket.get_state();

    // We can get away without updating the same value again: see detailed comment in `add_incident_edge`.
    if(state.edge_at(endpoint.side()) == cuttlefish::edge_encoding_t::N) // The incidence information has already been discarded.
        return true;

    // Discard the incidence information.
    state.update_edge_at(endpoint.side(), cuttlefish::edge_encoding_t::N);
    return hash_table.update(bucket);
}


template <uint16_t k>
bool Read_CdBG_Constructor<k>::add_path_cover_edge(const Edge<k>& e)
{
    // Fetch the hash table entry for the vertices associated to the endpoints.

    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket_u = hash_table[e.u().canonical()];
    State_Read_Space& st_u = bucket_u.get_state();
    if(st_u.edge_at(e.u().side()) != cuttlefish::edge_encoding_t::E)
        return false;
    
    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket_v = hash_table[e.v().canonical()];
    State_Read_Space& st_v = bucket_v.get_state();
    if(st_v.edge_at(e.v().side()) != cuttlefish::edge_encoding_t::E)
        return false;


    st_u.update_edge_at(e.u().side(), e.u().edge());
    st_v.update_edge_at(e.v().side(), e.v().edge());

    return hash_table.update_concurrent(bucket_u, bucket_v);
}



#endif

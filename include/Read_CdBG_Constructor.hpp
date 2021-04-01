
#ifndef READ_CDBG_CONSTRUCTOR_HPP
#define READ_CDBG_CONSTRUCTOR_HPP



#include "globals.hpp"
#include "Kmer_Hash_Table.hpp"
#include "State_Read_Space.hpp"
#include "Build_Params.hpp"
#include "Thread_Pool.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_SPMC_Iterator.hpp"


// A class to construct compacted read de Bruijn graphs.
template <uint16_t k>
class Read_CdBG_Constructor
{
    friend class Thread_Pool<k>;

private:

    const Build_Params params;  // Required parameters (wrapped inside).
    Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table; // Hash table for the vertices (canonical k-mers) of the graph.
    const Kmer_Container<k + 1> edge_container; // Wrapper container for the edge-database.
    Kmer_SPMC_Iterator<k + 1> edge_parser;  // Parser for the edges from the edge-database.

    // Members required to keep track of the total number of edges processed across different threads.
    mutable Spin_Lock lock;
    mutable uint64_t edges_processed = 0;


    // Distributes the DFA-states computation task to the worker threads in the thread pool `thread_pool`.
    void distribute_states_computation(Thread_Pool<k>& thread_pool);

    // Processes the edges provided to the thread with id `thread_id`, i.e. makes state-transitions for
    // the DFA as per the edges provided to that thread.
    void process_edges(uint16_t thread_id);

    // For the vertex `v`, adds information of the incidence of an `e_v`-encoded edge to its side `s_v`
    // — making the appropriate state transition for the DFA of `v`. Returns `false` iff an attempted
    // state transition failed.
    bool add_incident_edge(const Kmer<k>& v, cuttlefish::side_t s_v, cuttlefish::edge_encoding_t e_v);

    // For the vertex `v`, adds information of the incidence of a looping edge that connects its side
    // `s_u` to its side `s_v`, which may or may not be the same sides — making the appropriate state
    // transition for the DFA of `v`. Returns `false` iff an attempted state transition failed.
    bool add_loop(const Kmer<k>& u, cuttlefish::side_t s_u, cuttlefish::side_t s_v);

    // For the vertex `v`, adds information of the incidence of a loop that connects its different
    // sides — making the appropriate state transition for the DFA of `v`. Returns `false` iff an
    // attempted state transition failed.
    bool add_crossing_loop(const Kmer<k>& v);

    // For the vertex `v`, adds information of the incidence of a loop that connects its side `s_v`
    // to that side itself — making the appropriate state transition for the DFA of `v`. Returns
    // `false` iff an attempted state transition failed.
    bool add_one_sided_loop(const Kmer<k>& v, cuttlefish::side_t s_v);


public:

    // Consructs a read-CdBG builder object, with the required parameters wrapped in `params`, and uses
    // the Cuttlefish hash table `hash_table`.
    Read_CdBG_Constructor(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table);

    // Computes the states of the DFA in the de Bruijn graph.
    void compute_DFA_states();
};


template <uint16_t k>
inline bool Read_CdBG_Constructor<k>::add_incident_edge(const Kmer<k>& v, const cuttlefish::side_t s_v, const cuttlefish::edge_encoding_t e_v)
{
    // Fetch the hash table entry for the DFA of the vertex `v`.

    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket = hash_table[v];
    State_Read_Space& state = bucket.get_state();
    cuttlefish::edge_encoding_t e_curr = state.edge_at(s_v);

    // If we've already discarded the incidence information for this side, then a self-transition happens.
    if(e_curr == cuttlefish::edge_encoding_t::N)
        return true;    // Early return w/o updating the same value again is safe — see the note at the end of the method.

    
    const cuttlefish::edge_encoding_t e_old = e_curr;
    if(e_curr == cuttlefish::edge_encoding_t::E) // This side of the vertex is encountered for the first time.
        e_curr = e_v;
    else if(e_curr != e_v)  // This side has been visited earlier with a different edge — discard the incidence information.
        e_curr = cuttlefish::edge_encoding_t::N;

    
    // We can get away without updating the same value again, because — (1) even if this DFA's state changes
    // in the hash table by the time this method completes, making no updates at this point is theoretically
    // equivalent to returning instantaneously as soon as the hash table value had been read; and also (2) the
    // ordering of the edges processed does not matter in the algorithm.
    if(e_curr == e_old)
        return true;

    state.update_edge_at(s_v, e_curr);
    return hash_table.update(bucket);
}


template <uint16_t k>
inline bool Read_CdBG_Constructor<k>::add_loop(const Kmer<k>& v, const cuttlefish::side_t s_u, const cuttlefish::side_t s_v)
{
    return s_u == s_v ? add_one_sided_loop(v, s_u) : add_crossing_loop(v);
}


template <uint16_t k>
inline bool Read_CdBG_Constructor<k>::add_crossing_loop(const Kmer<k>& v)
{
    // Fetch the hash table entry for the DFA of the vertex `v`.

    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket = hash_table[v];
    State_Read_Space& state = bucket.get_state();
    const cuttlefish::edge_encoding_t e_front = state.edge_at(cuttlefish::side_t::front);
    const cuttlefish::edge_encoding_t e_back = state.edge_at(cuttlefish::side_t::back);

    const State_Read_Space state_old = state;

    if(e_front != cuttlefish::edge_encoding_t::N)
        state.update_edge_at(cuttlefish::side_t::front, cuttlefish::edge_encoding_t::N);
    
    if(e_back != cuttlefish::edge_encoding_t::N)
        state.update_edge_at(cuttlefish::side_t::back, cuttlefish::edge_encoding_t::N);

    // We can get away without updating the same value again: see detailed comment in `add_incident_edge`.
    return state == state_old ? true : hash_table.update(bucket);
}


template <uint16_t k>
inline bool Read_CdBG_Constructor<k>::add_one_sided_loop(const Kmer<k>& v, const cuttlefish::side_t s_v)
{
    // Fetch the hash table entry for the vertex `v`.

    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket = hash_table[v];
    State_Read_Space& state = bucket.get_state();
    const cuttlefish::edge_encoding_t e_v = state.edge_at(s_v);

    // We can get away without updating the same value again: see detailed comment in `add_incident_edge`.
    if(e_v == cuttlefish::edge_encoding_t::N)
        return true;
    
    state.update_edge_at(s_v, cuttlefish::edge_encoding_t::N);
    return hash_table.update(bucket);
}



#endif

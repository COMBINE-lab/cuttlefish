
#include "Read_CdBG_Constructor.hpp"


// The following methods are not used anymore with the current algorithm.
/*
template <uint16_t k>
bool Read_CdBG_Constructor<k>::add_incident_edge(const Endpoint<k>& endpoint, cuttlefish::edge_encoding_t& e_old, cuttlefish::edge_encoding_t& e_new)
{
    // Fetch the hash table entry for the vertex associated to the endpoint.

    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket = hash_table[endpoint.hash()];
    State_Read_Space& state = bucket.get_state();
    e_new = e_old = state.edge_at(endpoint.side());


    // If we've already discarded the incidence information for this side, then a self-transition happens.
    if(e_old == cuttlefish::edge_encoding_t::N)
        return true;    // Early return w/o updating the same value again is safe — see the note at the end of the method.

    if(e_old == cuttlefish::edge_encoding_t::E) // This side of the vertex is observed for the first time.
        e_new = endpoint.edge();
    else if(e_old != endpoint.edge())   // This side has been visited earlier, but with a different edge — discard the incidence information.
        e_new = cuttlefish::edge_encoding_t::N;

    
    // We can get away without updating the same value again, because — (1) even if this DFA's state changes
    // in the hash table by the time this method completes, making no updates at this point is theoretically
    // equivalent to returning instantaneously as soon as the hash table value had been read; and also (2) the
    // ordering of the edges processed does not matter in the algorithm.
    if(e_new == e_old)
        return true;

    state.update_edge_at(endpoint.side(), e_new);
    return hash_table.update(bucket);
}


template <uint16_t k>
bool Read_CdBG_Constructor<k>::add_crossing_loop(const Endpoint<k>& endpoint, cuttlefish::edge_encoding_t& e_front, cuttlefish::edge_encoding_t& e_back)
{
    // Fetch the hash table entry for the DFA of vertex associated to the endpoint.
    
    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket = hash_table[endpoint.hash()];
    State_Read_Space& state = bucket.get_state();
    e_front = state.edge_at(cuttlefish::side_t::front);
    e_back = state.edge_at(cuttlefish::side_t::back);

    const State_Read_Space state_old = state;

    if(e_front != cuttlefish::edge_encoding_t::N)   // Discard the front-incidence information, if not done already.
        state.update_edge_at(cuttlefish::side_t::front, cuttlefish::edge_encoding_t::N);
    
    if(e_back != cuttlefish::edge_encoding_t::N)    // Discard the back-incidence information, if not done already.
        state.update_edge_at(cuttlefish::side_t::back, cuttlefish::edge_encoding_t::N);

    // We can get away without updating the same value again: see detailed comment in `add_incident_edge`.
    return state == state_old ? true : hash_table.update(bucket);
}


template <uint16_t k>
bool Read_CdBG_Constructor<k>::add_one_sided_loop(const Endpoint<k>& endpoint, cuttlefish::edge_encoding_t& e_old)
{
    // Fetch the hash table entry for the vertex associated to the endpoint.

    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket = hash_table[endpoint.hash()];
    State_Read_Space& state = bucket.get_state();
    e_old = state.edge_at(endpoint.side());

    // We can get away without updating the same value again: see detailed comment in `add_incident_edge`.
    if(e_old == cuttlefish::edge_encoding_t::N) // The incidence information has already been discarded.
        return true;

    // Discard the incidence information.
    state.update_edge_at(endpoint.side(), cuttlefish::edge_encoding_t::N);
    return hash_table.update(bucket);
}


template <uint16_t k>
inline void Read_CdBG_Constructor<k>::propagate_discard(const Endpoint<k>& v_end, const cuttlefish::edge_encoding_t e)
{
    if(e != cuttlefish::edge_encoding_t::E && e != cuttlefish::edge_encoding_t::N)  // The incident edge is unique.
        discard_neighbor_side(v_end, e);
}


template <uint16_t k>
inline void Read_CdBG_Constructor<k>::propagate_discard(const Endpoint<k>& u_end, const Endpoint<k>& v_end, const cuttlefish::edge_encoding_t e)
{
    while(!discard_side(v_end));    // Discard the neighbor `v_end`.

    propagate_discard(u_end, e);    // Discard the other neighbor.
}


template <uint16_t k>
inline bool Read_CdBG_Constructor<k>::discard_side(const Endpoint<k>& v_end)
{
    // Fetch the hash table entry for the DFA of the vertex associated to the endpoint.

    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket = hash_table[v_end.hash()];
    State_Read_Space& state = bucket.get_state();
    const cuttlefish::edge_encoding_t e_curr = state.edge_at(v_end.side());

    // We can get away without updating the same value again: see detailed comment in `add_incident_edge`.
    if(e_curr == cuttlefish::edge_encoding_t::N)    // The incidende information has already been discarded.
        return true;

    // Discard the incidence information.
    state.update_edge_at(v_end.side(), cuttlefish::edge_encoding_t::N);
    return hash_table.update(bucket);
}


template <uint16_t k>
inline void Read_CdBG_Constructor<k>::discard_neighbor_side(const Endpoint<k>& v_end, const cuttlefish::edge_encoding_t e)
{
    const Endpoint<k> w = v_end.neighbor_endpoint(e, hash_table);   // Get the neighboring endpoint connected with `e`.
    
    while(!discard_side(w));    // Discard the incidence information off that neighbor.
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG_Constructor)
*/

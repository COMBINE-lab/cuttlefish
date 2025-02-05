
#ifndef DBG_UTILITIES_HPP
#define DBG_UTILITIES_HPP



#include "DNA_Utility.hpp"
#include "globals.hpp"

#include <algorithm>


// =============================================================================
namespace cuttlefish
{
    // Returns `true` iff the edge encoding `e` is fuzzy, i.e. a unique encoding
    // is not known for the corresponding edge(s).
    bool is_fuzzy_edge(const edge_encoding_t e);

    // Returns the opposite (or complement) side of the vertex-side `s`.
    side_t opposite_side(const side_t s);


    // Replaces the sequence `seq` in-place with its reverse complement.
    template <typename T_container_> void reverse_complement(T_container_& seq);
}


// TODO: consider lessening branches everywhere.


inline bool cuttlefish::is_fuzzy_edge(const edge_encoding_t e)
{
    return e == edge_encoding_t::N || e == edge_encoding_t::E;
}


inline cuttlefish::side_t cuttlefish::opposite_side(const side_t s)
{
    return s == side_t::back ? side_t::front : side_t::back;
}


template <typename T_container_>
inline void cuttlefish::reverse_complement(T_container_& seq)
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

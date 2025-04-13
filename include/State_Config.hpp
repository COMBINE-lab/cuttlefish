
#ifndef STATE_CONFIG_HPP
#define STATE_CONFIG_HPP



#include "DNA.hpp"
#include "Source_Hash.hpp"
#include "globals.hpp"
#include "utility.hpp"

#include <cstdint>
#include <climits>
#include <cassert>


namespace cuttlefish
{


// Class to count incident edges' ((k + 1)-mers') frequencies of vertices
// (k-mers).
// TODO: why is there no specialization for references? No advantages due to alignment-requirement issues?
class Edge_Frequency
{
private:

    static constexpr uint32_t max_f = 0b1111;   // Maximum supported frequency of a (k + 1)-mer.
    static uint8_t f_th;    // Edge-frequency threshold.

    uint32_t e_f;   // Packed-frequency of an associated vertex's incident edges.


    // Returns the frequency stored at offset `off`.
    uint32_t f_at(uint32_t off) const;


public:

    // Constructs an empty counter.
    Edge_Frequency(): e_f(0) {}

    // Sets the edge-frequency threshold to `f_th`.
    static void set_edge_threshold(uint8_t f_th);

    // Adds the edge-encoding `e` to the side `s` of a corresponding
    // vertex.
    template <side_t s> void add_edge(base_t e);

    // Returns the `Base`-encoding of the edge(s) passing frequency threshold
    // and incident to the side `s` of a vertex having this state.
    base_t edge_at(side_t s) const;

    // Returns the number of edges at side `s` of a corresponding vertex.
    uint32_t edge_count(side_t s) const;
};


// =============================================================================
// Neighborhood information of a vertex in a de Bruijn graph.
class Vertex_Neighborhood
{
private:

    Edge_Frequency e_f; // Frequency of the vertex's edges.

protected:

    // Constructs an empty neighborhood.
    Vertex_Neighborhood(): e_f() {}

public:

    // Adds the edge-encodings `front` and `back` to the associated sides of a
    // corresponding vertex.
    void update_edges(base_t front, base_t back);

    // Returns the `Base`-encoding of the edge(s) incident to the side `s` of a
    // vertex having this neighborhood.
    base_t edge_at(side_t s) const { return e_f.edge_at(s); }

    // Returns `true` iff some vertex having this neighborhood is branching
    // (i.e. has multiple incident edges) at its side `s`.
    bool is_branching_side(side_t s) const { return e_f.edge_count(s) > 1; }

    // Returns `true` iff some vertex having this neighborhood is empty (i.e.
    // no incident edges) at its side `s`.
    bool is_empty_side(side_t s) const { return e_f.edge_count(s) == 0; }

    // Returns `true` iff some vertex having this neighborhood is isolated off
    // the rest of the graph.
    bool is_isolated() const { return is_empty_side(side_t::back) && is_empty_side(side_t::front); }
};


// =============================================================================
// Status information of a vertex in a de Bruijn graph.
template <typename T_status_>
class Vertex_Status
{
private:

    static constexpr T_status_ visited = 0b0000'0001;   // Flag to denote a vertex as visited.
    static constexpr T_status_ discontinuity[2] = {0b0000'0010, 0b0000'0100};   // Flags to denote a vertex's sides as discontinuous.

protected:

    T_status_ status;   // Some status information of the vertex, bit-packed:
                            // whether it has discontinuous sides, whether it's been visited.


    // Constructs a clear status.
    Vertex_Status(): status(0) {}

    // Marks the associated vertex as discontinuous at side `s`, if `s` is a
    // valid side.
    void mark_discontinuous_optional(side_t s) { assert(as_int(s) <= 2); status |= (s == side_t::unspecified ? 0 : discontinuity[as_int(s)]); }

public:

    // Marks the associated vertex as visited.
    void mark_visited() { status |= visited; }

    // Marks the associated vertex as discontinuous at side `s`.
    void mark_discontinuous(side_t s) { assert(as_int(s) < 2); status |= discontinuity[as_int(s)]; }

    // Returns whether the associated vertex is visited.
    bool is_visited() const { return status & visited; }

    // Returns whether the associated vertex is discontinuous at side `s`.
    bool is_discontinuous(const side_t s) const { assert(as_int(s) < 2); return status & discontinuity[as_int(s)]; }

    // Returns whether the associated vertex has any discontinuous side.
    bool is_discontinuity() const;
};


// =============================================================================
// Class for a full state-configuration of a vertex in a de Bruijn graph: this
// is a configuration attached to vertices in subgraphs. `Colored_` denotes
// whether the state has color / annotation metadata associated to it.
template <bool Colored_>
class State_Config
{};


template <>
class State_Config<false>: public Vertex_Neighborhood, public Vertex_Status<uint8_t>
{
public:

    // Constructs an empty state.
    State_Config(): Vertex_Neighborhood(), Vertex_Status() {}

    // Adds the edge-encodings `front` and `back` to the associated sides of a
    // corresponding vertex, and marks the associated vertex as discontinuous
    // at sides `s_0` and `s_1`.
    void update(base_t front, base_t back, side_t s_0, side_t s_1)
    {
        update_edges(front, back);
        mark_discontinuous_optional(s_0), mark_discontinuous_optional(s_1);
    }
};


template <>
class State_Config<true> : public Vertex_Neighborhood, public Vertex_Status<uint32_t>
{
private:

    static constexpr uint32_t source_pos = 11;
    static constexpr uint32_t source_mask = 0x1FFFFF << source_pos; // Flag to extract last observed source's ID.

    uint64_t color_hash_;   // Hash of the associated color-set.

public:

    // Constructs an empty state.
    State_Config(): Vertex_Neighborhood(), Vertex_Status(), color_hash_(0) {}

    // Returns the hash of the associated color-set.
    uint64_t color_hash() const { return color_hash_; }

    // Adds the source ID `source` to the color-set of this state.
    void add_source(const source_id_t source)
    {
        assert(source <= (source_mask >> source_pos));
        const auto last_source = (status & source_mask) >> source_pos;
        if(source != last_source)   // Dealing with the problem of hashing multisets.
        {
            // assert(source >= last_source);  // Ensure sortedness of source-IDs.
            color_hash_ = hash_combine(color_hash_, source_hash(source));
            status = (status & ~source_mask) | (source << source_pos);
        }
    }
};


inline uint32_t Edge_Frequency::f_at(const uint32_t off) const
{
    assert(off % 4 == 0);

    return (e_f & (max_f << off)) >> off;
}


template <side_t s>
inline void Edge_Frequency::add_edge(const base_t e)
{
    static_assert(s != side_t::unspecified);
    assert(base_t::A <= e && e <= base_t::T);

    const uint32_t off = (s == side_t::front ? 0 : 16) + 4 * e;
    const uint32_t f_mask = max_f << off;
    const uint32_t f = (e_f & f_mask) >> off;
    e_f = (f < max_f ? (e_f & ~f_mask) | ((f + 1) << off) : e_f);
}


inline uint32_t Edge_Frequency::edge_count(const side_t s) const
{
    const uint32_t off = (s == side_t::front ? 0 : 16);
    const auto f_A = f_at(off + 0);
    const auto f_C = f_at(off + 4);
    const auto f_G = f_at(off + 8);
    const auto f_T = f_at(off + 12);

    return (f_A >= f_th) + (f_C >= f_th) + (f_G >= f_th) + (f_T >= f_th);
}


inline base_t Edge_Frequency::edge_at(const side_t s) const
{
    const uint32_t off = (s == side_t::front ? 0 : 16);
    const auto e_A = (f_at(off + 0) >= f_th);
    const auto e_C = (f_at(off + 4) >= f_th);
    const auto e_G = (f_at(off + 8) >= f_th);
    const auto e_T = (f_at(off + 12) >= f_th);
    const uint32_t e_c = e_A + e_C + e_G + e_T;

    switch(e_c)
    {
        case 0:
            return base_t::E;

        case 1:
            return static_cast<base_t>((e_A * base_t::A) + (e_C * base_t::C) + (e_G * base_t::G) + (e_T * base_t::T));

        default:
            return base_t::N;
    }
}


inline void Vertex_Neighborhood::update_edges(const base_t front, const base_t back)
{
    constexpr auto E = base_t::E;
    constexpr auto T = base_t::T;
    (void)T;

    assert(front == E || front <= T);
    if(front != E)
        e_f.add_edge<side_t::front>(front);

    assert(back == E || back <= T);
    if(back != E)
        e_f.add_edge<side_t::back>(back);
}


template <typename T_status_>
inline bool Vertex_Status<T_status_>::is_discontinuity() const
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wbitwise-instead-of-logical"

    return is_discontinuous(side_t::front) | is_discontinuous(side_t::back);

#pragma GCC diagnostic pop
}

}



#endif

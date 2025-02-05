
#ifndef DISCONTINUITY_GRAPH_CONTRACTOR_HPP
#define DISCONTINUITY_GRAPH_CONTRACTOR_HPP



#include "dBG_Contractor.hpp"
#include "Kmer.hpp"
#include "Kmer_Hasher.hpp"
#include "Discontinuity_Edge.hpp"
#include "Discontinuity_Graph.hpp"
#include "Concurrent_Hash_Table.hpp"
#include "utility.hpp"
#include "unordered_dense/unordered_dense.h"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <atomic>
#include <cassert>


class Data_Logistics;


namespace cuttlefish
{

// =============================================================================
// Contractor of discontinuity-graphs.
template <uint16_t k, bool Colored_>
class Discontinuity_Graph_Contractor
{
private:

    Discontinuity_Graph<k, Colored_>& G;    // The discontinuity-graph.

    typedef typename dBG_Contractor<k>::P_v_t P_v_t;
    P_v_t& P_v; // `P_v[j]` contains path-info for vertices in partition `j`—specifically, the meta-vertices.

    const std::string compressed_diagonal_path; // Path-prefix to the edges introduced in contracting diagonal blocks.

    class Other_End;
    Concurrent_Hash_Table<Kmer<k>, Other_End, Kmer_Hasher<k>> M;    // `M[v]` is the associated vertex to `v` at a given time.

    // TODO: remove `D_j` by adopting a more parallelization-amenable algorithm for diagonal contraction-expansion.
    std::vector<Discontinuity_Edge<k>> D_j; // Edges introduced in contracting a diagonal block. TODO: remove `D_j` by adding these edges to the diagonal block.
    std::vector<Padded<std::vector<Discontinuity_Edge<k>>>> D_c;    // `D_c[t]` contains the edges corresponding to compressed diagonal chains by worker `t`.
    Buffer<Discontinuity_Edge<k>> D_c_flat; // Flattened `D_c`.

    ankerl::unordered_dense::map<Kmer<k>, Other_End, Kmer_Hasher<k>> D; // `D[v]` is the associated vertex to `v` at a given time during diagonal compression.

    std::atomic_uint64_t phantom_count_;    // Number of phantom edges.
    std::atomic_uint64_t icc_count; // Number of ICCs.


    // Contracts the `[j, j]`'th edge-block. Buffer `buf` is used to read the
    // edges.
    void contract_diagonal_block(std::size_t j, Buffer<Discontinuity_Edge<k>>& buf);

    // Forms a meta-vertex in the contracted graph with the vertex `v` belonging
    // to the vertex-partition `part`. In the contracted graph, `v` has a `w_1`
    // weighted edge incident to its side `s_1` and a `w_2` weighted edge
    // incident to the other side. `is_cycle` denotes whether the meta-vertex
    // corresponds to a cycle.
    void form_meta_vertex(Kmer<k> v, std::size_t part, side_t s_1, weight_t w_1, weight_t w_2, bool is_cycle = false);

    // Forms a meta-vertex in the contracted graph with the vertex `v` belonging
    // to the vertex-partition `part`. In the contracted graph, `v` has a `w`-
    // weighted edge incident to its side `s`. `is_cycle` denotes whether the
    // meta-vertex corresponds to a cycle.
    void form_meta_vertex(Kmer<k> v, std::size_t part, side_t s, weight_t w, bool is_cycle = false);

    static constexpr auto now = std::chrono::high_resolution_clock::now;    // Current time-point in nanoseconds.

    // Returns the equivalent time-duration in seconds from `d` nanoseconds.
    static constexpr auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };


public:

    // Constructs a contractor for the discontinuity-graph `G`. `P_v[j]` is to
    // contain path-information for vertices at partition `j`. `logistics` is
    // the data logistics manager for the algorithm execution.
    Discontinuity_Graph_Contractor(Discontinuity_Graph<k, Colored_>& G, P_v_t& P_v, const Data_Logistics& logistics);

    // Contracts the discontinuity-graph.
    void contract();
};


// =============================================================================
// Other endpoint `v` associated to a current vertex `u` through an edge.
template <uint16_t k, bool Colored_>
class Discontinuity_Graph_Contractor<k, Colored_>::Other_End
{
private:

    Kmer<k> v_; // The other endpoint vertex `v`.
    side_t s_v_;    // Side of the endpoint `v` to which the associated edge is incident to.
    side_t s_u_;    // Side of the current vertex `u` to which the associated edge is incident to.
    bool is_phi_;   // Whether the endpoint is a ϕ vertex.
    weight_t w_;    // Weight of the associated edge.
    bool in_same_part_; // Whether the endpoints belong to the same partition.
    bool processed_; // Whether the endpoint has been processed, defined by the context.


public:

    // TODO: remove once we have our own hash-table.
    Other_End()
    {}

    // Constructs an endpoint with the vertex `v`, connected through its side
    // `s_v` to the current vertex's side `s_u`. `is_phi` should be `true` iff
    // `v` is the ϕ vertex. The connecting edge has weight `w`, and
    // `in_same_part` should be `true` iff the endpoints of the edge belong to
    // the same partition.
    Other_End(const Kmer<k>& v, const side_t s_v, const side_t s_u, const bool is_phi, const weight_t w, const std::size_t in_same_part, const bool processed = false):
          v_(v)
        , s_v_(s_v)
        , s_u_(s_u)
        , is_phi_(is_phi)
        , w_(w)
        , in_same_part_(in_same_part)
        , processed_(processed)
    {}

    // Mark the endpoint as processed, defined by the context.
    void process() { processed_ = true; }

    // Returns the endpoint vertex.
    auto v() const { return v_; }

    // Returns the side of the endpoint to which the associated edge is incident
    // to.
    auto s_v() const { return s_v_; }

    // Returns the side of the current vertex `u` to which the associated edge
    // is incident to.
    auto s_u() const { return s_u_; }

    // Returns whether the endpoint is a ϕ vertex.
    auto is_phi() const { return is_phi_; }

    // Returns the weight of the associated edge.
    auto w() const { return w_; }

    // Returns whether the endpoints belong to the same partition.
    auto in_same_part() const { return in_same_part_; }

    // Returns whether the endpoint has been processed, defined by the context.
    auto processed() const { return processed_; }
};


template <uint16_t k, bool Colored_>
inline void Discontinuity_Graph_Contractor<k, Colored_>::form_meta_vertex(const Kmer<k> v, const std::size_t part, const side_t s_1, const weight_t w_1, const weight_t w_2, const bool is_cycle)
{
    assert(w_1 > 0); assert(w_2 > 0);
    form_meta_vertex(v, part, side_t::front, s_1 == side_t::front ? w_1 : w_2, is_cycle);
}


template <uint16_t k, bool Colored_>
inline void Discontinuity_Graph_Contractor<k, Colored_>::form_meta_vertex(const Kmer<k> v, const std::size_t part, const side_t s, const weight_t w, const bool is_cycle)
{
    assert(w > 0);
    assert(part < P_v.size());
    P_v[part].unwrap().emplace(v, v, w, inv_side(s), is_cycle);   // The path-traversal enters `v` through its the side `s`.
}

}



#endif

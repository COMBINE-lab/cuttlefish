
#ifndef CONTRACTED_GRAPH_EXPANDER_HPP
#define CONTRACTED_GRAPH_EXPANDER_HPP



#include "dBG_Contractor.hpp"
#include "Discontinuity_Graph.hpp"
#include "Discontinuity_Edge.hpp"
#include "Path_Info.hpp"
#include "Kmer.hpp"
#include "Kmer_Hasher.hpp"
#include "Concurrent_Hash_Table.hpp"
#include "globals.hpp"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>


class Data_Logistics;


namespace cuttlefish
{

// =============================================================================
// Expander for contracted discontinuity-graphs.
template <uint16_t k, bool Colored_>
class Contracted_Graph_Expander
{
private:

    Discontinuity_Graph<k, Colored_>& G;    // The (augmented) discontinuity graph.

    typedef typename dBG_Contractor<k>::P_v_t P_v_t;
    typedef typename dBG_Contractor<k>::P_e_t P_e_t;
    P_v_t& P_v; // `P_v[i]` contains path-info for vertices in partition `i`.
    P_e_t& P_e; // `P_e[b]` contains path-info for edges in bucket `b`.

    const std::string compressed_diagonal_path; // Path-prefix to the edges introduced in contracting diagonal blocks.

    // TODO: remove `D_i` by adopting a more parallelization-amenable algorithm for diagonal contraction-expansion.
    Buffer<Discontinuity_Edge<k>> D_i;  // New edges introduced in contracted diagonal blocks.

    Concurrent_Hash_Table<Kmer<k>, Path_Info<k>, Kmer_Hasher<k, 0>> M;  // `M[v]` is the path-info for vertex `v`.


    // Expands the `[i, i]`'th (contracted) edge-block.
    void expand_diagonal_block(std::size_t i);

    // Infers a vertex v's path-info from that of vertex u's path-info `u_inf`
    // and returns it. The vertices are connected with an edge of weight `w`
    // through their sides `s_v` and `s_u`, respectively.
    Path_Info<k> infer(Path_Info<k> u_inf, side_t s_u, side_t s_v, weight_t w);

    // Computes the path-info of the edge `e` from its endpoints' path-info,
    // `u_inf` and `v_inf`, and adds the info to `e`'s path-info bucket.
    void add_edge_path_info(const Discontinuity_Edge<k>& e, Path_Info<k> u_inf, Path_Info<k> v_inf);

    // Computes the path-info of the diagonal edge `e` from its first endpoint
    // `u`'s path-info `u_inf`, and adds the info to `e`'s path-info bucket.
    void add_diagonal_edge_path_info(const Discontinuity_Edge<k>& e, Path_Info<k> u_inf);

    // Computes the path-info of the edge `e` of form `(ϕ, v)` from `v`'s path-
    // info, `v_inf`, and adds the info to `e`'s path-info bucket.
    void add_edge_path_info(const Discontinuity_Edge<k>& e, Path_Info<k> v_inf);

    // Collates worker local buffers in ID range `[beg, end)` from the
    // collection `source` into the global repository `dest`. Also clears the
    // local buffers.
    // template <typename T_s_, typename T_d_> static uint64_t collate_w_local_bufs(T_s_& source, std::size_t beg, size_t end, T_d_& dest);


    // TODO: replace w/ the general timer.
    static constexpr auto now = std::chrono::high_resolution_clock::now;    // Current time-point in nanoseconds.

    // Returns the equivalent time-duration in seconds from `d` nanoseconds.
    static constexpr auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };



public:

    // Constructs an expander for the contracted discontinuity-graph `G`.
    // `P_v[i]` is to contain path-information for vertices at partition `i`,
    // and `P_e[b]` is to contain path-information for edges at bucket `b`.
    // `logistics` is the data logistics manager for the algorithm execution.
    Contracted_Graph_Expander(Discontinuity_Graph<k, Colored_>& G, P_v_t& P_v, P_e_t& P_e, const Data_Logistics& logistics);

    // Expands the contracted discontinuity-graph.
    void expand();
};


template <uint16_t k, bool Colored_>
inline Path_Info<k> Contracted_Graph_Expander<k, Colored_>::infer(const Path_Info<k> u_inf, const side_t s_u, const side_t s_v, const weight_t w)
{
    assert(u_inf.r() > 0);

    // const auto p_v = u_inf.p(); // Path-ID.
    const auto r_v = (s_u == u_inf.o() ? u_inf.r() + w :    // Rank.
                                        (u_inf.r() > w ? u_inf.r() - w : 0));   // Trying to expand crossing a deleted edge from an ICC.
                                                                                // This works as no vertex can have a rank `0` in the model.
    const auto o_v = (s_u == u_inf.o() ? inv_side(s_v) : s_v);  // Orientation.
    // const auto is_cycle = u_inf.is_cycle();

    return Path_Info(u_inf.p(), r_v, o_v, u_inf.is_cycle());
}


template <uint16_t k, bool Colored_>
inline void Contracted_Graph_Expander<k, Colored_>::add_edge_path_info(const Discontinuity_Edge<k>& e, const Path_Info<k> u_inf, const Path_Info<k> v_inf)
{
    assert(e.w() == 1);
    assert(!e.x_is_phi() && !e.y_is_phi());
    assert(u_inf.p() == v_inf.p());

    // const auto p = u_inf.p();
    const auto r = std::min(u_inf.r(), v_inf.r());
    const auto o = (r == u_inf.r() ? e.o() : inv_side(e.o()));

    assert(e.b() > 0 && e.b() < P_e.size());
    P_e[e.b()].unwrap().emplace(e.b_idx(), u_inf.p(), r, o, u_inf.is_cycle());
}


template <uint16_t k, bool Colored_>
inline void Contracted_Graph_Expander<k, Colored_>::add_edge_path_info(const Discontinuity_Edge<k>& e, const Path_Info<k> v_inf)
{
    assert(e.w() == 1);
    assert(e.x_is_phi() && !e.y_is_phi());

    // const auto p = v_inf.p();
    const auto r = (v_inf.r() == 1 ? 0 : v_inf.r());
    const auto o = (r == 0 ? e.o() : inv_side(e.o()));

    assert(e.b() > 0 && e.b() < P_e.size());
    P_e[e.b()].unwrap().emplace(e.b_idx(), v_inf.p(), r, o, v_inf.is_cycle());
}


template <uint16_t k, bool Colored_>
inline void Contracted_Graph_Expander<k, Colored_>::add_diagonal_edge_path_info(const Discontinuity_Edge<k>& e, const Path_Info<k> u_inf)
{
    assert(e.w() == 1);
    assert(!e.x_is_phi() && !e.y_is_phi());

    // Edges in cycles belonging to diagonal blocks form a special case.
    // When the rank-1 vertex `v_1` in cycle `v_1, ..., v_p` propagates info
    // to `v_p` through their shared edge (the propagation cannot go the other
    // way due to the meta-vertex formation process for cycles), if `P(v_1) ≠
    // P(v_p)`, then `v_p` gets a "relative"-rank 0 from `v_1` (although
    // discarded), and their shared edge as a result gets ranked 0. Whereas for
    // the other case, i.e. when they are in the same partition, the rank of
    // the diagonal edges are computed in a different manner: the correct ranks
    // of `v_1` and `v_p` are already computed when the edge's rank is getting
    // computed. So the relative ranking capturing successive-ness disappears,
    // and needs to be introduced again.
    // Note that `(v_1, v_p)` need not necessarily be `(u, v)` in `e`. Hence,
    // `e` may get a rank `0` or `p + 1`, which does not matter in a cycle.
    const auto t = infer(u_inf, e.s_u(), e.s_v(), e.w());
    add_edge_path_info(e, u_inf, t);
}

}



#endif

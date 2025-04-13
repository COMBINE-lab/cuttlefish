
#ifndef DISCONTINUITY_EDGE_HPP
#define DISCONTINUITY_EDGE_HPP



#include "Kmer.hpp"
#include "globals.hpp"

#include <cstdint>
#include <iostream>


namespace cuttlefish
{

// =============================================================================
// An edge `e_{u, v} = ({(u, s_u), (v, s_v)}, w, e_b)` in a discontinuity-graph
// of `k`-mers.
template <uint16_t k>
class Discontinuity_Edge
{
private:

    // TODO: try packing booleans and sides.

    Kmer<k> u_; // An endpoint of the edge.
    Kmer<k> v_; // An endpoint of the edge.
    weight_t weight;    // Weight of the edge.
    uint16_t bucket_id; // ID of the bucket of the unitig corresponding to the edge.
    uni_idx_t b_idx_;   // Index of the corresponding unitig within its bucket.

    uint8_t mask;   // Bitmask for: each vertex's side to which the edge is incident to (u: 0, v: 1);
                    //              whether each vertex is ϕ (u: 2, v: 3);
                    //              and the exit-orientation of the corresponding literal unitig wrt the `(u, v)` orientation of the edge (6).


    static constexpr uint8_t side_u[2] = {0b0000'0000, 0b0000'0001};    // Flags to denote `u`'s side (`front: 0, back: 1`).
    static constexpr uint8_t side_v[2] = {0b0000'0000, 0b0000'0010};    // Flags to denote `v`'s side (`front: 0, back: 1`).
    static constexpr uint8_t phi_u[2]  = {0b0000'0000, 0b0000'0100};    // Flags to denote whether `u` is ϕ.
    static constexpr uint8_t phi_v[2]  = {0b0000'0000, 0b0000'1000};    // Flags to denote whether `v` is ϕ.
    static constexpr uint8_t unitig_o[2]  = {0b0000'0000, 0b0100'0000}; // Flags to denote the exit-orientation of the corresponding literal unitig wrt to the `(u, v)` orientation of the edge (`front: 0, back: 1`).


public:

    // Constructs a placeholder edge.
    Discontinuity_Edge(){}

    // Constructs an edge `{(u, s_u), (v, s_v)}` between the vertices `u` and
    // `v` that connects there sides `s_u` and `s_v` resp. It has weight `w`.
    // The locally-maximal unitig corresponding to this edge is stored in the
    // `b`'th bucket, at index `b_idx`. `u_is_phi` and `v_is_phi` denote whether
    // `u` and `v` are ϕ, respectively. `o` is the exit-orientation of the
    // corresponding literal unitig (if any) wrt to the `(u, v)` orientation.
    Discontinuity_Edge(const Kmer<k>& u, side_t s_u, const Kmer<k>& v, side_t s_v, weight_t w, uint16_t b, std::size_t b_idx, bool u_is_phi, bool v_is_phi, side_t o);

    // Returns the `u` endpoint of the edge.
    const Kmer<k>& u() const { return u_; }

    // Returns the `v` endpoint of the edge.
    const Kmer<k>& v() const { return v_; }

    // Returns the side of the `u` endpoint to which the edge is incident to.
    side_t s_u() const { return side_t(bool(mask & side_u[1])); }

    // Returns the side of the `v` endpoint to which the edge is incident to.
    side_t s_v() const { return side_t(bool(mask & side_v[1])); }

    // Returns the `u` endpoint of the edge.
    const Kmer<k>& x() const { return u(); }

    // Returns the `v` endpoint of the edge.
    const Kmer<k>& y() const { return v(); }

    // Returns the side of the `u` endpoint to which the edge is incident to.
    side_t s_x() const { return s_u(); }

    // Returns the side of the `v` endpoint to which the edge is incident to.
    side_t s_y() const { return s_v(); }

    // Returns the weight of the edge.
    weight_t w() const { return weight; }

    // Returns the ID of the bucket of this edge.
    uint16_t b() const { return bucket_id; }

    // Returns the index of the corresponding unitig within its bucket.
    std::size_t b_idx() const { return b_idx_; }

    // Returns whether `u` is the ϕ vertex.
    bool u_is_phi() const { return mask & phi_u[1]; }

    // Returns whether `v` is the ϕ vertex.
    bool v_is_phi() const { return mask & phi_v[1]; }

    // Returns whether `u` is the ϕ vertex.
    bool x_is_phi() const { return u_is_phi(); }

    // Returns whether `v` is the ϕ vertex.
    bool y_is_phi() const { return v_is_phi(); }

    // Returns exit-orientation of the corresponding literal unitig wrt the
    // `(u, v)` orientation of the edge.
    side_t o() const { return side_t(bool(mask & unitig_o[1])); }

    // (De)serializes the edge from / to the `cereal` archive `archive`.
    template <typename T_archive_> void serialize(T_archive_& archive);

    // Pretty-prints the edge `e` to the stream `os`.
    friend std::ostream& operator<<(std::ostream& os, const Discontinuity_Edge<k>& e)
    {
        return os   << "{"
                    << "(" << e.u() << ", " << (e.s_u() == side_t::front ? 'f' : 'b') << "), "
                    << "(" << e.v() << ", " << (e.s_v() == side_t::front ? 'f' : 'b') << ")"
                    << "}: " << e.w();
    }
};


template <uint16_t k>
inline Discontinuity_Edge<k>::Discontinuity_Edge(const Kmer<k>& u, const side_t s_u, const Kmer<k>& v, const side_t s_v, const weight_t w, const uint16_t b, const std::size_t b_idx, const bool u_is_phi, const bool v_is_phi, const side_t o):
      u_(u)
    , v_(v)
    , weight(w)
    , bucket_id(b)
    , b_idx_(b_idx)
    , mask( side_u[s_u == side_t::back] | side_v[s_v == side_t::back] |
            phi_u[u_is_phi] | phi_v[v_is_phi] |
            unitig_o[o == side_t::back])
{
    // Necessary condition for the status-mask to work.
    static_assert(as_int(side_t::front) == 0 && as_int(side_t::back) == 1);
}


template <uint16_t k>
template <typename  T_archive_>
inline void Discontinuity_Edge<k>::serialize(T_archive_& archive)
{
    archive(u_, v_, weight, bucket_id, b_idx_, mask);
}

}



#endif

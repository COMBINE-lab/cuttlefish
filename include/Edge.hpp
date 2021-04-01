
#ifndef EDGE_HPP
#define EDGE_HPP




#include "globals.hpp"

#include <cstdint>


// Class for an instance of a bidirected edge.
// NB: for some (k + 1)-mer `e`, `e` and `e_bar` denote the same bidirected edge `e_hat`;
// but these being different (k + 1)-mers, they are treated as different *instances* of the
// same edge. Semantically, the underlying edges are the same. This edge instance is in the
// tuple form `(u, s_\hat{u}, v, s_\hat{v})`.
template <uint16_t k>
class Edge
{
private:

    Kmer<k + 1> e_; // The edge (k + 1)-mer (need not be in canonical form).
    Kmer<k> u_; // One endpoint k-mer of this edge instance — source k-mer of the `e` form.
    Kmer<k> v_; // One endpoint k-mer of this edge instance — sink k-mer of the `e` form.
    Kmer<k> u_bar_; // Reverse complement of `u`.
    Kmer<k> v_bar_; // Reverse complement of `v`.
    const Kmer<k>* u_hat_ptr;   // Pointer to the canonical form of the k-mer `u`, i.e. ptr to `min(u, u_bar)`.
    const Kmer<k>* v_hat_ptr;   // Pointer to the canonical form of the k-mer `v`, i.e. ptr to `min(v, v_bar)`.
    cuttlefish::side_t s_u_hat_;    // The side of the vertex `u_hat` to which this edge instance is incident to.
    cuttlefish::side_t s_v_hat_;    // The side of the vertex `v_hat` to which this edge instance is incident to.


public:

    // Returns a mutable reference to the edge (k + 1)-mer.
    Kmer<k + 1>& e();

    // Configures the edge data. i.e. sets the relevant information of the
    // edge from the underlying (k + 1)-mer. Must be used whenever the edge
    // (k + 1)-mer (updatable using `e()`) is modified.
    void configure();

    // Returns `true` iff the edge is a loop.
    bool is_loop() const;

    // Returns the vertex (i.e. canonical k-mer) `u_hat` — which corresponds
    // to the source endpoint k-mer `u` of this edge instance `e`.
    const Kmer<k>& u_hat() const;

    // Returns the vertex (i.e. canonical k-mer) `v_hat` — which corresponds
    // to the sink endpoint k-mer `v` of this edge instance `e`.
    const Kmer<k>& v_hat() const;

    // Returns the side of the vertex (i.e. canonical k-mer) `u_hat` — which
    // corresponds to the source endpoint k-mer `u` of this edge instance `e`
    // — to which this edge instance `e` is incident to.
    cuttlefish::side_t s_u_hat() const;

    // Returns the side of the vertex (i.e. canonical k-mer) `v_hat` — which
    // corresponds to the sink endpoint k-mer `v` of this edge instance `e`
    // — to which this edge instance `e` is incident to.
    cuttlefish::side_t s_v_hat() const;

    // Returns the `DNA::Extended_Base` base-encoding of the underlying edge,
    // from the point-of-view of the vertex `u_hat`.
    cuttlefish::edge_encoding_t edge_encoding_u() const;

    // Returns the `DNA::Extended_Base` base-encoding of the underlying edge,
    // from the point-of-view of the vertex `v_hat`.
    cuttlefish::edge_encoding_t edge_encoding_v() const;
};


template <uint16_t k>
inline Kmer<k + 1>& Edge<k>::e()
{
    return e_;
}


template <uint16_t k>
inline void Edge<k>::configure()
{
    u_.from_prefix(e_),
    v_.from_suffix(e_);

    u_bar_.as_reverse_complement(u_),
    v_bar_.as_reverse_complement(v_);

    u_hat_ptr = Kmer<k>::canonical(u_, u_bar_),
    v_hat_ptr = Kmer<k>::canonical(v_, v_bar_);

    s_u_hat_ = (&u_ == u_hat_ptr ? cuttlefish::side_t::back : cuttlefish::side_t::front);
    s_v_hat_ = (&v_ == v_hat_ptr ? cuttlefish::side_t::front : cuttlefish::side_t::back);
}


template <uint16_t k>
inline bool Edge<k>::is_loop() const
{
    return *u_hat_ptr == *v_hat_ptr;
}


template <uint16_t k>
inline const Kmer<k>& Edge<k>::u_hat() const
{
    return *u_hat_ptr;
}


template <uint16_t k>
inline const Kmer<k>& Edge<k>::v_hat() const
{
    return *v_hat_ptr;
}


template <uint16_t k>
inline cuttlefish::side_t Edge<k>::s_u_hat() const
{
    return s_u_hat_;
}


template <uint16_t k>
inline cuttlefish::side_t Edge<k>::s_v_hat() const
{
    return s_v_hat_;
}


template <uint16_t k>
inline cuttlefish::edge_encoding_t Edge<k>::edge_encoding_u() const
{
    const DNA::Base base = (s_u_hat_ == cuttlefish::side_t::back ? e_.back() : DNA_Utility::complement(e_.back()));

    return DNA_Utility::map_extended_base(base);
}


template <uint16_t k>
inline cuttlefish::edge_encoding_t Edge<k>::edge_encoding_v() const
{
    const DNA::Base base = (s_v_hat_ == cuttlefish::side_t::front ? e_.front() : DNA_Utility::complement(e_.front()));

    return DNA_Utility::map_extended_base(base);
}



#endif

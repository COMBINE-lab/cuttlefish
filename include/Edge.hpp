
#ifndef EDGE_HPP
#define EDGE_HPP




#include "globals.hpp"
#include "Endpoint.hpp"

#include <cstdint>


template <uint16_t k, uint8_t BITS_PER_KEY> class Kmer_Hash_Table;


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
    Endpoint<k> u_; // One endpoint k-mer of this edge instance — source of the `e` form.
    Endpoint<k> v_; // One endpoint k-mer of this edge instance — sink of the `e` form.


public:

    // Returns a mutable reference to the edge (k + 1)-mer.
    Kmer<k + 1>& e();

    // Returns the source endpoint `u` of the edge instance.
    const Endpoint<k>& u() const;

    // Returns the sink endpoint `v` of the edge instance.
    const Endpoint<k>& v() const;

    // Configures the edge data, i.e. sets the relevant information of the
    // edge from the underlying (k + 1)-mer. Uses the hash table `hash` to
    // get the hash values of the endpoint vertices. Must be used whenever
    // the edge (k + 1)-mer (updatable using `e()`) is modified.
    void configure(const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Returns `true` iff the edge is a loop.
    bool is_loop() const;
};


template <uint16_t k>
inline Kmer<k + 1>& Edge<k>::e()
{
    return e_;
}


template <uint16_t k>
inline const Endpoint<k>& Edge<k>::u() const
{
    return u_;
}


template <uint16_t k>
inline const Endpoint<k>& Edge<k>::v() const
{
    return v_;
}


template <uint16_t k>
inline void Edge<k>::configure(const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
    u_.from_prefix(e_, hash),
    v_.from_suffix(e_, hash);
}


template <uint16_t k>
inline bool Edge<k>::is_loop() const
{
    return u_.canonical() == v_.canonical();
}



#endif

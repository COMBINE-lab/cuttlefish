
#ifndef ENDPOINT_HPP
#define ENDPOINT_HPP




#include "Kmer.hpp"
#include "globals.hpp"
#include "Kmer_Hash_Table.hpp"

#include <cstdint>


// A class denoting an endpoint of a bidirected edge instance.
template <uint16_t k>
class Endpoint
{
private:

    Kmer<k> kmer_;  // The endpoint k-mer spelled by the edge instance.
    Kmer<k> kmer_bar_;  // Reverse complement of the k-mer spelled by the edge instance.
    const Kmer<k>* kmer_hat_ptr;    // Pointer to the canonical form of the endpoint k-mer.
    cuttlefish::side_t s;   // The side of the endpoint vertex to which the edge instance is incident to.
    cuttlefish::edge_encoding_t e;  // The `DNA::Extended_Base` encoding of the edge instance incident to this endpoint.
    uint64_t h; // Hash value of the vertex, i.e. canonical k-mer.


    // Constructs an `Endpoint` object that appears in the form `kmer` in an edge instance, and
    // is the source (i.e. prefix) of that edge iff `is_source` is true — which decides the edge
    // incidence side to the corresponding vertex. Also gets the hash value of the vertex using
    // the hash table `hash`. The sole application of this constructor is get a specific side of
    // a vertex where the edge incidence information is to be discarded, hence no edge-encoding
    // is provided with the constructor, although the class has such a member.
    Endpoint(const Kmer<k>& kmer, bool is_source, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Returns the side of the associated vertex to which the edge instance corresponding to this
    // endpoint is incident to, if this endpoint is the source endpoint of the edge.
    cuttlefish::side_t exit_side() const;

    // Returns the side of the associated vertex to which the edge instance corresponding to this
    // endpoint is incident to, if this endpoint is the sink endpoint of the edge.
    cuttlefish::side_t entrance_side() const;

    // Returns the `DNA::Extended_Base` encoding of the edge `e` corresponding to this endpoint,
    // given the endpoint is the source endpoint of the edge.
    cuttlefish::edge_encoding_t exit_edge(const Kmer<k + 1>& e) const;

    // Returns the `DNA::Extended_Base` encoding of the edge `e` corresponding to this endpoint,
    // given the endpoint is the sink endpoint of the edge.
    cuttlefish::edge_encoding_t entrance_edge(const Kmer<k + 1>& e) const;


public:

    // Constructs an empty endpoint.
    Endpoint()
    {}

    // Configures the endpoint with the source (i.e. prefix) k-mer of the edge (k + 1)-mer `e`;
    // and uses the hash table `hash` to get the hash value of the vertex.
    void from_prefix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Configures the endpoint with the sink (i.e. suffix) k-mer of the edge (k + 1)-mer `e`;
    // and uses the hash table `hash` to get the hash value of the vertex.
    void from_suffix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Returns the neighboring endpoint of this endpoint that's connected with an edge encoded
    // with the code `e`, from the point-of-view of this endpoint. Uses the hash table `hash`
    // to get the hash value of the corresponding neighbor vertex.
    Endpoint neighbor_endpoint(cuttlefish::edge_encoding_t e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);
};


template <uint16_t k>
inline Endpoint<k>::Endpoint(const Kmer<k>& kmer, const bool is_source, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash):
    kmer_(kmer)
{
    kmer_bar_.as_reverse_complement(kmer_);
    kmer_hat_ptr = Kmer<k>::canonical(kmer_, kmer_bar_);

    s = (is_source ? exit_side() : entrance_side());

    h = hash(*kmer_hat_ptr);
}


template <uint16_t k>
inline void Endpoint<k>::from_prefix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
    kmer_.from_prefix(e);
    kmer_bar_.as_reverse_complement(kmer_);
    kmer_hat_ptr = Kmer<k>::canonical(kmer_, kmer_bar_);
    
    s = exit_side();
    this->e = exit_edge(e);

    h = hash(*kmer_hat_ptr);
}


template <uint16_t k>
inline void Endpoint<k>::from_suffix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
    kmer_.from_suffix(e);
    kmer_bar_.as_reverse_complement(kmer_);
    kmer_hat_ptr = Kmer<k>::canonical(kmer_, kmer_bar_);

    s = entrance_side();
    this->e = entrance_edge(e);

    h = hash(*kmer_hat_ptr);
}


template <uint16_t k>
inline cuttlefish::side_t Endpoint<k>::exit_side() const
{
    return &kmer_ == kmer_hat_ptr ? cuttlefish::side_t::back : cuttlefish::side_t::front;
}


template <uint16_t k>
inline cuttlefish::side_t Endpoint<k>::entrance_side() const
{
    return &kmer_ == kmer_hat_ptr ? cuttlefish::side_t::front : cuttlefish::side_t::back;
}


template <uint16_t k>
inline cuttlefish::edge_encoding_t Endpoint<k>::exit_edge(const Kmer<k + 1>& e) const
{
    return DNA_Utility::map_extended_base(s == cuttlefish::side_t::back ?
                                            e.back() : DNA_Utility::complement(e.back()));
}


template <uint16_t k>
inline cuttlefish::edge_encoding_t Endpoint<k>::entrance_edge(const Kmer<k + 1>& e) const
{
    return DNA_Utility::map_extended_base(s == cuttlefish::side_t::front ?
                                            e.front() : DNA_Utility::complement(e.front()));
}


template <uint16_t k>
inline Endpoint<k> Endpoint<k>::neighbor_endpoint(const cuttlefish::edge_encoding_t e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
    Kmer<k> kmer(*kmer_hat_ptr);

    if(s == cuttlefish::side_t::back)
    {
        kmer.roll_forward(e);
        return Endpoint<k>(kmer, false, hash);
    }
    
    kmer.roll_backward(e);
    return Endpoint<k>(kmer, true, hash);
}



#endif

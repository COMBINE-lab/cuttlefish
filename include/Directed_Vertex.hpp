
#ifndef DIRECTED_VERTEX_HPP
#define DIRECTED_VERTEX_HPP



#include "Kmer.hpp"
#include "globals.hpp"
#include "Kmer_Hash_Table.hpp"

#include <cstdint>


// A class denoting an instance of a vertex. It's "directed" in the sense that the k-mer
// observed for the vertex is in a particular orientation — although a vertex `v` has an
// unambiguous canonical k-mer `v_hat`, the vertex can be observed in two different k-mer
// forms: `v_hat` and `{v_hat}_bar` — the class keeps track of the particular k-mer form
// observed for the vertex instance.
template <uint16_t k>
class Directed_Vertex
{
private:

    Kmer<k> kmer_;  // The observed k-mer for the vertex.
    Kmer<k> kmer_bar_;  // Reverse complement of the k-mer observed for the vertex.
    const Kmer<k>* kmer_hat_ptr;    // Pointer to the canonical form of the k-mer associated to the vertex.
    uint64_t h; // Hash value of the vertex, i.e. hash of the canonical k-mer.

    // Initialize the data of the class once the observed k-mer `kmer_` is set.
    void init(const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);


public:

    // Constructs an empty vertex.
    Directed_Vertex()
    {}

    // Constructs a vertex observed for the k-mer `kmer`. Gets the hash value of the vertex using
    // the hash table `hash`.
    Directed_Vertex(const Kmer<k>& kmer, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Copy constructs the vertex from `rhs`.
    Directed_Vertex(const Directed_Vertex<k>& rhs);

    // Assigns the vertex `rhs` to this one, and returns a constant reference to this object.
    const Directed_Vertex<k>& operator=(const Directed_Vertex<k>& rhs);

    // Returns `true` iff the k-mer observed for the vertex is in its canonical form.
    bool in_canonical_form() const;

    // Configures the vertex with the source (i.e. prefix) k-mer of the edge (k + 1)-mer `e`;
    // and uses the hash table `hash` to get the hash value of the vertex.
    void from_prefix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Configures the vertex with the sink (i.e. suffix) k-mer of the edge (k + 1)-mer `e`;
    // and uses the hash table `hash` to get the hash value of the vertex.
    void from_suffix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Returns the observed k-mer for the vertex.
    const Kmer<k>& kmer() const;

    // Returns the reverse complement of the observed k-mer for the vertex.
    const Kmer<k>& kmer_bar() const;

    // Returns the canonical form of the vertex.
    const Kmer<k>& canonical() const;

    // Returns the hash value of the vertex.
    uint64_t hash() const;

    // Transforms this vertex to another by chopping off the first base from the associated
    // observed k-mer, and appending the nucleobase `b` to the end, i.e. effecitively
    // rolling the associated k-mer by one base "forward". The hash table `hash` is used
    // to get the hash value of the new vertex.
    void roll_forward(cuttlefish::base_t b, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Returns the side of the vertex which is to be the incidence side of some bidirected
    // edge instance if this vertex instance were to be the source vertex (i.e. prefix k-mer)
    // of that edge.
    cuttlefish::side_t exit_side() const;

    // Returns `true` iff this vertex and the vertex `v` are the same vertex, without the
    // directionality.
    bool is_same_vertex(const Directed_Vertex<k>& v) const;
};


template <uint16_t k>
inline void Directed_Vertex<k>::init(const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
    kmer_bar_.as_reverse_complement(kmer_);
    kmer_hat_ptr = Kmer<k>::canonical(kmer_, kmer_bar_);

    h = hash(*kmer_hat_ptr);
}


template <uint16_t k>
inline Directed_Vertex<k>::Directed_Vertex(const Kmer<k>& kmer, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash):
    kmer_(kmer)
{
    init(hash);
}


template <uint16_t k>
inline Directed_Vertex<k>::Directed_Vertex(const Directed_Vertex<k>& rhs):
    kmer_(rhs.kmer_),
    kmer_bar_(rhs.kmer_bar_),
    kmer_hat_ptr(Kmer<k>::canonical(kmer_, kmer_bar_)), // TODO: replace with pointer-check based assignment (check `operator=`).
    h(rhs.h)
{}


template <uint16_t k>
inline const Directed_Vertex<k>& Directed_Vertex<k>::operator=(const Directed_Vertex<k>& rhs)
{
    kmer_ = rhs.kmer_;
    kmer_bar_ = rhs.kmer_bar_;
    kmer_hat_ptr = (rhs.kmer_hat_ptr == &rhs.kmer_ ? &kmer_ : &kmer_bar_);
    h = rhs.h;

    return *this;
}


template <uint16_t k>
inline bool Directed_Vertex<k>::in_canonical_form() const
{
    return &kmer_ == kmer_hat_ptr;
}


template <uint16_t k>
inline void Directed_Vertex<k>::from_prefix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
    kmer_.from_prefix(e);
    init(hash);
}


template <uint16_t k>
inline void Directed_Vertex<k>::from_suffix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
    kmer_.from_suffix(e);
    init(hash);
}


template <uint16_t k>
inline const Kmer<k>& Directed_Vertex<k>::kmer() const
{
    return kmer_;
}


template <uint16_t k>
inline const Kmer<k>& Directed_Vertex<k>::kmer_bar() const
{
    return kmer_bar_;
}


template <uint16_t k>
inline const Kmer<k>& Directed_Vertex<k>::canonical() const
{
    return *kmer_hat_ptr;
}


template <uint16_t k>
inline uint64_t Directed_Vertex<k>::hash() const
{
    return h;
}


template <uint16_t k>
inline void Directed_Vertex<k>::roll_forward(const cuttlefish::base_t b, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
    kmer_.roll_to_next_kmer(b, kmer_bar_);
    kmer_hat_ptr = Kmer<k>::canonical(kmer_, kmer_bar_);

    h = hash(*kmer_hat_ptr);
}


template <uint16_t k>
inline cuttlefish::side_t Directed_Vertex<k>::exit_side() const
{
    return &kmer_ == kmer_hat_ptr ? cuttlefish::side_t::back : cuttlefish::side_t::front;
}


template <uint16_t k>
inline bool Directed_Vertex<k>::is_same_vertex(const Directed_Vertex<k>& v) const
{
    return hash() == v.hash();
}



#endif

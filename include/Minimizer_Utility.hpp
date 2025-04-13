
#ifndef MINIMIZER_UTILITY_HPP
#define MINIMIZER_UTILITY_HPP



#include "Kmer.hpp"
#include "globals.hpp"
#include "wyhash/wyhash.h"

#include <cstdint>
#include <cstddef>


template <typename T_seq_, uint16_t k, bool is_canonical_> class Minimizer_Iterator;


// =============================================================================
// A class containing various utility methods for k-mer minimizers.
class Minimizer_Utility
{

private:

    // Salt for `wyhash`.
    static constexpr uint64_t wy_salt[4] = {4167021922371662411llu, 7320285940802167691llu, 14307255741305819987llu, 10859488101230029397llu};

public:

    // Returns the hash value of the l-mer `lmer`. The seed-value `seed` is used
    // in hashing.
    static uint64_t hash(cuttlefish::minimizer_t lmer, uint64_t seed = 0);

    // Extracts the minimizer of the k-mer `kmer` into `min`, and its index into
    // `idx`.
    template <uint16_t k>
    static void get_minimizer(const Kmer<k>& kmer, const uint16_t l, cuttlefish::minimizer_t& min, std::size_t& idx);

    // Extracts the minimizer of the k-mer `kmer` into `min`.
    template <uint16_t k>
    static cuttlefish::minimizer_t canonical_minimizer(const Kmer<k>& kmer, const uint16_t l);
};


// A class to pack information regarding l-mers, to aid in computing
// l-minimizers of k-mers over sequences by `Minimizer_Iterator`.
class Lmer_Tuple
{
    friend class Minimizer_Utility;
    template <typename T_seq_, uint16_t k, bool is_canonical_> friend class Minimizer_Iterator;

    typedef cuttlefish::minimizer_t minimizer_t;

private:

    minimizer_t lmer;   // The l-mer.
    std::size_t index;  // Its index in the underlying sequence.
    uint64_t hash;      // Its hash value, determining the l-mer ordering.


public:

    Lmer_Tuple()
    {};

    // Constructs a tuple for an l-mer `lmer`, positioned at index `index` of
    // the underlying sequence, and having a hash value `hash`.
    Lmer_Tuple(minimizer_t lmer, std::size_t index, uint64_t hash);

    // Returns `true` iff this l-mer tuple is to be ordered as lesser to the
    // tuple `rhs`. The primary ordering is based on the tuples' hashes. If
    // equal, then the ordering is based on their literal form. If also equal,
    // then the tuple occurring earlier in the underlying sequence is lesser.
    bool operator<(const Lmer_Tuple& rhs) const;
};


inline Lmer_Tuple::Lmer_Tuple(const minimizer_t lmer, const std::size_t index, const uint64_t hash):
    lmer(lmer),
    index(index),
    hash(hash)
{}


inline bool Lmer_Tuple::operator<(const Lmer_Tuple& rhs) const
{
    if(hash != rhs.hash)
        return hash < rhs.hash;

    if(lmer != rhs.lmer)
        return lmer < rhs.lmer;

    return index < rhs.index;
}


inline uint64_t Minimizer_Utility::hash(const cuttlefish::minimizer_t lmer, const uint64_t seed)
{
#ifdef CF_DEVELOP_MODE
    return lmer;    // TODO: add as debug option for developer.
#endif

    // return XXH3_64bits_withSeed(&lmer, sizeof(lmer), seed);
    // TODO: the entropy of the minimizer-hash is getting diluted by considering all the leading 0 bits (64 - 2l bits, specifically.)
    // return wyhash(&lmer, (2 * l + 7) / 8, seed, wy_salt);
    return wyhash(&lmer, sizeof(lmer), seed, wy_salt);
}


template <uint16_t K>
inline void Minimizer_Utility::get_minimizer(const Kmer<K>& kmer, const uint16_t l, cuttlefish::minimizer_t& min, std::size_t& idx)
{
    const uint64_t* const kmer_data = kmer.data();

    const cuttlefish::minimizer_t last_lmer = kmer_data[0] & ((uint64_t(0b1) << (2 * l)) - 1);  // The last l-mer in `kmer`.
    Lmer_Tuple curr_lmer(last_lmer, K - l, hash(last_lmer));
    Lmer_Tuple min_lmer(curr_lmer);

    uint64_t base;
    for(std::size_t i = l; i < K; ++i)  // Backward scan of l-mers in `kmer`. Yes, backwards.
    {
        base = (kmer_data[i >> 5] >> (2 * (i & 31))) & uint64_t(0b11);
        curr_lmer.lmer = (curr_lmer.lmer >> 2) | (base << (2 * (l - 1)));
        curr_lmer.index--;
        curr_lmer.hash = hash(curr_lmer.lmer);

        if(curr_lmer < min_lmer)
            min_lmer = curr_lmer;
    }

    min = min_lmer.lmer, idx = min_lmer.index;
}


template <uint16_t k>
inline cuttlefish::minimizer_t Minimizer_Utility::canonical_minimizer(const Kmer<k>& kmer, const uint16_t l)
{
    cuttlefish::minimizer_t min_f, min_r;
    std::size_t min_idx;

    get_minimizer(kmer, l, min_f, min_idx);
    get_minimizer(kmer.reverse_complement(), l, min_r, min_idx);

    const auto h_f = hash(min_f), h_r = hash(min_r);
    return h_f != h_r ? (h_f < h_r ? min_f : min_r) : std::min(min_f, min_r);
}



#endif

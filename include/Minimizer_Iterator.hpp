
#ifndef MINIMIZER_ITERATOR_HPP
#define MINIMIZER_ITERATOR_HPP



#include "Minimizer_Utility.hpp"
#include "DNA_Utility.hpp"
#include "Kmer.hpp"
#include "globals.hpp"
#include "Fixed_Cap_Deque.hpp"
#include "xxHash/xxh3.h"

#include <cstdint>
#include <cassert>


// =============================================================================
// A class to iterate over the minimizers of the constituent `k`-mers of a given
// sequence of type `T_seq_`, computing the minimizer for each `k`-mer in
// amortized O(1) time. If `is_canonical_` is `true`, then canonical minimizers
// are computed, i.e. both strand-forms of the sequence is considered.
template <typename T_seq_, uint16_t k, bool is_canonical_ = false>
class Minimizer_Iterator
{
    typedef cuttlefish::minimizer_t minimizer_t;

private:

    T_seq_ seq; // The sequence on which to iterate over.
    std::size_t seq_len;    // Length of the given sequence.

    const uint16_t l;   // Size of the minimizers.

    const uint64_t seed;    // Seed for hashing l-mers.

    minimizer_t last_lmer;  // The last l-mer processed.
    minimizer_t last_lmer_bar;  // Reverse complement of the last l-mer processed.
    std::size_t last_lmer_idx;  // Index into the sequence of the last l-mer processed.
    const minimizer_t clear_MSN_mask;   // Bitmask to clear the most-significant nucleotide bits of l-mers.

    typedef cuttlefish::Fixed_Cap_Deque<Lmer_Tuple, 2 * k> deque_t;

    // Collection of l-mers that have already been seen, and cannot be ruled out
    // yet as candidate minimizers for the k-mers yet to be seen fully. `dq_f`
    // is for the forward-strand and `dq_r` is for the reverse-strand.
    deque_t dq_f, dq_r;

public:

    // Constructs a minimizer iterator to iterate over `l`-minimizers of
    // the `k`-mers of a given sequence in a streaming manner. The seed-value
    // `seed` is used in hashing the `l`-mers.
    Minimizer_Iterator(uint16_t l, uint64_t seed = 0);

    // Constructs a minimizer iterator to iterate over `l`-minimizers of
    // the `k`-mers of the sequence `seq`, of length `seq_len`, in a streaming
    // manner. The iterator sits at the first k-mer after the construction.
    // The seed-value `seed` is used in hashing the `l`-mers.
    Minimizer_Iterator(T_seq_ seq, std::size_t seq_len, uint16_t l, uint64_t seed = 0);

    // Resets the iterator to the sequence `seq` of length `seq_len`. The
    // iterator sits at the first k-mer after the construction.
    void reset(T_seq_ seq, std::size_t seq_len);

    // Moves the iterator to a next k-mer by extending it with the character
    // `ch`.
    void advance(char ch);

    // Moves the iterator to the next k-mer in the sequence. Returns `true` iff
    // the current k-mer is not the last k-mer in the sequence.
    bool operator++();

    // Puts the minimizer of the current k-mer into `minimizer`, and stores its
    // index in the sequence at `index`. The position information over the
    // sequence is maintained internally.
    void value_at(minimizer_t& minimizer, std::size_t& index) const;

    // Puts the minimizer of the current k-mer into `minimizer`, stores its
    // index in the sequence at `index`, and stores its 64-bit hash at `h`. The
    // position information over the sequence is maintained internally.
    void value_at(minimizer_t& minimizer, std::size_t& index, uint64_t& h) const;

    // Computes the `l`-minimizer of the `k`-length sequence `kmer` at `min`,
    // its hash at `h`, and index in `kmer` ar `idx`. The seed value `seed` is
    // used in hashing the `l`-mers.
    static void minimizer(T_seq_ kmer, uint16_t l, uint64_t seed, minimizer_t& min, uint64_t& h, std::size_t& idx);
};


template <typename T_seq_,  uint16_t k, bool is_canonical_>
inline Minimizer_Iterator<T_seq_, k, is_canonical_>::Minimizer_Iterator(const uint16_t l, const uint64_t seed):
      l(l)
    , seed(seed)
    , clear_MSN_mask(~(static_cast<minimizer_t>(0b11) << (2 * (l - 1))))
{
    assert(l <= k);
}


template <typename T_seq_, uint16_t k, bool is_canonical_>
inline Minimizer_Iterator<T_seq_, k, is_canonical_>::Minimizer_Iterator(T_seq_ const seq, const std::size_t seq_len, const uint16_t l, const uint64_t seed):
    Minimizer_Iterator(l, seed)
{
    reset(seq, seq_len);
}


template <typename T_seq_, uint16_t k, bool is_canonical_>
inline void Minimizer_Iterator<T_seq_, k, is_canonical_>::reset(T_seq_ const seq, const std::size_t seq_len)
{
    this->seq = seq;
    this->seq_len = seq_len;
    dq_f.clear(), dq_r.clear();

    last_lmer = 0;
    last_lmer_bar = 0;
    last_lmer_idx = 0;

    DNA::Base base;
    for(std::size_t idx = 0; idx < l; ++idx)
    {
        base = DNA_Utility::map_base(seq[idx]);
        assert(base != DNA::Base::N);
        last_lmer |= (static_cast<minimizer_t>(base) << (2 * (l - 1 - idx)));

        if constexpr(is_canonical_)
            last_lmer_bar |= (static_cast<minimizer_t>(DNA_Utility::complement(base)) << (2 * idx));
    }

    dq_f.emplace_back(last_lmer, last_lmer_idx, Minimizer_Utility::hash(last_lmer, seed));
    if constexpr(is_canonical_)
        dq_r.emplace_back(last_lmer_bar, last_lmer_idx, Minimizer_Utility::hash(last_lmer_bar, seed));

    while(last_lmer_idx + (l - 1) < static_cast<std::size_t>(k - 1))
        operator++();
}


template <typename T_seq_, uint16_t k, bool is_canonical_>
inline void Minimizer_Iterator<T_seq_, k, is_canonical_>::advance(const char ch)
{
    assert(DNA_Utility::is_DNA_base(ch));

    last_lmer_idx++;
    const DNA::Base base = DNA_Utility::map_base(ch);
    last_lmer = ((last_lmer & clear_MSN_mask) << 2) | base;
    if constexpr(is_canonical_)
        last_lmer_bar = (last_lmer_bar >> 2) | (static_cast<minimizer_t>(DNA_Utility::complement(base)) << (2 * (l - 1)));


    // End of the current k-mer, e = last_lmer_idx + l - 1; So start of the current k-mer: e - (k - 1).
    const std::size_t curr_kmer_idx = (last_lmer_idx + l - 1 >= k - 1lu ? last_lmer_idx + l - 1 - (k - 1) : 0);
    if(dq_f.front().index < curr_kmer_idx)    // This candidate l-mer falls out of the current k-mer.
        dq_f.pop_front();

    if constexpr(is_canonical_)
        if(dq_r.front().index < curr_kmer_idx)  // This candidate l-mer falls out of the current k-mer.
            dq_r.pop_front();


    const auto fix_dq =
        [](deque_t& dq, const Lmer_Tuple& lmer_tuple)
        {
            while(!dq.empty())
                if(dq.back() < lmer_tuple)
                    break;
                else
                    dq.pop_back();

            dq.push_back(lmer_tuple);
        };


    fix_dq(dq_f, Lmer_Tuple(last_lmer, last_lmer_idx, Minimizer_Utility::hash(last_lmer, seed)));
    if constexpr(is_canonical_)
        fix_dq(dq_r, Lmer_Tuple(last_lmer_bar, last_lmer_idx, Minimizer_Utility::hash(last_lmer_bar, seed)));
}


template <typename T_seq_, uint16_t k, bool is_canonical_>
inline bool Minimizer_Iterator<T_seq_, k, is_canonical_>::operator++()
{
    if(last_lmer_idx + l == seq_len)
        return false;

    advance(seq[last_lmer_idx + l]);
    return true;
}


template <typename T_seq_, uint16_t k, bool is_canonical_>
inline void Minimizer_Iterator<T_seq_, k, is_canonical_>::value_at(cuttlefish::minimizer_t& minimizer, std::size_t& index) const
{
    if constexpr(is_canonical_)
    {
        if(dq_f.front() < dq_r.front())
            minimizer = dq_f.front().lmer, index = dq_f.front().index;
        else
            minimizer = dq_r.front().lmer, index = dq_r.front().index;
    }
    else
        minimizer = dq_f.front().lmer, index = dq_f.front().index;
}


template <typename T_seq_, uint16_t k, bool is_canonical_>
inline void Minimizer_Iterator<T_seq_, k, is_canonical_>::value_at(cuttlefish::minimizer_t& minimizer, std::size_t& index, uint64_t& h) const
{
    if constexpr(is_canonical_)
    {
        if(dq_f.front() < dq_r.front())
            minimizer = dq_f.front().lmer, index = dq_f.front().index, h = dq_f.front().hash;
        else
            minimizer = dq_r.front().lmer, index = dq_r.front().index, h = dq_r.front().hash;
    }
    else
        minimizer = dq_f.front().lmer, index = dq_f.front().index, h = dq_f.front().hash;
}


template <typename T_seq_, uint16_t k, bool is_canonical_>
inline void Minimizer_Iterator<T_seq_, k, is_canonical_>::minimizer(T_seq_ kmer, uint16_t l, uint64_t seed, minimizer_t& min, uint64_t& h, std::size_t& idx)
{
    Minimizer_Iterator(kmer, k, l, seed).value_at(min, idx);
    h = Minimizer_Utility::hash(min);
}



// A class to iterate over the minimizers of the constituent `k`-mers of a given
// sequence, computing the minimizer for each `k`-mer in amortized O(1) time.
// Canonical minimizers are computed, i.e. both strand-forms of the sequence is
// considered.
template <uint16_t k>
class Min_Iterator
{
private:

    uint16_t l; // Minimizer size.

    uint64_t last_lmer; // The last l-mer processed from the sequence.
    uint64_t last_lmer_bar; // Reverse complement of the last l-mer processed.
    const uint64_t clear_MSN_mask;  // Bitmask to clear the most-significant nucleotide bits of l-mers.

    // The hashes of the l-mers in a k-mer window. The window sits at `H[0..k - l]` (in reversed order relative to
    // the l-mers' order for implementation concerns), and is conceptually broken into two sub-windows: a prefix and
    // a suffix one.
    uint64_t H[k];

    uint64_t M_pre[k];  // Suffix-mins of the prefix window.
    uint64_t M_suf; // Min of the suffix window.

    std::size_t pivot;  // Pivot-index, inclusive to the prefix window, breaking the window into the two sub-windows.


    // Returns the minimum of `x` and `y`.
    static constexpr auto min = [](const uint64_t x, const uint64_t y){ return x < y ? x : y; };

    // Returns the 64-bit equivalent of `val`.
    template <typename T_> static constexpr uint64_t as_u64(const T_ val) { return static_cast<uint64_t>(val); }

    // Returns the 64-bit hash value of the canonically equivalent l-mer pair
    // `(lmer, lmer_bar)`.
    static uint64_t lmer_hash(uint64_t lmer, uint64_t lmer_bar);

    // Moves the suffix window to the prefix window, and empties it—resetting
    // them to their default condition. Also recomputes the aggregate minimum
    // results of the windows.
    void reset_windows();


public:

    // Constructs a minimizer iterator to iterate over `l`-minimizers of
    // the `k`-mers of a given sequence in a streaming manner.
    Min_Iterator(uint16_t l);

    // Resets the iterator to the sequence `seq`. The iterator sits at the
    // first `k`-mer after the construction.
    void reset(const char* seq);

    // Returns the minimizer's 64-bit hash value.
    uint64_t hash() const;

    // Moves the iterator to a next k-mer extending it with the character `ch`.
    void advance(char ch);
};


template <uint16_t k>
inline uint64_t Min_Iterator<k>::lmer_hash(const uint64_t lmer, const uint64_t lmer_bar)
{
    constexpr uint64_t min_seed = 0;
    return min(Minimizer_Utility::hash(lmer, min_seed), Minimizer_Utility::hash(lmer_bar, min_seed));
}


template <uint16_t k>
inline Min_Iterator<k>::Min_Iterator(const uint16_t l):
      l(l)
    , clear_MSN_mask(~(as_u64(0b11) << (2 * (l - 1))))
{}


template <uint16_t k>
inline void Min_Iterator<k>::reset(const char* const seq)
{
    last_lmer = 0;
    last_lmer_bar = 0;

    for(std::size_t idx = 0; idx < l; ++idx)
    {
        const DNA::Base base = DNA_Utility::map_base(seq[idx]);
        assert(base != DNA::Base::N);
        last_lmer |= (as_u64(base) << (2 * (l - 1 - idx)));
        last_lmer_bar |= (as_u64(DNA_Utility::complement(base)) << (2 * idx));
    }

    pivot = k - l;
    H[pivot] = lmer_hash(last_lmer, last_lmer_bar);

    for(std::size_t idx = l; idx < k; ++idx)
    {
        const DNA::Base base = DNA_Utility::map_base(seq[idx]);
        assert(base != DNA::Base::N);
        last_lmer = ((last_lmer & clear_MSN_mask) << 2) | base;
        last_lmer_bar = (last_lmer_bar >> 2) | (as_u64(DNA_Utility::complement(base)) << (2 * (l - 1)));

        H[--pivot] = lmer_hash(last_lmer, last_lmer_bar);
    }


    reset_windows();
}


template <uint16_t k>
inline void Min_Iterator<k>::reset_windows()
{
    assert(pivot == 0);

    // TODO: try vectorization.

    const std::size_t ring_sz = k - l;
    M_pre[0] = H[0];
    for(std::size_t i = 1; i <= ring_sz; ++i)
        M_pre[i] = min(M_pre[i - 1], H[i]);

    M_suf = std::numeric_limits<uint64_t>::max();
    pivot = k - l;
}


template <uint16_t k>
inline uint64_t Min_Iterator<k>::hash() const
{
    return min(M_pre[pivot], M_suf);
}


template <uint16_t k>
inline void Min_Iterator<k>::advance(const char ch)
{
    assert(DNA_Utility::is_DNA_base(ch));

    const DNA::Base base = DNA_Utility::map_base(ch);
    last_lmer = ((last_lmer & clear_MSN_mask) << 2) | base;
    last_lmer_bar = (last_lmer_bar >> 2) | (as_u64(DNA_Utility::complement(base)) << (2 * (l - 1)));

    H[pivot] = lmer_hash(last_lmer, last_lmer_bar);
    M_suf = min(M_suf, H[pivot]);

    if(pivot > 0)
        pivot--;
    else
        reset_windows();
}



#endif


#ifndef KMER_HPP
#define KMER_HPP



#include "DNA_Utility.hpp"
#include "kmc_api/kmc_file.h"
#include "xxHash/xxh3.h"

#include <cstring>
#include <string>
#include <algorithm>


// Defining this macro states our intent that only odd k-values will be used for de Bruijn graph vertices.
// Hence, extraction of k-mers from (k + 1)-mers — vertices from edges — will only happen when k is odd.
#define ODD_K


template <uint16_t k>
class Kmer: public DNA_Utility
{
    // Make k-mers friend for (k + 1)-mer, so that de Bruijn graph vertices, i.e. k-mers,
    // may access private information (the raw data) from edges, i.e. (k + 1)-mers.
    friend class Kmer<k - 1>;

private:

    // Number of 64-bit integers required to compactly represent the underlying k-mer with 2-bits/base encoding.
    static constexpr const uint16_t NUM_INTS = (k + 31) / 32;

    // Bitmask used to clear the most significant DNA base character, i.e. the first base of the k-mer which is at the bits `2k-1 : 2k-2`.
    static constexpr const uint64_t CLEAR_MSN_MASK = ~(uint64_t(0b11) << (2 * ((k - 1) % 32)));

    // The underlying k-mer represented with 2-bit encoding, as a collection of 64-bit integers.
    // A k-mer `n_{k - 1} ... n_1 n_0` is stored in the array `kmer_data` such that, `kmer_data[0]`
    // stores the suffix `n_63 ... n_0`, then `kmer_data[1]` stores `n_127 ... n_64`, and so on.
    // That is, the suffix is aligned with a byte boundary.
    uint64_t kmer_data[NUM_INTS];


    // Left-shifts the collection of the bits at the `kmer_data` array by one base (2-bits).
    void left_shift();

    // Right-shifts the collection of the bits at the `kmer_data` array by one base (2-bits).
    void right_shift();

    // Left-shifts the collection of the bits at the `kmer_data` array by `B` bases (2B-bits).
    template <uint16_t B> void left_shift(char(*)[B != 0] = 0);
    template <uint16_t B> void left_shift(char(*)[B == 0] = 0);


public:

    // Default constructs the k-mer with a 0-value, equivalent to "AA...A".
    Kmer();

    // Constructs a k-mer from the provided characters at
    // `label[kmer_idx,...,kmer_idx + k - 1]`.
    Kmer(const char* label, size_t kmer_idx);

    // Constructs a k-mer from the provided string `label`.
    Kmer(const std::string& label);

    // Constructs a k-mer from the provided characters at
    // `label[kmer_idx,...,kmer_idx + k - 1]`.
    Kmer(const std::string& label, size_t kmer_idx);

    // Constructs a k-mer from `kmer_api`, a k-mer object built from KMC.
    Kmer(const CKmerAPI& kmer_api);

    // Copy constructs the k-mer from another k-mer `rhs`.
    Kmer(const Kmer<k>& rhs);

    // Copy assignment operator.
    Kmer<k>& operator=(const Kmer<k>& rhs);

    // Sets the value of the `k` parameter across the `Kmer` class.
    // Note: For this general k-mer class with `k` being fixed at compile-time,
    // this method effectively does nothing.
    static void set_k(uint16_t kmer_len);

    // Returns the k-parameter.
    constexpr static uint16_t get_k();

    // Returns a 64-bit hash value for the k-mer.
    uint64_t to_u64(uint64_t seed=0) const;

    // Gets the k-mer from the KMC api object `kmer_api`.
    void from_CKmerAPI(const CKmerAPI& kmer_api);

    // Gets the k-mer from its KMC raw-binary representation.
    void from_KMC_data(const uint64_t* kmc_data);

    // Gets the k-mer that is a prefix of the provided
    // (k + 1)-mer `k_plus_1_mer`.
    void from_prefix(const Kmer<k + 1>& k_plus_1_mer);

    // Gets the k-mer that is a suffix of the provided
    // (k + 1)-mer `k_plus_1_mer`.
    void from_suffix(const Kmer<k + 1>& k_plus_1_mer);

    // Returns the reverese complement of the k-mer.
    Kmer<k> reverse_complement() const;

    // Returns true iff the bitwise encoding of this k-mer is lesser to the
    // encoding of the other k-mer `rhs`.
    bool operator<(const Kmer<k>& rhs) const;

    // Returns true iff this k-mer is identical to the other k-mer `rhs`.
    bool operator==(const Kmer<k>& rhs) const;

    // Returns true iff this k-mer is not identical to the other k-mer `rhs`.
    bool operator!=(const Kmer<k>& rhs) const;

    // Returns `true` iff the k-mer is in the forward direction relative to
    // the other k-mer `kmer_hat`.
    bool in_forward(const Kmer<k>& kmer_hat) const;

    // Transforms this k-mer by chopping off the first base and
    // appending the next base `next_base` to the end, i.e.
    // rolls the k-mer by one base. Also sets the passed reverse
    // complement `rev_compl` of the k-mer accordingly.
    void roll_to_next_kmer(char next_base, Kmer<k>& rev_compl);

    // Returns the canonical version of the k-mer, comparing it to its
    // reverse complement `rev_compl`.
    Kmer canonical(const Kmer<k>& rev_compl) const;

    // Returns the canonical version of the k-mer.
    Kmer canonical() const;

    // Returns the string label of the k-mer.
    std::string string_label() const;

    // Prints the literal representation of the K-mer `kmer` to the
    // stream `ostream`.
    template <uint16_t K>
    friend std::ostream& operator<<(std::ostream& out, const Kmer<K>& kmer);
};


template <uint16_t k>
inline void Kmer<k>::left_shift()
{
    constexpr uint64_t mask_MSN = (uint64_t(0b11) << 62);

    for(uint16_t idx = NUM_INTS - 1; idx > 0; --idx)
        kmer_data[idx] = (kmer_data[idx] << 2) | ((kmer_data[idx - 1] & mask_MSN) >> 62);

    kmer_data[0] <<= 2;
}


template <uint16_t k>
inline void Kmer<k>::right_shift()
{
    constexpr uint64_t mask_LSN = uint64_t(0b11);

    for(uint16_t idx = 0; idx < NUM_INTS - 1; ++idx)
        kmer_data[idx] = (kmer_data[idx] >> 2) | ((kmer_data[idx + 1] & mask_LSN) << 62);

    kmer_data[NUM_INTS - 1] >>= 2;
}


template <uint16_t k>
template <uint16_t B>
inline void Kmer<k>::left_shift(char(*)[B != 0])
{
    static_assert(B < 32, "invalid bit-shift amount");

    constexpr uint16_t num_bit_shift = 2 * B;
    constexpr uint64_t mask_MSNs = ((static_cast<uint64_t>(1) << num_bit_shift) - 1) << (64 - num_bit_shift);

    for(uint16_t idx = NUM_INTS - 1; idx > 0; --idx)
        kmer_data[idx] = (kmer_data[idx] << num_bit_shift) | ((kmer_data[idx - 1] & mask_MSNs) >> (64 - num_bit_shift));

    kmer_data[0] <<= num_bit_shift;
}


template <uint16_t k>
template <uint16_t B>
inline void Kmer<k>::left_shift(char(*)[B == 0])
{}


template <uint16_t k>
inline uint64_t Kmer<k>::to_u64(uint64_t seed) const
{
    constexpr const uint16_t NUM_BYTES = (k + 3) / 4;
    return XXH3_64bits_withSeed(kmer_data, NUM_BYTES, seed);
}


template <uint16_t k>
inline Kmer<k>::Kmer():
    kmer_data() // Value-initializes the data array, i.e. zeroes it out.
{}


template <uint16_t k>
inline Kmer<k>::Kmer(const char* const label, const size_t kmer_idx):
    Kmer()
{
    for(size_t idx = kmer_idx; idx < kmer_idx + k; ++idx)
    {
        const DNA::Base base = map_base(label[idx]);

        left_shift();
        kmer_data[0] |= base;
    }
}


template <uint16_t k>
inline Kmer<k>::Kmer(const std::string& label):
    Kmer(label.c_str(), 0)
{}


template <uint16_t k>
inline Kmer<k>::Kmer(const std::string& label, const size_t kmer_idx):
    Kmer(label.c_str(), kmer_idx)
{}


template <uint16_t k>
inline Kmer<k>::Kmer(const CKmerAPI& kmer_api)
{
    from_CKmerAPI(kmer_api);
}


template <uint16_t k>
inline Kmer<k>::Kmer(const Kmer<k>& rhs)
{
    memcpy(kmer_data, rhs.kmer_data, NUM_INTS * sizeof(uint64_t));
    //for(uint16_t idx = 0; idx < NUM_INTS; ++idx)
    //    kmer_data[idx] = rhs.kmer_data[idx];
}


template <uint16_t k>
inline Kmer<k>& Kmer<k>::operator=(const Kmer<k>& rhs)
{
    memcpy(kmer_data, rhs.kmer_data, NUM_INTS * sizeof(uint64_t));
    //for(uint16_t idx = 0; idx < NUM_INTS; ++idx)
    //    kmer_data[idx] = rhs.kmer_data[idx];
    return *this;
}


template <uint16_t k>
inline void Kmer<k>::set_k(const uint16_t kmer_len)
{
    if(kmer_len != k)
    {
        std::cerr << "Expected k = " << k << ", recieved " << kmer_len << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
inline constexpr uint16_t Kmer<k>::get_k()
{
    return k;
}


template <uint16_t k>
inline void Kmer<k>::from_CKmerAPI(const CKmerAPI& kmer_api)
{
    kmer_api.to_u64<NUM_INTS>(kmer_data);
}


template <uint16_t k>
inline void Kmer<k>::from_KMC_data(const uint64_t* const kmc_data)
{
    // The endianness of the k-mer data array in the KMC database is in the opposite
    // order of the one that Cuttlefish uses. So the fetch gathers it in reverse.
    // Roughly, the KMC database stores a literal k-mer `b_{k - 1} ... b_1 b_0` such
    // that the prefix aligns with a byte boundary, i.e. `b_{k - 1} ... b_{k - 1 - 63}`
    // is stored as one 64-bit collection at its `kmer_data[0]` entry, and the rest
    // follows this alignment. This is why, the 64-bits corresponding to the substring
    // `b_63 ... b_0` might be shared between the two maximum indices of `kmer_data`;
    // so can be all the next disjoint 64-bit substrings. The remainder substring
    // (of < 64-bits) can be found from just the 0'th entry of `kmer_data`.
    
    constexpr uint8_t byte_alignment = (k % 4 != 0 ? 4 - (k % 4) : 0);
    constexpr uint32 offset = 62 - ((k - 1 + byte_alignment) & 31) * 2;
    //std::memcpy(kmer_data, kmc_data, NUM_INTS * sizeof(kmer_data[0]));
    
    if(offset)
    {
        for(int32 i = NUM_INTS - 1; i >= 1; --i)
        {
            kmer_data[NUM_INTS - 1 - i] = kmc_data[i] >> offset;
            kmer_data[NUM_INTS - 1 - i] += kmc_data[i - 1] << (64 - offset);
        }

        kmer_data[NUM_INTS - 1] = kmc_data[0] >> offset;
    }
    else
        for (int32 i = NUM_INTS - 1; i >= 0; --i)
            kmer_data[NUM_INTS - 1 - i] = kmc_data[i];
            
}


template <uint16_t k>
inline void Kmer<k>::from_prefix(const Kmer<k + 1>& k_plus_1_mer)
{
    // Note: `Kmer<k>` and `Kmer<k + 1>` always have the same number of words (i.e. `NUM_INTS`) for odd k-values.
    // The only time that they have different numbers of words is when `k` is a multiple of 32. In such cases,
    // a (k + 1)-mer contains just one base in its highest index word, and a k-mer's words are fully packed.

    std::memcpy(kmer_data, k_plus_1_mer.kmer_data, NUM_INTS * sizeof(uint64_t));
    right_shift();  // Clear the LSN of the (k + 1)-mer, from the k-mer.

    #ifndef ODD_K   // The following `if` conditional can only be `true` when `k` is a multiple of 32.
    constexpr uint16_t kp1_NUM_INTS = ((k + 1) + 31) / 32;
    if(kp1_NUM_INTS != NUM_INTS)    // Fetch the only base not copied from the (k + 1)-mer as the MSN for this k-mer.
        kmer_data[NUM_INTS - 1] |= (k_plus_1_mer.kmer_data[kp1_NUM_INTS - 1] << 62);
    #endif
}


template <uint16_t k>
inline void Kmer<k>::from_suffix(const Kmer<k + 1>& k_plus_1_mer)
{
    std::memcpy(kmer_data, k_plus_1_mer.kmer_data, NUM_INTS * sizeof(uint64_t));

    #ifndef ODD_K   // The following `if` conditional can only be `true` when `k` is a multiple of 32.
    constexpr uint16_t kp1_NUM_INTS = ((k + 1) + 31) / 32;
    if(kp1_NUM_INTS != NUM_INTS)    // The only base not copied from the (k + 1)-mer isn't required to be fetched — it will be cleared out anyways.
        return;
    #endif

    kmer_data[NUM_INTS - 1] &= Kmer<k + 1>::CLEAR_MSN_MASK; // Clear the MSN of the (k + 1)-mer from this k-mer.
}


template <uint16_t k>
inline Kmer<k> Kmer<k>::reverse_complement() const
{
    Kmer<k> kmer(*this);
    Kmer<k> rev_compl;
    
    constexpr uint64_t mask_LSN = uint64_t(0b11);

    for(uint16_t idx = 0; idx < k; ++idx)
    {
        rev_compl.left_shift();
        rev_compl.kmer_data[0] |= complement(DNA::Base(kmer.kmer_data[0] & mask_LSN));

        kmer.right_shift();
    }


    return rev_compl;
}


template <uint16_t k>
inline bool Kmer<k>::operator<(const Kmer<k>& rhs) const
{
    for(int16_t idx = NUM_INTS - 1; idx >= 0; --idx)
        if(kmer_data[idx] != rhs.kmer_data[idx])
            return kmer_data[idx] < rhs.kmer_data[idx];

    return false;
}


template <uint16_t k>
inline bool Kmer<k>::operator==(const Kmer<k>& rhs) const
{
    for(uint16_t idx = 0; idx < NUM_INTS; ++idx)
        if(kmer_data[idx] != rhs.kmer_data[idx])
            return false;

    return true;
}


template <uint16_t k>
inline bool Kmer<k>::operator!=(const Kmer<k>& rhs) const
{
    return !operator==(rhs);
}


template <uint16_t k>
inline bool Kmer<k>::in_forward(const Kmer<k>& kmer_hat) const
{
    return this->operator==(kmer_hat);
}


template <uint16_t k>
inline void Kmer<k>::roll_to_next_kmer(const char next_base, Kmer<k>& rev_compl)
{
    const DNA::Base mapped_base = map_base(next_base);

    // Logically, since a left shift moves the MSN out of the length `k` boundary, the clearing of the base
    // may seem redundant. But, the `to_u64` hashing method implementation works with bytes — not clearing
    // out this base breaks the consistency of the hashing.
    kmer_data[NUM_INTS - 1] &= CLEAR_MSN_MASK;
    left_shift();
    kmer_data[0] |= mapped_base;

    rev_compl.right_shift();
    rev_compl.kmer_data[NUM_INTS - 1] |= (uint64_t(complement(mapped_base)) << (2 * ((k - 1) & 31)));
}


template <uint16_t k>
inline Kmer<k> Kmer<k>::canonical(const Kmer<k>& rev_compl) const
{
    return std::min(*this, rev_compl);
}


template <uint16_t k>
inline Kmer<k> Kmer<k>::canonical() const
{
    return canonical(reverse_complement());
}


template <uint16_t k>
inline std::string Kmer<k>::string_label() const
{
    Kmer<k> kmer(*this);
    char* label = new char[k + 1];

    for(uint16_t idx = 0; idx < k; ++idx)
    {
        switch(kmer.kmer_data[0] & 0b11)
        {
        case DNA::A:
            label[idx] = 'A';
            break;
        
        case DNA::C:
            label[idx] = 'C';
            break;

        case DNA::G:
            label[idx] = 'G';
            break;

        case DNA::T:
            label[idx] = 'T';
            break;

        default:
            label[idx] = 'N';
        }


        kmer.right_shift();
    }


    std::reverse(label, label + k);
    std::string str_label(label, label + k);

    
    delete label;
    
    return str_label;
}


template <uint16_t k>
std::ostream& operator<<(std::ostream& out, const Kmer<k>& kmer)
{
    out << kmer.string_label();

    return out;
}



#endif

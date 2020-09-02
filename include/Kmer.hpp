
#ifndef KMER_HPP
#define KMER_HPP



#include "DNA_Utility.hpp"
#include "smhasher/MurmurHash3.h"
#include "kmc_api/kmc_file.h"
#include "xxhash/xxh3.h"

#include <string>
#include <algorithm>


template <uint16_t k>
class Kmer: public DNA_Utility
{
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

    // Returns the reverese complement of the k-mer.
    Kmer<k> reverse_complement() const;

    // Returns true iff the bitwise encoding of this k-mer is lesser to the
    // encoding of the other k-mer `rhs`.
    bool operator<(const Kmer<k>& rhs) const;

    // Returns true iff this k-mer is identical to the other k-mer `rhs`.
    bool operator==(const Kmer<k>& rhs) const;

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

    // For debugging purposes.
    template <uint16_t K>
    friend std::ostream& operator<<(std::ostream& out, const Kmer<k>& kmer);
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
inline uint64_t Kmer<k>::to_u64(uint64_t seed) const
{
    // The k-mer uses just one 64-bit integer, i.e. k <= 32.
    if(NUM_INTS == 1)
        return kmer_data[0];

    // The k-mer uses more than one 64-bit integer, i.e. k > 32.
    constexpr const uint16_t NUM_BYTES = (k + 3) / 4;
    return XXH3_64bits_withSeed(kmer_data, NUM_BYTES, seed);
    //uint64_t H[2];
    //MurmurHash3_x64_128(kmer_data, NUM_BYTES, seed, H);
    //return H[0] ^ H[1];
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
    /*
    for(uint16_t idx = 0; idx < NUM_INTS; ++idx)
        kmer_data[idx] = rhs.kmer_data[idx];
    */
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
inline bool Kmer<k>::in_forward(const Kmer<k>& kmer_hat) const
{
    return this->operator==(kmer_hat);
}


template <uint16_t k>
inline void Kmer<k>::roll_to_next_kmer(const char next_base, Kmer<k>& rev_compl)
{
    const DNA::Base mapped_base = map_base(next_base);

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

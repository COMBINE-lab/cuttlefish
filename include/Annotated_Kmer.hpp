
#ifndef ANNOTATED_KMER_HPP
#define ANNOTATED_KMER_HPP



#include "globals.hpp"
#include "Directed_Kmer.hpp"
#include "Kmer_Hash_Table.hpp"

#include <cstdint>
#include <cstddef>


// Complete k-mer information: the k-mer itself and its reverse complement, canonical
// form, direction, index in the corresponding sequence, and its state-class.
template <uint16_t k>
class Annotated_Kmer: public Directed_Kmer<k>
{
private:

    size_t idx_;
    cuttlefish::State_Class state_class_;


public:

    Annotated_Kmer()
    {}

    // Constructs an annotated k-mer with its complete information.
    // The template parameter `BITS_PER_KEY` is required to access
    // the hash table `hash`.
    template <uint8_t BITS_PER_KEY>
    Annotated_Kmer(const Kmer<k>& kmer, size_t kmer_idx, const Kmer_Hash_Table<k, BITS_PER_KEY>& hash);

    // Copy constructs the annotated k-mer from `rhs`.
    Annotated_Kmer(const Annotated_Kmer& rhs) = default;

    // Transforms this k-mer by chopping off the first base and
    // appending the next base `next_base` to the end, i.e. rolls
    // the k-mer by one base, sets all the relevant k-mer
    // information accordingly (k-mer state is set using the `hash`).
    template <uint8_t BITS_PER_KEY>
    void roll_to_next_kmer(char next_base, const Kmer_Hash_Table<k, BITS_PER_KEY>& hash);

    void operator=(const Annotated_Kmer<k>& rhs);

    // Returns the index of the k-mer.
    size_t idx() const;

    // Returns the state-class of the k-mer.
    cuttlefish::State_Class state_class() const;
};


template <uint16_t k>
template <uint8_t BITS_PER_KEY>
inline Annotated_Kmer<k>::Annotated_Kmer(const Kmer<k>& kmer, const size_t kmer_idx, const Kmer_Hash_Table<k, BITS_PER_KEY>& hash):
    Directed_Kmer<k>(kmer), idx_(kmer_idx), state_class_(hash[this->canonical_].state_class())
{}


template <uint16_t k>
template <uint8_t BITS_PER_KEY>
inline void Annotated_Kmer<k>::roll_to_next_kmer(const char next_base, const Kmer_Hash_Table<k, BITS_PER_KEY>& hash)
{
    Directed_Kmer<k>::roll_to_next_kmer(next_base);

    idx_++;
    state_class_ = hash[this->canonical_].state_class();
}


template <uint16_t k>
inline void Annotated_Kmer<k>::operator=(const Annotated_Kmer<k>& rhs)
{
    Directed_Kmer<k>::operator=(rhs);
    
    idx_ = rhs.idx_;
    state_class_ = rhs.state_class_;
}


template <uint16_t k>
inline size_t Annotated_Kmer<k>::idx() const
{
    return idx_;
}


template <uint16_t k>
inline cuttlefish::State_Class Annotated_Kmer<k>::state_class() const
{
    return state_class_;
}



#endif

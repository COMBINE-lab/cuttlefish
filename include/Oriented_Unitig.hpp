
#ifndef ORIENTED_UNITIG_HPP
#define ORIENTED_UNITIG_HPP



#include "globals.hpp"

#include <cstdint>
#include <cstddef>
#include <limits>


template <uint16_t k> class CdBG;


// An oriented unitig information: an id for the unitig, its direction,
// index of its first and last k-mers into the underlying reference sequence.
class Oriented_Unitig
{
    template <uint16_t k>
    friend class CdBG;

private:

    uint64_t unitig_id;
    cuttlefish::dir_t dir;
    size_t start_kmer_idx;  // Position, on the ref, of the last k-mer of the unitig occurrence
    size_t end_kmer_idx;    // Position, on the ref, of the first k-mer of the unitig occurrence

    constexpr static uint64_t INVALID_ID = std::numeric_limits<uint64_t>::max();


    Oriented_Unitig(uint64_t unitig_id, cuttlefish::dir_t dir, size_t start_kmer_idx, size_t end_kmer_idx);

    bool is_valid() const;

    size_t length(uint16_t k) const;

public:

    Oriented_Unitig();
};



inline Oriented_Unitig::Oriented_Unitig():
    unitig_id(INVALID_ID)
{}


inline Oriented_Unitig::Oriented_Unitig(const uint64_t unitig_id, const cuttlefish::dir_t dir, const size_t start_kmer_idx, const size_t end_kmer_idx):
    unitig_id(unitig_id), dir(dir), start_kmer_idx(start_kmer_idx), end_kmer_idx(end_kmer_idx)
{}


inline bool Oriented_Unitig::is_valid() const
{
    return unitig_id != INVALID_ID;
}


inline size_t Oriented_Unitig::length(const uint16_t k) const
{
    return end_kmer_idx - start_kmer_idx + k;
}



#endif


#ifndef ORIENTED_UNITIG_HPP
#define ORIENTED_UNITIG_HPP



#include "globals.hpp"

#include <limits>


class CdBG;


// An oriented unitig information: an id for the unitig, its direction,
// index of its first and last k-mers into the underlying reference sequence.
class Oriented_Unitig
{
    friend class CdBG;

private:

    // TODO: Add comments

    uint64_t unitig_id;
    cuttlefish::kmer_dir_t dir;
    size_t start_kmer_idx;
    size_t end_kmer_idx;

    const static uint64_t invalid_id = std::numeric_limits<uint64_t>::max();


    Oriented_Unitig(const uint64_t unitig_id, const cuttlefish::kmer_dir_t dir, const size_t start_kmer_idx, const size_t end_kmer_idx);

    bool is_valid() const;

public:

    Oriented_Unitig();
};



inline Oriented_Unitig::Oriented_Unitig(const uint64_t unitig_id, const cuttlefish::kmer_dir_t dir, const size_t start_kmer_idx, const size_t end_kmer_idx):
    unitig_id(unitig_id), dir(dir), start_kmer_idx(start_kmer_idx), end_kmer_idx(end_kmer_idx)
{}


inline Oriented_Unitig::Oriented_Unitig():
    unitig_id(invalid_id)
{}


inline bool Oriented_Unitig::is_valid() const
{
    return unitig_id != invalid_id;
}



#endif

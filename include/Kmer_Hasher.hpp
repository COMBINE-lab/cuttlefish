
#ifndef KMER_HASHER_HPP
#define KMER_HASHER_HPP



#include "globals.hpp"


template <uint16_t k>
class Kmer_Hasher
{
public:

    // Adopted from the BBHash library.
    // Ref: https://github.com/rizkg/BBHash/blob/48a854a378bce4e2fe4d4cd63bfe5e4f8755dc6e/BooPHF.h#L393
    uint64_t operator()(const Kmer<k>& key, uint64_t seed = 0xAAAAAAAA55555555ULL) const
    {
        return key.to_u64(seed);
        /*
        uint64_t hash = seed;
        const uint64_t key_u64 = key.to_u64();
        hash ^= (hash <<  7) ^  key_u64 * (hash >> 3) ^ (~((hash << 11) + (key_u64 ^ (hash >> 5))));
        hash = (~hash) + (hash << 21);
        hash = hash ^ (hash >> 24);
        hash = (hash + (hash << 3)) + (hash << 8);
        hash = hash ^ (hash >> 14);
        hash = (hash + (hash << 2)) + (hash << 4);
        hash = hash ^ (hash >> 28);
        hash = hash + (hash << 31);

        return hash;
        */
    }
};



#endif

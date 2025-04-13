
#ifndef SOURCE_HASH_HPP
#define SOURCE_HASH_HPP



#include "xxHash/xxhash.h"

#include <cstdint>
#include <cassert>


namespace cuttlefish
{

// Returns the 64-bit hash of the 21-bit source-ID `source`.
inline uint64_t source_hash(const uint32_t source)
{
    assert(source > 0 && source < (1 << 21));
    return XXH3_64bits(&source, 3);
}


// Combines the hashes `h_0` and `h_1` into one. Can be used incrementally. The
// order of the hashes matter.
inline uint64_t hash_combine(const uint64_t h_0, const uint64_t h_1)
{
    // Ref: https://www.boost.org/doc/libs/1_46_1/doc/html/hash/reference.html#boost.hash_combine
    // return h_0 ^ (h_1 + 0x9e3779b9 + (h_0 << 6) + (h_0 >> 2));
    return h_0 ^ h_1;
}

}



#endif

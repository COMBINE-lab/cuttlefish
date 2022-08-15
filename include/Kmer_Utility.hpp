
#ifndef KMER_UTILITY_HPP
#define KMER_UTILITY_HPP



#include "DNA_Utility.hpp"

#include <cstdint>


class Kmer_Utility
{
private:

    // Reverse complement (in the `DNA::Base` representation) of all possible bytes.
    static constexpr uint8_t REVERSE_COMPLEMENT_BYTE[256] =
    {
        255, 191, 127,  63, 239, 175, 111,  47, 223, 159,  95,  31, 207, 143,  79,  15,
        251, 187, 123,  59, 235, 171, 107,  43, 219, 155,  91,  27, 203, 139,  75,  11,
        247, 183, 119,  55, 231, 167, 103,  39, 215, 151,  87,  23, 199, 135,  71,   7,
        243, 179, 115,  51, 227, 163,  99,  35, 211, 147,  83,  19, 195, 131,  67,   3,
        254, 190, 126,  62, 238, 174, 110,  46, 222, 158,  94,  30, 206, 142,  78,  14,
        250, 186, 122,  58, 234, 170, 106,  42, 218, 154,  90,  26, 202, 138,  74,  10,
        246, 182, 118,  54, 230, 166, 102,  38, 214, 150,  86,  22, 198, 134,  70,   6,
        242, 178, 114,  50, 226, 162,  98,  34, 210, 146,  82,  18, 194, 130,  66,   2,
        253, 189, 125,  61, 237, 173, 109,  45, 221, 157,  93,  29, 205, 141,  77,  13,
        249, 185, 121,  57, 233, 169, 105,  41, 217, 153,  89,  25, 201, 137,  73,   9,
        245, 181, 117,  53, 229, 165, 101,  37, 213, 149,  85,  21, 197, 133,  69,   5,
        241, 177, 113,  49, 225, 161,  97,  33, 209, 145,  81,  17, 193, 129,  65,   1,
        252, 188, 124,  60, 236, 172, 108,  44, 220, 156,  92,  28, 204, 140,  76,  12,
        248, 184, 120,  56, 232, 168, 104,  40, 216, 152,  88,  24, 200, 136,  72,   8,
        244, 180, 116,  52, 228, 164, 100,  36, 212, 148,  84,  20, 196, 132,  68,   4,
        240, 176, 112,  48, 224, 160,  96,  32, 208, 144,  80,  16, 192, 128,  64,   0
    };


public:


    // Returns the reverse completement byte of the 4-mer `byte`;
    // both are to be in the `DNA::Base` representation.
    static uint8_t reverse_complement(const uint8_t byte)
    {
        return REVERSE_COMPLEMENT_BYTE[byte];
    }

    // Returns the binary encoding word of the literal k-mer `label`.
    template <uint16_t k>
    static uint64_t encode(const char* label);
};


template <uint16_t k>
inline uint64_t Kmer_Utility::encode(const char* const label)
{
    static_assert(0 < k && k <= 32, "invalid k-mer label length for machine word encoding");

    if constexpr(k > 1)
        return (static_cast<uint64_t>(DNA_Utility::map_base(*label)) << (2 * (k - 1))) | encode<k - 1>(label + 1);

    return static_cast<uint64_t>(DNA_Utility::map_base(*label));
}



#endif

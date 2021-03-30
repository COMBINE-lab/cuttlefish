
#ifndef DNA_UTILITY_HPP
#define DNA_UTILITY_HPP



#include "DNA.hpp"


class DNA_Utility
{
private:

    static constexpr DNA::Base A = DNA::A;
    static constexpr DNA::Base C = DNA::C;
    static constexpr DNA::Base G = DNA::G;
    static constexpr DNA::Base T = DNA::T;
    static constexpr DNA::Base N = DNA::N;

    // Mapped `DNA::Base` for the ASCII characters in the range [0, 127].
    static constexpr DNA::Base MAPPED_BASE[128] =
    {
        N, N, N, N, N, N, N, N, N, N,   // 0 - 9
        N, N, N, N, N, N, N, N, N, N,   // 10 - 19
        N, N, N, N, N, N, N, N, N, N,   // 20 - 29
        N, N, N, N, N, N, N, N, N, N,   // 30 - 39
        N, N, N, N, N, N, N, N, N, N,   // 40 - 49
        N, N, N, N, N, N, N, N, N, N,   // 50 - 59
        N, N, N, N, N, A, N, C, N, N,   // 60 - 69
        N, G, N, N, N, N, N, N, N, N,   // 70 - 79
        N, N, N, N, T, N, N, N, N, N,   // 80 - 89
        N, N, N, N, N, N, N, A, N, C,   // 90 - 99
        N, N, N, G, N, N, N, N, N, N,   // 100 - 109
        N, N, N, N, N, N, T, N, N, N,   // 110 - 119
        N, N, N, N, N, N, N, N          // 120 - 127
    };

    // Mapped complement `DNA::Base` for the ASCII characters in the range [0, 127].
    static constexpr DNA::Base COMPLEMENTED_BASE[5] =
    {
        T, G, C, A, N
    };

    // DNA-complement characters for the ASCII characters in the range [0, 127] 
    static constexpr char COMPLEMENTED_CHAR[128] = 
    {
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   // 0 - 9
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   // 10 - 19
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   // 20 - 29
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   // 30 - 39
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   // 40 - 49
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   // 50 - 59
        'N', 'N', 'N', 'N', 'N', 'T', 'N', 'G', 'N', 'N',   // 60 - 69
        'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',   // 70 - 79
        'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N',   // 80 - 89
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'T', 'N', 'G',   // 90 - 99
        'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N',   // 100 - 109
        'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N',   // 110 - 119
        'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'              // 120 - 127
    };

    // Booleans denoting a ASCII character is to be considered a placeholder base or not.
    static constexpr bool IS_PLACEHOLDER[128] =
    {
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,       // 0 - 9
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,       // 10 - 19
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,       // 20 - 29
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,       // 30 - 39
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,       // 40 - 49
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,       // 50 - 59
        1, 1, 1, 1, 1, 0, 1, 0, 1, 1,       // 60 - 69
        1, 0, 1, 1, 1, 1, 1, 1, 1, 1,       // 70 - 79
        1, 1, 1, 1, 0, 1, 1, 1, 1, 1,       // 80 - 89
        1, 1, 1, 1, 1, 1, 1, 0, 1, 0,       // 90 - 99
        1, 1, 1, 0, 1, 1, 1, 1, 1, 1,       // 100 - 109
        1, 1, 1, 1, 1, 1, 0, 1, 1, 1,       // 110 - 119
        1, 1, 1, 1, 1, 1, 1, 1              // 120 - 127
    };

    // TODO: Move these new k-mer specific (and not DNA-base specific) stuffs to a separate class.
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

    // Mapped `DNA::Extended_Base` for the corresponding `DNA::Base`, i.e.
    // a mapping from [0(A) — T(3)] to [1(A) — 4(T)].
    static constexpr DNA::Extended_Base MAPPED_EXTENDED_BASE[4] =
    {
        DNA::Extended_Base::A, DNA::Extended_Base::C, DNA::Extended_Base::G, DNA::Extended_Base::T
    };


public:

    // Returns the mapping integer value of the given character `base`.
    static DNA::Base map_base(const char base)
    {
        return MAPPED_BASE[uint8_t(base)];
    }

    // Returns the mapping integer value of the complement of `base`.
    static DNA::Base complement(const DNA::Base base)
    {
        return COMPLEMENTED_BASE[base];
    }
        
    // Returns the DNA-complement (upper-case) character of the character `base`.
    static char complement(const char base)
    {
        return COMPLEMENTED_CHAR[uint8_t(base)];
    }

    // Returns `true` iff the character `base` is a placeholder character.
    static bool is_placeholder(const char base)
    {
        return IS_PLACEHOLDER[uint8_t(base)];
    }

    // Returns the upper-case equivalent of the character `base`.
    static char upper(const char base)
    {
        return base <= 'T' ? base : (base - ('a' - 'A'));
    }

    // Returns the reverse completement byte of the 4-mer `byte`;
    // both are to be in the `DNA::Base` representation.
    static uint8_t reverse_complement(const uint8_t byte)
    {
        return REVERSE_COMPLEMENT_BYTE[byte];
    }

    // Returns the mapping `DNA::Extended_Base` representation of the
    // `DNA::Base` representation `base`.
    static DNA::Extended_Base map_extended_base(const DNA::Base base)
    {
        return MAPPED_EXTENDED_BASE[base];
    }
};



#endif

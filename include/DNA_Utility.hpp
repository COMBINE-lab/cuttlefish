
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
        return base == 'N' || base == 'n';
    }

    // Returns the upper-case equivalent of the character `base`.
    static char upper(const char base)
    {
        return base <= 'T' ? base : (base - ('a' - 'A'));
    }
};



#endif

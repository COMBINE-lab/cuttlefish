
#ifndef DNA_UTILITY_HPP
#define DNA_UTILITY_HPP



#include "DNA.hpp"

#include <cstdint>
#include <cstddef>


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

    // Mapped ASCII characters for the `DNA::Base` notations.
    static constexpr char MAPPED_CHAR[4] = 
    {
        'A', 'C', 'G', 'T'
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

    // Mapped `DNA::Extended_Base` for the corresponding `DNA::Base`, i.e.
    // a mapping from [0(A) — T(3)] to [1(A) — 4(T)].
    static constexpr DNA::Extended_Base MAPPED_EXTENDED_BASE[4] =
    {
        DNA::Extended_Base::A, DNA::Extended_Base::C, DNA::Extended_Base::G, DNA::Extended_Base::T
    };

    // Mapped `DNA::Base` for the corresponding `DNA::Extended_Base`, i.e.
    // a mapping from [1(A) — 4(3)] to [0(A) — 3(T)].
    static constexpr DNA::Base REVERSE_MAPPED_EXTENDED_BASE[5] =
    {
        DNA::Base::N, DNA::Base::A, DNA::Base::C, DNA::Base::G, DNA::Base::T
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

    // Returns the mapping character of the nucleobase `base`.
    static char map_char(const DNA::Base base)
    {
        return MAPPED_CHAR[static_cast<std::size_t>(base)];
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

    // Returns the mapping `DNA::Extended_Base` representation of the
    // `DNA::Base` representation `base`.
    static DNA::Extended_Base map_extended_base(const DNA::Base base)
    {
        return MAPPED_EXTENDED_BASE[base];
    }

    // Returns the mapping `DNA::Base` representation of the
    // `DNA::Extended_Base` representation `extended_base`.
    static DNA::Base map_base(const DNA::Extended_Base extended_base)
    {
        return REVERSE_MAPPED_EXTENDED_BASE[static_cast<std::size_t>(extended_base)];
    }
};



#endif

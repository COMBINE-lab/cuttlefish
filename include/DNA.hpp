
#ifndef DNA_HPP
#define DNA_HPP



#include <cstdint>


namespace DNA
{
    // A = 0, C = 1, G = 2, T = 3.
    // Note that, this is not possible to change this mapping w/o modifications to the
    // interfacing of our code with the KMC api. This mapping is essential for some
    // crucial performance hacks in the interfacing.
    // TODO: consider having it as `enum class`.
    enum Base: uint8_t
    {
        A = 0b000,  // 0
        C = 0b001,  // 1
        G = 0b010,  // 2
        T = 0b011,  // 3
        N = 0b100,  // 4 (ambiguous due to branching)
        E = 0b101   // 5 (empty)
    };


    // E = 0, A = 1, C = 2, G = 3, T = 4, N = 7;
    // O/P and — non-branching = 5, branching = 6.
    enum class Extended_Base: uint8_t
    {
        E = 0b000,  // 0
        A = 0b001,  // 1
        C = 0b010,  // 2
        G = 0b011,  // 3
        T = 0b100,  // 4
        N = 0b111,  // 7
        OP_non_branch = 0b101,  // 5
        OP_branching = 0b110,   // 6
    };
}



#endif
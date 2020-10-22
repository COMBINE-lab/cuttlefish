
#ifndef DNA_HPP
#define DNA_HPP



#include <cstdint>


namespace DNA
{
    // A = 0, C = 1, G = 2, T = 3.
    // Note that, this is not possible to change this mapping w/o modifications to the
    // interfacing of our code with the KMC api. This mapping is essential for some
    // crucial performance hacks in the interfacing.
    enum Base: uint8_t
    {
        A = 0b00,   // 0
        C = 0b01,   // 1
        G = 0b10,   // 2
        T = 0b11,   // 3
        N = 0b100   // 4
    };
}



#endif

#ifndef OUTPUT_FORMATS_HPP
#define OUTPUT_FORMATS_HPP



#include <cstdint>


namespace cuttlefish
{
    // Output format options for the algorithm.
    enum Output_Format: uint8_t
    {
        txt = 0,
        gfa1 = 1,
        gfa2 = 2,
        num_op_formats
    };
}



#endif

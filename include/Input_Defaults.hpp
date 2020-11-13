
#ifndef INPUT_DEFAULTS_HPP
#define INPUT_DEFAULTS_HPP



#include <Output_Format.hpp>


namespace cuttlefish
{
    // Default input options for the algorithm.
    namespace _default
    {
        constexpr char EMPTY[] = "";
        constexpr uint16_t K = 25;  // Set as per the KMC3 default.
        constexpr uint16_t THREAD_COUNT = 1;
        constexpr uint16_t OP_FORMAT = Output_Format::txt;
        constexpr char WORK_DIR[] = ".";
    }
}



#endif


#ifndef INPUT_DEFAULTS_HPP
#define INPUT_DEFAULTS_HPP



#include "Output_Format.hpp"


namespace cuttlefish
{
    // Default input options for the algorithm.
    namespace _default
    {
        constexpr char EMPTY[] = "";
        constexpr uint16_t K = 25;  // Set as per the KMC3 default.
        constexpr uint32_t CUTOFF_FREQ = 2; // Typical practice
        constexpr uint16_t THREAD_COUNT = 1;
        constexpr std::size_t MAX_MEMORY = 3;   // Set as per KMC3 stage 1 performance.
#ifdef CF_DEVELOP_MODE
        constexpr double GAMMA = 0;
#endif
        constexpr uint16_t OP_FORMAT = Output_Format::fa;
        constexpr char WORK_DIR[] = ".";
    }
}



#endif

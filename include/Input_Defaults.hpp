
#ifndef INPUT_DEFAULTS_HPP
#define INPUT_DEFAULTS_HPP



#include "Output_Format.hpp"


namespace cuttlefish
{
    // Default input options for the algorithm.
    namespace _default
    {
        constexpr char EMPTY[] = "";
        constexpr uint16_t K = 27;
        constexpr uint32_t CUTOFF_FREQ_READS = 2;   // Typical practice
        constexpr uint32_t CUTOFF_FREQ_REFS = 1;    // Typical assumption
        constexpr uint16_t THREAD_COUNT = 1;
        constexpr std::size_t MAX_MEMORY = 3;   // Set as per KMC3 stage 1 performance.
#ifdef CF_DEVELOP_MODE
        constexpr double GAMMA = 0;
#endif
        constexpr Output_Format OP_FORMAT = Output_Format::fa;
        constexpr char WORK_DIR[] = ".";
    }
}



#endif

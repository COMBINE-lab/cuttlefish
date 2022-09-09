
#ifndef INPUT_DEFAULTS_HPP
#define INPUT_DEFAULTS_HPP



#include "Output_Format.hpp"

#include <cstdint>
#include <cstddef>
#include <thread>


namespace cuttlefish
{
    // Default input options for the algorithm.
    namespace _default
    {
        constexpr char EMPTY[] = "";
        constexpr uint16_t K = 27;
        constexpr uint32_t CUTOFF_FREQ_READS = 2;   // Typical practice
        constexpr uint32_t CUTOFF_FREQ_REFS = 1;    // Typical assumption
        const uint16_t THREAD_COUNT = (std::thread::hardware_concurrency() ?
                                        (std::thread::hardware_concurrency()/ 4) : 8);  // A quarter of the total thread-count.
        constexpr std::size_t MAX_MEMORY = 3;   // Set as per KMC3 stage 1 performance.
#ifdef CF_DEVELOP_MODE
        constexpr double GAMMA = 0;
#endif
        constexpr Output_Format OP_FORMAT = Output_Format::fa;
        constexpr char WORK_DIR[] = ".";
    }
}



#endif

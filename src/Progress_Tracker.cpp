
#include "Progress_Tracker.hpp"


void Progress_Tracker::setup(uint64_t total_work_load, uint64_t work_chunk_threshold, const std::string& log_message)
{
    this->total_work_load = total_work_load;
    this->work_chunk_threshold = work_chunk_threshold;

    total_work_done = 0;
    percent_work_done = 0;

    this->log_message = log_message;

    std::cerr << "\n";
}

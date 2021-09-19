
#ifndef PROGRESS_TRACKER_HPP
#define PROGRESS_TRACKER_HPP



#include "Spin_Lock.hpp"

#include <string>
#include <cmath>
#include <iostream>


// A basic class to track and display progress for some work.
class Progress_Tracker
{
private:

    uint64_t total_work_load;  //   Total amount of work to be done over time.
    uint64_t work_chunk_threshold; // Granularity of the provided work chunk sizes that triggers tracking updates.

    uint64_t total_work_done;   // Amount of work done until now.
    uint16_t percent_work_done; // Percentage of the completed workload.
    std::string log_message;  // Message to display at the logs.

    Spin_Lock lock;  // Lock to ensure multiple threads can access the tracker safely.

public:

    // Sets up the tracker for some task with total size `total_work_load`; updates to the tracking are to
    // be triggered when some work-chunk of size at least `work_chunk_threshold` is provided to it. The
    // log message to be displayed over the course of tracking is `log_message`.
    void setup(uint64_t total_work_load, uint64_t work_chunk_threshold, const std::string& log_message);

    // Tracks progress made for a work-chunk of size `work_chunk_size`. If an update it made towards progress,
    // then the chunk-size is set to 0 to refresh it for the next cycle.
    // Note that, the chunk-size must be at least `work_chunk_threshold` for any updates to be made towards
    // the progress. All lesser sized chunk update requests are ignored. So, repeated invocation is suggested.
    void track_work(uint64_t& work_chunk_size);
};


inline void Progress_Tracker::track_work(uint64_t& work_chunk_size)
{
    if(work_chunk_size >= work_chunk_threshold)
    {
        lock.lock();
        
        total_work_done += work_chunk_size;

        const uint16_t new_percent = static_cast<uint16_t>(std::round((total_work_done * 100.0) / total_work_load));
        if(percent_work_done < new_percent)
        {
            percent_work_done = new_percent;
            std::cerr << "\r[" << log_message << "]\t" << percent_work_done << "%";
        }

        lock.unlock();


        work_chunk_size = 0;
    }
}



#endif

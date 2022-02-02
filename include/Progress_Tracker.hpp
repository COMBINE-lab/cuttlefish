
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

    // Tracks progress made for a work-chunk of size `work_chunk_size`. If the passed chunk-size is large
    // enough, then it is tracked and `true` is returned. All lesser sized chunk update requests are ignored
    // and `false` is returned. So, repeated invocation is suggested.
    bool track_work(uint64_t work_chunk_size);
};


inline bool Progress_Tracker::track_work(const uint64_t work_chunk_size)
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

        return true;
    }

    return false;
}



#endif

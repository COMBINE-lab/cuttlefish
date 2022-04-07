
#ifndef JOB_QUEUE_HPP
#define JOB_QUEUE_HPP



#include "Spin_Lock.hpp"

#include <cstdint>
#include <atomic>
#include <queue>


// A very basic single-producer single-consumer job queue — for jobs with
// id's and additional information of types  `T_id_` and `T_info_`, respectively.
template <typename T_id_, typename T_info_>
class Job_Queue
{

private:

    Spin_Lock job_lock; // Mutual-exclusion lock to avoid concurrent modification (post / fetch) of the job queue.
    std::queue<T_id_> job_id_queue;   // FIFO queue for the id (name) of the jobs.
    std::queue<T_info_> job_info_queue; // FIFO queue for additional information required for the jobs.

    std::atomic<bool> jobs_in_future;   // Whether more jobs are expected to be posted in future.
    std::atomic<uint64_t> jobs_posted;  // Total number of jobs that have been posted in the job queue.
    std::atomic<uint64_t> jobs_fetched; // Total number of jobs that have been fetched from the job queue.
    std::atomic<uint64_t> jobs_finished;    // Total number of jobs that have been completed from the job queue.


public:

    // Constructs an empty job queue.
    Job_Queue():
        jobs_in_future(true), jobs_posted(0), jobs_fetched(0), jobs_finished(0)
    {}

    // Posts a job to the queue with id `job_id` and additional information `job_info`.
    void post_job(const T_id_& job_id, const T_info_& job_info);

    // Fetches the next job from the queue — puts its id into `job_id` and additional information into `job_info`.
    void fetch_job(T_id_& job_id, T_info_& job_info);

    // Mark the last fetched job as completed.
    void finish_job();

    // Returns `true` iff there are unfinished jobs available in the queue.
    bool job_available() const;

    // Returns `true` iff there are jobs remaining to be finished, either present in the queue or expected in the future.
    bool jobs_remain() const;

    // Signals the queue that there are no jobs expected to be posted in the future.
    void signal_end();

    // Returns the sequence-id of the next job to complete (from the queue).
    uint64_t next_job_to_finish() const;

    // Returns the sequence-id of the next job to post (into the queue).
    uint64_t next_job_to_post() const;
};


template <typename T_id_, typename T_info_>
inline void Job_Queue<T_id_, T_info_>::post_job(const T_id_& job_id, const T_info_& job_info)
{
    job_lock.lock();

    job_id_queue.push(job_id);
    job_info_queue.push(job_info);

    job_lock.unlock();


    jobs_posted++;
}


template <typename T_id_, typename T_info_>
inline void Job_Queue<T_id_, T_info_>::fetch_job(T_id_& job_id, T_info_& job_info)
{
    if(job_id_queue.empty())
    {
        std::cerr << "Fetch attempted on an empty job queue. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    job_lock.lock();
    
    job_id = job_id_queue.front();
    job_info = job_info_queue.front();

    job_id_queue.pop();
    job_info_queue.pop();
    
    job_lock.unlock();


    jobs_fetched++;
}


template <typename T_id_, typename T_info_>
inline void Job_Queue<T_id_, T_info_>::finish_job()
{
    jobs_finished++;
}


template <typename T_id_, typename T_info_>
inline bool Job_Queue<T_id_, T_info_>::job_available() const
{
    return jobs_fetched < jobs_posted;
}


template <typename T_id_, typename T_info_>
inline bool Job_Queue<T_id_, T_info_>::jobs_remain() const
{
    return jobs_in_future || (jobs_finished < jobs_posted);
}


template <typename T_id_, typename T_info_>
inline void Job_Queue<T_id_, T_info_>::signal_end()
{
    jobs_in_future = false;
}


template <typename T_id_, typename T_info_>
inline uint64_t Job_Queue<T_id_, T_info_>::next_job_to_finish() const
{
    return jobs_finished + 1;
}


template <typename T_id_, typename T_info_>
inline uint64_t Job_Queue<T_id_, T_info_>::next_job_to_post() const
{
    return jobs_posted + 1;
}



#endif

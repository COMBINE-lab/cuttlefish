
#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP



#include "Task_Params.hpp"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <thread>


// A basic thread pool class to support avoidance of latency incurred with frequent
// construction and destruction of threads throughout the compaction algorithm.
template <uint16_t k>
class Thread_Pool
{
public:

    // Types of tasks supported by the thread pool.
    enum class Task_Type
    {
        classification,
        output_plain,
        output_gfa,
        output_gfa_reduced,
        compute_states_read_space,
        extract_unipaths_read_space,
    };


private:

    // Status of tasks provided to each thread in the pool.
    enum class Task_Status: uint8_t
    {
        pending,    // tasks yet to be provided;
        available,  // a task is available and waiting to be processed;
        no_more     // no tasks will be provided anymore.
    };
    

    // Number of threads in the pool.
    const uint16_t thread_count;
    
    // The de Bruijn graph object that this thread pool is operating with.
    void* const dBG;

    // The type of task that this thread pool will execute.
    const Task_Type task_type;

    // Collection of the task statuses of each thread.
    std::atomic<Task_Status>* const task_status;

    // The collection of the threads in the pool.
    std::vector<std::thread> thread_pool;

    // Collection of the task parameters for each thread.
    std::vector<Classification_Task_Params> classify_params;
    std::vector<Output_Task_Params> output_params;
    std::vector<Read_dBG_Compaction_Params> read_dBG_compaction_params;


    // Marks the thread number `thread_id` as busy with some task.
    void assign_task(uint16_t thread_id);

    // Marks the thread number `thread_id` as free to accept tasks.
    void free_thread(uint16_t thread_id);

    // Runs tasks with the thread number `thread_id` as long as new tasks are provided
    // to it. Halts at receiving a signal (through `task_status[thread_id]`) that no
    // more tasks will be provided.
    void task(uint16_t thread_id);


public:

    // Constructs a thread pool with `thread_count` number of threads to operate
    // on the de Brujin graph `dBG` for tasks of type `task_type`.
    Thread_Pool(uint16_t thread_count, void* dBG, Task_Type task_type);

    // Returns the id (number) of an idle thread from the pool.
    uint16_t get_idle_thread() const;

    // Waits until the thread number `thread_id` becomes idle.
    void get_thread(uint16_t thread_id) const;

    // Assigns a classification task to the thread number `thread_id` with the provided parameters.
    void assign_classification_task(uint16_t thread_id, const char* seq, size_t seq_len, size_t left_end, size_t right_end);

    // Assigns an outputting task to the thread number `thread_id` with the provided parameters.
    void assign_output_task(uint16_t thread_id, const char* seq, size_t seq_len, size_t left_end, size_t right_end);

    // Assigns a read-dBG compaction task, either DFA-states computation or maximal unitigs extraction,
    // to the thread number `thread_id`; the edges (i.e. (k + 1)-mers) or vertices (i.e. k-mers),
    // respectively, are parsed using `parser`.
    void assign_read_dBG_compaction_task(void* parser, uint16_t thread_id);

    // Waits until all the threads in the pool have completed their active tasks.
    void wait_completion() const;

    // Closes the thread pool.
    void close();
};



#endif

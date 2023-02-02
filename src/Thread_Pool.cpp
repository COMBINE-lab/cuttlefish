
#include "Thread_Pool.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "CdBG.hpp"
#include "Read_CdBG_Constructor.hpp"
#include "Read_CdBG_Extractor.hpp"

#include <iostream>


template <uint16_t k>
Thread_Pool<k>::Thread_Pool(const uint16_t thread_count, void* const dBG, const Task_Type task_type):
    thread_count(thread_count),
    dBG(dBG),
    task_type(task_type),
    task_status(new std::atomic<Task_Status>[thread_count])
{
    // Mark the status of the task for each thread as `pending`.
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
        task_status[t_id] = Task_Status::pending;

    
    // Resize the parameters collections.
    switch(task_type)
    {
    case Task_Type::classification:
        classify_params.resize(thread_count);
        break;
    
    case Task_Type::output_plain:
    case Task_Type::output_gfa:
    case Task_Type::output_gfa_reduced:
        output_params.resize(thread_count);
        break;

    case Task_Type::compute_states_read_space:
    case Task_Type::extract_unipaths_read_space:
        read_dBG_compaction_params.resize(thread_count);
        break;

    default:
        std::cerr << "Unrecognized task type encountered in thread pool. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    // Launch the threads.
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
        thread_pool.emplace_back(&Thread_Pool::task, this, t_id);
}


template <uint16_t k>
void Thread_Pool<k>::task(const uint16_t thread_id)
{
    while(true)
    {
        // Busy-wait for some task.
        while(task_status[thread_id] == Task_Status::pending);

        // No more tasks to come in on the future.
        if(task_status[thread_id] == Task_Status::no_more)
            return;


        // Some task is available for the thread number `thread_id`.
        if(task_status[thread_id] == Task_Status::available)
        {
            switch(task_type)
            {
            case Task_Type::classification:
                {
                    const Classification_Task_Params& params = classify_params[thread_id];
                    static_cast<CdBG<k>*>(dBG)->process_substring(params.seq, params.seq_len, params.left_end, params.right_end);
                }
                break;

            case Task_Type::output_plain:
                {
                    const Output_Task_Params& params = output_params[thread_id];
                    static_cast<CdBG<k>*>(dBG)->output_plain_off_substring(params.thread_id, params.seq, params.seq_len, params.left_end, params.right_end);
                }
                break;

            case Task_Type::output_gfa:
            case Task_Type::output_gfa_reduced:
                {
                    const Output_Task_Params& params = output_params[thread_id];
                    static_cast<CdBG<k>*>(dBG)->output_gfa_off_substring(params.thread_id, params.seq, params.seq_len, params.left_end, params.right_end);
                }
                break;
            
            case Task_Type::compute_states_read_space:
                {
                    const Read_dBG_Compaction_Params& params = read_dBG_compaction_params[thread_id];
                    static_cast<Read_CdBG_Constructor<k>*>(dBG)->
                        process_edges(static_cast<Kmer_SPMC_Iterator<k + 1>*>(params.parser), params.thread_id);
                }
                break;

            case Task_Type::extract_unipaths_read_space:
                {
                    const Read_dBG_Compaction_Params& params = read_dBG_compaction_params[thread_id];
                    static_cast<Read_CdBG_Extractor<k>*>(dBG)->
                        process_vertices(static_cast<Kmer_SPMC_Iterator<k>*>(params.parser), params.thread_id);
                }
                break;
            }


            free_thread(thread_id);
        }
    }
}


template <uint16_t k>
uint16_t Thread_Pool<k>::get_idle_thread() const
{
    int32_t idle_thread_id = -1;
    uint16_t t_id = 0;

    while(idle_thread_id == -1)
        if(task_status[t_id] == Task_Status::pending)
            idle_thread_id = t_id;
        else
            t_id = (t_id + 1) % thread_count;

    
    return idle_thread_id;
}


template <uint16_t k>
void Thread_Pool<k>::get_thread(const uint16_t thread_id) const
{
    while(task_status[thread_id] != Task_Status::pending);
}


template <uint16_t k>
void Thread_Pool<k>::assign_classification_task(const uint16_t thread_id, const char* const seq, const size_t seq_len, const size_t left_end, const size_t right_end)
{
    classify_params[thread_id] = Classification_Task_Params(seq, seq_len, left_end, right_end);

    assign_task(thread_id);
}


template <uint16_t k>
void Thread_Pool<k>::assign_output_task(const uint16_t thread_id, const char* const seq, const size_t seq_len, const size_t left_end, const size_t right_end)
{
    output_params[thread_id] = Output_Task_Params(thread_id, seq, seq_len, left_end, right_end);

    assign_task(thread_id);
}


template <uint16_t k>
void Thread_Pool<k>::assign_read_dBG_compaction_task(void* const parser, const uint16_t thread_id)
{
    read_dBG_compaction_params[thread_id] = Read_dBG_Compaction_Params(parser, thread_id);

    assign_task(thread_id);
}


template <uint16_t k>
void Thread_Pool<k>::assign_task(const uint16_t thread_id)
{
    if(task_status[thread_id] != Task_Status::pending)
    {
        std::cerr << "Expected thread " << thread_id << " to be idle while assigning a job, but found it unexpecedly busy. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    
    task_status[thread_id] = Task_Status::available;
}


template <uint16_t k>
void Thread_Pool<k>::free_thread(const uint16_t thread_id)
{
    if(task_status[thread_id] != Task_Status::available)
    {
        std::cerr << "Expected thread " << thread_id << " to be busy while assigning a job, but found it unexpecedly free. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    task_status[thread_id] = Task_Status::pending;
}


template <uint16_t k>
void Thread_Pool<k>::wait_completion() const
{
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
        while(task_status[t_id] == Task_Status::available);
}


template <uint16_t k>
void Thread_Pool<k>::close()
{
    // Wait for all the threads to finish.
    wait_completion();

    // Close all the threads.
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        // Signal each thread to stop running.
        task_status[t_id] = Task_Status::no_more;
        
        if(!thread_pool[t_id].joinable())
        {
            std::cerr << "Early termination of a worker thread encountered. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        thread_pool[t_id].join();
    }


    delete[] task_status;
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Thread_Pool)


#include "Read_CdBG_Constructor.hpp"


template <uint16_t k>
Read_CdBG_Constructor<k>::Read_CdBG_Constructor(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table):
    params(params),
    hash_table(hash_table)
{}


template <uint16_t k>
void Read_CdBG_Constructor<k>::compute_DFA_states()
{
    // Construct a thread pool.
    const uint16_t thread_count = params.thread_count();
    Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::compute_states_read_space);

    // Multi-threaded computation.
    distribute_states_computation(thread_pool);
    thread_pool.close();
}


template <uint16_t k>
void Read_CdBG_Constructor<k>::distribute_states_computation(Thread_Pool<k>& thread_pool)
{
    const uint16_t thread_count = params.thread_count();

    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        const uint16_t idle_thread_id = thread_pool.get_idle_thread();
        thread_pool.assign_compute_states_read_space_task(idle_thread_id);
    }
}


template <uint16_t k>
void Read_CdBG_Constructor<k>::process_edges(const uint16_t thread_id)
{
    (void)thread_id;
}



// Template instantiations for the required specializations.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG_Constructor)

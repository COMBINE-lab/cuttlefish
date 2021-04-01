
#include "Read_CdBG_Constructor.hpp"
#include "Edge.hpp"

#include "chrono"


template <uint16_t k>
Read_CdBG_Constructor<k>::Read_CdBG_Constructor(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table):
    params(params),
    hash_table(hash_table),
    edge_container(params.edge_db_path()),
    edge_parser(&edge_container, params.thread_count())
{
    std::cout << "Total number of distinct edges: " << edge_container.size() << ".\n";
}


template <uint16_t k>
void Read_CdBG_Constructor<k>::compute_DFA_states()
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Construct a thread pool.
    const uint16_t thread_count = params.thread_count();
    Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::compute_states_read_space);

    // Launch the reading (and parsing per demand) of the edges from disk.
    edge_parser.launch_production();

    // Launch (multi-threaded) computation of the states.
    distribute_states_computation(thread_pool);

    // Wait for the edges to be depleted from the database.
    edge_parser.seize_production();

    // Wait for the consumer threads to finish parsing and processing the edges.
    thread_pool.close();

    std::cout << "Number of processed egdes: " << edges_processed << "\n";


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done computing the DFA states. Time taken = " << elapsed_seconds << " seconds.\n";
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
    Edge<k> e;
    uint64_t edge_count = 0;

    while(edge_parser.tasks_expected(thread_id))
        if(edge_parser.value_at(thread_id, e.e()))
        {
            e.configure();  // A new edge (k + 1)-mer has been parsed; set the relevant k-mer and sides information.

            if(e.is_loop())
                while(!add_loop(e.u_hat(), e.s_u_hat(), e.s_v_hat()));
            else
            {
                while(!add_incident_edge(e.u_hat(), e.s_u_hat(), e.edge_encoding_u()));
                while(!add_incident_edge(e.v_hat(), e.s_v_hat(), e.edge_encoding_v()));
            }

            edge_count++;
        }

    lock.lock();
    std::cout << "Thread " << thread_id << " processed " << edge_count << " edges.\n";  // Temporary log. TODO: remove.
    edges_processed += edge_count;
    lock.unlock();
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG_Constructor)

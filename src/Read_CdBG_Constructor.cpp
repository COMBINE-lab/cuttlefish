
#include "Read_CdBG_Constructor.hpp"
#include "Edge.hpp"
#include "utility.hpp"

#include <cmath>
#include <chrono>


template <uint16_t k>
Read_CdBG_Constructor<k>::Read_CdBG_Constructor(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table):
    params(params),
    hash_table(hash_table)
{}


template <uint16_t k>
void Read_CdBG_Constructor<k>::compute_DFA_states(const std::string& edge_db_path)
{
    // std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    const Kmer_Container<k + 1> edge_container(edge_db_path);  // Wrapper container for the edge-database.
    Kmer_SPMC_Iterator<k + 1> edge_parser(&edge_container, params.thread_count());  // Parser for the edges from the edge-database.
    edge_count_ = edge_container.size();
    std::cout << "Total number of distinct edges: " << edge_count_ << ".\n";


    const std::string& buckets_file_path = params.buckets_file_path();
    if(!buckets_file_path.empty() && file_exists(buckets_file_path))    // The serialized hash table buckets, saved from some earlier execution, exists.
    {
        std::cout <<    "Found the hash table buckets at file " << buckets_file_path << ".\n"
                        "Loading the buckets.\n";
        hash_table.load_hash_buckets(buckets_file_path);
        std::cout << "Loaded the buckets into memory.\n";
    }
    else
    {
        // Construct a thread pool.
        const uint16_t thread_count = params.thread_count();
        Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::compute_states_read_space);

        // Launch the reading (and parsing per demand) of the edges from disk.
        edge_parser.launch_production();

        // Launch (multi-threaded) computation of the states.
        const uint64_t thread_load_percentile = static_cast<uint64_t>(std::round((edge_count_ / 100.0) / params.thread_count()));
        progress_tracker.setup(edge_count_, thread_load_percentile, "Computing DFA states");
        distribute_states_computation(&edge_parser, thread_pool);

        // Wait for the edges to be depleted from the database.
        edge_parser.seize_production();

        // Wait for the consumer threads to finish parsing and processing the edges.
        thread_pool.close();

        std::cout << "\nNumber of processed edges: " << edges_processed << "\n";
    }


    // std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    // double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    // std::cout << "Done computing the DFA states. Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
void Read_CdBG_Constructor<k>::distribute_states_computation(Kmer_SPMC_Iterator<k + 1>* const edge_parser, Thread_Pool<k>& thread_pool)
{
    const uint16_t thread_count = params.thread_count();

    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        const uint16_t idle_thread_id = thread_pool.get_idle_thread();
        thread_pool.assign_read_dBG_compaction_task(edge_parser, idle_thread_id);
    }
}


template <uint16_t k>
void Read_CdBG_Constructor<k>::process_edges(Kmer_SPMC_Iterator<k + 1>* const edge_parser, const uint16_t thread_id)
{
    // Data locations to be reused per each edge processed.
    Edge<k> e;  // For the edges to be processed one-by-one.
/*
    cuttlefish::edge_encoding_t e_front, e_back;    // Edges incident to the front and to the back of a vertex with a crossing loop.
    cuttlefish::edge_encoding_t e_u_old, e_u_new;   // Edges incident to some particular side of a vertex `u`, before and after the addition of a new edge.
    cuttlefish::edge_encoding_t e_v_old, e_v_new;   // Edges incident to some particular side of a vertex `v`, before and after the addition of a new edge.
*/

    uint64_t edge_count = 0;    // Number of edges processed by this thread.
    uint64_t progress = 0;  // Number of edges processed by the thread; is reset at reaching 1% of its approximate workload.

    while(edge_parser->tasks_expected(thread_id))
        if(edge_parser->value_at(thread_id, e.e()))
        {
            e.configure(hash_table);    // A new edge (k + 1)-mer has been parsed; set information for its two endpoints.

            if(e.is_loop())
                if(e.u().side() != e.v().side())    // It is a crossing loop.
                {
                    while(!add_crossing_loop(e.u()));
/*
                    while(!add_crossing_loop(e.u(), e_front, e_back));
                    
                    propagate_discard(e.u(), e.u().side() == cuttlefish::side_t::front ? e_front : e_back);
                    propagate_discard(e.v(), e.v().side() == cuttlefish::side_t::front ? e_front : e_back);
*/
                }
                else    // A one-sided loop.
                {
                    while(!add_one_sided_loop(e.u()));
/*
                    while(!add_one_sided_loop(e.u(), e_u_old));

                    propagate_discard(e.u(), e_u_old);
*/
                }
            else    // It connects two endpoints `u` and `v` of two distinct vertex.
            {
                while(!add_incident_edge(e.u()));
                while(!add_incident_edge(e.v()));
/*
                while(!add_incident_edge(e.u(), e_u_old, e_u_new));
                while(!add_incident_edge(e.v(), e_v_old, e_v_new));

                if(e_u_new == cuttlefish::edge_encoding_t::N)
                    propagate_discard(e.u(), e.v(), e_u_old);

                if(e_v_new == cuttlefish::edge_encoding_t::N)
                    propagate_discard(e.v(), e.u(), e_v_old);
*/
            }

            edge_count++;


            progress_tracker.track_work(++progress);
        }

    lock.lock();
    // std::cout << "Thread " << thread_id << " processed " << edge_count << " edges.\n";
    edges_processed += edge_count;
    lock.unlock();
}


template <uint16_t k>
uint64_t Read_CdBG_Constructor<k>::vertex_count() const
{
    return hash_table.size();
}


template <uint16_t k>
uint64_t Read_CdBG_Constructor<k>::edge_count() const
{
    return edge_count_;
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG_Constructor)


#include "Read_CdBG_Extractor.hpp"


template <uint16_t k>
Read_CdBG_Extractor<k>::Read_CdBG_Extractor(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table):
    params(params),
    hash_table(hash_table),
    vertex_container(params.vertex_db_path()),
    vertex_parser(&vertex_container, params.thread_count())
{
    std::cout << "Number of distinct vertices: " << vertex_container.size() << ".\n";
}


template <uint16_t k>
void Read_CdBG_Extractor<k>::extract_maximal_unitigs()
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Construct a thread pool.
    const uint16_t thread_count = params.thread_count();
    Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::extract_unipaths_read_space);

    // Launch the reading (and parsing per demand) of the vertices from disk.
    vertex_parser.launch_production();

    // Launch (multi-thread) extraction of the maximal unitigs.
    distribute_unipaths_extraction(thread_pool);

    // Wait for the vertices to be deplted from the database.
    vertex_parser.seize_production();

    // Wait for the consumer threads to finish parsing and processing edges.
    thread_pool.close();

    std::cout << "Number of processed vertices: " << vertices_processed << ".\n";


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done extracting the maximal unitigs. Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
void Read_CdBG_Extractor<k>::distribute_unipaths_extraction(Thread_Pool<k>& thread_pool)
{
    const uint16_t thread_count = params.thread_count();

    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        const uint16_t idle_thread_id = thread_pool.get_idle_thread();
        thread_pool.assign_read_dBG_compaction_task(idle_thread_id);
    }
}


template <uint16_t k>
void Read_CdBG_Extractor<k>::process_vertices(const uint16_t thread_id)
{
    Kmer<k> v;
    uint64_t vertex_count = 0;

    while(vertex_parser.tasks_expected(thread_id))
        if(vertex_parser.value_at(thread_id, v))
        {
            vertex_count++;
        }

    lock.lock();
    std::cout << "Thread " << thread_id << " processed " << vertex_count << " vertices.\n"; // TODO: remove.
    vertices_processed += vertex_count;
    lock.unlock();
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG_Extractor)

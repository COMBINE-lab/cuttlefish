
#include "Read_CdBG_Extractor.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "Character_Buffer.hpp"
#include "Thread_Pool.hpp"


template <uint16_t k>
Read_CdBG_Extractor<k>::Read_CdBG_Extractor(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table):
    params(params),
    hash_table(hash_table)
{}


template <uint16_t k>
void Read_CdBG_Extractor<k>::extract_maximal_unitigs(const std::string& vertex_db_path, const std::string& output_file_path)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Construct a thread pool.
    const uint16_t thread_count = params.thread_count();
    Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::extract_unipaths_read_space);

    // Launch the reading (and parsing per demand) of the vertices from disk.
    const Kmer_Container<k> vertex_container(vertex_db_path);  // Wrapper container for the vertex-database.
    Kmer_SPMC_Iterator<k> vertex_parser(&vertex_container, params.thread_count());  // Parser for the vertices from the vertex-database.
    std::cout << "Number of distinct vertices: " << vertex_container.size() << ".\n";

    vertex_parser.launch_production();

    // Clear the output file and initialize the output sink.
    clear_file(output_file_path);
    init_output_sink(output_file_path);

    // Launch (multi-threaded) extraction of the maximal unitigs.
    const uint64_t thread_load_percentile = static_cast<uint64_t>(std::round((vertex_count() / 100.0) / params.thread_count()));
    progress_tracker.setup(vertex_count(), thread_load_percentile, "Extracting maximal unitigs");
    distribute_unipaths_extraction(&vertex_parser, thread_pool);

    // Wait for the vertices to be depleted from the database.
    vertex_parser.seize_production();

    // Wait for the consumer threads to finish parsing and processing edges.
    thread_pool.close();

    // Close the output sink.
    close_output_sink();

    std::cout << "\nNumber of scanned vertices: " << vertices_scanned << ".\n";
    unipaths_meta_info_.print();


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Extracted the maximal unitigs. Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
void Read_CdBG_Extractor<k>::distribute_unipaths_extraction(Kmer_SPMC_Iterator<k>* const vertex_parser, Thread_Pool<k>& thread_pool)
{
    const uint16_t thread_count = params.thread_count();

    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        const uint16_t idle_thread_id = thread_pool.get_idle_thread();
        thread_pool.assign_read_dBG_compaction_task(vertex_parser, idle_thread_id);
    }
}


template <uint16_t k>
void Read_CdBG_Extractor<k>::process_vertices(Kmer_SPMC_Iterator<k>* const vertex_parser, const uint16_t thread_id)
{
    // Data structures to be reused per each vertex scanned.
    Kmer<k> v_hat;  // The vertex copy to be scanned one-by-one.
    Maximal_Unitig_Scratch<k> maximal_unitig;  // The scratch space to be used to construct the containing maximal unitig of `v_hat`.

    uint64_t vertex_count = 0;  // Number of vertices scanned by this thread.
    Unipaths_Meta_info<k> extracted_unipaths_info;  // Meta-information over the maximal unitigs extracted by this thread.
    uint64_t progress = 0;  // Number of vertices scanned by the thread; is reset at reaching 1% of its approximate workload.

    Character_Buffer<BUFF_SZ, sink_t> output_buffer(output_sink.sink());  // The output buffer for maximal unitigs.


    while(vertex_parser->tasks_expected(thread_id))
        if(vertex_parser->value_at(thread_id, v_hat))
        {
            if(extract_maximal_unitig(v_hat, maximal_unitig))
            {
                mark_maximal_unitig(maximal_unitig);

                extracted_unipaths_info.add_maximal_unitig(maximal_unitig);
                // output_buffer += maximal_unitig.fasta_rec();
                maximal_unitig.add_fasta_rec_to_buffer(output_buffer);
                
                if(progress_tracker.track_work(progress += maximal_unitig.size()))
                    progress = 0;
            }

            vertex_count++;
        }


    // Aggregate the meta-information over the extracted maximal unitigs and the thread-executions.
    lock.lock();

    vertices_scanned += vertex_count;
    unipaths_meta_info_.aggregate(extracted_unipaths_info);

    lock.unlock();
}


template <uint16_t k>
void Read_CdBG_Extractor<k>::init_output_sink(const std::string& output_file_path)
{
    output_sink.init_sink(output_file_path);
}


template <uint16_t k>
void Read_CdBG_Extractor<k>::close_output_sink()
{
    output_sink.close_sink();
}


template <uint16_t k>
const Build_Params& Read_CdBG_Extractor<k>::get_params() const
{
    return params;
}


template <uint16_t k>
const Unipaths_Meta_info<k>& Read_CdBG_Extractor<k>::unipaths_meta_info() const
{
    return unipaths_meta_info_;
}


template <uint16_t k>
uint64_t Read_CdBG_Extractor<k>::vertex_count() const
{
    return hash_table.size();
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG_Extractor)


#include "Read_CdBG_Extractor.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "FASTA_Record.hpp"
#include "Character_Buffer.hpp"
#include "utility.hpp"
#include "Thread_Pool.hpp"


template <uint16_t k>
Read_CdBG_Extractor<k>::Read_CdBG_Extractor(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table):
    params(params),
    hash_table(hash_table)
{}


template <uint16_t k>
void Read_CdBG_Extractor<k>::extract_maximal_unitigs(const std::string& vertex_db_path)
{
    // std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Construct a thread pool.
    const uint16_t thread_count = params.thread_count();
    Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::extract_unipaths_read_space);

    // Launch the reading (and parsing per demand) of the vertices from disk.
    const Kmer_Container<k> vertex_container(vertex_db_path);  // Wrapper container for the vertex-database.
    Kmer_SPMC_Iterator<k> vertex_parser(&vertex_container, params.thread_count());  // Parser for the vertices from the vertex-database.
    std::cout << "Number of distinct vertices: " << vertex_container.size() << ".\n";

    vertex_parser.launch_production();

    // Clear the output file and initialize the output sink.
    clear_file(params.output_file_path());
    init_output_sink();

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

    // Check for the existence of cycle(s).
    if(has_dcc())
        std::cout <<    "\nCycles disconnected from the rest of the graph are present."
                        " I.e. the cycles are graph components exclusively on their own.\n\n";


    // std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    // double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    // std::cout << "Done extracting the maximal unitigs. Time taken = " << elapsed_seconds << " seconds.\n";
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
void Read_CdBG_Extractor<k>::scan_vertices(Kmer_SPMC_Iterator<k>* const vertex_parser, const uint16_t thread_id)
{
    // Data structures to be reused per each vertex scanned.
    Kmer<k> v;  // The vertex copy to be scanned one-by-one.
    cuttlefish::side_t s_v; // The side of the vertex `v` that connects it to the maximal unitig `p` containing it, if `v` is flanking.
    State_Read_Space state; // State of the vertex `v`.
    uint64_t id;    // The unique ID of the maximal unitig `p`.
    std::vector<char> unipath;  // The extracted maximal unitig `p`.
    std::vector<uint64_t> path_hashes;  // Hash values of the vertices constituting the maximal unitig `p`.

    uint64_t vertex_count = 0;  // Number of vertices scanned by this thread.
    Unipaths_Meta_info<k> extracted_unipaths_info;  // Meta-information over the maximal unitigs extracted by this thread.
    uint64_t progress = 0;  // Number of vertices scanned by the thread; is reset at reaching 1% of its approximate workload.

    Character_Buffer<BUFF_SZ, sink_t> output_buffer(output_sink.sink());  // The output buffer for maximal unitigs.
    unipath.reserve(SEQ_SZ);

    const bool mark_unipaths = params.extract_cycles() || params.dcc_opt();
    if(mark_unipaths)
        path_hashes.reserve(BUFF_SZ);


    while(vertex_parser->tasks_expected(thread_id))
        if(vertex_parser->value_at(thread_id, v))
        {
            state = hash_table[v].state();
            
            if(!state.is_outputted() && is_flanking_state(state, s_v))
                if(extract_maximal_unitig(v, s_v, id, unipath, path_hashes))
                {
                    extracted_unipaths_info.add_maximal_unitig(unipath);

                    unipath.emplace_back('\n');
                    // output_buffer += unipath;
                    output_buffer += FASTA_Record<std::vector<char>>(id, unipath);
                    // unipath.clear();

                    if(mark_unipaths)
                        mark_path(path_hashes);
                }

            vertex_count++;


            progress_tracker.track_work(++progress);
        }


    // Aggregate the meta-information over the extracted maximal unitigs and the thread-executions.
    lock.lock();
    // std::cout << "Thread " << thread_id << " processed " << vertex_count << " vertices.\n";

    vertices_scanned += vertex_count;
    unipaths_meta_info_.aggregate(extracted_unipaths_info);

    lock.unlock();
}


template <uint16_t k>
bool Read_CdBG_Extractor<k>::extract_maximal_unitig(const Kmer<k>& v_hat, const cuttlefish::side_t s_v_hat, uint64_t& id, std::vector<char>& unipath, std::vector<uint64_t>& path_hashes)
{
    // Data structures to be reused per each vertex extension of the maximal unitig.
    cuttlefish::side_t s_v = s_v_hat;   // The side of the current vertex `v` through which to extend the maximal unitig, i.e. exit `v`.
    Directed_Vertex<k> v(s_v == cuttlefish::side_t::back ? v_hat : v_hat.reverse_complement(), hash_table); // Current vertex being added to the maximal unitig.
    State_Read_Space state = hash_table[v.hash()].state();  // State of the vertex `v`.
    cuttlefish::edge_encoding_t e_v;    // The next edge from `v` to include into the maximal unitig.
    cuttlefish::base_t b_ext;   // The nucleobase corresponding to the edge `e_v` and the exiting side `s_v` from `v` to add to the literal maximal unitig.
    const bool mark_unipaths = params.extract_cycles() || params.dcc_opt(); // Whether to mark the vertices present in the maximal unitigs.

    const Directed_Vertex<k> init_vertex(v);
    init_vertex.kmer().get_label(unipath);
    if(mark_unipaths)
    {
        path_hashes.clear();
        path_hashes.emplace_back(init_vertex.hash());
    }


    while(true)
    {
        if(state.is_outputted())    // The opposite end of the maximal unitig has been reached, and the unitig is found to have already been outputted.
            return false;

        if(is_flanking_side(state, s_v))
            break;


        e_v = state.edge_at(s_v);
        b_ext = (s_v == cuttlefish::side_t::back ? DNA_Utility::map_base(e_v) : DNA_Utility::complement(DNA_Utility::map_base(e_v)));

        v.roll_forward(b_ext, hash_table);
        s_v = v.exit_side();
        state = hash_table[v.hash()].state();
        
        unipath.emplace_back(Kmer<k>::map_char(b_ext));
        if(mark_unipaths)
            path_hashes.emplace_back(v.hash());
            // TODO: write-out to disk in case of the size crossing some threshold, and modify `mark_path` accordingly —
            // would prevent unwanted memory blow-up in presence of very large maximal unitigs.
    }

    const Directed_Vertex<k>& term_vertex = v;
    const bool in_canonical = (init_vertex.kmer() < term_vertex.kmer_bar());
    const Directed_Vertex<k>& sign_vertex = (in_canonical ? init_vertex : term_vertex);
    const Directed_Vertex<k>& cosign_vertex = (in_canonical ? term_vertex : init_vertex);

    // Mark the flanking vertices as outputted.
    if(!mark_flanking_vertices(sign_vertex, cosign_vertex))
        return false;

    if(!in_canonical)
        reverse_complement(unipath);

    id = sign_vertex.hash();

    return true;
}


template <uint16_t k>
void Read_CdBG_Extractor<k>::init_output_sink()
{
    output_sink.init_sink(params.output_file_path());
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


template <uint16_t k>
uint64_t Read_CdBG_Extractor<k>::unipaths_vertex_count() const
{
    return unipaths_meta_info_.kmer_count();
}


template <uint16_t k>
bool Read_CdBG_Extractor<k>::has_dcc() const
{
    return unipaths_vertex_count() != vertex_count();
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG_Extractor)

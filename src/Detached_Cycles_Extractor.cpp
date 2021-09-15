
#include "Read_CdBG_Extractor.hpp"
#include "Read_CdBG.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "FASTA_Record.hpp"
#include "Character_Buffer.hpp"
#include "Thread_Pool.hpp"
#include "dBG_Info.hpp"


template <uint16_t k>
void Read_CdBG_Extractor<k>::extract_detached_cycles(const std::string& vertex_db_path, const std::string& output_file_path, const dBG_Info<k>& dbg_info)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    if(!dbg_info.dcc_opt_performed())
    {
        std::cout << "Marking the vertices present in the extracted maximal unitigs.\n";
        mark_maximal_unitig_vertices(vertex_db_path);
        std::cout << "Done marking the vertices.\n";
    }

    std::cout << "Extracting the cycles.\n";
    extract_detached_chordless_cycles(vertex_db_path, output_file_path);

    std::cout <<    "\nNumber of detached chordless cycles: " << unipaths_meta_info_.dcc_count() << ".\n"
                    "Number of vertices in the cycles: " << unipaths_meta_info_.dcc_kmer_count() << ".\n";


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done extracting the cycles. Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
void Read_CdBG_Extractor<k>::mark_maximal_unitig_vertices(const std::string& vertex_db_path)
{
    // Construct a thread pool.
    const uint16_t thread_count = params.thread_count();
    Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::mark_unipath_vertices);

    // Launch the reading (and parsing per demand) of the vertices from disk.
    const Kmer_Container<k> vertex_container(vertex_db_path);   // Wrapper container for the vertex-database.
    Kmer_SPMC_Iterator<k> vertex_parser(&vertex_container, params.thread_count());  // Parser for the vertices from the vertex-database.
    std::cout << "Number of distinct vertices: " << vertex_container.size() << ".\n";

    vertex_parser.launch_production();


    // Launch (multi-threaded) marking of the vertices present in the maximal unitigs.
    const uint64_t thread_load_percentile = static_cast<uint64_t>(std::round((vertex_count() / 100.0) / params.thread_count()));
    progress_tracker.setup(vertex_count(), thread_load_percentile, "Marking maximal unitigs");
    distribute_unipaths_extraction(&vertex_parser, thread_pool);

    // Wait for the vertices to be depleted from the database.
    vertex_parser.seize_production();

    // Wait for the consumer threads to finish parsing and processing edges.
    thread_pool.close();

    std::cout << "Number of scanned vertices: " << vertices_scanned << ".\n";
    std::cout << "Number of marked vertices:  " << vertices_marked << ".\n";
}



template <uint16_t k>
void Read_CdBG_Extractor<k>::mark_maximal_unitig_vertices(Kmer_SPMC_Iterator<k>* const vertex_parser, const uint16_t thread_id)
{
    // Data structures to be reused per each vertex scanned.
    Kmer<k> v;  // The vertex copy to be scanned one-by-one.
    cuttlefish::side_t s_v; // The side of the vertex `v` that connects it to the maximal unitig containing it, if `v` is flanking.
    State_Read_Space state; // State of the vertex `v`.

    uint64_t vertex_count = 0;  // Number of vertices scanned by this thread.
    uint64_t marked_count = 0;  // Number of vertices marked as present in maximal unitigs by this thread.
    uint64_t progress = 0;  // Number of vertices scanned by the thread; is reset at reaching 1% of its approximate workload.

    while(vertex_parser->tasks_expected(thread_id))
        if(vertex_parser->value_at(thread_id, v))
        {
            state = hash_table[v].state();

            if(!state.is_outputted() && is_flanking_state(state, s_v))
                marked_count += mark_maximal_unitig(v, s_v);

            vertex_count++;


            progress_tracker.track_work(++progress);
        }


    // Aggregate the meta-information over the marked maximal unitigs and the thread-executions.
    lock.lock();

    std::cout << "Thread " << thread_id << " processed " << vertex_count << " vertices," // TODO: remove.
                " and marked " << marked_count << " vertices.\n";
    
    vertices_scanned += vertex_count;
    vertices_marked += marked_count;

    lock.unlock();
}


template <uint16_t k>
std::size_t Read_CdBG_Extractor<k>::mark_maximal_unitig(const Kmer<k>& v_hat, const cuttlefish::side_t s_v_hat)
{
    // Data structures to be reused per each vertex extension of the maximal unitig.
    cuttlefish::side_t s_v = s_v_hat;   // The side of the current vertex `v` through which to extend the maximal unitig, i.e. exit `v`.
    Directed_Vertex<k> v(s_v == cuttlefish::side_t::back ? v_hat : v_hat.reverse_complement(), hash_table); // Current vertex being added to the maximal unitig.
    State_Read_Space state = hash_table[v.hash()].state();  // State of the vertex `v`.
    cuttlefish::edge_encoding_t e_v;    // The next edge from `v` to include into the maximal unitig.
    cuttlefish::base_t b_ext;   // The nucleobase corresponding to the edge `e_v` and the exiting side `s_v` from `v` to add to the literal maximal unitig.

    std::size_t marked_count = 0;   // Number of vertices successfully marked by this thread.

    while(true)
    {
        if(!mark_vertex(v))
            break;

        marked_count++;

        if(is_flanking_side(state, s_v))
            break;
            
        e_v = state.edge_at(s_v);
        b_ext = (s_v == cuttlefish::side_t::back ? DNA_Utility::map_base(e_v) : DNA_Utility::complement(DNA_Utility::map_base(e_v)));

        v.roll_forward(b_ext, hash_table);
        s_v = v.exit_side();
        state = hash_table[v.hash()].state();
    }


    return marked_count;
}


template <uint16_t k>
void Read_CdBG_Extractor<k>::extract_detached_chordless_cycles(const std::string& vertex_db_path, const std::string& output_file_path)
{
     // Construct a thread pool.
    const uint16_t thread_count = params.thread_count();
    Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::extract_cycles);

    // Launch the reading (and parsing per demand) of the vertices from disk.
    const Kmer_Container<k> vertex_container(vertex_db_path);   // Wrapper container for the vertex-database.
    Kmer_SPMC_Iterator<k> vertex_parser(&vertex_container, params.thread_count());  // Parser for the vertices from the vertex-database.
    std::cout << "Number of distinct vertices: " << vertex_container.size() << ".\n";

    vertex_parser.launch_production();

    // Initialize the output sink.
    init_output_sink(output_file_path);

    // Launch (multi-threaded) marking of the vertices present in the maximal unitigs.
    const uint64_t thread_load_percentile = static_cast<uint64_t>(std::round((vertex_count() / 100.0) / params.thread_count()));
    progress_tracker.setup(vertex_count(), thread_load_percentile, "Extracting detached chordless cycles");
    distribute_unipaths_extraction(&vertex_parser, thread_pool);

    // Wait for the vertices to be depleted from the database.
    vertex_parser.seize_production();

    // Wait for the consumer threads to finish parsing and processing edges.
    thread_pool.close();

    // Close the output sink.
    close_output_sink();
}


template <uint16_t k>
void Read_CdBG_Extractor<k>::extract_detached_chordless_cycles(Kmer_SPMC_Iterator<k>* vertex_parser, const uint16_t thread_id)
{
    // Data structures to be reused per each vertex scanned.
    Kmer<k> v;  // The vertex copy to be scanned one-by-one.
    State_Read_Space state; // State of the vertex `v`.
    uint64_t id;    // The unique ID of the cycle.
    std::vector<char> cycle;  // The extracted cycle.
    std::size_t pivot;  // Index of the lexicographically lowest (canonical) k-mer in the cycle.

    uint64_t vertex_count = 0;  // Number of vertices scanned by this thread.
    uint64_t cycles_extracted = 0;  // Number of detached chordless cycles extracted by this thread.
    uint64_t cycle_vertices = 0;    // Number of vertices found to be in detached chordless cycles by this thread.
    uint64_t progress = 0;  // Number of vertices scanned by the thread; is reset at reaching 1% of its approximate workload.

    Character_Buffer<BUFF_SZ, sink_t> output_buffer(output_sink.sink());  // The output buffer for the cycles.
    cycle.reserve(SEQ_SZ);

    while(vertex_parser->tasks_expected(thread_id))
        if(vertex_parser->value_at(thread_id, v))
        {
            state = hash_table[v].state();

            if(!state.is_outputted())
                if(extract_cycle(v, id, cycle, pivot))
                {
                    cycles_extracted++;
                    cycle_vertices += cycle.size() - (k - 1);
                    unipaths_meta_info_.add_DCC(cycle);

                    // cycle.emplace_back('\n');
                    // output_buffer += FASTA_Record<std::vector<char>>(id, cycle);
                    output_buffer.rotate_append(FASTA_Record<std::vector<char>>(id, cycle), pivot);
                }

            vertex_count++;


            progress_tracker.track_work(++progress);
        }


    // Aggregate the meta-information over the marked maximal unitigs and the thread-executions.
    lock.lock();
    
    vertices_scanned += vertex_count;

    lock.unlock();
}


template <uint16_t k>
bool Read_CdBG_Extractor<k>::extract_cycle(const Kmer<k>& v_hat, uint64_t& id, std::vector<char>& cycle, std::size_t& pivot)
{
    // Data structures to be reused per each vertex extension of the cycle.
    cuttlefish::side_t s_v = cuttlefish::side_t::back;  // The side of the current vertex `v` through which to extend the cycle, i.e. exit `v`.
    Directed_Vertex<k> v(v_hat, hash_table);    // Current vertex being added to the cycle.
    State_Read_Space state = hash_table[v.hash()].state();  // State of the vertex `v`.
    cuttlefish::edge_encoding_t e_v;    // The next edge from `v` to include into the cycle.
    cuttlefish::base_t b_ext;   // The nucleobase corresponding to the edge `e_v` and the exiting side `s_v` from `v` to add to the literal cycle.
    std::size_t kmer_idx;   // Index of the vertex in the path being traversed.


    const Directed_Vertex<k> anchor(v);
    Directed_Vertex<k> sign_vertex(anchor);
    pivot = kmer_idx = 0;

    anchor.kmer().get_label(cycle);

    while(true)
    {
        if(state.is_outputted())    // The cycle is found to have already been outputted.
            return false;


        e_v = state.edge_at(s_v);
        b_ext = (s_v == cuttlefish::side_t::back ? DNA_Utility::map_base(e_v) : DNA_Utility::complement(DNA_Utility::map_base(e_v)));

        v.roll_forward(b_ext, hash_table);
        if(v.is_same_vertex(anchor))
            break;

        s_v = v.exit_side();
        state = hash_table[v.hash()].state();

        kmer_idx++;
        if(sign_vertex.canonical() > v.canonical())
            sign_vertex = v,
            pivot = kmer_idx;
        
        cycle.emplace_back(Kmer<k>::map_char(b_ext));
    }

    
    if(!mark_vertex(sign_vertex))
        return false;

    if(!sign_vertex.in_canonical_form())
    {
        reverse_complement(cycle);
        pivot = (cycle.size() - 1) - pivot - (k - 1);
    }

    id = sign_vertex.hash();

    return true;
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG_Extractor)

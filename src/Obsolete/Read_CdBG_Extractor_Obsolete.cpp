
#include "Read_CdBG_Extractor.hpp"
#include "Kmer_SPMC_Iterator.hpp"


// The following methods are not used anymore with the current algorithm.
/*
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

                    // unipath.emplace_back('\n');
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
        cuttlefish::reverse_complement(unipath);

    id = sign_vertex.hash();

    return true;
}


template <uint16_t k>
bool Read_CdBG_Extractor<k>::mark_flanking_vertices(const Directed_Vertex<k>& sign_vertex, const Directed_Vertex<k>& cosign_vertex)
{
    return mark_vertex(sign_vertex) && (sign_vertex.hash() == cosign_vertex.hash() || mark_vertex(cosign_vertex));
}


template <uint16_t k>
bool Read_CdBG_Extractor<k>::is_flanking_state(const State_Read_Space state, cuttlefish::side_t& unipath_side)
{
    if(is_flanking_side(state, cuttlefish::side_t::front))
    {
        unipath_side = cuttlefish::side_t::back;
        return true;
    }

    if(is_flanking_side(state, cuttlefish::side_t::back))
    {
        unipath_side = cuttlefish::side_t::front;
        return true;
    }

    return false;
}


template <uint16_t k>
bool Read_CdBG_Extractor<k>::is_flanking_side(const State_Read_Space state, const cuttlefish::side_t side)
{
    const cuttlefish::edge_encoding_t edge = state.edge_at(side);

    return edge == cuttlefish::edge_encoding_t::N || edge == cuttlefish::edge_encoding_t::E;
}


template <uint16_t k>
bool Read_CdBG_Extractor<k>::has_dcc() const
{
    return unipaths_vertex_count() != vertex_count();
}


template <uint16_t k>
uint64_t Read_CdBG_Extractor<k>::unipaths_vertex_count() const
{
    return unipaths_meta_info_.kmer_count();
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG_Extractor)
*/

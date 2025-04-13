
#include "Read_CdBG.hpp"
#include "kmer_Enumerator.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "kmer_Enumeration_Stats.hpp"
#include "Read_CdBG_Constructor.hpp"
#include "Read_CdBG_Extractor.hpp"
#include "Kmer_Index.hpp"
#include "kmc_runner.h"

#include <limits>


template <uint16_t k>
Read_CdBG<k>::Read_CdBG(const Build_Params& params):
    Read_CdBG(params, nullptr)
{}


template <uint16_t k>
Read_CdBG<k>::Read_CdBG(const Build_Params& params, Kmer_Index<k>* const kmer_idx):
    params(params),
    logistics(this->params),
    hash_table(nullptr),
    dbg_info(params.json_file_path()),
    kmer_idx(kmer_idx)
{}


template <uint16_t k>
Read_CdBG<k>::~Read_CdBG()
{
    if(hash_table != nullptr)
        hash_table->clear();

    dbg_info.dump_info();
}


template <uint16_t k>
void Read_CdBG<k>::construct()
{
    if(is_constructed())
    {
        std::cout << "\nThe compacted de Bruijn graph has been constructed earlier. Check " << dbg_info.file_path() << " for results.\n";
        return;
    }


    dbg_info.add_build_params(params);

    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

#ifdef CF_DEVELOP_MODE

    uint64_t edge_count;
    uint64_t vertex_count;

    if(params.edge_db_path().empty())
    {
        kmer_Enumeration_Stats<k + 1> edge_stats = enumerate_edges();
        kmer_Enumeration_Stats<k> vertex_stats = enumerate_vertices(edge_stats.max_memory());
        
        edge_count = edge_stats.counted_kmer_count();
        vertex_count = vertex_stats.counted_kmer_count();
    }
    else if(!params.vertex_db_path().empty())
    {
        edge_count = Kmer_Container<k + 1>::size(params.edge_db_path());
        vertex_count = Kmer_Container<k>::size(params.vertex_db_path());
    }
    else
    {
        std::cerr << "Vertex database must also be provided if edge database is passed. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    std::chrono::high_resolution_clock::time_point t_vertices = std::chrono::high_resolution_clock::now();
    std::cout << "Enumerated the edge and the vertex set of the graph. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_vertices - t_start).count() << " seconds.\n";
#else

    std::cout << "\nEnumerating the edges of the de Bruijn graph.\n";
    kmer_Enumeration_Stats<k + 1> edge_stats = enumerate_edges();
    edge_stats.log_stats();

    std::chrono::high_resolution_clock::time_point t_edges = std::chrono::high_resolution_clock::now();
    std::cout << "Enumerated the edge set of the graph. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_edges - t_start).count() << " seconds.\n";


    std::cout << "\nEnumerating the vertices of the de Bruijn graph.\n";
    kmer_Enumeration_Stats<k> vertex_stats = enumerate_vertices(edge_stats.max_memory());

    std::chrono::high_resolution_clock::time_point t_vertices = std::chrono::high_resolution_clock::now();
    std::cout << "Enumerated the vertex set of the graph. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_vertices - t_edges).count() << " seconds.\n";

    const uint64_t edge_count = edge_stats.counted_kmer_count();
    const uint64_t vertex_count = vertex_stats.counted_kmer_count();
#endif
    std::cout << "Number of edges:    " << edge_count << ".\n";
    std::cout << "Number of vertices: " << vertex_count << ".\n";


    std::cout << "\nConstructing the minimal perfect hash function (MPHF) over the vertex set.\n";
    construct_hash_table(vertex_count);

    std::chrono::high_resolution_clock::time_point t_mphf = std::chrono::high_resolution_clock::now();
    std::cout << "Constructed the minimal perfect hash function for the vertices. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_mphf - t_vertices).count() << " seconds.\n";


    std::cout << "\nComputing the DFA states.\n";
    compute_DFA_states();

#ifdef CF_DEVELOP_MODE
    if(params.edge_db_path().empty())
#endif
    Kmer_Container<k + 1>::remove(logistics.edge_db_path());
    
    std::chrono::high_resolution_clock::time_point t_dfa = std::chrono::high_resolution_clock::now();
    std::cout << "Computed the states of the automata. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_dfa - t_mphf).count() << " seconds.\n";


    std::cout << "\nExtracting " << (params.path_cover() ? "a maximal path cover" :  "the maximal unitigs") << ".\n";
    extract_maximal_unitigs();

#ifdef CF_DEVELOP_MODE
    if(params.vertex_db_path().empty())
#endif
    if(!params.save_vertices())
        Kmer_Container<k>::remove(logistics.vertex_db_path());

    std::chrono::high_resolution_clock::time_point t_extract = std::chrono::high_resolution_clock::now();
    std::cout << "Extracted the paths. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_extract - t_dfa).count() << " seconds.\n";


#ifndef CF_DEVELOP_MODE
    const double max_disk = static_cast<double>(max_disk_usage(edge_stats, vertex_stats)) / (1024.0 * 1024.0 * 1024.0);
    std::cout << "\nMaximum temporary disk-usage: " << max_disk << "GB.\n";
#endif
}


template <uint16_t k>
kmer_Enumeration_Stats<k + 1> Read_CdBG<k>::enumerate_edges() const
{
    const KMC::InputFileType ip_type = (params.is_read_graph() ? KMC::InputFileType::FASTQ : KMC::InputFileType::MULTILINE_FASTA);
    return kmer_Enumerator<k + 1>().enumerate(
        ip_type, logistics.input_paths_collection(), params.cutoff(), params.thread_count(),
        params.max_memory(), params.strict_memory(), params.strict_memory(), bits_per_vertex,
        logistics.working_dir_path(), logistics.edge_db_path());
}


template <uint16_t k>
kmer_Enumeration_Stats<k> Read_CdBG<k>::enumerate_vertices(const std::size_t max_memory) const
{
    return kmer_Enumerator<k>().enumerate(
        KMC::InputFileType::KMC, std::vector<std::string>(1, logistics.edge_db_path()), 1, params.thread_count(),
        max_memory, params.strict_memory(), false, bits_per_vertex,
        logistics.working_dir_path(), logistics.vertex_db_path());
}


template <uint16_t k>
void Read_CdBG<k>::construct_hash_table(const uint64_t vertex_count, const bool load)
{
    if(load)
    {
        hash_table = std::make_unique<Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>>(logistics.vertex_db_path(), vertex_count);
        hash_table->load(params);    
    }
    else
    {
        std::size_t max_memory = std::max(process_peak_memory(), params.max_memory() * 1024U * 1024U * 1024U);
        const std::size_t parser_memory = Kmer_SPMC_Iterator<k>::memory(params.thread_count());
        max_memory = (max_memory > parser_memory ? max_memory - parser_memory : 0);
        
        hash_table =
#ifdef CF_DEVELOP_MODE
                        std::make_unique<Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>>(logistics.vertex_db_path(), vertex_count, max_memory, params.gamma());
#else
                        (params.strict_memory() ?
                            std::make_unique<Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>>(logistics.vertex_db_path(), vertex_count, max_memory) :
                            std::make_unique<Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>>(logistics.vertex_db_path(), vertex_count, max_memory, std::numeric_limits<double>::max()));
#endif
        hash_table->construct(params.thread_count(), logistics.working_dir_path(), params.mph_file_path(), params.save_mph());
    }
}


template <uint16_t k>
void Read_CdBG<k>::compute_DFA_states()
{
    Read_CdBG_Constructor<k> cdBg_constructor(params, *hash_table);
    cdBg_constructor.compute_DFA_states(logistics.edge_db_path());

    dbg_info.add_basic_info(cdBg_constructor);
}


template <uint16_t k>
void Read_CdBG<k>::extract_maximal_unitigs()
{
    Read_CdBG_Extractor<k> cdBg_extractor(params, *hash_table, kmer_idx);

    cdBg_extractor.extract_maximal_unitigs(logistics.vertex_db_path(), logistics.output_file_path());
    dbg_info.add_unipaths_info(cdBg_extractor);
}


template <uint16_t k>
bool Read_CdBG<k>::is_constructed() const
{
    return file_exists(params.json_file_path());
}


template <uint16_t k>
std::size_t Read_CdBG<k>::max_disk_usage(const kmer_Enumeration_Stats<k + 1>& edge_stats, const kmer_Enumeration_Stats<k>& vertex_stats)
{
    const std::size_t at_edge_enum = std::max(edge_stats.temp_disk_usage(), edge_stats.db_size());
    const std::size_t at_vertex_enum = edge_stats.db_size() + std::max(vertex_stats.temp_disk_usage(), vertex_stats.db_size());

    const std::size_t max_disk = std::max(at_edge_enum, at_vertex_enum);
    return max_disk;
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG)

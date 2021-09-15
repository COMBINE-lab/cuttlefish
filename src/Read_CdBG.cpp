
#include "Read_CdBG.hpp"
#include "kmer_Enumerator.hpp"
#include "Read_CdBG_Constructor.hpp"
#include "Read_CdBG_Extractor.hpp"
#include "File_Extensions.hpp"
#include "utility.hpp"
#include "kmc_runner.h"

#include <iomanip>


template <uint16_t k>
Read_CdBG<k>::Read_CdBG(const Build_Params& params):
    params(params),
    dbg_info(params.json_file_path())
{}


template <uint16_t k>
void Read_CdBG<k>::construct()
{
    if(is_constructed(params) && (!dbg_info.has_dcc() || dbg_info.dcc_extracted() || !params.extract_cycles()))
    {
        std::cout << "\nThe compacted de Bruijn graph has already been constructed earlier.\n";
        if(dbg_info.has_dcc() && !dbg_info.dcc_extracted())
            std::cout << "There are Detached Chordless Cycles (DCC) present in the graph; run Cuttlefish with the `cycles` argument to extract those.\n";

        return;
    }

    dbg_info.add_build_params(params);

    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    std::cout << "\nEnumerating the edges of the de Bruijn graph.\n";
    const std::string edge_db_path = params.output_prefix() + cuttlefish::file_ext::edges_ext;
    kmer_Enumeration_Stats edge_stats = enumerate_edges(edge_db_path);

    std::chrono::high_resolution_clock::time_point t_edges = std::chrono::high_resolution_clock::now();
    std::cout << "Enumerated the edge set of the graph. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_edges - t_start).count() << " seconds.\n";


    std::cout << "\nEnumerating the vertices of the de Bruijn graph.\n";
    const std::string vertex_db_path = params.output_prefix() + cuttlefish::file_ext::vertices_ext;
    kmer_Enumeration_Stats vertex_stats = enumerate_vertices(edge_db_path, vertex_db_path, edge_stats.max_memory());

    std::chrono::high_resolution_clock::time_point t_vertices = std::chrono::high_resolution_clock::now();
    std::cout << "Enumerated the vertex set of the graph. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_vertices - t_edges).count() << " seconds.\n";

    std::cout << "Number of edges:    " << edge_stats.kmer_count() << ".\n";
    std::cout << "Number of vertices: " << vertex_stats.kmer_count() << ".\n";


    std::cout << "\nConstructing the minimal perfect hash function (MPHF) over the vertex set.\n";
    construct_hash_table(vertex_db_path, vertex_stats.kmer_count());

    std::chrono::high_resolution_clock::time_point t_mphf = std::chrono::high_resolution_clock::now();
    std::cout << "Constructed the minimal perfect hash function for the vertices. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_mphf - t_vertices).count() << " seconds.\n";


    std::cout << "\nComputing the DFA states.\n";
    compute_DFA_states(edge_db_path);

    if(!params.extract_cycles() && !params.dcc_opt())
        hash_table->save(params);
    
    std::chrono::high_resolution_clock::time_point t_dfa = std::chrono::high_resolution_clock::now();
    std::cout << "Computed the states of the automata. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_dfa - t_mphf).count() << " seconds.\n";


    std::cout << "\nExtracting the maximal unitigs.\n";
    extract_maximal_unitigs(vertex_db_path);

    if(!dbg_info.has_dcc() || dbg_info.dcc_extracted())
        hash_table->remove(params);

    std::chrono::high_resolution_clock::time_point t_extract = std::chrono::high_resolution_clock::now();
    std::cout << "Extracted the maximal unitigs. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_extract - t_dfa).count() << " seconds.\n";


    hash_table->clear();
    dbg_info.dump_info();
}


template <uint16_t k>
kmer_Enumeration_Stats Read_CdBG<k>::enumerate_edges(const std::string& edge_db_path)
{
    return kmer_Enumerator<k + 1>().enumerate(
        KMC::InputFileType::FASTQ, params.sequence_input().seqs(), params.cutoff(),
        params.thread_count(), params.max_memory(), params.strict_memory(), true,
        params.working_dir_path(), edge_db_path);
}


template <uint16_t k>
kmer_Enumeration_Stats Read_CdBG<k>::enumerate_vertices(const std::string& edge_db_path, const std::string& vertex_db_path, const std::size_t max_memory)
{
    return kmer_Enumerator<k>().enumerate(
        KMC::InputFileType::KMC, std::vector<std::string>(1, edge_db_path), 1,
        params.thread_count(), max_memory, params.strict_memory(), false,
        params.working_dir_path(), vertex_db_path);
}


template <uint16_t k>
void Read_CdBG<k>::construct_hash_table(const std::string& vertex_db_path, const uint64_t vertex_count)
{
    hash_table = std::make_unique<Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>>(vertex_db_path, vertex_count);
    hash_table->construct(params.thread_count(), params.working_dir_path(), params.mph_file_path());
}


template <uint16_t k>
void Read_CdBG<k>::compute_DFA_states(const std::string& edge_db_path)
{
    Read_CdBG_Constructor<k> cdBg_constructor(params, *hash_table);
    cdBg_constructor.compute_DFA_states(edge_db_path);

    dbg_info.add_basic_info(cdBg_constructor);
}


template <uint16_t k>
void Read_CdBG<k>::extract_maximal_unitigs(const std::string& vertex_db_path)
{
    Read_CdBG_Extractor<k> cdBg_extractor(params, *hash_table);
    const std::string temp_output_path = params.working_dir_path() + filename(params.output_prefix()) + cuttlefish::file_ext::temp;
    const std::string output_file_path = params.output_file_path();

    if(!is_constructed(params))
    {
        cdBg_extractor.extract_maximal_unitigs(vertex_db_path, temp_output_path);
        
        dbg_info.add_unipaths_info(cdBg_extractor);

        if(cdBg_extractor.has_dcc())
        {
            if(params.extract_cycles())
            {
                cdBg_extractor.extract_detached_cycles(vertex_db_path, temp_output_path, dbg_info);

                dbg_info.add_DCC_info(cdBg_extractor);
            }
            else if(params.dcc_opt())
                hash_table->save(params);
        }

        move_file(temp_output_path, output_file_path);
    }
    else if(params.extract_cycles())
    {
        if(dbg_info.has_dcc())
        {
            if(!dbg_info.dcc_extracted())
            {
                cdBg_extractor.extract_detached_cycles(vertex_db_path, output_file_path, dbg_info);

                dbg_info.add_DCC_info(cdBg_extractor);
            }
            else
                std::cout << "\nThe DCCs (Detached Chordless Cycles) have already been extracted earlier.\n";
        }
        else
            std::cout << "\nThe de Bruijn graph has no DCCs (Detached Chordless Cycles).\n";
    }
    else
        std::cout << "\nNothing to do.\n";
}


template <uint16_t k>
bool Read_CdBG<k>::is_constructed(const Build_Params& params)
{
    return file_exists(params.json_file_path());
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG)

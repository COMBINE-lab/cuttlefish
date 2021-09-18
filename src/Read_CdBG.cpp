
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
    hash_table(nullptr),
    dbg_info(params.json_file_path())
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
        std::cout << "\nThe compacted de Bruijn graph has been constructed earlier.\n";
        extract_DCCs(params.output_file_path(), true);
        return;
    }


    dbg_info.add_build_params(params);

    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    std::cout << "\nEnumerating the edges of the de Bruijn graph.\n";
    kmer_Enumeration_Stats edge_stats = enumerate_edges();

    std::chrono::high_resolution_clock::time_point t_edges = std::chrono::high_resolution_clock::now();
    std::cout << "Enumerated the edge set of the graph. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_edges - t_start).count() << " seconds.\n";


    std::cout << "\nEnumerating the vertices of the de Bruijn graph.\n";
    kmer_Enumeration_Stats vertex_stats = enumerate_vertices(edge_stats.max_memory());

    std::chrono::high_resolution_clock::time_point t_vertices = std::chrono::high_resolution_clock::now();
    std::cout << "Enumerated the vertex set of the graph. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_vertices - t_edges).count() << " seconds.\n";

    std::cout << "Number of edges:    " << edge_stats.kmer_count() << ".\n";
    std::cout << "Number of vertices: " << vertex_stats.kmer_count() << ".\n";


    std::cout << "\nConstructing the minimal perfect hash function (MPHF) over the vertex set.\n";
    construct_hash_table(vertex_stats.kmer_count());

    std::chrono::high_resolution_clock::time_point t_mphf = std::chrono::high_resolution_clock::now();
    std::cout << "Constructed the minimal perfect hash function for the vertices. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_mphf - t_vertices).count() << " seconds.\n";


    std::cout << "\nComputing the DFA states.\n";
    compute_DFA_states();

    Kmer_Container<k + 1>::remove(edge_db_path());
    if(!params.extract_cycles() && !params.dcc_opt())
        hash_table->save(params);
    
    std::chrono::high_resolution_clock::time_point t_dfa = std::chrono::high_resolution_clock::now();
    std::cout << "Computed the states of the automata. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_dfa - t_mphf).count() << " seconds.\n";


    std::cout << "\nExtracting the maximal unitigs.\n";
    extract_maximal_unitigs();

    if(!dbg_info.has_dcc() || dbg_info.dcc_extracted()) // Either there are no DCCs, or the DCCs have already been extracted in this run.
    {
        // Kmer_Container<k>::remove(vertex_db_path());
        hash_table->remove(params);
    }

    std::chrono::high_resolution_clock::time_point t_extract = std::chrono::high_resolution_clock::now();
    std::cout << "Extracted the maximal unitigs. Time taken = " << std::chrono::duration_cast<std::chrono::duration<double>>(t_extract - t_dfa).count() << " seconds.\n";
}


template <uint16_t k>
kmer_Enumeration_Stats Read_CdBG<k>::enumerate_edges() const
{
    return kmer_Enumerator<k + 1>().enumerate(
        KMC::InputFileType::FASTQ, params.sequence_input().seqs(), params.cutoff(),
        params.thread_count(), params.max_memory(), params.strict_memory(), true,
        params.working_dir_path(), edge_db_path());
}


template <uint16_t k>
kmer_Enumeration_Stats Read_CdBG<k>::enumerate_vertices(const std::size_t max_memory) const
{
    return kmer_Enumerator<k>().enumerate(
        KMC::InputFileType::KMC, std::vector<std::string>(1, edge_db_path()), 1,
        params.thread_count(), max_memory, params.strict_memory(), false,
        params.working_dir_path(), vertex_db_path());
}


template <uint16_t k>
void Read_CdBG<k>::construct_hash_table(const uint64_t vertex_count, const bool load)
{
    hash_table = std::make_unique<Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>>(vertex_db_path(), vertex_count);
    load ?  hash_table->load(params) :
            hash_table->construct(params.thread_count(), params.working_dir_path(), params.mph_file_path());
}


template <uint16_t k>
void Read_CdBG<k>::compute_DFA_states()
{
    Read_CdBG_Constructor<k> cdBg_constructor(params, *hash_table);
    cdBg_constructor.compute_DFA_states(edge_db_path());

    dbg_info.add_basic_info(cdBg_constructor);
}


template <uint16_t k>
void Read_CdBG<k>::extract_maximal_unitigs()
{
    Read_CdBG_Extractor<k> cdBg_extractor(params, *hash_table);
    const std::string temp_output_path = params.working_dir_path() + filename(params.output_prefix()) + cuttlefish::file_ext::temp;
    const std::string output_file_path = params.output_file_path();

    cdBg_extractor.extract_maximal_unitigs(vertex_db_path(), temp_output_path);
    dbg_info.add_unipaths_info(cdBg_extractor);

    if(!extract_DCCs(temp_output_path) && params.dcc_opt())
        hash_table->save(params);

    move_file(temp_output_path, output_file_path);
}


template <uint16_t k>
bool Read_CdBG<k>::extract_DCCs(const std::string& output_file_path, const bool rerun)
{
    if(!extract_DCCs_this_run())
        return !dbg_info.has_dcc();

    if(rerun)
    {
        if(!DCC_data_structs_exist())
        {
            std::cout <<    "The data structure(s) required for the cycles extraction have been removed.\n"
                            "Please re-run Cuttlefish with the originial parameters to recover those.\n";
            return false;
        }

        construct_hash_table(Kmer_Container<k>::size(vertex_db_path()), true);
    }

    
    Read_CdBG_Extractor<k> cdBg_extractor(params, *hash_table);
    cdBg_extractor.extract_detached_cycles(vertex_db_path(), output_file_path, dbg_info);

    dbg_info.add_DCC_info(cdBg_extractor);

    return true;
}


template <uint16_t k>
bool Read_CdBG<k>::extract_DCCs_this_run() const
{
    if(!dbg_info.has_dcc())
    {
        std::cout << "The graph does not contain any detached chordless cycles.\n";
        return false;
    }

    if(dbg_info.dcc_extracted())
    {
        std::cout << "The detached chordless cycles have been extracted earlier.\n";
        return false;
    }

    if(!params.extract_cycles())
    {
        std::cout <<    "There are Detached Chordless Cycles (DCC) present in the graph.\n"
                        "Run Cuttlefish with the `cycles` argument to extract those.\n";
        return false;
    }


    return true;
}


template <uint16_t k>
bool Read_CdBG<k>::DCC_data_structs_exist() const
{
    const std::string vertex_db_path = params.output_prefix() + cuttlefish::file_ext::vertices_ext;
    const std::string mph_path = params.mph_file_path();
    const std::string buckets_path = params.buckets_file_path();

    return Kmer_Container<k>::exists(vertex_db_path) && file_exists(mph_path) && file_exists(buckets_path);
}


template <uint16_t k>
bool Read_CdBG<k>::is_constructed() const
{
    return file_exists(params.json_file_path());
}


template <uint16_t k>
const std::string Read_CdBG<k>::edge_db_path() const
{
    return params.output_prefix() + cuttlefish::file_ext::edges_ext;
}


template <uint16_t k>
const std::string Read_CdBG<k>::vertex_db_path() const
{
    return params.output_prefix() + cuttlefish::file_ext::vertices_ext;
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG)


#include "Read_CdBG.hpp"
#include "Read_CdBG_Extractor.hpp"


// The following methods are not used anymore with the current algorithm.
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



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG)

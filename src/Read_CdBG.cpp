
#include "Read_CdBG.hpp"
#include "Read_CdBG_Constructor.hpp"
#include "Read_CdBG_Extractor.hpp"

#include <iomanip>


template <uint16_t k>
Read_CdBG<k>::Read_CdBG(const Build_Params& params):
    params(params),
    dbg_info(params.json_file_path())
{}


template <uint16_t k>
void Read_CdBG<k>::construct()
{
    if(is_constructed(params) && (!dbg_info.has_dcc() || dbg_info.dcc_extracted()))
    {
        std::cout << "\nThe compacted de Bruijn graph has already been completely constructed earlier.\n";
        return;
    }


    dbg_info.add_build_params(params);


    std::cout << "\nConstructing the minimal perfect hash function (MPHF) over the vertex set.\n";
    hash_table = std::make_unique<Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>>(params.vertex_db_path());
    hash_table->construct(params.thread_count(), params.working_dir_path(), params.mph_file_path());

    std::cout << "\nComputing the DFA states.\n";
    compute_DFA_states();

    if(!params.extract_cycles() && !params.dcc_opt())
        hash_table->save(params);

    std::cout << "\nExtracting the maximal unitigs.\n";
    extract_maximal_unitigs();

    if(!dbg_info.has_dcc() || dbg_info.dcc_extracted())
        hash_table->remove(params);


    hash_table->clear();
    dbg_info.dump_info();
}


template <uint16_t k>
void Read_CdBG<k>::compute_DFA_states()
{
    Read_CdBG_Constructor<k> cdBg_constructor(params, *hash_table);
    cdBg_constructor.compute_DFA_states();

    dbg_info.add_basic_info(cdBg_constructor);
}


template <uint16_t k>
void Read_CdBG<k>::extract_maximal_unitigs()
{
    Read_CdBG_Extractor<k> cdBg_extractor(params, *hash_table);
    if(!is_constructed(params))
    {
        cdBg_extractor.extract_maximal_unitigs();
        
        dbg_info.add_unipaths_info(cdBg_extractor);

        if(cdBg_extractor.has_dcc())
        {
            if(params.extract_cycles())
            {
                cdBg_extractor.extract_detached_cycles(dbg_info);

                dbg_info.add_DCC_info(cdBg_extractor);
            }
            else if(params.dcc_opt())
                hash_table->save(params);
        }
    }
    else if(params.extract_cycles())
    {
        if(dbg_info.has_dcc())
        {
            if(!dbg_info.dcc_extracted())
            {
                cdBg_extractor.extract_detached_cycles(dbg_info);

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

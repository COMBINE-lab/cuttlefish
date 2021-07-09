
#include "Read_CdBG.hpp"
#include "Read_CdBG_Constructor.hpp"
#include "Read_CdBG_Extractor.hpp"

#include <iomanip>


template <uint16_t k>
Read_CdBG<k>::Read_CdBG(const Build_Params& params):
    params(params),
    hash_table(params.vertex_db_path()),
    dbg_info(params.json_file_path().empty() ? (params.output_file_path() + ".json") : params.json_file_path())
{}


template <uint16_t k>
void Read_CdBG<k>::construct()
{
    std::cout << "\nConstructing the minimal perfect hash function (MPHF) over the vertex set.\n";
    hash_table.construct(params.thread_count(), params.working_dir_path(), params.mph_file_path());


    std::cout << "\nComputing the DFA states.\n";
    Read_CdBG_Constructor<k> cdBg_constructor(params, hash_table);
    cdBg_constructor.compute_DFA_states();

    dbg_info.add_basic_info(cdBg_constructor);


    std::cout << (!params.extract_cycles() ?
                    "\nExtracting the maximal unitigs.\n": "\nExtracting the detached chordless cycles.\n");
    Read_CdBG_Extractor<k> cdBg_extractor(params, hash_table);
    !params.extract_cycles() ?
        (cdBg_extractor.extract_maximal_unitigs(), dbg_info.add_unipaths_info(cdBg_extractor)):
        (cdBg_extractor.extract_detached_cycles(), dbg_info.add_DCC_info(cdBg_extractor));


    hash_table.clear();

    dbg_info.dump_info();
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG)


#include "dBG_Info.hpp"
#include "Read_CdBG_Constructor.hpp"
#include "Read_CdBG_Extractor.hpp"
#include "Build_Params.hpp"
#include "utility.hpp"

#include <iomanip>
#include <fstream>


template <uint16_t k>
dBG_Info<k>::dBG_Info(const std::string& file_path):
    file_path(file_path)
{
    if(file_exists(file_path))
        load_from_file();
}


template <uint16_t k>
void dBG_Info<k>::load_from_file()
{
    std::ifstream input(file_path.c_str());
        
    input >> dBg_info;

    if(input.fail())
    {
        std::cerr << "Error loading JSON object from file " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    input.close();
}


template <uint16_t k>
void dBG_Info<k>::add_basic_info(const Read_CdBG_Constructor<k>& cdbg_constructor)
{
    dBg_info[basic_field]["vertex count"] = cdbg_constructor.vertex_count();
    dBg_info[basic_field]["edge count"] = cdbg_constructor.edge_count();
}


template <uint16_t k>
void dBG_Info<k>::add_unipaths_info(const Read_CdBG_Extractor<k>& cdbg_extractor)
{
    const Unipaths_Meta_info<k>& unipaths_info = cdbg_extractor.unipaths_meta_info();

    dBg_info[contigs_field]["maximal unitig count"] = unipaths_info.unipath_count();
    dBg_info[contigs_field]["vertex count in the maximal unitigs"] = unipaths_info.kmer_count();
    dBg_info[contigs_field]["shortest maximal unitig length"] = unipaths_info.min_len();
    dBg_info[contigs_field]["longest maximal unitig length"] = unipaths_info.max_len();
    dBg_info[contigs_field]["sum maximal unitig length"] = unipaths_info.sum_len();
    dBg_info[contigs_field]["avg. maximal unitig length"] = unipaths_info.avg_len();
    dBg_info[contigs_field]["_comment"] = "lengths are in bases";
}


template <uint16_t k>
void dBG_Info<k>::add_DCC_info(const Read_CdBG_Extractor<k>& cdbg_extractor)
{
    const Unipaths_Meta_info<k>& unipaths_info = cdbg_extractor.unipaths_meta_info();

    dBg_info[dcc_field]["DCCs extracted?"] = true;
    dBg_info[dcc_field]["DCC count"] = unipaths_info.dcc_count();
    dBg_info[dcc_field]["vertex count in the DCCs"] = unipaths_info.dcc_kmer_count();
    dBg_info[dcc_field]["sum DCC length (in bases)"] = unipaths_info.dcc_sum_len();
}


template <uint16_t k>
void dBG_Info<k>::add_build_params(const Build_Params& params)
{
    dBg_info[params_field]["input"] = concat_strings(params.sequence_input().seqs());
    dBg_info[params_field]["k"] = params.k();
    dBg_info[params_field]["output prefix"] = params.output_prefix();
}


template <uint16_t k>
void dBG_Info<k>::dump_info() const
{
    std::ofstream output(file_path.c_str());
    output << std::setw(4) << dBg_info << "\n"; // Pretty-print the JSON wrapper with overloaded `std::setw`.

    if(output.fail())
    {
        std::cerr << "Error writing to the information file " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    output.close();

    std::cout << "\nStructural information for the de Bruijn graph is written to " << file_path << ".\n";
}


template <uint16_t k>
bool dBG_Info<k>::has_dcc() const
{
    return dBg_info[dcc_field]["DCCs present?"];
}


template <uint16_t k>
bool dBG_Info<k>::dcc_opt_performed() const
{
    return dBg_info[dcc_field]["DCC optimization performed?"];
}


template <uint16_t k>
bool dBG_Info<k>::dcc_extracted() const
{
    return dBg_info[dcc_field]["DCCs extracted?"];
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, dBG_Info)

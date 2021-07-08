
#include "dBG_Info.hpp"
#include "Read_CdBG_Constructor.hpp"
#include "Read_CdBG_Extractor.hpp"

#include <iomanip>


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
void dBG_Info<k>::dump_info(const std::string& file_path) const
{
    std::ofstream output(file_path.c_str());
    output << std::setw(4) << dBg_info << "\n"; // Pretty-print the JSON wrapper with overloaded `std::setw`.

    if(output.fail())
    {
        std::cerr << "Error writing to the information file " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    output.close();
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, dBG_Info)

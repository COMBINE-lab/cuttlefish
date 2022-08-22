
#ifndef DBG_INFO_HPP
#define DBG_INFO_HPP



#include "nlohmann/json.hpp"

#include <cstdint>
#include <string>


// Forward declarations.
template <uint16_t k> class Read_CdBG_Constructor;
template <uint16_t k> class Read_CdBG_Extractor;
template <uint16_t k> class CdBG;
template <uint16_t k> class Unipaths_Meta_info;
class Build_Params;


// A class to wrap the structural information of a de Bruijn graph and some execution
// information of Cuttlefish over it.
template <uint16_t k>
class dBG_Info
{
private:

    nlohmann::ordered_json dBg_info;    // A JSON object wrapping all the information.

    const std::string file_path_;   // Path to the disk-file to store the JSON object.

    static constexpr const char* basic_field = "basic info";    // Category header for basic graph information.
    static constexpr const char* contigs_field = "contigs info";    // Category header for information about the contigs (maximal unitigs).
    static constexpr const char* short_seqs_field = "short seqs";   // Category header for information about sequences shorter than length `k`.
    static constexpr const char* dcc_field = "detached chordless cycles (DCC) info";  // Category header for information about the DCCs.
    static constexpr const char* params_field = "parameters info"; // Category header for the graph build parameters.


    // Loads the JSON file from disk, if the corresponding file exists.
    void load_from_file();

    // Adds information about maximal unitigs tracked in `unipaths_info`.
    void add_unipaths_info(const Unipaths_Meta_info<k>& unipaths_info);


public:

    // Constructs a `dBG_Info` object that would correspond to the file at
    // path `file_path`.
    dBG_Info(const std::string& file_path);

    // Returns the path to the disk-file storing the corresponding JSON object.
    std::string file_path() const;

    // Adds build parameters information of the Cuttlefish algorithm from `params`.
    void add_build_params(const Build_Params& params);

    // Adds basic graph structural information from `cdbg_constructor`.
    void add_basic_info(const Read_CdBG_Constructor<k>& cdbg_constructor);

    // Adds basic graph structural information from `cdbg`.
    void add_basic_info(const CdBG<k>& cdbg);

    // Adds information about the extracted maximal unitigs from `cdbg_extractor`.
    void add_unipaths_info(const Read_CdBG_Extractor<k>& cdbg_extractor);

    // Adds information about the extracted maximal unitigs from `cdbg`.
    void add_unipaths_info(const CdBG<k>& cdbg);

    // Adds information about the references shorter than length k.
    void add_short_seqs_info(const std::vector<std::pair<std::string, std::size_t>>& short_seqs);

    // Writes the JSON object to its corresponding disk-file.
    void dump_info() const;
};



#endif

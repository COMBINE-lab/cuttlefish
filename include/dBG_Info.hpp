
#ifndef DBG_INFO_HPP
#define DBG_INFO_HPP



#include "globals.hpp"

#include "nlohmann/json.hpp"


// Forward declarations.
template <uint16_t k> class Read_CdBG_Constructor;
template <uint16_t k> class Read_CdBG_Extractor;


// A class to wrap the structural information of a de Bruijn graph and some execution
// information of Cuttlefish over it.
template <uint16_t k>
class dBG_Info
{
private:

    nlohmann::ordered_json dBg_info;    // A JSON object wrapping all the information.

    static constexpr const char* basic_field = "basic info";    // Category header for basic graph information.
    static constexpr const char* contigs_field = "contigs info";    // Category header for information about the contigs (maximal unitigs).


public:

    // Adds basic graph structural information from `cdbg_constructor`.
    void add_basic_info(const Read_CdBG_Constructor<k>& cdbg_constructor);

    // Adds information about the extracted maximal unitigs from `cdbg_extractor`.
    void add_unipaths_info(const Read_CdBG_Extractor<k>& cdbg_extractor);

    // Writes the information to a file at path `file_path`.
    void dump_info(const std::string& file_path) const;
};



#endif


#include "Kmer_Index_Utility.hpp"
#include "File_Extensions.hpp"

#include <cassert>


const std::string Kmer_Index_Utility::index_file_path(const std::string& idx_pref)
{
    return idx_pref + cuttlefish::file_ext::idx_file_ext;
}


uint16_t Kmer_Index_Utility::kmer_len(const std::string& idx_path)
{
    std::ifstream idx_file(idx_path, std::ios::in | std::ios::binary);

    uint16_t k;
    idx_file.read(reinterpret_cast<char*>(&k), sizeof(k));
    assert(idx_file.gcount() == sizeof(k));

    if(!idx_file)
    {
        std::cerr << "Error reading from the index file at " << idx_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    idx_file.close();

    return k;
}


uint16_t Kmer_Index_Utility::minimizer_len(const std::string& idx_path)
{
    std::ifstream idx_file(idx_path, std::ios::in | std::ios::binary);

    idx_file.seekg(sizeof(uint16_t), std::ios::beg);    // Skip the `k`-value.

    uint16_t l;
    idx_file.read(reinterpret_cast<char*>(&l), sizeof(l));
    assert(idx_file.gcount() == sizeof(l));

    if(!idx_file)
    {
        std::cerr << "Error reading from the index file at " << idx_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    idx_file.close();

    return l;
}

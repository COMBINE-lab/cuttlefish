
#ifndef VALIDATION_PARAMS_HPP
#define VALIDATION_PARAMS_HPP



#include <string>


class Validation_Params
{
private:

    const std::string ref_file_path_;   // Path to the file containing the reference.
    const uint16_t k_;  // The k-parameter of the compacted edge-centric de Bruijn graph.
    const std::string kmc_db_path_;  // Prefix of the KMC database of the k-mer set of the reference.
    const std::string cdbg_file_path_;   // Path to the file containing the maximal unitigs.
    const uint16_t thread_count_;   // Number of threads to work with.
    const std::string& mph_file_path_;   // Optional path to file storing an MPH over the k-mer set.


public:

    // Constructs a parameters wrapper object with the self-explanatory parameters.
    Validation_Params(  const std::string& ref_file_path,
                        uint16_t k,
                        const std::string& kmc_db_path,
                        const std::string& cdbg_file_path,
                        uint16_t thread_count,
                        const std::string& mph_file_path):
        ref_file_path_(ref_file_path),
        k_(k),
        kmc_db_path_(kmc_db_path),
        cdbg_file_path_(cdbg_file_path),
        thread_count_(thread_count),
        mph_file_path_(mph_file_path)
    {}


    // Returns the path to the reference file.
    const std::string& ref_file_path() const
    {
        return ref_file_path_;
    }


    // Returns the k-parameter.
    uint16_t k() const
    {
        return k_;
    }


    // Returns the path to the KMC database.
    const std::string& kmc_db_path() const
    {
        return kmc_db_path_;
    }


    // Returns the path to the algorithm output file containing the maximal unitigs.
    const std::string& cdbg_file_path() const
    {
        return cdbg_file_path_;
    }


    // Returns the number of threads to use.
    uint16_t thread_count() const
    {
        return thread_count_;
    }


    // Returns the path to the optional MPH file.
    const std::string& mph_file_path() const
    {
        return mph_file_path_;
    }
};



#endif


#ifndef VALIDATION_PARAMS_HPP
#define VALIDATION_PARAMS_HPP



#include "Seq_Input.hpp"

#include <string>
#include <thread>


class Validation_Params
{
private:

    const Seq_Input reference_input_;   // Collection of the input references.
    const uint16_t k_;  // The k-parameter of the compacted edge-centric de Bruijn graph.
    const std::string kmc_db_path_;  // Prefix of the KMC database of the k-mer set of the reference.
    const std::string cdbg_file_path_;   // Path to the file containing the maximal unitigs.
    const uint16_t thread_count_;   // Number of threads to work with.
    const std::string& working_dir_path_;    // Path to the working directory (for temporary files).
    const std::string& mph_file_path_;   // Optional path to file storing an MPH over the k-mer set.


public:

    // Constructs a parameters wrapper object with the self-explanatory parameters.
    Validation_Params(  const std::vector<std::string>& ref_paths,
                        const std::vector<std::string>& list_paths,
                        const std::vector<std::string>& dir_paths,
                        uint16_t k,
                        const std::string& kmc_db_path,
                        const std::string& cdbg_file_path,
                        uint16_t thread_count,
                        const std::string& working_dir_path,
                        const std::string& mph_file_path):
        reference_input_(ref_paths, list_paths, dir_paths),
        k_(k),
        kmc_db_path_(kmc_db_path),
        cdbg_file_path_(cdbg_file_path),
        thread_count_(thread_count),
        working_dir_path_(working_dir_path),
        mph_file_path_(mph_file_path)
    {}


    // Returns the reference input collections.
    const Seq_Input& reference_input() const
    {
        return reference_input_;
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


    // Returns the path to the working directory (for temporary files).
    const std::string& working_dir_path() const
    {
        return working_dir_path_;
    }


    // Returns the path to the optional MPH file.
    const std::string& mph_file_path() const
    {
        return mph_file_path_;
    }


    // Returns `true` iff the parameters selections are valid.
    bool is_valid() const;
};


inline bool Validation_Params::is_valid() const
{
    // Even `k` values are not consistent with the theory.
    // Also, `k` needs to be in the range `[1, MAX_K]`.
    if((k_ & 1) == 0 || (k_ > cuttlefish::MAX_K))
    {
        std::cout << "The k-mer length (k) needs to be odd and within " << cuttlefish::MAX_K << ".\n";
        return false;
    }


    // Discard unsupported thread counts.
    const auto num_threads = std::thread::hardware_concurrency();
    if(num_threads > 0 && thread_count_ > num_threads)
    {
        std::cout << "At most " << num_threads << " concurrent threads are supported at the machine.\n";
        return false;
    }


    return true;
}



#endif

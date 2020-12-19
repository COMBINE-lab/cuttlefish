
#ifndef BUILD_PARAMS_HPP
#define BUILD_PARAMS_HPP



#include "Reference_Input.hpp"
#include "Output_Format.hpp"

#include <string>
#include <vector>
#include <thread>


class Build_Params
{
private:

    const Reference_Input reference_input_; // Collection of the input references.
    const uint16_t k_;   // The k parameter for the edge-centric de Bruijn graph to be compacted.
    const std::string kmc_db_path_; // Path to the KMC database containing the k-mer set.
    const uint16_t thread_count_;    // Number of threads to work with.
    const std::string& output_file_path_;   // Path to the output file.
    const cuttlefish::Output_Format output_format_;   // Output format (0: txt, 1: GFAv1, 2: GFAv2).
    const std::string& working_dir_path_;    // Path to the working directory (for temporary files).
    const bool remove_kmc_db_;  // Option to remove the KMC database, once no longer required.
    const std::string& mph_file_path_;   // Optional path to file storing an MPH over the k-mer set.
    const std::string& buckets_file_path_;  // Optional path to file storing the hash table buckets for the k-mer set.


public:

    // Constructs a parameters wrapper object with the self-explanatory parameters.
    Build_Params(   const std::vector<std::string>& ref_paths,
                    const std::vector<std::string>& list_paths,
                    const std::vector<std::string>& dir_paths,
                    const uint16_t k,
                    const std::string& kmc_db_path,
                    const uint16_t thread_count,
                    const std::string& output_file_path,
                    const uint8_t output_format,
                    const std::string& working_dir_path,
                    const bool remove_kmc_db,
                    const std::string& mph_file_path,
                    const std::string& buckets_file_path):
        reference_input_(ref_paths, list_paths, dir_paths),
        k_(k),
        kmc_db_path_(kmc_db_path),
        thread_count_(thread_count),
        output_file_path_(output_file_path),
        output_format_(cuttlefish::Output_Format(output_format)),
        working_dir_path_(working_dir_path),
        remove_kmc_db_(remove_kmc_db),
        mph_file_path_(mph_file_path),
        buckets_file_path_(buckets_file_path)
    {}


    // Returns the reference input collections.
    const Reference_Input& reference_input() const
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


    // Returns the number of threads to use.
    uint16_t thread_count() const
    {
        return thread_count_;
    }


    // Returns the path to the output file.
    const std::string& output_file_path() const
    {
        return output_file_path_;
    }


    // Returns the output format.
    cuttlefish::Output_Format output_format() const
    {
        return output_format_;
    }


    // Returns the working directory (for temporary files).
    const std::string& working_dir_path() const
    {
        return working_dir_path_;
    }


    // Returns the boolean flag for removing the KMC database.
    bool remove_kmc_db() const
    {
        return remove_kmc_db_;
    }


    // Returns the path to the optional MPH file.
    const std::string& mph_file_path() const
    {
        return mph_file_path_;
    }


    // Returns the path to the optional file storing the hash table buckets.
    const std::string& buckets_file_path() const
    {
        return buckets_file_path_;
    }


    // Returns `true` iff the parameters selections are valid.
    bool is_valid() const;
};


inline bool Build_Params::is_valid() const
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


    // Discard invalid output formats.
    if(output_format_ >= cuttlefish::num_op_formats)
    {
        std::cout << "Invalid output file format.\n";
        return false;
    }


    return true;
}



#endif


#ifndef BUILD_PARAMS_HPP
#define BUILD_PARAMS_HPP



#include <string>


class Build_Params
{
private:

    const std::string input_file_path_;   // Path to the file containing the reference or list of references.
    const bool is_list_;    // Whether the input file corresponds to list of multiple references, or is a reference by itself.
    const uint16_t k_;   // The k parameter for the edge-centric de Bruijn graph to be compacted.
    const std::string kmc_db_path_; // Path to the KMC database containing the k-mer set.
    const uint16_t thread_count_;    // Number of threads to work with.
    const std::string& output_file_path_;   // Path to the output file.
    const uint8_t output_format_;   // Output format (0: txt, 1: GFAv1, 2: GFAv2).
    const std::string& working_dir_path_;    // Path to the working directory (for temporary files).
    const std::string& mph_file_path_;   // Optional path to file storing an MPH over the k-mer set.


public:

    // Constructs a parameters wrapper object with the self-explanatory parameters.
    Build_Params(   const std::string& input_file_path,
                    const bool is_list,
                    uint16_t k,
                    const std::string& kmc_db_path,
                    uint16_t thread_count,
                    const std::string& output_file_path,
                    uint8_t output_format,
                    const std::string& working_dir_path,
                    const std::string& mph_file_path):
        input_file_path_(input_file_path),
        is_list_(is_list),
        k_(k),
        kmc_db_path_(kmc_db_path),
        thread_count_(thread_count),
        output_file_path_(output_file_path),
        output_format_(output_format),
        working_dir_path_(working_dir_path),
        mph_file_path_(mph_file_path)
    {}


    // Returns the path to the reference file.
    const std::string& input_file_path() const
    {
        return input_file_path_;
    }


    // Returns whether the input file corresponds to list of multiple references, or is a reference by itself.
    bool is_list() const
    {
        return is_list_;
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
    uint8_t output_format() const
    {
        return output_format_;
    }


    // Returns the working directory (for temporary files).
    const std::string& working_dir_path() const
    {
        return working_dir_path_;
    }


    // Returns the path to the optional MPH file.
    const std::string& mph_file_path() const
    {
        return mph_file_path_;
    }
};



#endif


#ifndef BUILD_PARAMS_HPP
#define BUILD_PARAMS_HPP



#include "globals.hpp"
#include "Reference_Input.hpp"
#include "Output_Format.hpp"

#include <string>
#include <vector>
#include <thread>
#include <iostream>


class Build_Params
{
private:

    const bool is_read_graph_;  // Whether to build a read- or a reference-compacted de Bruijn graph.
    const Reference_Input reference_input_; // Collection of the input references.
    const uint16_t k_;   // The k parameter for the edge-centric de Bruijn graph to be compacted.
    const std::string vertex_db_path_;  // Path to the KMC database containing the vertices (canonical k-mers).
    const std::string edge_db_path_;    // Path to the KMC database containing the edges (canonical (k + 1)-mers).
    const uint16_t thread_count_;    // Number of threads to work with.
    const std::string& output_file_path_;   // Path to the output file.
    const cuttlefish::Output_Format output_format_;   // Output format (0: txt, 1: GFAv1, 2: GFAv2).
    const std::string& working_dir_path_;    // Path to the working directory (for temporary files).
    const bool remove_kmc_db_;  // Option to remove the KMC database, once no longer required.
    const std::string& mph_file_path_;   // Optional path to file storing an MPH over the k-mer set.
    const std::string& buckets_file_path_;  // Optional path to file storing the hash table buckets for the k-mer set.
    const bool extract_cycles_; // Option to extract detached chordless cycles from the de Bruijn graph after compaction.


public:

    // Constructs a parameters wrapper object with the self-explanatory parameters.
    Build_Params(   const bool is_read_graph,
                    const std::vector<std::string>& ref_paths,
                    const std::vector<std::string>& list_paths,
                    const std::vector<std::string>& dir_paths,
                    const uint16_t k,
                    const std::string& vertex_db_path,
                    const std::string& edge_db_path,
                    const uint16_t thread_count,
                    const std::string& output_file_path,
                    const uint8_t output_format,
                    const std::string& working_dir_path,
                    const bool remove_kmc_db,
                    const std::string& mph_file_path,
                    const std::string& buckets_file_path,
                    const bool extract_cycles):
        is_read_graph_(is_read_graph),
        reference_input_(ref_paths, list_paths, dir_paths),
        k_(k),
        vertex_db_path_(vertex_db_path),
        edge_db_path_(edge_db_path),
        thread_count_(thread_count),
        output_file_path_(output_file_path),
        output_format_(cuttlefish::Output_Format(output_format)),
        working_dir_path_(working_dir_path),
        remove_kmc_db_(remove_kmc_db),
        mph_file_path_(mph_file_path),
        buckets_file_path_(buckets_file_path),
        extract_cycles_(extract_cycles)
    {}


    // Returns the boolean flag to whether to build a read- or a reference-compacted de Bruijn graph.
    bool is_read_graph() const
    {
        return is_read_graph_;
    }


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


    // Returns the path to the vertex database.
    const std::string& vertex_db_path() const
    {
        return vertex_db_path_;
    }


    // Returns the path to the edge database.
    const std::string& edge_db_path() const
    {
        return edge_db_path_;
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


    // Returns whether the option of extracting detached chordless cycles is specified.
    bool extract_cycles() const
    {
        return extract_cycles_;
    }


    // Returns `true` iff the parameters selections are valid.
    bool is_valid() const;
};



#endif


#ifndef BUILD_PARAMS_HPP
#define BUILD_PARAMS_HPP



#include "globals.hpp"
#include "Seq_Input.hpp"
#include "Output_Format.hpp"
#include "File_Extensions.hpp"

#include <string>
#include <vector>
#include <thread>
#include <iostream>


class Build_Params
{
private:

    const bool is_read_graph_;  // Whether to build a compacted read de Bruijn graph or not.
    const bool is_ref_graph_;   // Whether to build a compacted reference de Bruijn graph or not.
    const Seq_Input seq_input_; // Collection of the input sequences.
    const uint16_t k_;   // The k parameter for the edge-centric de Bruijn graph to be compacted.
    const uint32_t cutoff_; // Frequency cutoff for the (k + 1)-mers (for short-read set input).
    const std::string vertex_db_path_;  // Path to the KMC database containing the vertices (canonical k-mers).
    const std::string edge_db_path_;    // Path to the KMC database containing the edges (canonical (k + 1)-mers).
    const uint16_t thread_count_;    // Number of threads to work with.
    const std::size_t max_memory_;  // Soft maximum memory limit (in GB).
    const bool strict_memory_;  // Whether strict memory limit restriction is specifiied.
    const std::string output_file_path_;    // Path to the output file.
    const cuttlefish::Output_Format output_format_;   // Output format (0: txt, 1: GFAv1, 2: GFAv2).
    const std::string working_dir_path_;    // Path to the working directory (for temporary files).
    const bool path_cover_; // Whether to extract a maximal path cover of the de Bruijn graph.
    const std::string mph_file_path_;   // Optional path to file storing an MPH over the k-mer set.
    const std::string buckets_file_path_;   // Optional path to file storing the hash table buckets for the k-mer set.
    const bool save_vertices_;  // Option to save the vertex set of the de Bruijn graph (in KMC database format).
#ifdef CF_DEVELOP_MODE
    const double gamma_;    // The gamma parameter for the BBHash MPHF.
#endif


public:

    // Constructs a parameters wrapper object with the self-explanatory parameters.
    Build_Params(   const bool is_read_graph,
                    const bool is_ref_graph,
                    const std::vector<std::string>& seq_paths,
                    const std::vector<std::string>& list_paths,
                    const std::vector<std::string>& dir_paths,
                    const uint16_t k,
                    const uint32_t cutoff,
                    const std::string& vertex_db_path,
                    const std::string& edge_db_path,
                    const uint16_t thread_count,
                    const std::size_t max_memory,
                    const bool strict_memory,
                    const std::string& output_file_path,
                    const uint8_t output_format,
                    const std::string& working_dir_path,
                    const bool path_cover,
                    const std::string& mph_file_path,
                    const std::string& buckets_file_path,
                    const bool save_vertices
#ifdef CF_DEVELOP_MODE
                    , const double gamma
#endif
                    );


    // Returns the boolean flag to whether to build a compacted read de Bruijn graph or not.
    bool is_read_graph() const
    {
        return is_read_graph_;
    }


    // Returns the boolean flag to whether to build a compacted reference de Bruijn graph or not.
    bool is_ref_graph() const
    {
        return is_ref_graph_;
    }


    // Returns the sequence input collection.
    const Seq_Input& sequence_input() const
    {
        return seq_input_;
    }


    // Returns the k-parameter.
    uint16_t k() const
    {
        return k_;
    }


    // Returns the frequency cutoff for the (k + 1)-mers (for short-reads set input).
    uint32_t cutoff() const
    {
        return cutoff_;
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


    // Returns the soft maximum memory limit (in GB).
    std::size_t max_memory() const
    {
        return max_memory_;
    }


    // Returns whether strict memory limit restriction is specifiied.
    bool strict_memory() const
    {
        return strict_memory_;
    }


    // Returns the path prefix for all outputs of the algorithm.
    const std::string output_prefix() const
    {
        return output_file_path_;
    }


    // Returns the path to the output file.
    const std::string output_file_path() const
    {
        return (is_read_graph() || is_ref_graph()) ? (output_file_path_ + cuttlefish::file_ext::unipaths_ext) : output_file_path_;
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


    // Returns whether to extract a maximal path cover of the de Bruijn graph.
    bool path_cover() const
    {
        return path_cover_;
    }


    // Returns the path to the optional MPH file.
    const std::string mph_file_path() const
    {
        return (is_read_graph() || is_ref_graph()) ? (output_file_path_ + cuttlefish::file_ext::hash_ext) : mph_file_path_;
    }


    // Returns the path to the optional file storing the hash table buckets.
    const std::string buckets_file_path() const
    {
        return (is_read_graph() || is_ref_graph()) ? (output_file_path_ + cuttlefish::file_ext::buckets_ext) : buckets_file_path_;
    }


    // Returns whether the option to save the vertex set of the de Bruijn graph (in KMC database format) is specified or not.
    bool save_vertices() const
    {
        return save_vertices_;
    }


    // Returns the path to the optional file storing meta-information about the graph and cuttlefish executions.
    const std::string json_file_path() const
    {
        return output_file_path_ + cuttlefish::file_ext::json_ext;
    }


#ifdef CF_DEVELOP_MODE
    // Returns the gamma parameter for the BBHash MPHF.
    double gamma() const
    {
        return gamma_;
    }
#endif


    // Returns `true` iff the parameters selections are valid.
    bool is_valid() const;
};



#endif

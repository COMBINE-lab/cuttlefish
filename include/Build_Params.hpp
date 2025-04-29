
#ifndef BUILD_PARAMS_HPP
#define BUILD_PARAMS_HPP



#include "Seq_Input.hpp"
#include "Output_Format.hpp"
#include "File_Extensions.hpp"
#include "Input_Defaults.hpp"

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <optional>


class Build_Params
{
private:

    const bool is_read_graph_;  // Whether to build a compacted read de Bruijn graph or not.
    const bool is_ref_graph_;   // Whether to build a compacted reference de Bruijn graph or not.
    const Seq_Input seq_input_; // Collection of the input sequences.
    const uint16_t k_;   // The k parameter for the edge-centric de Bruijn graph to be compacted.
    const std::optional<uint32_t> cutoff_;  // Frequency cutoff for the (k + 1)-mers.
    const bool color_;  // Whether to color the compacted graph or not.
    const std::size_t vertex_part_count_;   // Number of vertex-partitions in the discontinuity graph; needs to be a power of 2.
    const std::size_t lmtig_bucket_count_;  // Number of buckets storing literal locally-maximal unitigs.
    const std::size_t gmtig_bucket_count_;  // Number of buckets storing literal globally-maximal unitigs.
    const uint16_t thread_count_;    // Number of threads to work with.
    const bool idx_;    // Whether to construct a k-mer index of the de Bruijn graph.
    const uint16_t min_len_;    // Length of the l-minimizers used in the k-mer index.
    const std::string output_file_path_;    // Path to the output file.
    const std::optional<cuttlefish::Output_Format> output_format_;  // Output format (0: FASTA, 1: GFAv1, 2: GFAv2, 3: GFA-reduced).
    const bool track_short_seqs_;   // Whether to track input sequences shorter than `k` bases.
    const bool poly_n_stretch_; // Whether to include tiles in GFA-reduced output that track the polyN stretches in the input.
    const std::string working_dir_path_;    // Path to the working directory (for temporary files).
    const bool path_cover_; // Whether to extract a maximal path cover of the de Bruijn graph.


    // Returns the extension of the output file, depending on the output format requested.
    const std::string output_file_ext() const;

public:

    // Constructs a parameters wrapper object with the self-explanatory parameters.
    Build_Params(   bool is_read_graph,
                    bool is_ref_graph,
                    const std::optional<std::vector<std::string>>& seq_paths,
                    const std::optional<std::vector<std::string>>& list_paths,
                    const std::optional<std::vector<std::string>>& dir_paths,
                    uint16_t k,
                    std::optional<uint32_t> cutoff,
                    bool color,
                    std::size_t vertex_part_count,
                    std::size_t lmtig_bucket_count,
                    std::size_t gmtig_bucket_count,
                    uint16_t thread_count,
                    const bool idx,
                    const uint16_t l,
                    const std::string& output_file_path,
                    std::optional<cuttlefish::Output_Format> output_format,
                    bool track_short_seqs,
                    bool poly_n_stretch,
                    const std::string& working_dir_path,
                    bool path_cover
                    );

    // Returns the boolean flag to whether to build a compacted read de Bruijn graph or not.
    auto is_read_graph() const { return is_read_graph_; }

    // Returns the boolean flag to whether to build a compacted reference de Bruijn graph or not.
    auto is_ref_graph() const { return is_ref_graph_; }

    // Returns the sequence input collection.
    const auto& sequence_input() const { return seq_input_; }

    // Returns the k-parameter.
    auto k() const { return k_; }

    // Returns whether to color the compacted graph or not.
    auto color() const { return color_; }

    // Returns the frequency cutoff for the (k + 1)-mers (for short-reads set input).
    auto cutoff() const { return cutoff_.value_or(is_read_graph() ? cuttlefish::_default::CUTOFF_FREQ_READS : cuttlefish::_default::CUTOFF_FREQ_REFS); }

    // Returns the number of vertex-partitions in the discontinuity graph.
    auto vertex_part_count() const { return vertex_part_count_; }

    // Returns the number of buckets storing literal locally-maximal unitigs.
    auto lmtig_bucket_count() const { return lmtig_bucket_count_; }

    // Returns the number of buckets storing literal globally-maximal unitigs.
    auto gmtig_bucket_count() const { return gmtig_bucket_count_; }

    // Returns the number of threads to use.
    auto thread_count() const { return thread_count_; }

    // Returns whether to construct a k-mer index of the de Bruijn graph.
    auto idx() const { return idx_; }

    // Returns the length of the l-minimizers used in the k-mer index.
    auto min_len() const { return min_len_; }

    // Returns the path prefix for all outputs of the algorithm.
    auto output_prefix() const { return output_file_path_; }

    // Returns the path to the output file.
    auto output_file_path() const { return output_file_path_ + output_file_ext(); }

    // Returns the output format.
    auto output_format() const { return output_format_.value_or(cuttlefish::_default::OP_FORMAT); }

    // Returns whether to track input sequences shorter than `k` bases.
    auto track_short_seqs() const { return track_short_seqs_; }

    // Returns whether to include tiles in GFA-reduced output that track the polyN stretches in the input.
    auto poly_n_stretch() const { return poly_n_stretch_; }

    // Returns the path to the output segment-file for the GFA-reduced format.
    auto segment_file_path() const { return output_file_path_ + cuttlefish::file_ext::seg_ext; }

    // Returns the path to the output sequence-file for the GFA-reduced format.
    auto sequence_file_path() const { return output_file_path_ + cuttlefish::file_ext::seq_ext; }

    // Returns the working directory (for temporary files).
    auto working_dir_path() const { return working_dir_path_; }

    // Returns whether to extract a maximal path cover of the de Bruijn graph.
    auto path_cover() const { return path_cover_; }

    // Returns the path to the optional file storing meta-information about the graph and cuttlefish executions.
    auto json_file_path() const { return output_file_path_ + cuttlefish::file_ext::json_ext; }

    // Returns `true` iff the parameters selections are valid.
    bool is_valid() const;
};



#endif

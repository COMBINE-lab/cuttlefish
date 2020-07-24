
#ifndef CDBG_HPP
#define CDBG_HPP



#include "globals.hpp"
#include "Vertex.hpp"
#include "State.hpp"
#include "Kmer.hpp"
#include "Kmer_Hash_Table.hpp"
#include "Annotated_Kmer.hpp"
#include "Oriented_Unitig.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/async.h"
#include "spdlog/sinks/basic_file_sink.h"

#include <sstream>

class CdBG
{
private:

    const std::string ref_file;   // Name of the file containing the reference.
    const uint16_t k;   // The k parameter for the edge-centric de Bruijn graph to be compacted.
    Kmer_Hash_Table Vertices;   // The hash table for the vertices (canonical k-mers) of the de Bruijn graph.
    
    // A running counter kept to track the number of sequence being processed.
    uint32_t seq_count = 0;

    // `output_buffer[t_id]` holds output lines yet to be written to the disk from thread number `t_id`.
    std::vector<std::stringstream> output_buffer;

    // `buffer_size[t_id]` holds the count of lines currently stored at the buffer of thread number `t_id`.
    std::vector<uint64_t> buffer_size;

    // `path_output[t_id]` and `overlap_output[t_id]` are the output streams for the paths
    // and the overlaps between the links in the paths respectively, produced from the
    // underlying sequence, by the thread number `t_id`.
    std::vector<std::ofstream> path_output, overlap_output;

    // After all the threads finish the parallel GFA outputting, `first_unitig[t_id]`,
    // `second_unitig[t_id]`, and `last_unitig[t_id]` contain the first unitig, the
    // second unitig, and the last unitig respectively, seen in their entirety by the
    // thread number `t_id`.
    std::vector<Oriented_Unitig> first_unitig, second_unitig, last_unitig;

    // The GFA header lines.
    const static std::string GFA1_HEADER;
    const static std::string GFA2_HEADER;
    
    // The prefixes for the names of the temprorary files used to store the thread-specific
    // paths and overlaps.
    // TODO: Apply some one-time randomness into the file names to avoid possible name conflicts in the file system.
    const static std::string PATH_OUTPUT_PREFIX;
    const static std::string OVERLAP_OUTPUT_PREFIX;


    // Classifies the vertices into different types (or, classes), using up-to
    // `thread_count` number of threads.
    void classify_vertices(const uint16_t thread_count);

    // Processes classification of the valid k-mers present at the sequence `seq`
    // (of length `seq_len`) that have their starting indices between (inclusive)
    // `left_end` and `right_end`.
    void process_substring(const char* const seq, const size_t seq_len, const size_t left_end, const size_t right_end);

    // Returns the index of the first valid k-mer, i.e. the first k-mer without
    // the placeholder nucleotide 'N', in the index range `[left_end, right_end]`
    // of the sequence `seq`. If no such k-mer is found, returns the first invalid
    // index after its assigned range, i.e. `right_end + 1`.
    size_t search_valid_kmer(const char* const seq, const size_t left_end, const size_t right_end);

    // Processes classification for the canonical versions of the valid k-mers of
    // the sequence `seq` (of length `seq_len`) that are present at its contiguous
    // subsequence starting from the index `start_idx` and have their starting
    // indices at most up-to the index `right_end`. The processing goes up-to
    // either the ending of the assigned range itself (i.e. `right_end`), or to the
    // last k-mer before the first encountered placeholder nucleotide 'N', whichever
    // comes first. Also, returns the non-inclusive point of termination of the
    // processed subsequence, i.e. the index following the end of it.
    size_t process_contiguous_subseq(const char* const seq, const size_t seq_len, const size_t right_end, const size_t start_idx);

    // Process classification for the canonical version `kmer_hat` of some k-mer
    // in the sequence that is isolated, i.e. does not have any adjacent k-mers.
    // Returns `false` iff an attempted state transition for the k-mer failed.
    bool process_isolated_kmer(const cuttlefish::kmer_t& kmer_hat);

    // Processes classification (partially) for the canonical version `kmer_hat` of
    // the first k-mer of some sequence, where the k-mer is encountered in the
    // direction `dir`, the canonical version of the next k-mer in the sequence is
    // `next_kmer_hat`, and the nucletiode succeeding the first k-mer is `next_nucl`.
    // Returns `false` iff an attempted state transition for the k-mer failed.
    bool process_leftmost_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::dir_t dir, const cuttlefish::kmer_t& next_kmer_hat, const cuttlefish::nucleotide_t next_nucl);

    // Processes classification (partially) for the canonical version `kmer_hat` of
    // the last k-mer of some sequence, where the k-mer is encountered in the
    // direction `dir`, and the nucletiode preceding the last k-mer is `prev_nucl`.
    // Returns `false` iff an attempted state transition for the k-mer failed.
    bool process_rightmost_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::dir_t dir, const cuttlefish::nucleotide_t prev_nucl);

    // Processes classification (partially) for the canonical version `kmer_hat` of
    // some internal k-mer of some sequence, where the k-mer is encountered in the
    // direction `dir`, the canoninal version of the next k-mer in the sequence is
    // `next_kmer_hat`, the nucletiode preceding the k-mer is `prev_nucl`, and the
    // nucletiode succeeding the k-mer is `next_nucl`.
    // Returns `false` iff an attempted state transition for the k-mer failed.
    bool process_internal_kmer(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::dir_t dir, const cuttlefish::kmer_t& next_kmer_hat, const cuttlefish::nucleotide_t prev_nucl, const cuttlefish::nucleotide_t next_nucl);

    // Returns a Boolean denoting whether the canonical k-mer `kmer_hat` forms a
    // self loop with the canonical k-mer `next_kmer_hat` in the sequence. This
    // should only be used with the canonical versions of two adjacent k-mers.
    bool is_self_loop(const cuttlefish::kmer_t& kmer_hat, const cuttlefish::kmer_t& next_kmer_hat) const;

    // Returns the plain DNA-complement character of the provided `nucleotide` character.
    static cuttlefish::nucleotide_t complement(const cuttlefish::nucleotide_t nucleotide);

    // Outputs all the distinct maximal unitigs of the compacted de Bruijn graph
    // (in canonical form) to a file named `output_file`, using up-to `thread_count`
    // number of threads.
    void output_maximal_unitigs(const std::string& output_file, const uint16_t thread_count);

    // Writes the maximal unitigs at the sequence `seq` (of length `seq_len`) that
    // have their starting indices between (inclusive) `left_end` and `right_end`,
    // to the stream `output`.
    void output_off_substring(const uint64_t thread_id, const char* const seq, const size_t seq_len, const size_t left_end, const size_t right_end, cuttlefish::logger_t output);

    // Outputs the distinct maximal unitigs of the sequence `seq` (of length
    // `seq_len`) to the stream `output`, that are present at its contiguous
    // subsequence starting from the index `start_idx`, going up-to either
    // the ending of the maximal unitig containing the index `right_end`, or
    // up-to the first encountered placeholder nucleotide 'N'. Also, returns
    // the non-inclusive point of termination of the processed subsequence,
    // i.e. the index following the end of it.
    size_t output_maximal_unitigs(const uint64_t thread_id, const char* const seq, const size_t seq_len, const size_t right_end, const size_t start_idx, cuttlefish::logger_t output);

    // Returns a Boolean denoting whether a k-mer with state `state` traversed in
    // the direction `dir` starts a maximal unitig, where `prev_kmer_state` and
    // `prev_kmer_dir` are the state and the direction of the previous k-mer in
    // the sequence, respectively.
    bool is_unipath_start(const cuttlefish::Vertex_Class vertex_class, const cuttlefish::dir_t dir, const cuttlefish::Vertex_Class prev_kmer_class, const cuttlefish::dir_t prev_kmer_dir) const;

    // Returns a Boolean denoting whether a k-mer with state `state` traversed in
    // the direction `dir` ends a maximal unitig, where `next_kmer_state` and
    // `next_kmer_dir` are the state and the direction of the next k-mer in the
    // sequence, respectively.
    bool is_unipath_end(const cuttlefish::Vertex_Class vertex_class, const cuttlefish::dir_t dir, const cuttlefish::Vertex_Class next_kmer_class, const cuttlefish::dir_t next_kmer_dir) const;

    // Outputs the unitig at the k-mer range between the annotated k-mers
    // `start_kmer` and `end_kmer` of the sequence `seq` (if the unitig had not
    // been output already), to the stream `output`.
    void output_unitig(const uint64_t thread_id, const char* const ref, const Annotated_Kmer& start_kmer, const Annotated_Kmer& end_kmer, cuttlefish::logger_t output);
    
    // Writes the path in the sequence `seq` with its starting and ending k-mers
    // located at the indices `start_kmer_idx` and `end_kmer_idx` respectively,
    // to the stream `output`. If `in_forward` is true, then the string spelled
    // by the path is written; otherwise its reverse complement is written.
    // Note that, the output operation appends a newline at the end.
    void write_path(const uint64_t thread_id, const char* const seq, const size_t start_kmer_idx, const size_t end_kmer_idx, const bool in_forward);

    // Increases the buffer size for this thread, i.e. `buffer_size[thread_id]`
    // by `fill_amount`. If the resulting buffer size overflows a predefined constant
    // size (TODO: refactor the fixed 128 to a const field), then the buffer content
    // at `output_buffer[thread_id]` are dumped into the stream `output` and the buffer
    // is emptied.
    void fill_buffer(const uint64_t thread_id, const uint64_t fill_amount, cuttlefish::logger_t output);

    // Writes the string `str` to the output object `output`.
    static void write(cuttlefish::logger_t output, const std::string& str);
    
    // Flushes the output buffers (one for each thread) to the stream `output`.
    void flush_buffers(const uint16_t thread_count, cuttlefish::logger_t output);

    // Outputs the distinct maximal unitigs (in canonical form) of the compacted de
    // Bruijn graph in GFA format (version `gfa_v`), to the file named `gfa_file_name`,
    // using up-to `thread_count` number of threads.
    void output_maximal_unitigs_gfa(const std::string& gfa_file_name, const uint8_t gfa_v, const uint16_t thread_count);

    // Resets the path output streams (depending on GFA version `gfa_v`) for each
    // thread. Needs to be invoked before processing each new sequence.
    void reset_path_streams(const uint8_t gfa_v, const uint16_t thread_count);

    // Writes the maximal unitigs (in the GFA version `gfa_v`) from the sequence `seq`
    // (of length `seq_len`) that have their starting indices between (inclusive)
    // `left_end` and `right_end`, to the stream `output`.
    void output_gfa_off_substring(const uint64_t thread_id, const char* const seq, const size_t seq_len, const size_t left_end, const size_t right_end, const uint8_t gfa_v, cuttlefish::logger_t output);

    // Outputs the distinct maximal unitigs (in GFA version `gfa_v`) of the sequence `seq`
    // (of length `seq_len`) to the stream `output`, that are present at its contiguous
    // subsequence starting from the index `start_idx`, going up-to either the ending
    // of the maximal unitig containing the index `right_end`, or up-to the first
    // encountered placeholder nucleotide 'N'. Also, returns the non-inclusive point of
    // termination of the processed subsequence, i.e. the index following the end of it.
    size_t output_maximal_unitigs_gfa(const uint64_t thread_id, const char* const seq, const size_t seq_len, const size_t right_end, const size_t start_idx, const uint8_t gfa_v, cuttlefish::logger_t output);

    // Outputs the unitig (in GFA version `gfa_v`) at the k-mer range between the annotated
    // k-mers `start_kmer` and `end_kmer` of the sequence `seq` (if the unitig had not been
    // output already), to the stream `output`.
    void output_unitig_gfa(const uint64_t thread_id, const char* const ref, const Annotated_Kmer& start_kmer, const Annotated_Kmer& end_kmer, const uint8_t gfa_v, cuttlefish::logger_t output);

    // Writes the GFA header record (for version `gfa_v`) to the stream `output`.
    void write_gfa_header(const uint8_t gfa_v, std::ofstream& output) const;

    // Writes the GFA segment of the sequence `seq` having its starting and ending k-mers
    // located at the indices `start_kmer_idx` and `end_kmer_idx` respectively, to the
    // stream `output`, in the GFA format (of version `gfa_v`). The GFA segment is named
    // as `segment_name`. If `dir` is `cuttlefish::FWD`, then the string spelled by the
    // path is written; otherwise its reverse complement is written.
    // Note that, the output operation appends a newline at the end.
    void write_gfa_segment(const uint64_t thread_id, const char* const seq, const uint64_t segment_name, const size_t start_kmer_idx, const size_t end_kmer_idx, const cuttlefish::dir_t dir, const uint8_t gfa_v, cuttlefish::logger_t output);

    // Writes a GFA connection (link, edge, or gap depending on GFA version `gfa_v`) between
    // the oriented unitigs `left_unitig` and `right_unitig`, to the stream `output`.
    void write_gfa_connection(const uint64_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig, const uint8_t gfa_v, cuttlefish::logger_t output);

    // Writes a GFA1 link between the oriented unitigs `left_unitig` and `right_unitig`,
    // to the stream `output`.
    void write_gfa_link(const uint64_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig, cuttlefish::logger_t output);

    // Writes a GFA2 edge between the oriented unitigs `left_unitig` and `right_unitig`,
    // to the stream `output`.
    void write_gfa_edge(const uint64_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig, cuttlefish::logger_t output);

    // Writes a GFA2 gap between the oriented unitigs `left_unitig` and `right_unitig`,
    // to the stream `output`.
    void write_gfa_gap(const uint64_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig, cuttlefish::logger_t output);

    // Appends a link between the oriented unitigs `left_unitig` and `right_unitig` to
    // the path and the overlap output streams of the thread number `thread_id`.
    void append_link_to_path(const uint64_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig);

    // Appends an edge between the oriented unitigs `left_unitig` and `right_unitig` to
    // the path output stream of the thread number `thread_id`.
    void append_edge_to_path(const uint64_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig);

    // Writes the connections (links, edges, or gaps) present between unitigs processed
    // by different threads, for GFA version `gfa_v`.
    void write_inter_thread_connections(const uint16_t thread_count, const uint8_t gfa_v, cuttlefish::logger_t output);

    // Searches for the very first connection (link, edge, or gap) present at the underlying
    // sequence being processed, by scanning through the `first_unitig` and `second_unitig`
    // entries for the `thread_count` number of threads. Puts the left and the right unitigs
    // producing that connectiom into `left_unitig` and `right_unitig` respectively.
    void search_first_connection(const uint16_t thread_count, Oriented_Unitig& left_unitig, Oriented_Unitig& right_unitig) const;

    // Writes a GFA1 path that completely tiles the underlying sequence being processed, at
    // the end of the output file named `gfa_file_name`. It basicaly stiches together the
    // path and overlap outputs produced by the `thread_count` number of threads.
    void write_gfa_path(const uint16_t thread_count, const std::string& gfa_file_name);

    // Writes a GFA2 path that completely tiles the underlying sequence being processed, at
    // the end of the output file named `gfa_file_name`. It basicaly stiches together the
    // path outputs produced by the `thread_count` number of threads.
    void write_gfa_ordered_group(const uint16_t thread_count, const std::string& gfa_file_name);

    // Removes the temporary files used for the thread-specific path output streams (depending
    // on the GFA version `gfa_v`) from the disk.
    void remove_temp_files(const uint16_t thread_count, const uint8_t gfa_v) const;

    // Prints the distribution of the vertex classes for the canonical k-mers present
    // at the database named `kmc_file_name`.
    // For debugging purposes.
    void print_vertex_class_dist(const std::string& kmc_file_name) const;


public:

    CdBG(const std::string& ref_file, const uint16_t k);

    // Constructs the compacted de Bruijn graph with its distinct k-mers collection
    // being present in the KMC database with prefix `kmc_file_name` using
    // `thread_count` threads, and outputs the maximal unitigs into the file named
    // `output_file_name`.
    // TODO: Move `kmc_file_name` as one of the class fields.
    void construct(const std::string& kmc_file_name, const std::string& bbhash_file_name, const uint16_t thread_count, const std::string& output_file_name);
};



inline cuttlefish::nucleotide_t CdBG::complement(const cuttlefish::nucleotide_t nucleotide)
{
    switch (nucleotide)
    {
    case 'A':
        return 'T';

    case 'C':
        return 'G';

    case 'G':
        return 'C';

    case 'T':
        return 'A';
    
    default:
        // Placeholder rule to handle `N` nucleotides.
        // TODO: Need to make an informed rule for this.
        
        std::cerr << "Invalid nucleotide " << nucleotide << " encountered. Aborting.";
        std::exit(EXIT_FAILURE);
    }
}



#endif

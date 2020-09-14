
#ifndef CDBG_HPP
#define CDBG_HPP



#include "globals.hpp"
#include "Kmer_Hash_Table.hpp"
#include "Annotated_Kmer.hpp"
#include "Oriented_Unitig.hpp"
#include "Build_Params.hpp"
#include "Thread_Pool.hpp"
#include "spdlog/async_logger.h"

#include <string>


// De Bruijn graph class to support the compaction algorithm.
template <uint16_t k>
class CdBG
{
    friend class Thread_Pool<k>;

private:

    const Build_Params params;    // Required parameters wrapped in one object.
    Kmer_Hash_Table<k> Vertices;   // The hash table for the vertices (canonical k-mers) of the de Bruijn graph.

    // Minimum size of a partition to be processed by one thread.
    static constexpr uint16_t PARTITION_SIZE_THRESHOLD = 1;

    // `output_buffer[t_id]` holds output content yet to be written to the disk from thread number `t_id`.
    std::vector<std::string> output_buffer;

    // `path_buffer[t_id]` and `overlap_buffer[t_id]` (applicable for GFA1) holds path and overlap
    // output content yet to be written to the disk from the thread number `t_id`.
    std::vector<std::string> path_buffer, overlap_buffer;

    // Capacities for the memory pre-allocation of each output buffer, and the threshold buffer size
    // that triggers a disk-flush.
    static constexpr size_t BUFFER_THRESHOLD = 100 * 1024;  // 100 KB.
    static constexpr size_t BUFFER_CAPACITY = 1.1 * BUFFER_THRESHOLD;   // 110% of the buffer threshold.

    // `spdlog`'s queue size, i.e. the maximum number of log units that it can contain before a flush.
    static constexpr size_t ASYNC_LOG_QUEUE_SZ = 1024;

    // Number of backing worker threads for `spdlog`, i.e. the threads that actually make the writes to sink.
    static constexpr uint16_t ASYNC_LOG_N_THREADS = 1;

    // `spdlog` thread pool for outputting unitigs (or GFA segments and connections).
    std::shared_ptr<spdlog::details::thread_pool> tp_output;

    // `spdlog` thread pool for outputting GFA paths.
    std::shared_ptr<spdlog::details::thread_pool> tp_path;

    // The asynchronous output logger (for GFA segments and connections).
    cuttlefish::logger_t output;

    // Copies of the asynchronous logger `output` for each thread.
    std::vector<cuttlefish::logger_t> output_;

    // `path_output_[t_id]` and `overlap_output_[t_id]` are the output loggers for the paths
    // and the overlaps between the links in the paths respectively, produced from the
    // underlying sequence, by the thread number `t_id`.
    std::vector<cuttlefish::logger_t> path_output_, overlap_output_;

    // After all the threads finish the parallel GFA outputting, `first_unitig[t_id]`,
    // `second_unitig[t_id]`, and `last_unitig[t_id]` contain the first unitig, the
    // second unitig, and the last unitig respectively, seen in their entirety by the
    // thread number `t_id`.
    std::vector<Oriented_Unitig> first_unitig, second_unitig, last_unitig;

    // The GFA header lines.
    const static std::string GFA1_HEADER, GFA2_HEADER;
    
    // The prefixes for the names of the temprorary files used to store the thread-specific
    // paths and overlaps.
    std::string path_file_prefix = "cuttlefish-path-output-";
    std::string overlap_file_prefix = "cuttlefish-overlap-output-";
    static constexpr size_t TEMP_FILE_PREFIX_LEN = 10;

    // Debug
    std::vector<double> seg_write_time;
    std::vector<double> link_write_time;
    std::vector<double> buff_flush_time;
    std::vector<double> path_write_time;
    std::vector<double> path_flush_time;
    double path_concat_time = 0;
    double logger_flush_time = 0;


    /* Build methods */
    
    // Classifies the vertices into different types (or, classes).
    void classify_vertices();

    // Distributes the classification task for the sequence `seq` of length
    // `seq_len` to the thread pool `thread_pool`.
    void distribute_classification(const char* seq, size_t seq_len, Thread_Pool<k>& thread_pool);

    // Processes classification of the valid k-mers present at the sequence `seq`
    // (of length `seq_len`) that have their starting indices between (inclusive)
    // `left_end` and `right_end`.
    void process_substring(const char* seq, size_t seq_len, size_t left_end, size_t right_end);

    // Returns the index of the first valid k-mer, i.e. the first k-mer without
    // a placeholder base, in the index range `[left_end, right_end]` of the
    // sequence `seq`. If no such k-mer is found, returns the first invalid
    // index after its assigned range, i.e. `right_end + 1`.
    size_t search_valid_kmer(const char* seq, size_t left_end, size_t right_end) const;

    // Processes classification for the canonical versions of the valid k-mers of
    // the sequence `seq` (of length `seq_len`) that are present at its contiguous
    // subsequence starting from the index `start_idx` and have their starting
    // indices at most up-to the index `right_end`. The processing goes up-to
    // either the ending of the assigned range itself (i.e. `right_end`), or to the
    // last k-mer before the first encountered placeholder base, whichever
    // comes first. Also, returns the non-inclusive point of termination of the
    // processed subsequence, i.e. the index following the end of it.
    size_t process_contiguous_subseq(const char* seq, size_t seq_len, size_t right_end, size_t start_idx);

    // Process classification for the canonical version `kmer_hat` of some k-mer
    // in the sequence that is isolated, i.e. does not have any adjacent k-mers.
    // Returns `false` iff an attempted state transition for the k-mer failed.
    bool process_isolated_kmer(const Kmer<k>& kmer_hat);

    // Processes classification (partially) for the canonical version `kmer_hat` of
    // the first k-mer of some sequence, where the k-mer is encountered in the
    // direction `dir`, the canonical version of the next k-mer in the sequence is
    // `next_kmer_hat`, and the base character succeeding the first k-mer is `next_char`.
    // Returns `false` iff an attempted state transition for the k-mer failed.
    bool process_leftmost_kmer(const Kmer<k>& kmer_hat, cuttlefish::dir_t dir, const Kmer<k>& next_kmer_hat, char next_char);

    // Processes classification (partially) for the canonical version `kmer_hat` of
    // the last k-mer of some sequence, where the k-mer is encountered in the
    // direction `dir`, and the base character preceding the last k-mer is `prev_char`.
    // Returns `false` iff an attempted state transition for the k-mer failed.
    bool process_rightmost_kmer(const Kmer<k>& kmer_hat, cuttlefish::dir_t dir, char prev_char);

    // Processes classification (partially) for the canonical version `kmer_hat` of
    // some internal k-mer of some sequence, where the k-mer is encountered in the
    // direction `dir`, the canoninal version of the next k-mer in the sequence is
    // `next_kmer_hat`, the base character preceding the k-mer is `prev_char`, and
    // the base character succeeding the k-mer is `next_char`.
    // Returns `false` iff an attempted state transition for the k-mer failed.
    bool process_internal_kmer(const Kmer<k>& kmer_hat, cuttlefish::dir_t dir, const Kmer<k>& next_kmer_hat, char prev_char, char next_char);

    // Returns a Boolean denoting whether the canonical k-mer `kmer_hat` forms a
    // self loop with the canonical k-mer `next_kmer_hat` in the sequence. This
    // should only be used with the canonical versions of two adjacent k-mers.
    bool is_self_loop(const Kmer<k>& kmer_hat, const Kmer<k>& next_kmer_hat) const;


    /* Writer methods */

    // Outputs the compacted de Bruijn graph.
    void output_maximal_unitigs();

    // Outputs all the distinct maximal unitigs of the compacted de Bruijn graph
    // (in canonical form) in plain text format.
    void output_maximal_unitigs_plain();

    // Distributes the outputting task of the maximal unitigs in plain format for
    // the sequence `seq` of length `seq_len` to the thread pool `thread_pool`.
    void distribute_output_plain(const char* seq, size_t seq_len, Thread_Pool<k>& thread_pool);

    // Outputs the distinct maximal unitigs (in canonical form) of the compacted de
    // Bruijn graph in GFA format.
    void output_maximal_unitigs_gfa();

    // Distributes the outputting task of the maximal unitigs in GFA format for
    // the sequence `seq` of length `seq_len` to the thread pool `thread_pool`.
    void distribute_output_gfa(const char* seq, size_t seq_len, Thread_Pool<k>& thread_pool);

    // Clears the output file content.
    void clear_output_file() const;

    // Sets a unique prefix for the temporary files to be used during GFA output.
    void set_temp_file_prefixes(const std::string& working_dir);

    // Initializes the output logger of each thread.
    void init_output_loggers();

    // Resets the path output streams (depending on the GFA version) for each
    // thread. Needs to be invoked before processing each new sequence.
    void reset_path_loggers();

    // Allocates memory for the output buffer of each thread.
    void allocate_output_buffers();

    // Allocate memory for the path (and overlap for GFA1) buffer of each thread.
    void allocate_path_buffers();

    // Writes the maximal unitigs at the sequence `seq` (of length `seq_len`) that
    // have their starting indices between (inclusive) `left_end` and `right_end`.
    // The process is executed by the thread number `thread_id`.
    void output_plain_off_substring(uint16_t thread_id, const char* seq, size_t seq_len, size_t left_end, size_t right_end);

    // Outputs the distinct maximal unitigs of the sequence `seq` (of length
    // `seq_len`), that are present at its contiguous subsequence starting
    // from the index `start_idx`, going up-to either the ending of the maximal
    // unitig containing the index `right_end`, or up-to the first encountered
    // placeholder base. Also, returns the non-inclusive point of termination of
    // the processed subsequence, i.e. the index following the end of it. The
    // process is executed by the thread number `thread_id`.
    size_t output_maximal_unitigs_plain(uint16_t thread_id, const char* seq, size_t seq_len, size_t right_end, size_t start_idx);

    // Returns a Boolean denoting whether a k-mer with state `state` traversed in
    // the direction `dir` starts a maximal unitig, where `prev_kmer_state` and
    // `prev_kmer_dir` are the state and the direction of the previous k-mer in
    // the sequence, respectively.
    bool is_unipath_start(cuttlefish::Vertex_Class vertex_class, cuttlefish::dir_t dir, cuttlefish::Vertex_Class prev_kmer_class, cuttlefish::dir_t prev_kmer_dir) const;

    // Returns a Boolean denoting whether a k-mer with state `state` traversed in
    // the direction `dir` ends a maximal unitig, where `next_kmer_state` and
    // `next_kmer_dir` are the state and the direction of the next k-mer in the
    // sequence, respectively.
    bool is_unipath_end(cuttlefish::Vertex_Class vertex_class, cuttlefish::dir_t dir, cuttlefish::Vertex_Class next_kmer_class, cuttlefish::dir_t next_kmer_dir) const;

    // Outputs the unitig at the k-mer range between the annotated k-mers
    // `start_kmer` and `end_kmer` of the sequence `seq` (if the unitig had not
    // been output already). The process is executed by the thread number `thread_id`.
    void output_plain_unitig(uint16_t thread_id, const char* seq, const Annotated_Kmer<k>& start_kmer, const Annotated_Kmer<k>& end_kmer);
    
    // Writes the path in the sequence `seq` with its starting and ending k-mers
    // located at the indices `start_kmer_idx` and `end_kmer_idx` respectively to
    // the output buffer of the thread number `thread_id`, putting into the logger
    // of the thread, if necessary. If `dir` is `FWD`, then the string spelled by the
    // path is written; otherwise its reverse complement is written.
    // Note that, the output operation appends a newline at the end.
    void write_path(uint16_t thread_id, const char* seq, size_t start_kmer_idx, size_t end_kmer_idx, cuttlefish::dir_t dir);

    // Writes the maximal unitigs from the sequence `seq` (of length `seq_len`) that
    // have their starting indices between (inclusive) `left_end` and `right_end`.
    // The process is executed by the thread number `thread_id`.
    void output_gfa_off_substring(uint16_t thread_id, const char* seq, size_t seq_len, size_t left_end, size_t right_end);

    // Outputs the distinct maximal unitig of the sequence `seq` (of length `seq_len`)
    // that are present at its contiguous subsequence starting from the index `start_idx`,
    // going up-to either the ending of the maximal unitig containing the index `right_end`,
    // or up-to the first encountered placeholder base. Also, returns the non-inclusive
    // point of termination of the processedsubsequence, i.e. the index following the end of
    // it. The process is executed by the thread number `thread_id`.
    size_t output_maximal_unitigs_gfa(uint16_t thread_id, const char* seq, size_t seq_len, size_t right_end, size_t start_idx);

    // Resets the first, the second, and the last unitigs seen for each thread.
    void reset_extreme_unitigs();

    // Outputs the unitig at the k-mer range between the annotated k-mers `start_kmer` and
    // `end_kmer` of the sequence `seq` (if the unitig had not been output already). The
    // process is executed by the thread number `thread_id`.
    void output_gfa_unitig(uint16_t thread_id, const char* ref, const Annotated_Kmer<k>& start_kmer, const Annotated_Kmer<k>& end_kmer);

    // Writes the GFA header record to output.
    void write_gfa_header() const;

    // Writes the GFA segment of the sequence `seq` having its starting and ending k-mers
    // located at the indices `start_kmer_idx` and `end_kmer_idx` respectively, in the GFA
    // format (of version `gfa_v`). The GFA segment is named as `segment_name`. If `dir` is
    // `cuttlefish::FWD`, then the string spelled by the path is written; otherwise its
    // reverse complement is written. The process is executed by the thread number `thread_id`.
    // Note that, the output operation appends a newline at the end.
    void write_gfa_segment(uint16_t thread_id, const char* seq, uint64_t segment_name, size_t start_kmer_idx, size_t end_kmer_idx, cuttlefish::dir_t dir);

    // Writes a GFA connection (link, edge, or gap depending on GFA version `gfa_v`) between
    // the oriented unitigs `left_unitig` and `right_unitig`. The process is executed by the
    // thread number `thread_id`.
    void write_gfa_connection(uint16_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig);

    // Writes a GFA1 link between the oriented unitigs `left_unitig` and `right_unitig`.
    // The process is executed by the thread number `thread_id`.
    void write_gfa_link(uint16_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig);

    // Writes a GFA2 edge between the oriented unitigs `left_unitig` and `right_unitig`.
    // The process is executed by the thread number `thread_id`.
    void write_gfa_edge(uint16_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig);

    // Writes a GFA2 gap between the oriented unitigs `left_unitig` and `right_unitig`.
    // The process is executed by the thread number `thread_id`.
    void write_gfa_gap(uint16_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig);

    // Appends a link between the oriented unitigs `left_unitig` and `right_unitig` to
    // the path and the overlap output streams of the thread number `thread_id`.
    void append_link_to_path(uint16_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig);

    // Appends an edge between the oriented unitigs `left_unitig` and `right_unitig` to
    // the path output stream of the thread number `thread_id`.
    void append_edge_to_path(uint16_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig);

    // Writes the connections (links, edges, or gaps) present between unitigs processed
    // by different threads.
    void write_inter_thread_connections();

    // Searches for the very first connection (link, edge, or gap) present at the underlying
    // sequence being processed, by scanning through the `first_unitig` and `second_unitig`
    // entries for the threads. Puts the left and the right unitigs producing that connection
    // into `left_unitig` and `right_unitig` respectively.
    void search_first_connection(Oriented_Unitig& left_unitig, Oriented_Unitig& right_unitig) const;

    // Writes a GFA1 path that completely tiles the underlying sequence being processed, at
    // the end of the output file. It basically stiches together the path and overlap outputs
    // produced by the threads. The name of the path is written as provided in `path_name`.
    void write_gfa_path(const std::string& path_name);

    // Writes a GFA2 path (formally referred to as "ordered groups") that completely tiles
    // the underlying sequence being processed, at the end of the output file. It basically
    // stiches together the path produced by the threads. The id of the group is written as
    // provided in `path_id`.
    void write_gfa_ordered_group(const std::string& path_id);

    // Ensures that the string `buf` has enough free space to append a log of length
    // `log_len` at its end without overflowing its capacity by flushing its content
    // to the logger `log` if necessary. The request is non-binding in the sense that
    // if the capacity of the buffer `str` is smaller than `log_len`, then this method
    // does not ensure enough buffer space.
    static void ensure_buffer_space(std::string& buf, size_t log_len, const cuttlefish::logger_t& log);

    // Writes the string `str` to the logger `log`, and empties `str`.
    static void flush_buffer(std::string& str, const cuttlefish::logger_t& log);

    // Puts the content of the string `str` to the logger `log`.
    static void write(const std::string& str, const cuttlefish::logger_t& log);

    // Checks the output buffer for the thread number `thread_id`. If the buffer
    // size exceeds `BUFFER_THRESHOLD`, then the buffer content is put into the
    // corresponding logger of the thread, and the buffer is emptied.
    void check_output_buffer(uint16_t thread_id);

    // Checks the path buffer (and overlap buffer if using GFA1). If the buffer
    // size exceeds `BUFFER_THRESHOLD`, then the buffer content is put into the
    // corresponding logger of the thread, and the buffer is emptied.
    void check_path_buffer(uint16_t thread_id);
    
    // Flushes the output buffers (one for each thread).
    void flush_output_buffers();

    // Flushes the output buffers (one for each thread).
    void flush_path_buffers();

    // Flushes (non-blocking) all the loggers (output, path, and overlap (GFA1-specific)).
    void flush_loggers();   // TODO: Make it const after removing timing profiles.

    // Closes (shuts down) all the loggers, with required flushing as necessary.
    void close_loggers();

    // Removes the temporary files used for the thread-specific path output streams
    // (depending on the GFA version) from the disk.
    void remove_temp_files() const;

    // Prints the distribution of the vertex classes for the canonical k-mers present
    // at the database named `kmc_file_name`.
    // For debugging purposes.
    void print_vertex_class_dist() const;


public:

    // Constructs a `CdBG` object with the parameters wrapped at `params`.
    CdBG(const Build_Params& params);

    // Constructs the compacted de Bruijn graph using up-to `thread_count` threads, and
    // outputs the maximal unitigs into the file named `output_file_name`.
    void construct();
};



#endif

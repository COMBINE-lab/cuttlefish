
#ifndef CDBG_HPP
#define CDBG_HPP



#include "globals.hpp"
#include "Oriented_Unitig.hpp"
#include "Build_Params.hpp"
#include "Data_Logistics.hpp"
#include "Unipaths_Meta_info.hpp"
#include "dBG_Info.hpp"
#include "spdlog/sinks/basic_file_sink.h"

#include <cstdint>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>


template <uint16_t k, uint8_t BITS_PER_KEY> class Kmer_Hash_Table;
template <uint16_t k> class Directed_Kmer;
template <uint16_t k> class Annotated_Kmer;
template <uint16_t k> class kmer_Enumeration_Stats;
template <uint16_t k> class Thread_Pool;
template <typename T_id_, typename T_info_> class Job_Queue;


// De Bruijn graph class to support the compaction algorithm.
template <uint16_t k>
class CdBG
{
    friend class Thread_Pool<k>;

private:

    const Build_Params params;    // Required parameters wrapped in one object.
    const Data_Logistics logistics; // Data logistics manager for the algorithm execution.
    std::unique_ptr<Kmer_Hash_Table<k, cuttlefish::BITS_PER_REF_KMER>> hash_table;  // Hash table for the vertices (canonical k-mers) of the graph.

    Unipaths_Meta_info<k> unipaths_meta_info_;  // Meta-information over the extracted maximal unitigs.
    std::vector<Unipaths_Meta_info<k>> unipaths_info_local; // Meta-information over the extracted maximal unitigs per thread.

    dBG_Info<k> dbg_info;   // Wrapper object for structural information of the graph.

    static constexpr double bits_per_vertex = 8.71; // Expected number of bits required per vertex by Cuttlefish 2.
    static constexpr std::size_t parser_memory = 256 * 1024U * 1024U;   // An empirical estimation of the memory used by the sequence parser. 256 MB.

    // Minimum size of a partition to be processed by one thread.
    static constexpr uint16_t PARTITION_SIZE_THRESHOLD = 1;

    // `output_buffer[t_id]` holds output content yet to be written to the disk from thread number `t_id`.
    std::vector<std::string> output_buffer;

    // `path_buffer[t_id]` and `overlap_buffer[t_id]` (applicable for GFA1) holds path and overlap
    // output content yet to be written to the disk from the thread number `t_id`.
    std::vector<std::string> path_buffer, overlap_buffer;

    // `link_added[t_id]` is `true` iff at least one link has been added to the output for thread id `t_id`.
    std::vector<uint64_t> link_added;

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

    // Hold information about references too short for processing
    std::vector<std::pair<std::string, size_t>> short_refs;
    
    // Debug
    std::vector<double> seg_write_time;
    std::vector<double> link_write_time;
    std::vector<double> buff_flush_time;
    std::vector<double> path_write_time;
    std::vector<double> path_flush_time;
    double path_concat_time = 0;
    double logger_flush_time = 0;


    /* Build methods */

    // Returns `true` iff the compacted de Bruijn graph to be built from the parameters
    // collection `params` had been constructed in an earlier execution.
    // NB: only the existence of the output meta-info file is checked for this purpose.
    bool is_constructed() const;

    // Enumerates the vertices of the de Bruijn graph and returns summary statistics of the
    // enumearation.
    kmer_Enumeration_Stats<k> enumerate_vertices() const;

    // Constructs the Cuttlefish hash table for the `vertex_count` vertices of the graph.
    void construct_hash_table(uint64_t vertex_count);

    // TODO: rename the "classify" methods with appropriate terminology that are consistent with the theory.
    
    // Classifies the vertices into different types (or, classes).
    void classify_vertices(std::vector<std::pair<std::string, size_t>>& short_refs_info);

    // Returns the maximum temporary disk-usage incurred by some execution of the algorithm,
    // that has its vertices-enumeration stats in `vertex_stats`.
    static std::size_t max_disk_usage(const kmer_Enumeration_Stats<k>& vertex_stats);

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

    // Processes classification for the directed version `kmer` of some k-mer
    // in the sequence that is isolated, i.e. does not have any adjacent k-mers.
    // Returns `false` iff an attempted state transition for the k-mer failed.
    bool process_isolated_kmer(const Directed_Kmer<k>& kmer);

    // Processes classification (partially) for the directed version `kmer` of
    // the first k-mer in some sequence, where the directed version of the next
    // k-mer in the sequence is `next_kmer`, and the base character succeeding
    // the first k-mer is `next_char`. Returns `false` iff an attempted state
    // transition for the k-mer failed.
    bool process_leftmost_kmer(const Directed_Kmer<k>& kmer, const Directed_Kmer<k>& next_kmer, char next_char);

    // Processes classification (partially) for the directed version `kmer` of
    // the last k-mer in some sequence, where the base character preceding the
    // last k-mer is `prev_char`. Returns `false` iff an attempted state transition
    // for the k-mer failed.
    bool process_rightmost_kmer(const Directed_Kmer<k>& kmer, char prev_char);

    // Processes classification (partially) for the directed version `kmer` of
    // some internal k-mer in some sequence, where the directed version of the
    // next k-mer in the sequence is `next_kmer`, the base character preceding
    // the k-mer is `prev_char`, and the base character succeeding the k-mer is
    // `next_char`. Returns `false` iff an attempted state transition for the
    // k-mer failed.
    bool process_internal_kmer(const Directed_Kmer<k>& kmer, const Directed_Kmer<k>& next_kmer, char prev_char, char next_char);

    // Returns a Boolean denoting whether the canonical k-mer `kmer_hat` forms a
    // self loop with the canonical k-mer `next_kmer_hat` in the sequence. This
    // should only be used with the canonical versions of two adjacent k-mers.
    bool is_self_loop(const Kmer<k>& kmer_hat, const Kmer<k>& next_kmer_hat) const;

    // Processes (partial) classification for the directed version `kmer` of
    // some k-mer in a sequence that forms a self-loop with its next k-mer in
    // the sequence. The directed version of the next k-mer is `next_kmer`, and
    // the base character preceding the k-mer is `prev_char`. Returns `false`
    // iff an attempted state transition for the k-mer failed.
    bool process_loop(const Directed_Kmer<k>& kmer, const Directed_Kmer<k>& next_kmer, char prev_char = 0);


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

    // Outputs the distinct maximal unitigs (in canonical form) of the compacted de
    // Bruijn graph in a GFA-reduced format.
    void output_maximal_unitigs_gfa_reduced();

    // Distributes the outputting task of the maximal unitigs in a GFA-reduced
    // format for the sequence `seq` of length `seq_len` to the thread pool `thread_pool`.
    void distribute_output_gfa_reduced(const char* seq, size_t seq_len, Thread_Pool<k>& thread_pool);

    // Clears the output file content.
    void clear_output_file() const;

    // Sets a unique prefix for the temporary files to be used during GFA output.
    void set_temp_file_prefixes(const std::string& working_dir);

    // Initializes the output logger of each thread.
    void init_output_loggers();

    // Resets the path output streams (depending on the GFA version) for each
    // thread. Needs to be invoked before processing each new sequence. The
    // thread-specific output file names are appended with a `file_id` if a
    // non-zero value is provided for such.
    void reset_path_loggers(uint64_t file_id = 0);

    // Returns the name of the thread-specific path-output file for the thread
    // number `thread_id`. The file names have a suffix `file_id` if a non-zero
    // value is provided for such.
    std::string path_file_name(const uint16_t thread_id, const uint64_t file_id = 0) const;

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

    // Returns a Boolean denoting whether a k-mer with state-class `state_class`
    // traversed in the direction `dir` initiates a maximal unitig traversal, where
    // `prev_kmer_class` and `prev_kmer_dir` are the state-class and the direction
    // of the previous k-mer in the walk, respectively.
    bool is_unipath_start(cuttlefish::State_Class state_class, cuttlefish::dir_t dir, cuttlefish::State_Class prev_kmer_class, cuttlefish::dir_t prev_kmer_dir) const;

    // Returns a Boolean denoting whether a k-mer with state-class `state_class`
    // traversed in the direction `dir` terminates a maximal unitig traversal,
    // where `next_kmer_class` and `next_kmer_dir` are the state-class and the
    // direction of the next k-mer in the walk, respectively.
    bool is_unipath_end(cuttlefish::State_Class state_class, cuttlefish::dir_t dir, cuttlefish::State_Class next_kmer_class, cuttlefish::dir_t next_kmer_dir) const;

    // Outputs the unitig at the k-mer range between the annotated k-mers
    // `start_kmer` and `end_kmer` of the sequence `seq` (if the unitig had not
    // been output already). The process is executed by the thread number `thread_id`.
    void output_plain_unitig(uint16_t thread_id, const char* seq, const Annotated_Kmer<k>& start_kmer, const Annotated_Kmer<k>& end_kmer);
    
    // Writes the path in the sequence `seq` with its starting and ending k-mers
    // located at the indices `start_kmer_idx` and `end_kmer_idx` respectively to
    // the output buffer of the thread number `thread_id`, putting into the logger
    // of the thread, if necessary. The unitig is named as `unitig_id`. If `dir` is
    // `FWD`, then the string spelled by the path is written; otherwise its reverse
    // complement is written. Note that, the output operation appends a newline at the end.
    void write_path(uint16_t thread_id, const char* seq, const uint64_t unitig_id, size_t start_kmer_idx, size_t end_kmer_idx, cuttlefish::dir_t dir);

    // Writes the maximal unitigs from the sequence `seq` (of length `seq_len`) that
    // have their starting indices between (inclusive) `left_end` and `right_end`.
    // The process is executed by the thread number `thread_id`.
    void output_gfa_off_substring(uint16_t thread_id, const char* seq, size_t seq_len, size_t left_end, size_t right_end);

    // Outputs the distinct maximal unitig of the sequence `seq` (of length `seq_len`)
    // that are present at its contiguous subsequence starting from the index `start_idx`,
    // going up-to either the ending of the maximal unitig containing the index `right_end`,
    // or up-to the first encountered placeholder base. Also, returns the non-inclusive
    // point of termination of the processed subsequence, i.e. the index following the end of
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

    // Writes the maximal unitig of the sequence `seq` having its starting and ending k-mers
    // located at the indices `start_kmer_idx` and `end_kmer_idx` respectively. The unitig is
    // named as `segment_name`. If `dir` is `cuttlefish::FWD`, then the string spelled by the
    // path is written; otherwise its reverse complement is written. The process is executed
    // by the thread number `thread_id`. Note that, the output operation appends a newline at
    // the end.
    void write_segment(uint16_t thread_id, const char* seq, uint64_t segment_name, size_t start_kmer_idx, size_t end_kmer_idx, cuttlefish::dir_t dir);

    // Writes the order of the maximal unitigs that completely tiles the underlying sequence
    // being processed, at the end of the sequence-tiling file. It basically stiches together
    // the disjoint tilings produced by the threads. The id of the group is written as `path_id`.
    void write_sequence_tiling(const std::string& path_id);

    // Writes the orders of the maximal unitigs that completely tile the input sequences, at the
    // sequence-tiling file (in the same order as the input sequences). It basically stiches
    // together the disjoint tilings produced by the threads for each input sequence. A producer
    // thread puts a job into the queue `job_queue` after the completion of writings (of disjoint
    // tilings) into the thread-specific files for an input sequence; a consumer thread (executing
    // this method) fetches the jobs from `job_queue`, concatenates the tilings into the output
    // file, and deletes the tiling files.
    void write_sequence_tiling(Job_Queue<std::string, Oriented_Unitig>& job_queue);

    // Ensures that the string `buf` has enough free space to append a log of length
    // `log_len` at its end without overflowing its capacity by flushing its content
    // to the logger `log` if necessary. The request is non-binding in the sense that
    // if the capacity of the buffer `buf` is smaller than `log_len`, then this method
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

    // Flushes the path buffers (one for each thread).
    void flush_path_buffers();

    // Flushes (non-blocking) all the loggers (output, path, and overlap (GFA1-specific)).
    void flush_loggers();   // TODO: Make it const after removing timing profiles.

    // Flushes (non-blocking) the output (segments and connections) logger.
    void flush_output_logger();

    // Flushes (non-blocking) GFA path-output specific loggers.
    void flush_path_loggers();

    // Closes (shuts down) all the loggers, with required flushing as necessary.
    void close_loggers();

    // Closes the path-output specific loggers, with required flushing as necessary.
    void close_path_loggers();

    // Removes the temporary files used for the thread-specific path output streams
    // (depending on the GFA version) from the disk. The thread-specific output file
    // names have a suffix `file_id` if a non-zero value is provided.
    void remove_temp_files(uint64_t file_id = 0) const;


public:

    // Constructs a `CdBG` object with the parameters required for the construction of the
    // compacted representation of the underlying reference de Bruijn graph wrapped in `params`.
    CdBG(const Build_Params& params);

    // Destructs the compacted graph builder object, freeing its hash table and dumping the
    // graph information to disk.
    ~CdBG();

    // Constructs the compacted reference de Bruijn graph, employing the parameters received
    // with the object-constructor.
    void construct();

    // Returns a wrapper over the meta-information of the extracted unitigs.
    const Unipaths_Meta_info<k>& unipaths_meta_info() const;

    // Returns the number of distinct vertices in the underlying graph.
    uint64_t vertex_count() const;
};



#endif

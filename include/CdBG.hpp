
#ifndef CDBG_HPP
#define CDBG_HPP



#include "globals.hpp"
#include "Kmer_Hash_Table.hpp"
#include "Annotated_Kmer.hpp"
#include "Oriented_Unitig.hpp"
#include "Build_Params.hpp"
#include "Thread_Pool.hpp"

#include <sstream>


// De Bruijn graph class to support the compaction algorithm.
template <uint16_t k>
class CdBG
{
    friend class Thread_Pool<k>;

private:

    const Build_Params params;    // Required parameters wrapped in one object.
    Kmer_Hash_Table<k> Vertices;   // The hash table for the vertices (canonical k-mers) of the de Bruijn graph.
    
    uint32_t seq_count = 0; // Running counter to track the number of sequence being processed.

    // Minimum size of a partition to be processed by one thread.
    static constexpr uint16_t PARTITION_SIZE_THRESHOLD = 1;

    // `output_buffer[t_id]` holds output lines yet to be written to the disk from thread number `t_id`.
    std::vector<std::stringstream> output_buffer;

    // `buffer_size[t_id]` holds the count of lines currently stored at the buffer of thread number `t_id`.
    std::vector<uint64_t> buffer_size;

    // Maximum output buffer size that triggers a flush. Currently, the size is measured with line counts. 
    // TODO: Do better.
    const static uint64_t MAX_BUFF_SIZE = 128;

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
    static std::string PATH_OUTPUT_PREFIX;
    static std::string OVERLAP_OUTPUT_PREFIX;
    constexpr static size_t TEMP_FILE_PREFIX_LEN = 10;


    // Sets a unique prefix for the temporary files to be used during GFA output.
    static void set_temp_file_prefixes(const std::string& working_dir);

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

    // Outputs the compacted de Bruijn graph.
    void output_maximal_unitigs();

    // Outputs all the distinct maximal unitigs of the compacted de Bruijn graph
    // (in canonical form) in a plain text format.
    void output_maximal_unitigs_plain();

    // Distributes the outputting task of the maximal unitigs in plain format for
    // the sequence `seq` of length `seq_len` to the thread pool `thread_pool`.
    void distribute_output_plain(const char* seq, size_t seq_len, cuttlefish::logger_t output, Thread_Pool<k>& thread_pool);

    // Outputs the distinct maximal unitigs (in canonical form) of the compacted de
    // Bruijn graph in GFA format.
    void output_maximal_unitigs_gfa();

    // Distributes the outputting task of the maximal unitigs in GFA format for
    // the sequence `seq` of length `seq_len` to the thread pool `thread_pool`.
    void distribute_output_gfa(const char* seq, size_t seq_len, cuttlefish::logger_t output, Thread_Pool<k>& thread_pool);

    // Writes the maximal unitigs at the sequence `seq` (of length `seq_len`) that
    // have their starting indices between (inclusive) `left_end` and `right_end`,
    // to the stream `output`.
    void output_plain_off_substring(uint16_t thread_id, const char* seq, size_t seq_len, size_t left_end, size_t right_end, cuttlefish::logger_t output);

    // Outputs the distinct maximal unitigs of the sequence `seq` (of length
    // `seq_len`) to the stream `output`, that are present at its contiguous
    // subsequence starting from the index `start_idx`, going up-to either
    // the ending of the maximal unitig containing the index `right_end`, or
    // up-to the first encountered placeholder base. Also, returns the
    // non-inclusive point of termination of the processed subsequence, i.e.
    // the index following the end of it.
    size_t output_maximal_unitigs_plain(uint16_t thread_id, const char* seq, size_t seq_len, size_t right_end, size_t start_idx, cuttlefish::logger_t output);

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
    // been output already), to the stream `output`.
    void output_plain_unitig(uint16_t thread_id, const char* seq, const Annotated_Kmer<k>& start_kmer, const Annotated_Kmer<k>& end_kmer, cuttlefish::logger_t output);
    
    // Writes the path in the sequence `seq` with its starting and ending k-mers
    // located at the indices `start_kmer_idx` and `end_kmer_idx` respectively to
    // the output buffer of ghe thread number `thread_id`, flushing to the stream
    // `output` if necessary. If `dir` is `FWD`, then the string spelled by the
    // path is written; otherwise its reverse complement is written.
    // Note that, the output operation appends a newline at the end.
    void write_path(uint16_t thread_id, const char* seq, size_t start_kmer_idx, size_t end_kmer_idx, cuttlefish::dir_t dir, cuttlefish::logger_t output);

    // Increases the buffer size for this thread, i.e. `buffer_size[thread_id]`
    // by `fill_amount`. If the resulting buffer size overflows `MAX_BUFF_SIZE`,
    // then the buffer content at `output_buffer[thread_id]` are dumped into the
    // stream `output` and the buffer is emptied.
    void fill_buffer(uint16_t thread_id, uint64_t fill_amount, cuttlefish::logger_t output);

    // Writes the string `str` to the output object `output`.
    static void write(cuttlefish::logger_t output, const std::string& str);
    
    // Flushes the output buffers (one for each thread) to the stream `output`.
    void flush_buffers(cuttlefish::logger_t output);

    // Resets the path output streams (depending on the GFA version) for each
    // thread. Needs to be invoked before processing each new sequence.
    void reset_path_streams();

    // Writes the maximal unitigs from the sequence `seq` (of length `seq_len`) that
    // have their starting indices between (inclusive) `left_end` and `right_end`,
    // to the stream `output`.
    void output_gfa_off_substring(uint16_t thread_id, const char* seq, size_t seq_len, size_t left_end, size_t right_end, cuttlefish::logger_t output);

    // Outputs the distinct maximal unitig of the sequence `seq` (of length `seq_len`)
    // to the stream `output`, that are present at its contiguous subsequence starting
    // from the index `start_idx`, going up-to either the ending of the maximal unitig
    // containing the index `right_end`, or up-to the first encountered placeholder
    // base. Also, returns the non-inclusive point of termination of the processed
    // subsequence, i.e. the index following the end of it.
    size_t output_maximal_unitigs_gfa(uint16_t thread_id, const char* seq, size_t seq_len, size_t right_end, size_t start_idx, cuttlefish::logger_t output);

    // Outputs the unitig at the k-mer range between the annotated k-mers `start_kmer` and
    // `end_kmer` of the sequence `seq` (if the unitig had not been output already), to the
    // stream `output`.
    void output_gfa_unitig(uint16_t thread_id, const char* ref, const Annotated_Kmer<k>& start_kmer, const Annotated_Kmer<k>& end_kmer, cuttlefish::logger_t output);

    // Writes the GFA header record to the stream `output`.
    void write_gfa_header(std::ofstream& output) const;

    // Writes the GFA segment of the sequence `seq` having its starting and ending k-mers
    // located at the indices `start_kmer_idx` and `end_kmer_idx` respectively, to the
    // stream `output`, in the GFA format (of version `gfa_v`). The GFA segment is named
    // as `segment_name`. If `dir` is `cuttlefish::FWD`, then the string spelled by the
    // path is written; otherwise its reverse complement is written.
    // Note that, the output operation appends a newline at the end.
    void write_gfa_segment(uint16_t thread_id, const char* seq, uint64_t segment_name, size_t start_kmer_idx, size_t end_kmer_idx, cuttlefish::dir_t dir, cuttlefish::logger_t output);

    // Writes a GFA connection (link, edge, or gap depending on GFA version `gfa_v`) between
    // the oriented unitigs `left_unitig` and `right_unitig`, to the stream `output`.
    void write_gfa_connection(uint16_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig, cuttlefish::logger_t output);

    // Writes a GFA1 link between the oriented unitigs `left_unitig` and `right_unitig`,
    // to the stream `output`.
    void write_gfa_link(uint16_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig, cuttlefish::logger_t output);

    // Writes a GFA2 edge between the oriented unitigs `left_unitig` and `right_unitig`,
    // to the stream `output`.
    void write_gfa_edge(uint16_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig, cuttlefish::logger_t output);

    // Writes a GFA2 gap between the oriented unitigs `left_unitig` and `right_unitig`,
    // to the stream `output`.
    void write_gfa_gap(uint16_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig, cuttlefish::logger_t output);

    // Appends a link between the oriented unitigs `left_unitig` and `right_unitig` to
    // the path and the overlap output streams of the thread number `thread_id`.
    void append_link_to_path(uint16_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig);

    // Appends an edge between the oriented unitigs `left_unitig` and `right_unitig` to
    // the path output stream of the thread number `thread_id`.
    void append_edge_to_path(uint16_t thread_id, const Oriented_Unitig& left_unitig, const Oriented_Unitig& right_unitig);

    // Writes the connections (links, edges, or gaps) present between unitigs processed
    // by different threads.
    void write_inter_thread_connections(cuttlefish::logger_t output);

    // Searches for the very first connection (link, edge, or gap) present at the underlying
    // sequence being processed, by scanning through the `first_unitig` and `second_unitig`
    // entries for the threads. Puts the left and the right unitigs producing that connection
    // into `left_unitig` and `right_unitig` respectively.
    void search_first_connection(Oriented_Unitig& left_unitig, Oriented_Unitig& right_unitig) const;

    // Writes a GFA1 path that completely tiles the underlying sequence being processed, at
    // the end of the output file. It basically stiches together the path and overlap outputs
    // produced by the threads.
    void write_gfa_path();

    // Writes a GFA1 path (formally referred to as "ordered groups") that completely tiles
    // the underlying sequence being processed, at the end of the output file. It basically
    // stiches together the path and overlap outputs produced by the threads.
    void write_gfa_ordered_group();

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

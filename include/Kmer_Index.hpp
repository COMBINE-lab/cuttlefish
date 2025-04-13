
#ifndef KMER_INDEX_HPP
#define KMER_INDEX_HPP



#include "Kmer_Index_Utility.hpp"
#include "Kmer.hpp"
#include "DNA_Utility.hpp"
#include "Kmer_Utility.hpp"
#include "Spin_Lock.hpp"
#include "Minimizer_Iterator.hpp"
#include "Minimizer_Utility.hpp"
#include "Kmer_Hasher.hpp"
#include "File_Extensions.hpp"
#include "Build_Params.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "key-value-collator/Key_Value_Collator.hpp"
#include "compact_vector/compact_vector.hpp"
#include "BBHash/BooPHF.h"
#include "elias_fano/sequence.hpp"

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <utility>
#include <fstream>


// =============================================================================


// TODO: use `assert`s aggressively.


// A class to index the k-mers of provided de Bruijn graph path-sequences,
// potentially from many producer threads, based on the k-mers' minimizers.
template <uint16_t k>
class Kmer_Index : private Kmer_Index_Utility
{
    template <uint16_t, uint16_t> friend class Index_Validator;

public:

    class Kmer_Alignment;

private:

    const uint16_t l_;  // Size of the l-minimizers.

    const std::string output_pref;  // Prefix of the output path for the index.
    const std::string working_dir;  // Path to the working directory for the index construction.

    const uint16_t producer_count;  // Number of producer threads supplying the paths to the indexer.
    const uint16_t worker_count;    // Number of worker threads for various multi-threaded tasks.

    const std::optional<Build_Params> params;   // Build parameters wrapped inside.

    const bool retain;  // Whether to retain the index in memory after construction.

    typedef cuttlefish::minimizer_t minimizer_t;    // Minimizers can be represented using 64-bit integers.
    typedef key_value_collator::Identity_Functor<minimizer_t> minimizer_collate_hasher_t;   // The minimizer hasher class for the minimizer-instance collator.
    typedef key_value_collator::Key_Value_Collator<minimizer_t, std::size_t, minimizer_collate_hasher_t> min_collator_t;    // Type of the minimizer-instance collator.
    typedef compact::vector<uint8_t, 2> path_vector_t;      // Type of the bitvector storing the path sequences.
    typedef compact::vector<std::size_t> path_end_vector_t; // Type of the bitvector storing the path endpoints.
    typedef compact::ts_vector<std::size_t> min_vector_t;   // Type of the bitvectors storing the minimizer instance-counts and -indices.
    // typedef std::vector<char> path_vector_t; // For testing.

    path_vector_t paths;    // The concatenated path-sequences.

    std::vector<std::size_t> path_ends_vec; // Vector for the ending indices, possibly with many unused bits, of the paths in the concatenated sequence.
    elias_fano::sequence<true> path_ends;   // Elias-Fano encoded ending indices of the paths in the concatenated sequence.

    std::size_t path_count_;    // Number of paths deposited to the sequence.
    std::size_t sum_paths_len_; // Sum length of all the path-sequences.
    uint64_t num_instances_;    // Number of minimizer instances in the paths.
    uint64_t min_count_;    // Number of unique minimizers in the paths.
    uint64_t max_inst_count_;   // Maximum count of instances for some minimizer.

    std::vector<path_vector_t> producer_path_buf;   // Separate buffer for each producer, to contain their deposited paths.
    std::vector<std::vector<std::size_t>> producer_path_end_buf;    // Separate buffer for each producer, to contain the (exclusive) indices of the path endpoints in their deposited paths.

    std::vector<min_collator_t::buf_t*> producer_min_inst_buf;  // Separate buffer for each producer, to pass the minimizer-instances of their deposited paths to the collator.
    typedef min_collator_t::key_val_pair_t min_inst_t;  // Type of the minimizer-instances when passed to the collator through buffers.
    min_collator_t min_collator;    // Key-value collator to collate the minimizer-instances produced from different producers.

    typedef boomphf::SingleHashFunctor<minimizer_t> minimizer_hasher_t; // The seeded hasher class for minimizers.
                                                                        // TODO: placeholder for now; think it through.
    typedef boomphf::mphf<minimizer_t, minimizer_hasher_t, false> minimizer_mphf_t; // The minimizer-MPHF type.
    std::optional<const minimizer_mphf_t> min_mphf; // MPHF of the minimizers.

    std::optional<min_vector_t> min_inst_count_bv;  // Count of instances per each unique minimizer. Only used in index construction and not a part of the index.
    elias_fano::sequence<false> min_inst_count; // Elias-Fano encoded counts of instances per each unique minimizer.

    std::optional<min_vector_t> min_offset; // Offsets of the instances for the unique minimizers, laid flat all together.

    uint64_t overflow_min_count_;   // Number of unique minimizers with instance-counts exceeding the preset threshold.
    uint64_t overflow_kmer_count_;  // Number of k-mers associated to the overflowing minimizers.

    typedef boomphf::mphf<Kmer<k>, Kmer_Hasher<k>, false> kmer_mphf_t;  // The MPHF type for the overflown k-mers.
    std::optional<const kmer_mphf_t> kmer_mphf; // MPHF of the overflown k-mers.

    std::optional<min_vector_t> overflow_kmer_map;  // Mapping of each overflown k-mer to its corresponding minimizer-instance's index (into its block).

    mutable std::ofstream serialize_stream; // Serialization stream for the index.

    std::size_t curr_token; // Number of tokens generated for the producers so far.

    mutable Spin_Lock lock; // Mutually-exclusive access lock for different producers.

    static constexpr char OVERFLOW_KMER[] = ".overflow.kmers";  // Extension of the temporary file for the k-mers corresponding to the overflowing minimizers.
    static constexpr char OVERFLOW_MIN_INST_IDX[] = ".overflow.offset"; // Extension of the temporary file of the overflowing minimizer-instances' relative indices.


    // Saves the configuration constants of the index, such as the k-mer and
    // the minimizer lengths.
    void save_config() const;

    // Flushes the buffers of the producer with ID `producer_id`.
    void flush(std::size_t producer_id);

    // Constructs the MPHF `min_mphf` over the unique minimizers found.
    void construct_minimizer_mphf();

    // Counts the number of instances for each distinct minimizer.
    void count_minimizer_instances();

    // Enumerates the offsets of the instances for each distinct minimizer.
    void get_minimizer_offsets();

    // Constructs the minimizer-overflow index.
    void construct_overflow_index();

    // Collects the k-mers (and their minimizer-instances' indices into the
    // instance-blocks) corresponding to the minimizers that each have their
    // instance counts exceeding the preset threshold.
    void collect_overflown_kmers();

    // Gathers the k-mers (and their minimizer-instances' indices into the
    // instance-blocks) corresponding to the minimizers within the ID range
    // `[low, high)` that each have their instance counts exceeding the preset
    // threshold. The k-mers and the instance-indices are output to `kmer_op`
    // and `inst_idx_op`, respectively, in a resource-locked manner. Puts the
    // number of such minimizers and their associated k-mers into `num_min` and
    // `num_kmer`, respectively.
    void collect_overflown_kmers(std::size_t low, std::size_t high, std::ofstream& kmer_op, std::ofstream& inst_idx_op, std::size_t& num_min, std::size_t& num_kmer) const;

    // Constructs the MPHF `kmer_mphf` over the overflown k-mers found.
    void construct_overflow_kmer_mphf();

    // Maps the overflown k-mers to their corresponding minimizer-instances'
    // indices (into the instance-blocks).
    void map_overflown_kmers();

    // Looks up the minimizer `min` in the MPHF and returns its hash value + 1.
    uint64_t hash(minimizer_t min) const;

    // Extracts the k-mer situated at the index `idx` of the concatenated paths
    // sequence into `kmer`.
    void get_kmer(std::size_t idx, Kmer<k>& kmer) const;

    // Closes the deposit stream incoming from the producers and flushes the
    // remaining content to disk.
    void close_deposit_stream();

    // Removes the temporary files used in the index construction and closes the
    // output stream of the index.
    void close_output() const;

    // Tries to align the k-mer `kmer` to the concatenated paths sequence at the
    // index `idx`. Returns `true` iff the alignment succeeds.
    bool align(const Kmer<k>& kmer, std::size_t idx) const;

    // Tries to align the k-mer `kmer` to the concatenated paths sequence such
    // that the `l_`-mer at the index `kmer_min_idx` of `kmer` docks at the
    // index `min_idx` of the concatenated paths sequence. The alignment must be
    // contained in a path—`kmer` is to be present entirely in a single path.
    // Iff the alignment succeeds, the resultant information is stored into
    // `alignment` and `true` is returned. Otherwise, `false` is returned.
    bool align_contained(const Kmer<k>& kmer, std::size_t kmer_min_idx, std::size_t min_idx, Kmer_Alignment& alignment) const;

    uint16_t l() const { return l_; }   // Returns the size of the l-minimizers.
    const std::string min_instance_path_pref() const { return working_dir + filename(output_pref) + cuttlefish::file_ext::min_inst_file_ext; }  // Returns the path-prefix for the working files of the minimizer-instance collator.
    const std::string overflow_kmers_path() const { return working_dir + filename(output_pref) + OVERFLOW_KMER; }   // Returns the file-path for the k-mers corresponding to the overflowing minimizers.
    const std::string overflow_min_insts_path() const { return working_dir + filename(output_pref) + OVERFLOW_MIN_INST_IDX; }   // Returns the file-path for the overflowing minimizer-instances' relative index.

public:

    // Constructs a k-mer indexer that will index sequences produced from at
    // most `producer_count` producers, based on the k-mers' `l`-minimizers.
    // Retains the index into memory after construction if `retain` is `true`.
    // The index is stored at the path prefix `output_pref`, and `working_dir`
    // is used to store temporary files during the construction. An optional
    // build-parameters pack `params` can be provided, if compacted de Bruijn
    // graph construction is to be done from this object.
    Kmer_Index(uint16_t l, uint16_t producer_count, bool retain, const std::string& output_pref, const std::string& working_dir, std::optional<Build_Params> params = std::optional<Build_Params>());

    // Loads the k-mer index stored at path `idx_path`.
    Kmer_Index(const std::string& idx_path);

    // Constructs a k-mer indexer object with the parameters required for the
    // index construction wrapped in `params`.
    Kmer_Index(const Build_Params& params);

    // Constructs an index over the underlying de Bruijn graph's k-mers.
    void construct();

    class Producer_Token;

    // Returns a unique token object.
    const Producer_Token get_token();

    // Deposits the sequence `seq` of length `len` to the index, from a producer
    // having the token `token`.
    void deposit(const Producer_Token &token, const char* seq, std::size_t len);

    // Indexes the deposited path-sequences. Should only be invoked after all the
    // sequences to be indexed have been deposited.
    void index();

    std::size_t path_count() const { return path_count_; }  // Returns the number of paths deposited to the sequence.
    std::size_t sum_paths_len() const { return sum_paths_len_; }    // Returns the sum length of all the path-sequences.
    uint64_t num_instances() const { return num_instances_; }   // Returns the number of minimizer instances in the paths.
    uint64_t min_count() const { return min_count_; }   // Returns the number of unique minimizers in the paths.
    uint64_t max_inst_count() const { return max_inst_count_; } // Returns the maximum count of instances for some minimizer.

    // Returns the number of k-mers in the index.
    std::size_t size() const;

    // Returns the size of the path (in k-mers) having sequence-ID `path_id`.
    std::size_t path_size(std::size_t path_id) const;

    // Returns the prefix-sum of path-sizes upto but not including the path with
    // sequence-ID `path_id`.
    std::size_t prefix_sum_path_size(std::size_t path_id) const;

    // Extracts the k-mer situated at the index `idx` of the path with ID
    // `path_id` into `kmer`.
    void get_kmer(std::size_t path_id, std::size_t idx, Kmer<k>& kmer) const;

    // Queries the k-mer `kmer` (in its literal form) into the index. If found,
    // stores alignment information (into the index) of the k-mer into `result`:
    // the ID of the k-mer's containing path, the ID of the k-mer itself in the
    // k-mer ordering of the index, and the k-mer's ID within its containing
    // path. All IDs are the sequence-orders of the corresponding entities.
    // Returns `true` iff the k-mer is found.
    bool query(const Kmer<k>& kmer, Kmer_Alignment& result) const;
};


// A token class to provide distinct token objects to different producers, to
// distinguish them for the indexer.
template <uint16_t k>
class Kmer_Index<k>::Producer_Token
{
    friend class Kmer_Index;

private:

    std::size_t id;

    Producer_Token(const std::size_t id): id(id) {}

    std::size_t get_id() const { return id; }
};


template <uint16_t k>
inline void Kmer_Index<k>::deposit(const Producer_Token& token, const char* const seq, const std::size_t len)
{
    if(len < k)
        return;

    const std::size_t id = token.get_id();          // The producer's id.
    if(producer_min_inst_buf[id] == nullptr)
        producer_min_inst_buf[id] = &(min_collator.get_buffer());

    auto& path_buf = producer_path_buf[id];         // The producer's concatenated-paths sequence buffer.
    auto& path_end_buf = producer_path_end_buf[id]; // The producer's path-endpoints buffer.
    auto& min_inst_buf = *producer_min_inst_buf[id];    // The producer's minimizer-instance buffer.

    // Get the minimizers and the path sequence.

    Minimizer_Iterator<const char*, k> minimizer_iterator(seq, len, l_);    // To iterate over each minimizer in `seq`.
    minimizer_t minimizer;                                  // The minimizer itself.
    std::size_t min_idx;                                    // Index of the minimizer within `seq`.
    std::size_t last_min_idx = len;                         // To track minimizer shifts.

    const std::size_t rel_offset = path_buf.size(); // Relative offset of each minimizer in this producer's path buffer.
    for(std::size_t i = 0; i + (k - 1) < len; ++i)
    {
        minimizer_iterator.value_at(minimizer, min_idx);
        if(min_idx != last_min_idx)
            min_inst_buf.emplace_back(minimizer, rel_offset + min_idx),
            last_min_idx = min_idx;

        path_buf.push_back(DNA_Utility::map_base(seq[i]));
        // path_buf.push_back(seq[i]); // For testing.

        ++minimizer_iterator;
    }

    // Get the (k - 1)-length tail of the last k-mer.
    for(std::size_t i = len - (k - 1); i < len; ++i)
        path_buf.push_back(DNA_Utility::map_base(seq[i]));
        // path_buf.push_back(seq[i]); // For testing

    // Get the ending index (exclusive) of this path in the concatenated paths.
    path_end_buf.emplace_back(path_buf.size());

    const std::size_t buf_size = ((path_buf.size() * 2) / 8) + (path_end_buf.size() * sizeof(std::size_t)) + (min_inst_buf.size() * sizeof(min_inst_t));    // Total buffer size (in bytes) of the producer.
    if(buf_size >= buf_sz_th)
        flush(id);
}


template <uint16_t k>
inline void Kmer_Index<k>::flush(const std::size_t producer_id)
{
    auto& path_buf = producer_path_buf[producer_id];            // The producer's concatenated paths sequence buffer.
    auto& path_end_buf = producer_path_end_buf[producer_id];    // The producer's path-endpoints buffer.
    auto& min_inst_buf = *producer_min_inst_buf[producer_id];   // The producer's minimizer-instance buffer.

    // Dump the producer-specific path buffer and the endpoints buffer to the global concatenated-paths and -endpoints.
    lock.lock();

    const std::size_t offset_shift = paths.size();

    // TODO: think of faster copying, possibly with chunks.
    for(auto p = path_buf.cbegin(); p != path_buf.cend(); ++p)
        paths.push_back(*p);

    // Shift the producer-specific relative indices of the path-endpoints to their absolute indices into the concatenated paths.
    std::for_each(path_end_buf.cbegin(), path_end_buf.cend(),
                    [offset_shift, &p_ends = path_ends_vec](const auto& p) { p_ends.emplace_back(p + offset_shift); }
                );

    lock.unlock();

    path_buf.clear();
    path_end_buf.clear();


    // Shift the producer-specific relative indices of the minimizers to their absolute indices into the concatenated paths.
    std::for_each(min_inst_buf.begin(), min_inst_buf.end(),
                    [offset_shift](auto& p){ p.second += offset_shift; }
                );

    // Dump the producer-specific minimizer information to disk.
    min_collator.return_buffer(min_inst_buf);
    producer_min_inst_buf[producer_id] = nullptr;
}


template <uint16_t k>
inline uint64_t Kmer_Index<k>::hash(const cuttlefish::minimizer_t min) const
{
    return min_mphf->lookup(min) + 1;
}


template <uint16_t k>
inline std::size_t Kmer_Index<k>::size() const
{
    return sum_paths_len_ - path_count_ * (k - 1);
}


template <uint16_t k>
inline std::size_t Kmer_Index<k>::path_size(const std::size_t path_id) const
{
    return (path_id == 0 ? path_ends[path_id] : (path_ends[path_id] - path_ends[path_id - 1])) - (k - 1);
}


template <uint16_t k>
inline std::size_t Kmer_Index<k>::prefix_sum_path_size(const std::size_t path_id) const
{
    return path_id == 0 ? 0 : path_ends[path_id - 1] - path_id * (k - 1);
}


template <uint16_t k>
inline void Kmer_Index<k>::get_kmer(const std::size_t path_id, const std::size_t idx, Kmer<k>& kmer) const
{
    assert(path_id < path_count_);

    const std::size_t base_idx = (path_id > 0 ? path_ends[path_id - 1] : 0);
    assert(base_idx + k <= path_ends[path_id]);

    get_kmer(base_idx + idx, kmer);
}


template <uint16_t k>
__attribute__((optimize("unroll-loops")))
inline void Kmer_Index<k>::get_kmer(const std::size_t idx, Kmer<k>& kmer) const
{
    constexpr std::size_t packed_word_count = k / 32;
    uint64_t* const kmer_data = kmer.data();

    // Note: the endianness of the DNA-bases in the path vector is in opposite orientation to Cuttlefish k-mers.

    // Get the partially packed word, i.e. the highest indexed one in `kmer`.
    if constexpr(k % 32)
        kmer_data[packed_word_count] = Kmer_Utility::base_reverse<(k % 32)>(paths.get_int<uint64_t, k % 32>(idx));

    // Get the completely packed words.
    for(std::size_t word_idx = 0; word_idx < packed_word_count; ++word_idx)
        kmer_data[packed_word_count - 1 - word_idx] = Kmer_Utility::base_reverse<32>(paths.get_int<uint64_t, 32>(idx + (k % 32) + word_idx * 32));
}


// A class to pack a k-mer alignment result from `Kmer_Index`.
template <uint16_t k>
class Kmer_Index<k>::Kmer_Alignment
{
    friend class Kmer_Index<k>;

private:

    std::size_t path_id_;   // The ID of the path containing the query k-mer.
    std::size_t kmer_id_;   // The ID of the query k-mer itself in the k-mer ordering of the index.
    std::size_t kmer_id_in_path_;   // The ID of the query k-mer within its containing path.

public:

    std::size_t path_id() const { return path_id_; }    // Returns the ID of the path containing the query k-mer.
    std::size_t kmer_id() const { return kmer_id_; }    // Returns the ID of the query k-mer itself in the k-mer ordering of the index.
    std::size_t kmer_id_in_path() const { return kmer_id_in_path_; }    // Returns the ID of the query k-mer within its containing path.
};


template <uint16_t k>
inline bool Kmer_Index<k>::query(const Kmer<k>& kmer, Kmer_Alignment& result) const
{
    minimizer_t kmer_min;   // The minimizer of `kmer`.
    std::size_t kmer_min_idx;   // The index of the minimizer in `kmer`.
    Minimizer_Utility::get_minimizer(kmer, l_, kmer_min, kmer_min_idx);

    const uint64_t h = hash(kmer_min) - 1;
    if(h >= min_count_) // The k-mer's minimizer is absent in the index.
        return false;

    const auto& m_offset = *min_offset;

    const std::size_t idx_begin = (h > 0 ? min_inst_count[h - 1] : 0);  // Offset of the block of indices of the minimizer's instances.
    const std::size_t idx_end = min_inst_count[h];  // End-offset of the instance block (exclusive).

    const std::size_t inst_count = idx_end - idx_begin; // Number of instances of this minimizer.
    if(inst_count >= overflow_threshold)
    {
        const uint64_t kmer_hash = kmer_mphf->lookup(kmer);
        if(kmer_hash >= overflow_kmer_count_)
            return false;

        return align_contained(kmer, kmer_min_idx, m_offset[idx_begin + (*overflow_kmer_map)[kmer_hash]], result);
    }

    // Try to align the k-mer to each instance of this minimizer.
    for(std::size_t i = idx_begin; i < idx_end; ++i)
        if(align_contained(kmer, kmer_min_idx, m_offset[i], result))
            return true;

    return false;
}


template <uint16_t k>
inline bool Kmer_Index<k>::align_contained(const Kmer<k>& kmer, const std::size_t kmer_min_idx, const std::size_t min_idx, Kmer_Alignment& alignment) const
{
    // Docking the k-mer at this instance would make it fall off off either end of the concatenated paths.
    if(min_idx < kmer_min_idx || (min_idx + (k - kmer_min_idx)) > sum_paths_len_)
        return false;

    // The k-mer doesn't align with the concatenated paths sequence when docked at this instance.
    if(!align(kmer, min_idx - kmer_min_idx))
        return false;

    const int64_t l = path_ends.prev_leq(min_idx);  // ID of the path preceding this path.
    const std::size_t l_end = (l < 0 ? 0 : path_ends[l]);   // Index of the left end of the path containing this instance.
    if(min_idx - kmer_min_idx < l_end)  // Alignment starting position of the k-mer exceeds the left end.
        return false;

    const int64_t r = l + 1;    // ID of this path.
    const std::size_t r_end = path_ends[r]; // Index of the right end (exclusive) of the path containing this instance.
    if(min_idx + (k - kmer_min_idx) > r_end)    // Alignment ending position of the k-mer exceeds the right end.
        return false;

    alignment.path_id_ = r;
    alignment.kmer_id_ = min_idx - kmer_min_idx - r * (k - 1);
    alignment.kmer_id_in_path_ = min_idx - kmer_min_idx - l_end;

    return true;
}


template <uint16_t k>
__attribute__((optimize("unroll-loops")))
inline bool Kmer_Index<k>::align(const Kmer<k>& kmer, const std::size_t idx) const
{
    constexpr std::size_t packed_word_count = k / 32;
    const uint64_t* const kmer_data = kmer.data();

    // Note: the endianness of the DNA-bases in the path vector is in opposite orientation to Cuttlefish k-mers.

    // Align the completely-packed words, i.e. except for possibly the highest-indexed word.
    for(std::size_t word_num = 0; word_num < packed_word_count; ++word_num)
        if(paths.get_int<uint64_t, 32>((idx + k) - word_num * 32 - 32) != Kmer_Utility::base_reverse<32>(kmer_data[word_num]))
            return false;

    // Align the (only) partially-packed word.
    if constexpr(k & 31)
        if(paths.get_int<uint64_t, k & 31>(idx) != Kmer_Utility::base_reverse<k & 31>(kmer_data[packed_word_count]))
            return false;

    return true;
}



#endif

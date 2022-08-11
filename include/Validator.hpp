
#ifndef VALIDATOR_HPP
#define VALIDATOR_HPP



#include "globals.hpp"
#include "Kmer_Hasher.hpp"
#include "Validation_Params.hpp"
#include "BBHash/BooPHF.h"
#include "spdlog/sinks/stdout_color_sinks.h"

#include <cstddef>
#include <cstdint>
#include <vector>
#include <string>


template <uint16_t k>
class Validator
{
    typedef boomphf::mphf<Kmer<k>, Kmer_Hasher<k>> mphf_t;    // The MPH function type.

private:

    const Validation_Params params; // Required parameters wrapped in one object.
    mphf_t* mph = NULL; // Minimal perfect hash function over the set of canonical k-mers of the reference.

    constexpr static size_t PROGRESS_GRAIN_SIZE = 1000000;  // 1M
    
    // The gamma factor for the BBHash algorithm. Lowest bits/elem is achieved with gamma = 1,
    // higher values lead to larger mphf but faster construction/query.
    constexpr static double GAMMA_FACTOR = 2.0;
    
    std::vector<std::string> U; // Collection of the maximal unitig strings produced by the compaction algorithm.
    
    // For some k-mer `kmer` with (canonical) hash value `h`, the entry `unitig_id[h]` contains an
    // index `id`, i.e. `id = unitig_id[h]`, such that the unitig present at `U[id]` has `kmer` as
    // one of its flanking k-mers (both if the unitig length is `k`).
    std::vector<size_t> unitig_id;

    // An unitig, when encountered through one of its canonical k-mers, must be traversed in exactly
    // one of these directions.
    enum class Unitig_Dir: uint8_t
    {
        none = 0,   // Direction for "invalid" unitigs, i.e. unitigs absent at the algorithm result;
        fwd = 1,    // Unitig must be traversed in the canonical form;
        bwd = 2,    // Unitig must be traversed in the non-canonical form;
        either = 3  // Unitig can be traversed in either of the forms.
    };

    // For some k-mer `kmer` with (canonical) hash value `h`, the entry `unitig_dir[h]` contains
    // the direction of an unitig to traverse upon when starting the traversal with the k-mer `kmer`.
    // Obviously, `kmer` must be a flanking k-mer of that unitig, i.e. the unitig is at `U [ unitig_id [h] ]`.
    std::vector<Unitig_Dir> unitig_dir;

    typedef std::shared_ptr<spdlog::logger> logger_t;

    // Console logger to display log messages.
    logger_t console;


    // Builds the minimal perfect hash function `mph` or loads it from disk.
    void build_mph_function();

    // Loads the unitigs from the algorithm output file into the collection `U`, and
    // builds the tables `unitig_id` and `unitig_dir`.
    void build_unitig_tables();

    // Traverses the sequence `seq` of length `seq_len` to check that the sequence can be
    // spelled out by the unitigs at `U`. Returns `true` iff the spelling is successful.
    void walk_sequence(const char* seq, size_t seq_len, bool& result) const;

    // Returns the index of the first valid k-mer, i.e. the first k-mer without a placeholder
    // base, of the sequence `seq` (of length `seq_len`), searching onwards from the
    // index `start_idx`. If no such k-mer is found, returns `seq_len`.
    size_t search_valid_kmer(const char* seq, size_t seq_len, size_t start_idx) const;

    // Traverses the sequence `seq` (of length `seq_len`) partially, starting from the index
    // `start_idx` such that, the traversal spells out an unitig from the unitigs collection
    // `U`. If the spelling fails, then returns `std::numeric_limits<size_t>::max()`. Otherwise,
    // returns the index of the immediately following k-mer, which might be invalid; and if
    // valid, that k-mer is the starting point of another unitig traversal.
    size_t walk_first_unitig(const char* seq, size_t seq_len, size_t start_idx) const;

    // Traverses the sequence `seq` (of length `seq_len`) partially, starting from the index
    // `start_idx` such that, the traversal spells out the string `unitig` in the direction `dir`.
    // Returns `true` iff the spelling is successful.
    bool walk_unitig(const char* seq, size_t seq_len, size_t start_idx, const std::string& unitig, Unitig_Dir dir) const;

    // Traverses the sequence `seq` (of length `seq_len`) partially, starting from the index
    // `start_idx` such that, the traversal spells out the string `unitig` in the forward or
    // the backward direction, based on whether `in_forward` is true or false. Returns `true`
    // iff the spelling is successful.
    bool walk_unitig(const char* seq, size_t seq_len, size_t start_idx, const std::string& unitig, bool in_forward) const;

    // Returns validation result of the set of k-mers present at the supposed maximal unitigs
    // produced by the algorithm. Checks if the set of k-mers present at the supposed maximal
    // unitigs at the file `cdbg_file_name` is the same as the set of k-mers at the KMC k-mers
    // database `kmc_db_name`. The boolean result is not guaranteed to hold true always as a
    // minimal perfect hash function is used to hash the k-mers, and it may hash extra k-mers
    // from the results to valid hash values. But the implementation logic ensures that the
    // validation result holds true with very high probability. For this result to hold false,
    // the number of "extra" k-mers encountered in the result must be equal to the number of
    // "missing" k-mers from the results; and also, those "extra" k-mers must hash, without
    // collisions, to exactly the hash values of the "missing" k-mers.
    void validate_kmer_set(bool& result) const;

    // Returns the validation result of the completeness of the coverage of the reference sequence by
    // the resulting unitigs of the compaction algorithm. I.e., walks the reference and checks if it
    // can be spelled out completely with unitigs at collection `U`, using up-to `thread_count` threads.
    void validate_sequence_completion(bool& result);

    // Clears the hash table used.
    void clear();


public:

    // Constructs a `CdBG` object with the parameters wrapped at `params`, and uses `console` for logging.
    Validator(const Validation_Params& params, logger_t console);

    // Constructs a `CdBG` object with the parameters wrapped at `params`.
    Validator(const Validation_Params& params);

    // ~Validator();

    // Performs validation of the uniqueness and completeness of the k-mer set present at the
    // unitigs produced by the compaction algorithm, and also the validation of the complete
    // coverage of the reference by those unitigs. Returns `true` iff the validation succeeds.
    bool validate();
};



#endif

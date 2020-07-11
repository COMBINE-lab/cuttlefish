
#ifndef VALIDATOR_HPP
#define VALIDATOR_HPP



#include "globals.hpp"

#include <string>
#include <vector>
#include <iostream>


class Validator
{
private:

    const std::string ref_file_name; // Name of the file containing the reference.
    const uint16_t k;   // The k-parameter of the compacted edge-centric de Bruijn graph.
    const std::string kmc_db_name;  // Prefix of the KMC database for the k-mer set of the reference.
    const std::string cdbg_file_name;   // File containing the maximal unitigs.
    cuttlefish::mphf_t* mph;    // MPH function over the set of canonical k-mers of the reference.
    
    std::vector<std::string> U; // Collection of the maximal unitig strings produced by the compaction algorithm.
    
    // For some k-mer `kmer` with (canonical) hash value `h`, the entry `unitig_id[h]` contains an
    // index `id`, i.e. `id = unitig_id[h]`, such that the unitig present at `U[id]` has `kmer` as
    // one of its flanking k-mers (both if the unitig length is `k`).
    std::vector<size_t> unitig_id;

    // Possible directions of the produced unitigs.
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


    // static cuttlefish::nucleotide_t complement(const cuttlefish::nucleotide_t nucleotide);

    // Builds the minimal perfect hash function `mph` at the disk-file `bbhash_file_name`
    // (or loads from it), using `thread_count` number of threads.
    void build_mph_function(const std::string& bbhash_file_name, const uint16_t thread_count);

    // Loads the unitigs from the file `cdbg_file_name` into the collection `U`, and
    // builds the tables `unitig_id` and `unitig_dir`.
    void build_unitig_tables();

    // Traverses the sequence `seq` of length `seq_len` to check that the sequence can be
    // spelled out by the unitigs. Returns `true` iff the spelling is successful.
    bool walk_sequence(const char* seq, const size_t seq_len);

    // Returns the index of the first valid k-mer, i.e. the first k-mer without the placeholder
    // nucleotide 'N', of the sequence `seq` (of length `seq_len`), searching onwards from the
    // index `start_idx`. If no such k-mer is found, returns `seq_len`.
    size_t search_valid_kmer(const char* seq, const size_t seq_len, const size_t start_idx) const;

    // Traverses the sequence `seq` (of length `seq_len`) partially, starting from the index
    // `star_idx` such that, the traversal spells out an unitig produced from the compaction
    // algorithm. If the spelling fails, then returns `std::numeric_limits<size_t>::max()`.
    // Otherwise, returns the index of the immediately following k-mer, which might be invalid;
    // and if valid, that k-mer is the starting point of another unitig traversal.
    size_t walk_first_unitig(const char* seq, const size_t seq_len, const size_t start_idx) const;

    // Traverses the sequence `seq` (of length `seq_len`) partially, starting from the index
    // `start_idx` such that, the traversal spells out the unitig `unitig` in the direction `dir`.
    // Returns `true` iff the spelling is successful.
    bool walk_unitig(const char* seq, const size_t seq_len, const size_t start_idx, const std::string& unitig, const Unitig_Dir dir) const;

    // Traverses the sequence `seq` (of length `seq_len`) partially, starting from the index
    // `start_idx` such that, the traversal spells out the unitig `unitig` in the forward or
    // the backward direction, based on whether `in_forward` is true or false. Returns `true`
    // iff the spelling is successful.
    bool walk_unitig(const char* seq, const size_t seq_len, const size_t start_idx, const std::string& unitig, const bool in_forward) const;

    // Returns validation result of the set of k-mers present at the supposed maximal unitigs
    // produced by the algorithm. Checks if the set of k-mers present at the supposed maximal
    // unitigs at the file `cdbg_file_name` is the same as the set of k-mers at the KMC k-mers
    // database `kmc_db_name`. The boolean result is not guaranteed to hold true always as a
    // minimal perfect hash function is used to hash the k-mers, and it may hash extra k-mers
    // from the results to valid hash values. But the implementation logic ensures that the
    // validation result holds true with very high probability. For this result to hold false,
    // the number of "extra" k-mers encountered in the result must be equal to the number of
    // "missing" k-mers from the results, and also, those "extra" k-mers must hash without
    // collision to the hash values of the "missing" k-mers.
    bool validate_kmer_set() const;

    // Returns the validation result of the completeness of the coverage of the reference sequence
    // by the resulting unitigs of the compaction algorithm. I.e., walks the reference and checks
    // if it can be spelled out completely with the result unitigs.
    bool validate_sequence_completion();


public:

    // Constructs a validator object to validate the resultant unitigs produced at the file
    // `cdbg_file_name` by the compaction algorithm on the reference file `ref_file_name`,
    // for k-mers of length `k`. The KMC database of the k-mers are stored at `kmc_db_name`. 
    Validator(const std::string& ref_file_name, const uint16_t k, const std::string& kmc_db_name, const std::string& cdbg_file_name);

    // TODO: Add comment.
    bool validate(const std::string& bbhash_file_name, const uint16_t thread_count);
};



#endif

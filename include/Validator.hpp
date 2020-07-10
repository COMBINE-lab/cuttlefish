
#ifndef VALIDATOR_HPP
#define VALIDATOR_HPP



#include "globals.hpp"

#include <string>


class Validator
{
private:

    const std::string ref_file_name; // Name of the file containing the reference.
    const uint16_t k;   // The k-parameter of the compacted edge-centric de Bruijn graph.
    const std::string kmc_db_name;  // Prefix of the KMC database for the k-mer set of the reference.
    const std::string cdbg_file_name;   // File containing the maximal unitigs.
    cuttlefish::mphf_t* mph;    // MPH function over the set of canonical k-mers of the reference.


    // Builds the minimal perfect hash function `mph` at
    // the disk-file `bbhash_file_name` (or loads from it),
    // using `thread_count` number of threads.
    void build_mph_function(const std::string& bbhash_file_name, const uint16_t thread_count);

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


public:

    // 
    Validator(const std::string& ref_file_name, const uint16_t k, const std::string& kmc_db_name, const std::string& cdbg_file_name);

    // 
    bool validate(const std::string& bbhash_file_name, const uint16_t thread_count);
};



#endif


#ifndef KMER_CONTAINER_HPP
#define KMER_CONTAINER_HPP



#include "kmc_api/kmc_file.h"

#include <cstdint>
#include <cstddef>
#include <string>


template <uint16_t k> class Kmer_Iterator;
template <uint16_t k> class Kmer_Buffered_Iterator;
template <uint16_t k> class Kmer_SPMC_Iterator;

// Wrapper class for KMC databases on disk.
template <uint16_t k>
class Kmer_Container
{
    typedef Kmer_Iterator<k> iterator;
    typedef Kmer_Buffered_Iterator<k> buf_iterator;
    typedef Kmer_SPMC_Iterator<k> spmc_iterator;


private:

    const std::string kmc_file_path;  // Name of the KMC database.
    CKMCFileInfo kmer_database_info; // Meta-information of the database.


public:

    // Constructs a wrapper k-mer container over the KMC database named
    // `kmc_file_path`. Some metadata are loaded into memory.
    Kmer_Container(const std::string& kmc_file_path);

    // Returns the path to the KMC database on disk.
    const std::string& container_location() const;

    // Returns the length of the k-mers present in the underlying k-mer database.
    uint32_t kmer_length() const;

    // Returns the number of k-mers present in the underlying k-mer database.
    uint64_t size() const;

    // Returns the number of k-mers present in the k-mer database with path prefix `kmc_db_path`.
    static uint64_t size(const std::string& kmc_db_path);

    // Returns `true` iff a KMC database is present at the path prefix `kmc_db_prefix`.
    static bool exists(const std::string& kmc_db_prefix);

    // Removes the KMC database at path `kmc_db_prefix` from disk.
    static void remove(const std::string& kmc_db_prefix);

    // Returns the size in bytes of the k-mer database with path prefix `kmc_db_prefix`.
    static std::size_t database_size(const std::string& kmc_db_prefix);    

    // Returns an iterator pointing to the beginning of the underlying k-mer
    // database.
    // buf_iterator buf_begin() const;

    // Returns an iterator pointing to the ending (exclusive) of the underlying
    // k-mer database.
    // buf_iterator buf_end() const;

    // Returns an iterator pointing to the beginning of the underlying k-mer
    // database.
    // iterator begin() const;

    // Returns an iterator pointing to the ending (exclusive) of the underlying
    // k-mer database.
    // iterator end() const;

    // Returns an SPMC iterator pointing to the beginning of the underlying k-mer
    // database, that can support `consumer_count` consumers.
    spmc_iterator spmc_begin(size_t consumer_count) const;

    // Returns an SPMC iterator pointing to the ending (exclusive) of the underlying
    // k-mer database, that can support `consumer_count` consumers.
    spmc_iterator spmc_end(size_t consumer_count) const;
};



#endif

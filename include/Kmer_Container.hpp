
#ifndef KMER_CONTAINER_HPP
#define KMER_CONTAINER_HPP



#include "globals.hpp"
#include "kmc_api/kmc_file.h"


template <uint16_t k> class Kmer_Iterator;


// Wrapper class for KMC databases on disk.
template <uint16_t k>
class Kmer_Container
{
    typedef Kmer_Iterator<k> iterator;


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

    // Returns an iterator pointing to the beginning of the underlying k-mer
    // database.
    iterator begin() const;

    // Returns an iterator pointing to the ending (exclusive) of the underlying
    // k-mer database.
    iterator end() const;
};



#endif

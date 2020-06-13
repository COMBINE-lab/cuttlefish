
#ifndef KMER_CONTAINER_HPP
#define KMER_CONTAINER_HPP



#include "kmc_api/kmc_file.h"


class Kmer_Iterator;


// Wrapper class for KMC databases on disk.
class Kmer_Container
{
    friend class Kmer_Iterator;

    typedef Kmer_Iterator iterator;


private:

    std::string kmc_file_name;  // Name of the KMC database.
    CKMCFileInfo kmer_database_info; // Meta-information of the database.



public:

    // Constructs a wrapper k-mer container over the KMC database named
    // `kmc_file_name`. Some metadata are loaded into memory.
    Kmer_Container(const std::string& kmc_file_name);

    // Returns the length of the k-mers present in the underlying k-mer database.
    uint32_t kmer_length() const;

    // Returns the number of k-mers present in the underlying k-mer database.
    uint64_t size() const;

    // Returns an iterator pointing to the beginning of the underlying k-mer
    // database.
    iterator begin();

    // Returns an iterator pointing to the ending (exclusive) of the underlying
    // k-mer database.
    iterator end();
};



#endif

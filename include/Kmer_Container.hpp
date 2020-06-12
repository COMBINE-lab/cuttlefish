
#ifndef KMER_CONTAINER_HPP
#define KMER_CONTAINER_HPP



#include "kmc_api/kmc_file.h"


class Kmer_Iterator;


// Wrapper class for KMC databases on disk.
class Kmer_Container
{
    friend class Kmer_Iterator;


private:

    CKMCFile kmer_database; // The KMC database.



public:

    // Constructs a k-mer container over the KMC database named `kmc_file_name`.
    // The entire k-mer database is not loaded into memory, some metadata are.
    Kmer_Container(const std::string& kmc_file_name);

    // Returns the length of the k-mers present in the underlying database.
    uint32_t kmer_length() const;

    // Returns the number of k-mers present in the underlying k-mer database.
    uint64_t size() const;

    // Returns an iterator pointing to the beginning of the underlying k-mer
    // database.
    Kmer_Iterator begin();

    // Returns an iterator pointing to the ending (exclusive) of the underlying
    // k-mer database.
    Kmer_Iterator end();

    // Closes the underlying k-mer database.
    void close();
};



#endif

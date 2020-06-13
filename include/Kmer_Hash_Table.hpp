
#ifndef KMER_HASH_TABLE_HPP
#define KMER_HASH_TABLE_HPP


#include "globals.hpp"
#include "Kmer.hpp"
#include "Vertex_Encoding.hpp"
#include "BBHash/BooPHF.h"
#include "Kmer_Hasher.hpp"

#include <vector>


class Kmer_Hash_Table
{
    // Tells BBHash to use the custom hash functor `Kmer_Hasher` for the input keys (`kmer_t`). 
    typedef boomphf::mphf<cuttlefish::kmer_t, Kmer_Hasher> boophf_t;    // The MPH function type.

private:

    // Lowest bits/elem is achieved with gamma = 1, higher values lead to larger mphf but faster construction/query.
    constexpr static double gamma_factor = 1.0;
    
    // The MPH function.
    boophf_t* mph = NULL;

    // The values (`Vertex_Encoding`) collection for the hash table;
    // keys (`kmer_t`) are passed to the MPHF, and the resulting function-value is used as index in the values table.
    std::vector<Vertex_Encoding> hash_table;


public:
    Kmer_Hash_Table()
    {}

    // Constructs a minimal perfect hash function (specifically, the BBHash) for
    // the collection of k-mers present at the KMC database named `kmc_file_name`.
    void construct(const std::string& kmc_file_name, const uint16_t thread_count);

    // Returns a mutable reference to the value (in the hash-table) for the
    // key `kmer`.
    Vertex_Encoding& operator [](const cuttlefish::kmer_t& kmer);

    // Returns the value (in the hash-table) for the key `kmer`.
    const Vertex_Encoding operator [](const cuttlefish::kmer_t& kmer) const;

    // Clears the hash-table. Do not invoke on an unused object.
    void clear();
};



inline Vertex_Encoding& Kmer_Hash_Table::operator [](const cuttlefish::kmer_t& kmer)
{
    return hash_table[mph->lookup(kmer)];
}


inline const Vertex_Encoding Kmer_Hash_Table::operator [](const cuttlefish::kmer_t& kmer) const
{
    return hash_table[mph->lookup(kmer)];
}


#endif

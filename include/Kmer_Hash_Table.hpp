
#ifndef KMER_HASH_TABLE_HPP
#define KMER_HASH_TABLE_HPP


#include "globals.hpp"
#include "Kmer.hpp"
#include "Vertex_Encoding.hpp"
#include "BBHash/BooPHF.h"

#include <map>
#include <vector>


// Tells BBHash to use included hash function working on u_int64_t input keys. 
typedef boomphf::SingleHashFunctor<u_int64_t> hasher_t;

// Bloom-filter based minimal perfect hash function type.
typedef boomphf::mphf<u_int64_t, hasher_t> boophf_t;


class Kmer_Hash_Table
{
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

    // Constructs a minimal perfect hash function for the collection of k-mers
    // (in 64-bits format) present in the file named `kmers_file`, that has
    // `kmer_count` number of k-mers.
    void construct(const std::string& kmers_file, const uint64_t kmer_count);

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
    return hash_table[mph->lookup(kmer.int_label())];
}


inline const Vertex_Encoding Kmer_Hash_Table::operator [](const cuttlefish::kmer_t& kmer) const
{
    return hash_table[mph->lookup(kmer.int_label())];
}


#endif

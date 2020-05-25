
#ifndef KMER_HASHTABLE_HPP
#define KMER_HASHTABLE_HPP


#include "globals.hpp"
#include "Kmer.hpp"
#include "Vertex_Encoding.hpp"

#include <map>


class Kmer_HashTable
{
private:
    std::map<cuttlefish::kmer_t, Vertex_Encoding> hash;


public:
    Kmer_HashTable()
    {}

    // Returns true iff the provided k-mer `kmer` has an entry in the hash-table.
    bool is_present(const cuttlefish::kmer_t& kmer) const;

    // Returns a mutable reference to the value (in the hash-table) for the key
    // `kmer`.
    Vertex_Encoding& operator [](const cuttlefish::kmer_t& kmer);

    const Vertex_Encoding operator [](const cuttlefish::kmer_t& kmer) const;

    // Clears the hash-table.
    void clear();

    // For debugging purposes.
    void print_hash_table() const;

    // For debugging purposes.
    void print_stats() const;
};


inline bool Kmer_HashTable::is_present(const cuttlefish::kmer_t& kmer) const
{
    return hash.find(kmer) != hash.end();
}


inline Vertex_Encoding& Kmer_HashTable::operator [](const cuttlefish::kmer_t& kmer)
{
    return hash[kmer];
}


inline const Vertex_Encoding Kmer_HashTable::operator [](const cuttlefish::kmer_t& kmer) const
{
    return hash.find(kmer) -> second;
}


#endif

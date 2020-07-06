
#ifndef KMER_HASH_TABLE_HPP
#define KMER_HASH_TABLE_HPP



#include "globals.hpp"
#include "Kmer.hpp"
#include "Vertex_Encoding.hpp"
#include "BBHash/BooPHF.h"
#include "Kmer_Hasher.hpp"
#include "compact_vector/compact_vector.hpp"
#include "Kmer_Hash_Entry_API.hpp"
#include "SpinLock/SpinLock.hpp"


class Kmer_Hash_Table
{
    // Tells BBHash to use the custom hash functor `Kmer_Hasher` for the input keys (`kmer_t`). 
    typedef boomphf::mphf<cuttlefish::kmer_t, Kmer_Hasher> boophf_t;    // The MPH function type.

private:

    // Lowest bits/elem is achieved with gamma = 1, higher values lead to larger mphf but faster construction/query.
    constexpr static double gamma_factor = 2.0;
    constexpr static const uint64_t num_chunks{65536};

    // The MPH function.
    boophf_t* mph = NULL;

    // The values (`Vertex_Encoding`) collection for the hash table;
    // keys (`kmer_t`) are passed to the MPHF, and the resulting function-value is used as index in the values table.
    cuttlefish::bitvector_t hash_table;

    std::array<SpinLock, num_chunks> locks_;

public:

    Kmer_Hash_Table()
    {}

    // Constructs a minimal perfect hash function (specifically, the BBHash) for
    // the collection of k-mers present at the KMC database named `kmc_file_name`.
    void construct(const std::string& kmc_file_name, const uint16_t thread_count);

    // Returns an API to the entry (in the hash table) for the key `kmer`. The API
    // wraps the hash table position and the vertex encoding value at that position.
    Kmer_Hash_Entry_API operator[](const cuttlefish::kmer_t& kmer);

    // Returns the value (in the hash-table) for the key `kmer`.
    const Vertex_Encoding operator[](const cuttlefish::kmer_t& kmer) const;

    // Attempts to update the entry (in the hash-table) for the API object according
    // to its wrapped vertex encoding values, and returns true or false as per success
    // status. If the corresponding hash table position now contains a different
    // vertex encoding than the one that had been read earlier, then the update fails.
    bool update(Kmer_Hash_Entry_API& api);

    // Clears the hash-table. Do not invoke on an unused object.
    void clear();
};



inline Kmer_Hash_Entry_API Kmer_Hash_Table::operator[](const cuttlefish::kmer_t& kmer)
{
    auto v = mph->lookup(kmer);
    uint64_t lidx = v / num_chunks; 
    locks_[lidx].lock();
    auto r = Kmer_Hash_Entry_API(hash_table[v]);
    locks_[lidx].unlock();
    return r;
}


inline const Vertex_Encoding Kmer_Hash_Table::operator[](const cuttlefish::kmer_t& kmer) const
{
    // NOTE: this makes the `const` a lie.  Should be a better solution here.
    auto v = mph->lookup(kmer);
    uint64_t lidx = v / num_chunks; 
    auto* tp = const_cast<Kmer_Hash_Table*>(this);
    const_cast<decltype(tp->locks_[lidx])>(tp->locks_[lidx]).lock();
    auto ve = Vertex_Encoding(hash_table[v]);
    const_cast<decltype(tp->locks_[lidx])>(tp->locks_[lidx]).unlock();
    return ve;
}


inline bool Kmer_Hash_Table::update(Kmer_Hash_Entry_API& api)
{
    auto it = &(api.bv_entry);
    uint64_t lidx = (std::distance(hash_table.begin(), it)) / num_chunks;
    locks_[lidx].lock();
    bool success = (api.bv_entry == api.get_read_encoding());
    if (success) {
        api.bv_entry = api.get_current_encoding();
    }
    locks_[lidx].unlock();
    return success;
}



#endif

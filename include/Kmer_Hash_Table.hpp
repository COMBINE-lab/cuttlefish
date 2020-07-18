
#ifndef KMER_HASH_TABLE_HPP
#define KMER_HASH_TABLE_HPP



#include "globals.hpp"
#include "Kmer.hpp"
#include "State.hpp"
#include "Kmer_Container.hpp"
#include "BBHash/BooPHF.h"
#include "Kmer_Hasher.hpp"
#include "compact_vector/compact_vector.hpp"
#include "Kmer_Hash_Entry_API.hpp"
#include "SpinLock/SpinLock.hpp"


class Kmer_Hash_Table
{
private:

    // Lowest bits/elem is achieved with gamma = 1, higher values lead to larger mphf but faster construction/query.
    constexpr static double gamma_factor = 2.0;
    constexpr static const uint64_t num_chunks{65536};  // TODO: Comment.
    uint64_t chunk_size;    // TODO: Comment.

    // The MPH function.
    cuttlefish::mphf_t* mph = NULL;

    // The values (`State`) collection for the hash table;
    // keys (`kmer_t`) are passed to the MPHF, and the resulting function-value is used as index in the values table.
    cuttlefish::bitvector_t hash_table;

    std::array<SpinLock, num_chunks> locks_;    // TODO: Comment.


    // Builds the minimal perfect hash function `mph` over the set of
    // k-mers present at the KMC database container `kmer_container`,
    // with `bbhash_file_name` being the file to use for BBHash build
    // using `thread_count` number of threads.
    void build_mph_function(const Kmer_Container& kmer_container, const std::string& bbhash_file_name, const uint16_t thread_count);

public:

    // TODO: Make everything private and add `CdBG` as friend.

    Kmer_Hash_Table()
    {}

    // Constructs a minimal perfect hash function (specifically, the BBHash) for
    // the collection of k-mers present at the KMC database named `kmc_file_name`.
    void construct(const std::string& kmc_file_name, const std::string& bbhash_file_name, const uint16_t thread_count);

    // Returns an API to the entry (in the hash table) for the key `kmer`. The API
    // wraps the hash table position and the state value at that position.
    Kmer_Hash_Entry_API operator[](const cuttlefish::kmer_t& kmer);

    // Returns the value (in the hash-table) for the key `kmer`.
    const State operator[](const cuttlefish::kmer_t& kmer) const;

    // Attempts to update the entry (in the hash-table) for the API object according
    // to its wrapped state values, and returns true or false as per success
    // status. If the corresponding hash table position now contains a different
    // state than the one that had been read earlier, then the update fails.
    bool update(Kmer_Hash_Entry_API& api);

    // Clears the hash-table. Do not invoke on an unused object.
    void clear();
};



inline Kmer_Hash_Entry_API Kmer_Hash_Table::operator[](const cuttlefish::kmer_t& kmer)
{
    auto v = mph->lookup(kmer);
    uint64_t lidx = v / chunk_size; 
    locks_[lidx].lock();
    auto r = Kmer_Hash_Entry_API(hash_table[v]);
    locks_[lidx].unlock();
    return r;
}


inline const State Kmer_Hash_Table::operator[](const cuttlefish::kmer_t& kmer) const
{
    // NOTE: this makes the `const` a lie.  Should be a better solution here.
    auto v = mph->lookup(kmer);
    uint64_t lidx = v / chunk_size; 
    auto* tp = const_cast<Kmer_Hash_Table*>(this);
    const_cast<decltype(tp->locks_[lidx])>(tp->locks_[lidx]).lock();
    auto ve = State(hash_table[v]);
    const_cast<decltype(tp->locks_[lidx])>(tp->locks_[lidx]).unlock();
    return ve;
}


inline bool Kmer_Hash_Table::update(Kmer_Hash_Entry_API& api)
{
    auto it = &(api.bv_entry);
    uint64_t lidx = (std::distance(hash_table.begin(), it)) / chunk_size;
    locks_[lidx].lock();
    bool success = (api.bv_entry == api.get_read_state());
    if (success) {
        api.bv_entry = api.get_current_state();
    }
    locks_[lidx].unlock();
    return success;
}



#endif

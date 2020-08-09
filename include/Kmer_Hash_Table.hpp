
#ifndef KMER_HASH_TABLE_HPP
#define KMER_HASH_TABLE_HPP



#include "globals.hpp"
#include "State.hpp"
#include "Kmer_Container.hpp"
#include "BBHash/BooPHF.h"
#include "Kmer_Hasher.hpp"
#include "compact_vector/compact_vector.hpp"
#include "Kmer_Hash_Entry_API.hpp"
#include "SpinLock/SpinLock.hpp"


class CdBG;


class Kmer_Hash_Table
{
    friend class CdBG;

private:

    // Lowest bits/elem is achieved with gamma = 1, higher values lead to larger mphf but faster construction/query.
    constexpr static double gamma_factor = 2.0;

    // The MPH function.
    cuttlefish::mphf_t* mph = NULL;

    // The values (`State`) collection for the hash table;
    // keys (`kmer_t`) are passed to the MPHF, and the resulting function-value is used as index in the values table.
    cuttlefish::bitvector_t hash_table;

    // Number of locks for thread-safe access to the bitvector `hash_table`.
    constexpr static const uint64_t lock_count{65536};

    // Number of contiguous entries of the bitvector that each lock is assigned to.
    uint64_t lock_range_size;

    // The locks to maintain mutually exclusive access for threads, to the bitvector `hash_table`.
    std::array<SpinLock, lock_count> locks_;


    // Builds the minimal perfect hash function `mph` over the set of
    // k-mers present at the KMC database container `kmer_container`,
    // with `mph_file_path` being the file to use for BBHash build
    // using `thread_count` number of threads.
    void build_mph_function(const Kmer_Container& kmer_container, const uint16_t thread_count, const std::string& mph_file_path);

    // Loads an MPH function stored at the file named `file_path` into `mph`.
    void load_mph_function(const std::string& file_path);

    // Saves the MPH function `mph` into a file named `file_path`.
    void save_mph_function(const std::string& file_path) const;

    // Returns the id / number of the bucket in the hash table that is
    // supposed to store value items for the key `kmer`.
    uint64_t bucket_id(const cuttlefish::kmer_t& kmer) const;

    // Returns an API to the entry (in the hash table) for a k-mer hashing
    // to the bucket number `bucket_id` of the hash table. The API wraps
    // the hash table position and the state value at that position.
    Kmer_Hash_Entry_API operator[](const uint64_t bucket_id);

public:

    Kmer_Hash_Table()
    {}

    // Constructs a minimal perfect hash function (specifically, the BBHash) for
    // the collection of k-mers present at the KMC database at path `kmc_db_path`,
    // using up-to `thread_count` number of threads. If a non-empty path is passed
    // with `mph_file_path`, either an MPH is loaded from there (instead of building
    // from scratch), or the newly built MPH is saved there.
    void construct(const std::string& kmc_db_path, const uint16_t thread_count, const std::string& mph_file_path);

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



inline uint64_t Kmer_Hash_Table::bucket_id(const cuttlefish::kmer_t& kmer) const
{
    return mph->lookup(kmer);
}


inline Kmer_Hash_Entry_API Kmer_Hash_Table::operator[](const uint64_t bucket_id)
{
    uint64_t lidx = bucket_id / lock_range_size; 
    locks_[lidx].lock();
    auto r = Kmer_Hash_Entry_API(hash_table[bucket_id]);
    locks_[lidx].unlock();
    return r;
}


inline Kmer_Hash_Entry_API Kmer_Hash_Table::operator[](const cuttlefish::kmer_t& kmer)
{
    return operator[](mph->lookup(kmer));
}


inline const State Kmer_Hash_Table::operator[](const cuttlefish::kmer_t& kmer) const
{
    // NOTE: this makes the `const` a lie.  Should be a better solution here.
    auto v = mph->lookup(kmer);
    uint64_t lidx = v / lock_range_size; 
    auto* tp = const_cast<Kmer_Hash_Table*>(this);
    const_cast<decltype(tp->locks_[lidx])>(tp->locks_[lidx]).lock();
    auto ve = State(hash_table[v]);
    const_cast<decltype(tp->locks_[lidx])>(tp->locks_[lidx]).unlock();
    return ve;
}


inline bool Kmer_Hash_Table::update(Kmer_Hash_Entry_API& api)
{
    auto it = &(api.bv_entry);
    uint64_t lidx = (std::distance(hash_table.begin(), it)) / lock_range_size;
    locks_[lidx].lock();
    bool success = (api.bv_entry == api.get_read_state());
    if (success) {
        api.bv_entry = api.get_current_state();
    }
    locks_[lidx].unlock();
    return success;
}



#endif

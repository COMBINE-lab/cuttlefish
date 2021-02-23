
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


template <uint16_t k> class CdBG;


template <uint16_t k, uint8_t BITS_PER_KEY>
class Kmer_Hash_Table
{
    friend class CdBG<k>;

    typedef boomphf::mphf<Kmer<k>, Kmer_Hasher<k>> mphf_t;    // The MPH function type.

    typedef compact::ts_vector<cuttlefish::state_code_t, BITS_PER_KEY, uint64_t, std::allocator<uint64_t>> bitvector_t;

private:

    // Lowest bits/elem is achieved with gamma = 1, higher values lead to larger mphf but faster construction/query.
    constexpr static double GAMMA_FACTOR = 2.0;

    // The MPH function.
    mphf_t* mph = NULL;

    // The buckets collection (raw `State` representations) for the hash table structure.
    // Keys (`Kmer<k>`) are passed to the MPHF, and the resulting function-value is used as index into the buckets table.
    bitvector_t hash_table;

    // Number of locks for mutually exclusive access for threads to the same indices into the bitvector `hash_table`.
    // TODO: increase locks and check note at the end about the false `const` issue.
    constexpr static uint64_t lock_count{65536};

    // Number of contiguous entries of the bitvector that each lock is assigned to.
    // TODO: try making it `const`.
    uint64_t lock_range_size;

    // The locks to maintain mutually exclusive access for threads to the same indices into the bitvector `hash_table`.
    std::array<SpinLock, lock_count> locks_;


    // Builds the minimal perfect hash function `mph` over the set of
    // k-mers present at the KMC database container `kmer_container`,
    // with `mph_file_path` being the file to use for BBHash build
    // using `thread_count` number of threads. Uses the directory
    // at `working_dir_path` to store temporary files.
    void build_mph_function(const Kmer_Container<k>& kmer_container, uint16_t thread_count, const std::string& working_dir_path, const std::string& mph_file_path);

    // Loads an MPH function from the file at `file_path` into `mph`.
    void load_mph_function(const std::string& file_path);

    // Saves the MPH function `mph` into a file at `file_path`.
    void save_mph_function(const std::string& file_path) const;

    // Saves the hash table buckets `hash_table` into a file at `file_path`.
    void save_hash_buckets(const std::string& file_path) const;

    // Loads the hash table buckets `hash_table` from the file at `file_path`.
    void load_hash_buckets(const std::string& file_path);

    // Returns the id / number of the bucket in the hash table that is
    // supposed to store value items for the key `kmer`.
    uint64_t bucket_id(const Kmer<k>& kmer) const;

    // Returns an API to the entry (in the hash table) for a k-mer hashing
    // to the bucket number `bucket_id` of the hash table. The API wraps
    // the hash table position and the state value at that position.
    Kmer_Hash_Entry_API<BITS_PER_KEY> operator[](uint64_t bucket_id);


public:

    // Constructs a minimal perfect hash function (specifically, the BBHash) for
    // the collection of k-mers present at the KMC database at path `kmc_db_path`,
    // using up-to `thread_count` number of threads. If a non-empty path is passed
    // with `mph_file_path`, either an MPH is loaded from there (instead of building
    // from scratch), or the newly built MPH is saved there.
    void construct(const std::string& kmc_db_path, uint16_t thread_count, const std::string& working_dir_path, const std::string& mph_file_path);

    // Returns an API to the entry (in the hash table) for the key `kmer`. The API
    // wraps the hash table position and the state value at that position.
    Kmer_Hash_Entry_API<BITS_PER_KEY> operator[](const Kmer<k>& kmer);

    // Returns the value (in the hash-table) for the key `kmer`.
    State operator[](const Kmer<k>& kmer) const;

    // Attempts to update the entry (in the hash-table) for the API object according
    // to its wrapped state values, and returns `true` or `false` as per success
    // status. If the corresponding hash table position now contains a different
    // state than the one that had been read earlier, then the update fails.
    bool update(Kmer_Hash_Entry_API<BITS_PER_KEY>& api);

    // Clears the hash-table. Do not invoke on an unused object.
    void clear();
};


template <uint16_t k, uint8_t BITS_PER_KEY>
inline uint64_t Kmer_Hash_Table<k, BITS_PER_KEY>::bucket_id(const Kmer<k>& kmer) const
{
    return mph->lookup(kmer);
}


template <uint16_t k, uint8_t BITS_PER_KEY>
inline Kmer_Hash_Entry_API<BITS_PER_KEY> Kmer_Hash_Table<k, BITS_PER_KEY>::operator[](const uint64_t bucket_id)
{
    uint64_t lidx = bucket_id / lock_range_size; 
    locks_[lidx].lock();
    auto r = Kmer_Hash_Entry_API<BITS_PER_KEY>(hash_table[bucket_id]);
    locks_[lidx].unlock();
    return r;
}


template <uint16_t k, uint8_t BITS_PER_KEY>
inline Kmer_Hash_Entry_API<BITS_PER_KEY> Kmer_Hash_Table<k, BITS_PER_KEY>::operator[](const Kmer<k>& kmer)
{
    return operator[](bucket_id(kmer));
}


template <uint16_t k, uint8_t BITS_PER_KEY>
inline State Kmer_Hash_Table<k, BITS_PER_KEY>::operator[](const Kmer<k>& kmer) const
{
    // NOTE: this makes the `const` a lie.  Should be a better solution here.
    // TODO: Design a sparse-locks collection class, moving the locks array there. Have a pointer to `Sparse_Lock` in this class.
    auto v = mph->lookup(kmer);
    uint64_t lidx = v / lock_range_size; 
    auto* tp = const_cast<Kmer_Hash_Table*>(this);
    const_cast<decltype(tp->locks_[lidx])>(tp->locks_[lidx]).lock();
    auto ve = State(hash_table[v]);
    const_cast<decltype(tp->locks_[lidx])>(tp->locks_[lidx]).unlock();
    return ve;
}


template <uint16_t k, uint8_t BITS_PER_KEY>
inline bool Kmer_Hash_Table<k, BITS_PER_KEY>::update(Kmer_Hash_Entry_API<BITS_PER_KEY>& api)
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

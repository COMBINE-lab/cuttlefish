
#ifndef KMER_HASH_TABLE_HPP
#define KMER_HASH_TABLE_HPP



#include "globals.hpp"
#include "Kmer.hpp"
#include "Kmer_Hasher.hpp"
#include "State.hpp"
#include "Kmer_Hash_Entry_API.hpp"
#include "Spin_Lock.hpp"
#include "Sparse_Lock.hpp"
#include "BBHash/BooPHF.h"
#include "compact_vector/compact_vector.hpp"

#include <cstdint>
#include <cstddef>
#include <string>


class Build_Params;


template <uint16_t k, uint8_t BITS_PER_KEY>
class Kmer_Hash_Table
{
    typedef boomphf::mphf<Kmer<k>, Kmer_Hasher<k>> mphf_t;    // The MPH function type.

    typedef compact::ts_vector<cuttlefish::state_code_t, BITS_PER_KEY, uint64_t, std::allocator<uint64_t>> bitvector_t;

private:

    // The minimum gamma-value that we require for BBHash.
    static constexpr double gamma_min = 2.0;

    // The maximum gamma-value that we may use with BBHash.
    static constexpr double gamma_max = 10.0;

    // The minimum bits per hash key we require for BBHash.
    static constexpr double min_bits_per_hash_key = 3.71;

    // Empiricial bits-per-key requirement for each gamma in the range (0, 10].
    static constexpr double bits_per_gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                3.06, 3.07, 3.11, 3.16, 3.22, 3.29, 3.36, 3.44, 3.53, 3.62,
                                                3.71, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.61,
                                                4.71, 4.82, 4.92, 5.03, 5.13, 5.24, 5.35, 5.45, 5.56, 5.67,
                                                5.78, 5.89, 6.00, 6.10, 6.21, 6.32, 6.43, 6.54, 6.65, 6.76,
                                                6.87, 6.98, 7.09, 7.20, 7.31, 7.42, 7.53, 7.64, 7.75, 7.86,
                                                7.97, 8.08, 8.20, 8.31, 8.42, 8.53, 8.64, 8.75, 8.86, 8.97,
                                                9.08, 9.20, 9.31, 9.42, 9.53, 9.64, 9.75, 9.86, 9.98, 10.09,
                                                10.20, 10.31, 10.42, 10.53, 10.64, 10.76, 10.87, 10.98, 11.09, 11.20,
                                                11.31, 11.43, 11.54, 11.65, 11.76, 11.87, 11.99, 12.10, 12.21, 12.32,
                                                12.43};

    // The resolution of gamma that we support.
    static constexpr double gamma_resolution = 0.1;

    // The gamma parameter of the BBHash function.
    // Lowest bits/elem is achieved with gamma = 1, higher values lead to larger mphf but faster construction/query.
    double gamma;

    // Path to the underlying k-mer database, over which the hash table is constructed.
    const std::string kmc_db_path;

    // Number of keys (`Kmer<k>`s) in the hash table.
    const uint64_t kmer_count;

    // The MPH function.
    const mphf_t* mph = nullptr;

    // The buckets collection (raw `State` representations) for the hash table structure.
    // Keys (`Kmer<k>`) are passed to the MPHF, and the resulting function-value is used as index into the buckets table.
    bitvector_t hash_table;

    // Number of locks for mutually exclusive access for threads to the same indices into the bitvector `hash_table`.
    // TODO: increase locks and check note at the end about the false `const` issue.
    constexpr static uint64_t lock_count{65536};

    // The locks to maintain mutually exclusive access for threads to the same indices into the bitvector `hash_table`.
    mutable Sparse_Lock<Spin_Lock> sparse_lock;

    
    // Sets the `gamma` parameter of the hash function to the maximum amount so that the
    // hash table does not incur more than `max_memory` bytes of space.
    void set_gamma(std::size_t max_memory);

    // Builds the minimal perfect hash function `mph` over the set of
    // k-mers present at the KMC database container `kmer_container`,
    // using `thread_count` number of threads. Uses the directory
    // at `working_dir_path` to store temporary files. If the MPHF is
    // found present at the file `mph_file_path`, then it is loaded
    // instead.
    void build_mph_function(uint16_t thread_count, const std::string& working_dir_path, const std::string& mph_file_path);

    // Loads an MPH function from the file at `file_path` into `mph`.
    void load_mph_function(const std::string& file_path);

    // Saves the MPH function `mph` into a file at `file_path`.
    void save_mph_function(const std::string& file_path) const;


public:

    // Constructs a k-mer hash table where the table is to be built over the k-mer
    // database with path prefix `kmer_db_path`.
    Kmer_Hash_Table(const std::string& kmer_db_path);

    // Constructs a k-mer hash table where the table is to be built over the k-mer
    // database having path prefix `kmer_db_path` and `kmer_count` distinct k-mers.
    Kmer_Hash_Table(const std::string& kmc_db_path, uint64_t kmer_count);

    // Constructs a k-mer hash table where the table is to be built over the k-mer
    // database having path prefix `kmer_db_path` and `kmer_count` distinct k-mers.
    // The hash table may use at most `max_memory` bytes of memory.
    Kmer_Hash_Table(const std::string& kmc_db_path, uint64_t kmer_count, std::size_t max_memory);

    // Constructs a k-mer hash table where the table is to be built over the k-mer
    // database having path prefix `kmer_db_path` and `kmer_count` distinct k-mers.
    // The gamma factor of the BBHash MPHF of the table is attempted to be set to
    // `gamma`, if it is non-zero. Otherwise, it is set such that the the hash
    // table may use at most `max_memory` bytes of memory.
    Kmer_Hash_Table(const std::string& kmc_db_path, uint64_t kmer_count, std::size_t max_memory, double gamma);

    // Constructs a minimal perfect hash function (specifically, the BBHash) for
    // the collection of k-mers present at the KMC database at path `kmc_db_path`,
    // using up-to `thread_count` number of threads. The existence of an MPHF is
    // checked at the path `mph_file_path`—if found, it is loaded from the file.
    // If `save_mph` is specified, then the MPHF is saved into the file `mph_file_path`.
    void construct(uint16_t thread_count, const std::string& working_dir_path, const std::string& mph_file_path, const bool save_mph = false);

    // Returns the id / number of the bucket in the hash table that is
    // supposed to store value items for the key `kmer`.
    uint64_t bucket_id(const Kmer<k>& kmer) const;

    // Returns the hash value of the k-mer `kmer`.
    uint64_t operator()(const Kmer<k>& kmer) const;

    // Returns an API to the entry (in the hash table) for a k-mer hashing
    // to the bucket number `bucket_id` of the hash table. The API wraps
    // the hash table position and the state value at that position.
    Kmer_Hash_Entry_API<BITS_PER_KEY> operator[](uint64_t bucket_id);

    // Returns an API to the entry (in the hash table) for the key `kmer`. The API
    // wraps the hash table position and the state value at that position.
    Kmer_Hash_Entry_API<BITS_PER_KEY> operator[](const Kmer<k>& kmer);

    // Returns the value (in the hash-table) for the key `kmer`.
    const State operator[](const Kmer<k>& kmer) const;

    // Returns an API to the entry (in the hash table) for the key `kmer`. The API
    // wraps the hash table position and the state value at that position.
    Kmer_Hash_Entry_API<BITS_PER_KEY> at(const Kmer<k>& kmer);

    // Returns an API to the entry (in the hash table) for a k-mer hashing
    // to the bucket number `bucket_id` of the hash table. The API wraps
    // the hash table position and the state value at that position.
    Kmer_Hash_Entry_API<BITS_PER_KEY> at(uint64_t bucket_id);

    // Attempts to update the entry (in the hash-table) for the API object according
    // to its wrapped state values, and returns `true` or `false` as per success
    // status. If the corresponding hash table position now contains a different
    // state than the one that had been read earlier, then the update fails.
    bool update(Kmer_Hash_Entry_API<BITS_PER_KEY>& api);

    // Updates the state-entry in the hash-table that's at the bucket with ID
    // `bucket_id` with the state-value `state`.
    void update(uint64_t bucket_id, const State_Read_Space& state);

    // Transforms the state-entry in the hash-table that's at the bucket with ID
    // `bucket_id` through the function `transform`.
    void update(uint64_t bucket_id, cuttlefish::state_code_t (*transform)(cuttlefish::state_code_t));
    
    // Attempts to update the hash table entries for the API objects `api_1` and
    // `api_2` concurrently, i.e. both the updates need to happen in a tied manner
    // — both successful or failing. Returns `true` iff the updates succeed. If
    // either of the table positions contains a different state than the one
    // expected by the API objects, then the concurrent update fails.
    bool update_concurrent(Kmer_Hash_Entry_API<BITS_PER_KEY>& api_1, Kmer_Hash_Entry_API<BITS_PER_KEY>& api_2);

    // Returns the number of keys in the hash table.
    uint64_t size() const;

    // Clears the hash-table. Do not invoke on an unused object.
    void clear();

    // Saves the hash table buckets `hash_table` into a file at `file_path`.
    void save_hash_buckets(const std::string& file_path) const;

    // Loads the hash table buckets `hash_table` from the file at `file_path`.
    void load_hash_buckets(const std::string& file_path);

    // Saves the hash table (i.e. the hash function and the buckets) into file
    // paths determined from the parameters collection `params`.
    void save(const Build_Params& params) const;

    // Loads the hash table from disk files, determined from the parameters
    // collection `params`.
    void load(const Build_Params& params);

    // Removes the hash table files (if exists) from disk, with the file paths
    // being determined from the parameters collection `params`.
    void remove(const Build_Params& params) const;

    // Destructs the hash table.
    ~Kmer_Hash_Table();
};


template <uint16_t k, uint8_t BITS_PER_KEY>
inline uint64_t Kmer_Hash_Table<k, BITS_PER_KEY>::bucket_id(const Kmer<k>& kmer) const
{
    return mph->lookup(kmer);
}


template <uint16_t k, uint8_t BITS_PER_KEY>
inline uint64_t Kmer_Hash_Table<k, BITS_PER_KEY>::operator()(const Kmer<k>& kmer) const
{
    return bucket_id(kmer);
}


template <uint16_t k, uint8_t BITS_PER_KEY>
inline Kmer_Hash_Entry_API<BITS_PER_KEY> Kmer_Hash_Table<k, BITS_PER_KEY>::operator[](const uint64_t bucket_id)
{
    sparse_lock.lock(bucket_id);
    const Kmer_Hash_Entry_API<BITS_PER_KEY> r(hash_table[bucket_id]);
    sparse_lock.unlock(bucket_id);
    
    return r;
}


template <uint16_t k, uint8_t BITS_PER_KEY>
inline Kmer_Hash_Entry_API<BITS_PER_KEY> Kmer_Hash_Table<k, BITS_PER_KEY>::operator[](const Kmer<k>& kmer)
{
    return operator[](bucket_id(kmer));
}


template <uint16_t k, uint8_t BITS_PER_KEY>
inline const State Kmer_Hash_Table<k, BITS_PER_KEY>::operator[](const Kmer<k>& kmer) const
{
    const uint64_t bucket = bucket_id(kmer);

    sparse_lock.lock(bucket);
    const State state(hash_table[bucket]);
    sparse_lock.unlock(bucket);

    return state;
}


template <uint16_t k, uint8_t BITS_PER_KEY>
inline Kmer_Hash_Entry_API<BITS_PER_KEY> Kmer_Hash_Table<k, BITS_PER_KEY>::at(const Kmer<k>& kmer)
{
    return operator[](kmer);
}


template <uint16_t k, uint8_t BITS_PER_KEY>
inline Kmer_Hash_Entry_API<BITS_PER_KEY> Kmer_Hash_Table<k, BITS_PER_KEY>::at(const uint64_t bucket_id)
{
    return operator[](bucket_id);
}


template <uint16_t k, uint8_t BITS_PER_KEY>
inline bool Kmer_Hash_Table<k, BITS_PER_KEY>::update(Kmer_Hash_Entry_API<BITS_PER_KEY>& api)
{
    // const auto it = &(api.bv_entry);
    // const uint64_t lidx = (std::distance(hash_table.begin(), it)) / lock_range_size;
    // locks_[lidx].lock();
    // const bool success = (api.bv_entry == api.get_read_state());
    // if (success) {
    //     api.bv_entry = api.get_current_state();
    // }
    // locks_[lidx].unlock();
    // return success;

    const uint64_t bucket = std::distance(hash_table.begin(), &(api.bv_entry));

    sparse_lock.lock(bucket);
    const bool success = (api.bv_entry == api.get_read_state());
    if(success)
        api.bv_entry = api.get_current_state();
    sparse_lock.unlock(bucket);
    
    return success;
}


template <uint16_t k, uint8_t BITS_PER_KEY>
inline void Kmer_Hash_Table<k, BITS_PER_KEY>::update(const uint64_t bucket_id, const State_Read_Space& state)
{
    sparse_lock.lock(bucket_id);
    hash_table[bucket_id] = state.get_state();
    sparse_lock.unlock(bucket_id);
}


template <uint16_t k, uint8_t BITS_PER_KEY>
inline void Kmer_Hash_Table<k, BITS_PER_KEY>::update(const uint64_t bucket_id, cuttlefish::state_code_t (* const transform)(cuttlefish::state_code_t))
{
    sparse_lock.lock(bucket_id);
    hash_table[bucket_id] = transform(hash_table[bucket_id]);
    sparse_lock.unlock(bucket_id);
}


template <uint16_t k, uint8_t BITS_PER_KEY>
inline bool Kmer_Hash_Table<k, BITS_PER_KEY>::update_concurrent(Kmer_Hash_Entry_API<BITS_PER_KEY>& api_1, Kmer_Hash_Entry_API<BITS_PER_KEY>& api_2)
{
    Kmer_Hash_Entry_API<BITS_PER_KEY>* api_l = &api_1;
    Kmer_Hash_Entry_API<BITS_PER_KEY>* api_r = &api_2;
    uint64_t bucket_l = std::distance(hash_table.begin(), &(api_1.bv_entry));
    uint64_t bucket_r = std::distance(hash_table.begin(), &(api_2.bv_entry));

    // Resolution for potential deadlocks.
    if(bucket_l > bucket_r)
        std::swap(api_l, api_r),
        std::swap(bucket_l, bucket_r);


    sparse_lock.lock(bucket_l);
    bool success = (api_l->bv_entry == api_l->get_read_state());
    if(success)
    {
        sparse_lock.lock_if_different(bucket_l, bucket_r);

        success = (api_r->bv_entry == api_r->get_read_state());
        if(success)
            api_l->bv_entry = api_l->get_current_state(),
            api_r->bv_entry = api_r->get_current_state();
        
        sparse_lock.unlock_if_different(bucket_l, bucket_r);
    }
    sparse_lock.unlock(bucket_l);

    return success;
}


template <uint16_t k, uint8_t BITS_PER_KEY>
inline uint64_t Kmer_Hash_Table<k, BITS_PER_KEY>::size() const
{
    return kmer_count;
}



#endif

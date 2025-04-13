
#ifndef CONCURRENT_HASH_TABLE_HPP
#define CONCURRENT_HASH_TABLE_HPP



#include "Spin_Lock.hpp"
#include "utility.hpp"
#include "xxHash/xxhash.h"
#include "parlay/parallel.h"

#include <cstddef>
#include <vector>
#include <cstring>
#include <numeric>
#include <cstdlib>
#include <algorithm>
#include <cmath>

namespace cuttlefish
{

// =============================================================================
template <typename T_key_, typename T_val_, typename T_hasher_>
class Concurrent_Hash_Table
{
    class Iterator;
    friend class Iterator;

    static constexpr bool atomic_read = (sizeof(T_key_) <= 8);  // Whether reads are atomic. TODO: do it based on arch.

private:

    struct Key_Val_Pair
    {
        T_key_ key;
        T_val_ val;
    };

    static constexpr double lf_default = 0.75;   // Default maximum load-factor supported.

    T_key_ empty_key_;  // The empty key; currently it's set to all 1-bits.

    const T_hasher_ hash;   // The hasher object.

    const std::size_t capacity_;    // True capacity of the table; adjusted to be a power of 2.
    const std::size_t idx_wrapper_mask; // Bitmask to wrap indexing into the table.

    Key_Val_Pair* const T;  // The flat table of key-value collection.
    // std::size_t size_;  // Number of elements in the table. Probably not usable cheaply with concurrency.

    mutable std::vector<Spin_Lock> lock;    // Locks for atomic storing of key-value entries.
    // TODO: experiment moving the spin-lock to the key-value structure for seemingly better cache-behavior.

    // Maps the hash value `h` to an index into the table.
    std::size_t hash_to_idx(std::size_t h) const { return h & idx_wrapper_mask; }

    // Returns the next (wrapped) index for `i`.
    std::size_t next_index(std::size_t i) const { return hash_to_idx(i + 1); }

    // Atomically compare-and-swaps the memory location `ptr`'s value to
    // `new_key` from `old_key`. Returns `true` iff this succeeds.
    static bool CAS(T_key_* ptr, T_key_ old_key, T_key_ new_key);

    // Returns the memory-equivalent value of `val` in type `T_to`.
    template <typename T_to_, typename T_from_>
    static T_to_ pun_type(const T_from_ val) { return *reinterpret_cast<const T_to_*>(&val); }

    // Returns a 64-bit signature of the key-set of the hash table if
    // `hash_key_set_` is true. Otherwise returns a 64-bit signature of the
    // value-collection of the table.
    template <bool hash_key_set_> uint64_t signature() const;


public:

    // Constructs a concurrent hash table to support upto `max_n` elements, with
    // a maximum load-factor of `load_factor`. The object `hasher` is used to
    // hash the keys in the table.
    Concurrent_Hash_Table(std::size_t max_n, double load_factor = lf_default, T_hasher_ hasher = T_hasher_());

    ~Concurrent_Hash_Table() { deallocate(T); }

    // T_key_ empty_key() const { return empty_key_; }

    // Returns the capacity of the hash table.
    std::size_t capacity() const { return capacity_; }

    // Clears the hash table.
    // TODO: consider whether a generic empty-key might be required in our particular use-cases ever or not.
    void clear();

    // Inserts the key `key` with value `val` into the table. Returns `false` if
    // the key already exists in the table. Otherwise returns `true` iff the
    // insertion succeeds, i.e. free space was found for the insertion.
    bool insert(T_key_ key, T_val_ val);

    // Inserts the key `key` with value `val` into the table. Returns `false` if
    // the key already exists in the table; and in that case, the address of the
    // existing value for the key is stored in `val_add`. Otherwise returns
    // `true` iff the insertion succeeds, i.e. free space was found for the
    // insertion.
    bool insert(T_key_ key, T_val_ val, T_val_*& val_add);

    // Inserts the key `key` with value `val` into the table. Returns `false` if
    // the key already exists in the table; and in that case, the existing value
    // associated to `key` is overwritten with `val`. Otherwise returns `true`
    // iff the insertion succeeds, i.e. free space was found for the insertion.
    bool insert_or_assign(T_key_ key, T_val_ val);

    // Searches for `key` in the table and returns the address of the value
    // associated to it iff it is found. Returns `nullptr` otherwise.
    const T_val_* find(T_key_ key) const;

    // Searches for `key` in the table and returns the address of the value
    // associated to it iff it is found. Returns `nullptr` otherwise.
    T_val_* find(T_key_ key);

    // Searches for `key` in the table and returns `true` iff it is found. If
    // found, the associated value is stored in `val`. `val` remains unchanged
    // otherwise.
    bool find(T_key_ key, T_val_& val) const;

    // Returns an iterator for the key-value pairs in the table.
    Iterator iterator() { return Iterator(*this, 1, 0); }

    // Returns an iterator for the key-value pairs in the table that belongs to
    // a group of `it_count` iterators and has an ID `it_id` in the group.
    Iterator iterator(const std::size_t it_count, const std::size_t it_id) { return Iterator(*this, it_count, it_id); }

    // Returns a 64-bit signature of the key-set of the hash table.
    uint64_t signature() const { return signature<true>(); };

    // Returns a 64-bit signature of the values in the hash table.
    uint64_t signature_vals() const { return signature<false>(); };
};


// =============================================================================
template <typename T_key_, typename T_val_, typename T_hasher_>
class Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::Iterator
{
    friend class Concurrent_Hash_Table;

private:

    Concurrent_Hash_Table& M;   // The hash table to iterate on.
    std::size_t idx;    // Current slot-index the iterator is in.
    std::size_t end;    // End-index of the slot-range for the iterator.


    // Constructs an iterator for the hash table `M`, that would belong to a
    // group of `it_count` iterators and has an ID `it_id` in the group.
    Iterator(Concurrent_Hash_Table& M, std::size_t it_count, std::size_t it_id);


public:

    // Moves the iterator to the next key in the table. Iff some key is found
    // within its remaining range, the key and the associated value are put at
    // `key` and `val` respectively, and returns true.
    bool next(T_key_& key, T_val_& val);
};


template <typename T_key_, typename T_val_, typename T_hasher_>
inline Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::Concurrent_Hash_Table(const std::size_t max_n, const double load_factor, const T_hasher_ hasher):
      empty_key_()
    , hash(hasher)
    , capacity_(ceil_pow_2(static_cast<std::size_t>(std::ceil(max_n / load_factor))))
    , idx_wrapper_mask(capacity_ - 1)
    , T(allocate<Key_Val_Pair>(capacity_))
    , lock(capacity_)
{
    std::memset(reinterpret_cast<char*>(&empty_key_), -1, sizeof(empty_key_));

    clear();
}

template <typename T_key_, typename T_val_, typename T_hasher_>
inline void Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::clear()
{
    // TODO: try version-stamping.
    // Straightforward way.
    // parlay::parallel_for(0, capacity_, [&](std::size_t idx){ T[idx].key = empty_key_; });

    const auto byte_count = capacity_ * sizeof(Key_Val_Pair);
    const auto cache_line_count = byte_count / L1_CACHE_LINE_SIZE;
    const auto lines_per_w = cache_line_count / parlay::num_workers();
    const auto bytes_per_w = lines_per_w * L1_CACHE_LINE_SIZE;

    const auto clear_segment = [&](const std::size_t w_id)
    {
        const auto bytes_to_clear = (w_id < parlay::num_workers() - 1 ?
                                        bytes_per_w : byte_count - bytes_per_w * w_id);

        std::memset(reinterpret_cast<char*>(T) + bytes_per_w * w_id, -1, bytes_to_clear);
    };

    parlay::parallel_for(0, parlay::num_workers(), clear_segment, 1);
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline bool Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::insert(const T_key_ key, const T_val_ val)
{
#ifndef NDEBUG
    std::size_t tried_slots = 0;
#endif

    bool success = false;

    for(std::size_t i = hash_to_idx(hash(key)); ; i = next_index(i))
    {
        assert(++tried_slots <= capacity_);

    if constexpr(atomic_read)
    {
        if(T[i].key == empty_key_)  // TODO: check atomic-read / partial-read guarantees.
        {
            lock[i].lock();
            if(T[i].key == empty_key_)
                T[i].key = key, T[i].val = val,
                success = true;
            lock[i].unlock();

            if(success)
                return true;
        }

        if(T[i].key == key) // TODO: check atomic-read / partial-read guarantees.
            return false;
    }
    else
    {
        lock[i].lock();

        if(T[i].key == empty_key_)
        {
            T[i].key = key, T[i].val = val;

            lock[i].unlock();
            return true;
        }

        if(T[i].key == key)
        {
            lock[i].unlock();
            return false;
        }

        lock[i].unlock();
    }
    }

    return false;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline bool Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::insert(const T_key_ key, const T_val_ val, T_val_*& val_add)
{
#ifndef NDEBUG
    std::size_t tried_slots = 0;
#endif

    bool success = false;

    for(std::size_t i = hash_to_idx(hash(key)); ; i = next_index(i))
    {
        assert(++tried_slots <= capacity_);

    if constexpr(atomic_read)
    {
        if(T[i].key == empty_key_)  // TODO: check atomic-read / partial-read guarantees.
        {
            lock[i].lock();
            if(T[i].key == empty_key_)
                T[i].key = key, T[i].val = val,
                success = true;
            lock[i].unlock();

            if(success)
                return true;
        }

        if(T[i].key == key) // TODO: check atomic-read / partial-read guarantees.
        {
            lock[i].lock(); // Some other worker can be manipulating the value.
            val_add = &T[i].val;
            lock[i].unlock();

            return false;
        }
    }
    else
    {
        lock[i].lock();

        if(T[i].key == empty_key_)
        {
            T[i].key = key, T[i].val = val;

            lock[i].unlock();
            return true;
        }

        if(T[i].key == key)
        {
            val_add = &T[i].val;

            lock[i].unlock();
            return false;
        }

        lock[i].unlock();
    }
    }

    return false;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline bool Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::insert_or_assign(const T_key_ key, const T_val_ val)
{
#ifndef NDEBUG
    std::size_t tried_slots = 0;
#endif

    bool success = false;

    for(std::size_t i = hash_to_idx(hash(key)); ; i = next_index(i))
    {
        assert(++tried_slots <= capacity_);

    if constexpr(atomic_read)
    {
        if(T[i].key == empty_key_)  // TODO: check atomic-read / partial-read guarantees.
        {
            lock[i].lock();
            if(T[i].key == empty_key_)
                T[i].key = key, T[i].val = val,
                success = true;
            lock[i].unlock();

            if(success)
                return true;
        }

        if(T[i].key == key) // TODO: check atomic-read / partial-read guarantees.
        {
            lock[i].lock();
            T[i].val = val;
            lock[i].unlock();

            return false;
        }
    }
    else
    {
        lock[i].lock();

        if(T[i].key == empty_key_)
        {
            T[i].key = key, T[i].val = val;

            lock[i].unlock();
            return true;
        }

        if(T[i].key == key)
        {
            T[i].val = val;

            lock[i].unlock();
            return false;
        }

        lock[i].unlock();
    }
    }

    return false;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline T_val_* Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::find(const T_key_ key)
{
#ifndef NDEBUG
    std::size_t tried_slots = 0;
#endif

    T_val_* val_add = nullptr;
    for(std::size_t i = hash_to_idx(hash(key)); ; i = next_index(i))
    if constexpr(atomic_read)
    {
        if(T[i].key == key) // TODO: check atomic-read / partial-read guarantees.
        {
            lock[i].lock(); // To ensure that some other thread is not updating this val atm—a stable value is guaranteed.
                            // This only works correctly because in our use-case of this table, a key is accessed only twice.
            val_add = &T[i].val;
            lock[i].unlock();
            break;
        }
        else if(T[i].key == empty_key_) // TODO: check atomic-read / partial-read guarantees.
            break;
        else
            assert(++tried_slots <= capacity_);
    }
    else
    {
        lock[i].lock();

        if(T[i].key == key)
        {
            val_add = &T[i].val;    // This only works correctly because in our use-case of this table, a key is accessed only twice.
            lock[i].unlock();
            break;
        }
        else if(T[i].key == empty_key_)
        {
            lock[i].unlock();
            break;
        }
        else
            assert(++tried_slots <= capacity_);

        lock[i].unlock();
    }

    return val_add;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline const T_val_* Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::find(const T_key_ key) const
{
    return const_cast<Concurrent_Hash_Table*>(this)->find(key);
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline bool Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::find(const T_key_ key, T_val_& val) const
{
    const auto val_add = find(key);
    if(val_add != nullptr)
    {
        val = *val_add;
        return true;
    }

    return false;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline bool Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::CAS(T_key_* const ptr, const T_key_ old_key, const T_key_ new_key)
{
    constexpr auto key_bytes = sizeof(T_key_);
    static_assert(  key_bytes == 1 || key_bytes == 2 || key_bytes == 4 || key_bytes == 8 || key_bytes == 16,
                    "Unsupported-sized keys for CAS in hash-table.");

    if constexpr(sizeof(T_key_) == 1)
        return __sync_bool_compare_and_swap(reinterpret_cast<uint8_t*>(ptr), pun_type<uint8_t>(old_key), pun_type<uint8_t>(new_key));

    if constexpr(sizeof(T_key_) == 2)
        return __sync_bool_compare_and_swap(reinterpret_cast<uint16_t*>(ptr), pun_type<uint16_t>(old_key), pun_type<uint16_t>(new_key));

    if constexpr(sizeof(T_key_) == 4)
        return __sync_bool_compare_and_swap(reinterpret_cast<uint32_t*>(ptr), pun_type<uint32_t>(old_key), pun_type<uint32_t>(new_key));

    if constexpr(sizeof(T_key_) == 8)
        return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t*>(ptr), pun_type<uint64_t>(old_key), pun_type<uint64_t>(new_key));

    if constexpr(sizeof(T_key_) == 16)
        return __sync_bool_compare_and_swap(reinterpret_cast<__uint128_t*>(ptr), pun_type<__uint128_t>(old_key), pun_type<__uint128_t>(new_key));

    return false;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
template <bool hash_key_set_>
inline uint64_t Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::signature() const
{
    std::vector<Padded<uint64_t>> sign(parlay::num_workers(), 0);

    const auto hash =
        [&](const std::size_t idx)
        {
            if(T[idx].key != empty_key_)
            {
                if constexpr(hash_key_set_)
                    sign[parlay::worker_id()].unwrap() ^= XXH3_64bits(&T[idx].key, sizeof(T[idx].key));
                else
                    sign[parlay::worker_id()].unwrap() ^= XXH3_64bits(&T[idx].val, sizeof(T[idx].val));
            }
        };

    parlay::parallel_for(0, capacity_, hash, capacity_ / parlay::num_workers());

    return std::accumulate(sign.cbegin(), sign.cend(), 0lu, [](const uint64_t r, const Padded<uint64_t>& p_data){ return r ^ p_data.unwrap(); });
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::Iterator::Iterator(Concurrent_Hash_Table& M, std::size_t it_count, std::size_t it_id):
      M(M)
{
    const auto range_sz = (M.capacity_ + it_count - 1) / it_count;
    idx = it_id * range_sz;
    end = std::min((it_id + 1) * range_sz, M.capacity_);
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline bool Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::Iterator::next(T_key_& key, T_val_& val)
{
    for(; idx < end; idx++)
        if(M.T[idx].key != M.empty_key_)
        {
            key = M.T[idx].key, val = M.T[idx].val;
            idx++;
            return true;
        }

    return false;
}

}



#endif

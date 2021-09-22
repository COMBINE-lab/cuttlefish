
#ifndef SPARSE_LOCK_HPP
#define SPARSE_LOCK_HPP



#include <vector>
#include <cstdint>
#include <cmath>


// A collection of locks, of type `T_Lock`.
// Intended to be used when a set of sparsely distributed locks over some index range is required.
template <typename T_Lock>
class Sparse_Lock
{
private:

    // Each lock is assigned a power-of-two number of entries to guard.

    const size_t num_entries;           // Number of entries to guard.
    const uint8_t lg_per_lock_range;    // Base-2 log of the number of entries assigned to each lock.
    const size_t per_lock_range;        // Number of contiguous entries (indices) that each lock is assigned to.
    const size_t num_locks;             // Number of locks in the collection.
    std::vector<T_Lock> lock_;          // The collection of locks.


public:

    // Constructs a sparse-lock collection consisting of `lock_count` locks, for `range_size` number of entries.
    Sparse_Lock(size_t range_size, size_t lock_count);

    // Acquires lock for the entry with index `idx`.
    void lock(size_t idx);

    // Releases lock for the entry with index `idx`.
    void unlock(size_t idx);

    // Acquires lock for the entry with index `curr_idx` iff the corresponding lock for the index `prev_idx`
    // is a different lock.
    void lock_if_different(std::size_t prev_idx, std::size_t curr_idx);

    // Releases lock for the entry with index `curr_idx` iff the corresponding lock for the index `prev_idx`
    // is a different lock.
    void unlock_if_different(std::size_t prev_idx, std::size_t curr_idx);
};


template <typename T_Lock>
inline Sparse_Lock<T_Lock>::Sparse_Lock(const size_t range_size, const size_t lock_count):
    num_entries(range_size),
    lg_per_lock_range(static_cast<uint8_t>(std::floor(std::log2((num_entries + lock_count - 1) / lock_count)))),
    per_lock_range(static_cast<size_t>(1) << lg_per_lock_range),
    num_locks((num_entries + per_lock_range - 1) / per_lock_range),
    lock_(num_locks)
{}


template <typename T_Lock>
inline void Sparse_Lock<T_Lock>::lock(const size_t idx)
{
    lock_[idx >> lg_per_lock_range].lock();
}


template <typename T_Lock>
inline void Sparse_Lock<T_Lock>::unlock(const size_t idx)
{
    lock_[idx >> lg_per_lock_range].unlock();
}


template <typename T_Lock>
inline void Sparse_Lock<T_Lock>::lock_if_different(const std::size_t prev_idx, const std::size_t curr_idx)
{
    if((curr_idx >> lg_per_lock_range) != (prev_idx >> lg_per_lock_range))
        lock_[curr_idx >> lg_per_lock_range].lock();
}


template <typename T_Lock>
inline void Sparse_Lock<T_Lock>::unlock_if_different(const std::size_t prev_idx, const std::size_t curr_idx)
{
    if((curr_idx >> lg_per_lock_range) != (prev_idx >> lg_per_lock_range))
        lock_[curr_idx >> lg_per_lock_range].unlock();
}



#endif

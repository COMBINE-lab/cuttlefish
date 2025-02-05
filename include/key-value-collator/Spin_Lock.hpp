
#ifndef SPIN_LOCK_HPP
#define SPIN_LOCK_HPP



#include <atomic>


namespace key_value_collator
{


// A lightweight lock-free mutex class.
// It is based on `std::atomic_flag`, which is guaranteed to be a lock-free atomic construct .
// Reference: https://en.cppreference.com/w/cpp/atomic/atomic_flag
class Spin_Lock
{
private:

    std::atomic_flag lock_ = ATOMIC_FLAG_INIT;


public:

    // Acquires the lock for mutually-exlcusive access to it.
    void lock();

    // Releases the lock, giving up the exclusive access to it.
    void unlock();
};


inline void Spin_Lock::lock()
{
    // Due to the memory access order `memory_order_acquire`, no reads or writes in the current thread can be
    // reordered before this load of the variable `lock_` (enforced by the compiler and the processor) —
    // ensuring that memory-access instructions after a `lock` invokation stays after it.

    while(lock_.test_and_set(std::memory_order_acquire))
        ;// while(lock_.test(std::memory_order_relaxed));   // C++20 optimization to avoid the redundant stores from the spinning `test_and_set`.
}


inline void Spin_Lock::unlock()
{
    // Due to the memory access order `memory_order_release`, no reads or writes in the current thread can be
    // reordered after this store of the variable `lock_` (enforced by the compiler and the processor) —
    // ensuring that memory-access instructions before an `unlock` invokation stays before it.

    lock_.clear(std::memory_order_release);
}


}



#endif

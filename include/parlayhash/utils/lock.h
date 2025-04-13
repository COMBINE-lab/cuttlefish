#include <atomic>
#include <cstddef>
#include <vector>

#ifndef PARLAYLOCK_H_
#define PARLAYLOCK_H_

namespace parlay {

// creates 2^16 lock slots.
// locks.try_lock(i, f) will hash i to the h(i) % 2^16th lock.
// If the lock is not taken then f is run and the try_lock returns the
// boolean result of f then releasing the lock.   Otherwise it returns false.
struct lock_set {
private:
  using lck = std::atomic<bool>;
  const int bucket_bits = 16;
  const size_t mask = ((1ul) << bucket_bits) - 1;
  std::vector<lck> locks;

  static inline uint64_t hash64(uint64_t x) {
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    x = x ^ (x >> 31);
    return x;
  }
public:
  template <typename F>
  bool try_lock(long i, F f) {
    bool old = false;
    bool result = false;
    lck& x = locks[hash64(i) & mask];
    if (x.compare_exchange_strong(old, true)) {
      result = f();
      x = false;
    }
    return result;
  }
  lock_set() : locks(std::vector<lck>(1ul << bucket_bits)) {
    std::fill(locks.begin(), locks.end(), false);
  }
};

  extern inline lock_set& get_locks() {
    static lock_set locks;
    return locks;
  }

      template <typename F>
  auto try_loop(const F& f, int delay = 200, const int max_multiplier = 20) {
    int multiplier = 1;
    int cnt = 0;
    while (true)  {
      if (cnt++ == 10000000000ul/(delay*max_multiplier)) {
  	std::cerr << "problably in an infinite retry loop" << std::endl;
  	abort(); 
      }
      auto r = f();
      if (r.has_value()) return *r;
      multiplier = std::min(2*multiplier, max_multiplier);
      for (volatile int i=0; i < delay * multiplier; i++);
    }
  }

}

#endif // PARLAYLOCK_H_

// ***************************
// Epoch-based memory reclamation
// Supports:
//     epoch::with_epoch(F f),
// which runs f within an epoch, as well as:
//     epoch::New<T>(args...)
//     epoch::Retire(T* a)   -- delays destruction and free
//     epoch::Delete(T* a)   -- destructs and frees immediately
// Retire delays destruction and free until no operation that was in a
// with_epoch at the time it was run is still within the with_epoch.
//
// All operations take constant time overhead (beyond the cost of the
// system malloc and free).
//
// Designed to work with C++ threads, or compatible threading
// libraries.  In particular it uses thread_local variables, and no two
// concurrent processes can share the same instance of the  variable.
//
// When NDEBUG is not set, the operations check for memory corruption
// of the bytes immediately before and after the object, and check
// for double retires/deletes.  Also:
//     epoch::check_ptr(T* a)
// will check that an object allocated using epoch::New(..) has not
// been corrupted.
//
// Supports undoing retires.  This can be useful in transactional
// system in which an operation aborts, and any retires done during
// the operations have to be undone.  In particular Retire returns a
// pointer to a boolean.  Running
//    epoch::undo_retire(bool* x)
// will undo the retire.  Must be run in same with_epoch as the retire
// was run, otherwise it is too late to undo.  If you don't want
// to undo retires, you can ignore this feature.
//
// New<T>, Retire and Delete use a shared pool for the retired lists,
// which, although not very large, is not cleared until program
// termination.  A private pool can be created with
//     epoch::memory_pool<T> a;
// which then supports a->New(args...), a->Retire(T*) and
// a->Delete(T*).  On destruction of "a", all elements of the retired
// lists will be destructed and freed.
//
// Achieves constant times overhead by incrementally taking steps.
// In particular every Retire takes at most a constant number of
// incremental steps towards updating the epoch and clearing the
// retired lists.
//
// Developed as part of parlay project at CMU, initially for flock then
// used for verlib, and parlayhash.
// Current dependence on parlay is just for parlay::my_thread_id() and
// parlay::num_thread_ids() which are from "parlay/thread_specific.h".
// ***************************

#include <atomic>
#include <cstdlib>
#include <functional>
#include <list>
#include <ostream>
#include <vector>
#include <type_traits>
#include <utility>
// Needed for parlay::my_thread_id of parlay::num_thread_ids
#include "parlay/thread_specific.h"

#ifndef PARLAY_EPOCH_H_
#define PARLAY_EPOCH_H_

#ifndef NDEBUG
// Checks for corruption of bytes before and after allocated structures, as well as double frees.
// Requires some extra memory to pad the front and back of a structure.
#define EpochMemCheck 1
#endif
//#define EpochMemCheck 1

#define USE_STEPPING 1
//#define USE_UNDO 1

#ifdef USE_PARLAY_ALLOC
#include "parlay/alloc.h"
#endif

// ***************************
// epoch structure
// ***************************

namespace epoch {

  namespace internal {

  inline int worker_id() {return parlay::my_thread_id(); }
  inline int num_workers() {return parlay::num_thread_ids();}
  constexpr int max_num_workers = 1024;

struct alignas(64) epoch_s {

  // functions to run when epoch is incremented
  std::vector<std::function<void()>> before_epoch_hooks;
  std::vector<std::function<void()>> after_epoch_hooks;

  struct alignas(64) announce_slot {
    std::atomic<long> last;
    announce_slot() : last(-1l) {}
  };

  std::vector<announce_slot> announcements;
  std::atomic<long> current_epoch;
  epoch_s() :
    announcements(std::vector<announce_slot>(max_num_workers)),
    current_epoch(0),
    epoch_state(0) {}

  long get_current() {
    return current_epoch.load();
  }

  long get_my_epoch() {
    return announcements[worker_id()].last;
  }

  void set_my_epoch(long e) {
    announcements[worker_id()].last = e;
  }

  int announce() {
    size_t id = worker_id();
    assert(id < max_num_workers);
    while (true) {
      long current_e = get_current();
      long tmp = current_e;
      // apparently an exchange is faster than a store (write and fence)
      announcements[id].last.exchange(tmp, std::memory_order_seq_cst);
      if (get_current() == current_e) return id;
    }
  }

  void unannounce(size_t id) {
    announcements[id].last.store(-1l, std::memory_order_release);
  }

  // top 16 bits are used for the process id, and the bottom 48 for
  // the epoch number
  using state = size_t;
  std::atomic<state> epoch_state;

  // Attempts to takes num_steps checking the announcement array to
  // see that all slots are up-to-date with the current epoch.  Once
  // they are, the epoch is updated.  Designed to deamortize the cost
  // of sweeping the announcement array--every thread only does
  // constant work.
  state update_epoch_steps(state prev_state, int num_steps) {
    state current_state = epoch_state.load();
    if (prev_state != current_state)
      return current_state;
    size_t i = current_state >> 48;
    size_t current_e = ((1ul << 48) - 1) & current_state;
    size_t workers = num_workers();
    if (i == workers) {
      for (const auto h : before_epoch_hooks) h();
      long tmp = current_e;
      if (current_epoch.load() == current_e &&
	  current_epoch.compare_exchange_strong(tmp, current_e+1)) {
	for (const auto h : after_epoch_hooks) h();
      }
      state new_state = current_e + 1;
      epoch_state.compare_exchange_strong(current_state, new_state);
      return epoch_state.load();
    }
    size_t j;
    for (j = i ; j < i + num_steps && j < workers; j++)
      if ((announcements[j].last != -1l) && announcements[j].last < current_e)
	return current_state;
    state new_state = (j << 48 | current_e);
    if (epoch_state.compare_exchange_strong(current_state, new_state))
      return new_state;
    return current_state;
  }

  // this version does the full sweep
  void update_epoch() {
    long current_e = get_current();

    // check if everyone is done with earlier epochs
    int workers;
    do {
      workers = num_workers();
      if (workers > max_num_workers) {
	std::cerr << "number of threads: " << workers
		  << ", greater than max_num_threads: " << max_num_workers << std::endl;
	abort();
      }
      for (int i=0; i < workers; i++)
	if ((announcements[i].last != -1l) && announcements[i].last < current_e)
	  return;
    } while (num_workers() != workers); // this is unlikely to loop

    // if so then increment current epoch
    for (const auto h : before_epoch_hooks) h();
    if (current_epoch.compare_exchange_strong(current_e, current_e+1)) {
      for (const auto h : after_epoch_hooks) h();
    }
  }
};

  // Juat one epoch structure shared by all
  extern inline epoch_s& get_epoch() {
    static epoch_s epoch;
    return epoch;
  }

// ***************************
// type specific memory pools
// ***************************

template <typename T>
struct alignas(64) memory_pool {
private:

  struct list_entry {
    T* ptr;
#ifdef USE_UNDO
    bool keep_;
    bool keep() {return keep_;}
    list_entry() : keep_(false) {}
    list_entry(T* ptr) : ptr(ptr), keep_(false) {}
#else
    bool keep() {return false;}
#endif
  };

  // each thread keeps one of these
  struct alignas(256) old_current {
    std::list<list_entry> old;  // linked list of retired items from previous epoch
    std::list<list_entry> current; // linked list of retired items from current epoch
    std::list<list_entry> reserve;  // linked list of items that could be destructed, but delayed so they can be reused
    long epoch; // epoch on last retire, updated on a retire
    long retire_count; // number of retires so far, reset on updating the epoch
    long alloc_count;
    epoch_s::state e_state;
    old_current() : e_state(0), epoch(0), retire_count(0), alloc_count(0) {}
  };

  std::vector<old_current> pools;
    
  // wrapper used so can pad for the memory checked version
  struct wrapper {
#ifdef EpochMemCheck    
    long pad;
    std::atomic<long> head;
    T value;
    std::atomic<long> tail;
#else
    T value;
#endif
  };

  // values used to check for corruption or double delete
  static constexpr long default_val = 10;
  static constexpr long deleted_val = 55;

  // given a pointer to value in a wrapper, return a pointer to the wrapper.
  wrapper* wrapper_from_value(T* p) {
    size_t offset = ((char*) &((wrapper*) p)->value) - ((char*) p);
    return (wrapper*) (((char*) p) - offset);
  }

  // destructs entries on a list
  void clear_list(std::list<list_entry>& lst) {
    for (list_entry& x : lst)
      if (!x.keep()) {
	x.ptr->~T();
	free_wrapper(wrapper_from_value(x.ptr));
      }
    lst.clear();
  }

  void advance_epoch(int i, old_current& pid) {
#ifndef USE_UNDO
    int delay = 1;
#else
    int delay = 2;
#endif
    if (pid.epoch + delay < get_epoch().get_current()) {
      pid.reserve.splice(pid.reserve.end(), pid.old);
      pid.old = std::move(pid.current);
      pid.epoch = get_epoch().get_current();
    }
    // a heuristic
#ifdef USE_STEPPING
    long update_threshold = 10;
#else
    long update_threshold = 10 * num_workers();
#endif
    if (++pid.retire_count == update_threshold) {
      pid.retire_count = 0;
#ifdef USE_STEPPING
      pid.e_state = get_epoch().update_epoch_steps(pid.e_state, 8);
#else
      get_epoch().update_epoch();
#endif
    }
  }

#ifdef USE_PARLAY_ALLOC
  using Allocator = parlay::type_allocator<wrapper>;
#endif

  void check_wrapper_on_destruct(wrapper* x) {
#ifdef EpochMemCheck
    // check nothing is corrupted or double deleted
    if (x->head != default_val || x->tail != default_val) {
      if (x->head == deleted_val) std::cerr << "double free" << std::endl;
      else if (x->head != default_val)  std::cerr << "corrupted head" << x->head << std::endl;
      if (x->tail != default_val) std::cerr << "corrupted tail: " << x->tail << std::endl;
      abort();
    }
    x->head = deleted_val;
#endif
  }

  void set_wrapper_on_construct(wrapper* x) {
#ifdef EpochMemCheck
    x->pad = x->head = x->tail = default_val;
#endif
  }

  void free_wrapper(wrapper* x) {
    check_wrapper_on_destruct(x);
#ifdef USE_PARLAY_ALLOC
    return Allocator::freex);
#else
    return std::free(x);
#endif
  }
  
  wrapper* allocate_wrapper() {
    auto &pid = pools[worker_id()];
    if (!pid.reserve.empty()) {
      list_entry x = pid.reserve.front();
      pid.reserve.pop_front();
      if (!x.keep()) {
	x.ptr->~T();
	wrapper* w = wrapper_from_value(x.ptr);
	check_wrapper_on_destruct(w);
	set_wrapper_on_construct(w);
	return w;
      }
    }
#ifdef USE_PARLAY_ALLOC
    wrapper* w = Allocator::alloc();
#else
    wrapper* w = (wrapper*) std::malloc(sizeof(wrapper));
#endif
    set_wrapper_on_construct(w);
    return w;
  }

 public:
  memory_pool() {
    long workers = max_num_workers;
    pools = std::vector<old_current>(workers);
    for (int i = 0; i < workers; i++) {
      pools[i].retire_count = 0;
    }
  }

  memory_pool(const memory_pool&) = delete;
  ~memory_pool() { clear(); }

  // for backwards compatibility
  void acquire(T* p) { }

  template <typename ... Args>
  T* New(Args... args) {
    wrapper* x = allocate_wrapper();
    T* newv = &x->value;
    new (newv) T(args...);
    return newv;
  }

  // f is a function that initializes a new object before it is shared
  template <typename F, typename ... Args>
  T* New_Init(F f, Args... args) {
    T* x = New(args...);
    f(x);
    return x;
  }

  // retire and return a pointer if want to undo the retire
#ifdef USE_UNDO
  bool* Retire(T* p) {
#else
  void Retire(T* p) {
#endif
    auto i = worker_id();
    auto &pid = pools[i];
    if (pid.reserve.size() > 500) {
      list_entry x = pid.reserve.front();
      if (!x.keep()) {
	x.ptr->~T();
	free_wrapper(wrapper_from_value(x.ptr));
      }
      pid.reserve.pop_front();
    }
    advance_epoch(i, pid);
    pid.current.push_back(list_entry{p});
#ifdef USE_UNDO
    return &pid.current.back().keep_;
#endif
  }

  // destructs and frees the object immediately
  void Delete(T* p) {
    p->~T();
    free_wrapper(wrapper_from_value(p));
  }

  bool check_ptr(T* ptr, bool silent=false) {
#ifdef EpochMemCheck
    if (ptr == nullptr) return true;
    wrapper* x = wrapper_from_value(ptr);
    if (!silent) {
      if (x->pad != default_val) std::cerr << "memory_pool, check: pad word corrupted" << x->pad << std::endl;
      if (x->head != default_val) std::cerr << "memory_pool, check: head word corrupted" << x->head << std::endl;
      if (x->tail != default_val) std::cerr << "memory_pool, check: tail word corrupted: " << x->tail << std::endl;
    }
    return (x->pad == default_val && x->head == default_val && x->tail == default_val);
#endif
    return true;
  }

  // Clears all the lists, to be used on termination, or could be use
  // at a quiescent point when noone is reading any retired items.
  void clear() {
    // for (int i=0; i < num_workers(); i++)
    //   std::cout << i << ": " << pools[1].old.size() << ", "
    // 		<< pools[i].current.size() << ", "
    // 		<< pools[i].reserve.size() << std::endl;
    get_epoch().update_epoch();
    for (int i=0; i < num_workers(); i++) {
      clear_list(pools[i].old);
      clear_list(pools[i].current);
      clear_list(pools[i].reserve);
    }
    //Allocator::print_stats();
  }

  void stats() {}
};
  
template <typename T>
struct alignas(64) retire_pool {
private:

  struct list_entry {
    char data[sizeof(T)];
  };

  // each thread keeps one of these
  struct alignas(256) old_current {
    std::list<list_entry> old;  // linked list of retired items from previous epoch
    std::list<list_entry> current; // linked list of retired items from current epoch
    long epoch; // epoch on last retire, updated on a retire
    long retire_count; // number of retires so far, reset on updating the epoch
    epoch_s::state e_state;
    old_current() : e_state(0), epoch(0), retire_count(0) {}
  };

  std::vector<old_current> pools;

  // destructs entries on a list
  void clear_list(std::list<list_entry>& lst) {
    for (list_entry& x : lst)
      ((T*) (&(x.data)))->~T();
    lst.clear();
  }

  void advance_epoch(int i, old_current& pid) {
    if (pid.epoch + 1 < get_epoch().get_current()) {
      clear_list(pid.old);
      pid.old = std::move(pid.current);
      pid.epoch = get_epoch().get_current();
    }
#ifdef USE_STEPPING
    long update_threshold = 10;
#else
    long update_threshold = 10 * num_workers();
#endif
    if (++pid.retire_count == update_threshold) {
      pid.retire_count = 0;
#ifdef USE_STEPPING
      pid.e_state = get_epoch().update_epoch_steps(pid.e_state, 8);
#else
      get_epoch().update_epoch();
#endif
    }
  }

 public:
  retire_pool() {
    long workers = max_num_workers;
    pools = std::vector<old_current>(workers);
    for (int i = 0; i < workers; i++) 
      pools[i].retire_count = 0;
  }

  retire_pool(const retire_pool&) = delete;
  ~retire_pool() { clear(); }

  void Retire(T* p) {
    auto i = worker_id();
    auto &pid = pools[i];
    advance_epoch(i, pid);
    list_entry x;
    strncpy(x.data, (char*) p, sizeof(T));
    pid.current.push_back(x);
  }

  // Clears all the lists, to be used on termination, or could be use
  // at a quiescent point when noone is reading any retired items.
  void clear() {
    get_epoch().update_epoch();
    for (int i=0; i < num_workers(); i++) {
      clear_list(pools[i].old);
      clear_list(pools[i].current);
    }
  }

  void stats() {}
};

} // namespace internal

// ***************************
// The public interface
// ***************************
  
  // x should point to the skip field of a link
  inline void undo_retire(bool* x) { *x = true;}
  
  template <typename T>
  using memory_pool = internal::memory_pool<T>;

  template <typename T>
  extern inline memory_pool<T>& get_default_pool() {
    static memory_pool<T> pool;
    return pool;
  }

  template <typename T>
  using retire_pool = internal::retire_pool<T>;

  template <typename T>
  extern inline retire_pool<T>& get_default_retire_pool() {
    static retire_pool<T> pool;
    return pool;
  }

  template <typename T, typename ... Args>
  static T* New(Args... args) {
    return get_default_pool<T>().New(std::forward<Args>(args)...);}

  template <typename T>
  static void Delete(T* p) {get_default_pool<T>().Delete(p);}

  template <typename T>
#ifdef USE_UNDO
  static bool* Retire(T* p) {return get_default_pool<T>().Retire(p);}
#else
  void Retire(T* p) {return get_default_pool<T>().Retire(p);}
#endif
    
  template <typename T>
  static bool check_ptr(T* p, bool silent=false) {
    return get_default_pool<T>().check_ptr(p, silent);}

  template <typename T>
  static void clear() {get_default_pool<T>().clear();}

  //template <typename T>
  //static void stats() {get_default_pool<T>().stats();}

  template <typename Thunk>
  auto with_epoch(Thunk f) {
    int id = internal::get_epoch().announce();
    if constexpr (std::is_void_v<std::invoke_result_t<Thunk>>) {
      f();
      internal::get_epoch().unannounce(id);
    } else {
      auto v = f();
      internal::get_epoch().unannounce(id);
      return v;
    }
  }

} // end namespace epoch

#endif //PARLAY_EPOCH_H_

#ifdef USE_PARLAY
#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <parlay/delayed.h>
namespace parlay {
#define PARLAY_USE_STD_ALLOC 1

  using scheduler_type = internal::scheduler_type;

  template <typename F>
  long tabulate_reduce(long n, const F& f) {
    return parlay::reduce(parlay::delayed::tabulate(n, [&] (size_t i) {
	     return f(i);}));
  }
}
#else
namespace parlay {

  struct scheduler_type {
    scheduler_type(int num_procs) {}
  };

  template <typename F>
  long tabulate_reduce(long n, const F& f) {
    long r = 0;
    for (long i=0; i < n; i++)
      r += f(i);
    return r;
  }
  
  template <typename F>
  void parallel_for(long n, const F& f) {
    for (long i=0; i < n; i++) f(i);
  }
}
#endif

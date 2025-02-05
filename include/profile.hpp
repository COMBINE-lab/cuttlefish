
#ifndef PROFILE_HPP
#define PROFILE_HPP



#include <string>
#include <functional>


// Executes `f()` with optional profiling at record-file `tag`.
// Empty `__VA_ARGS__` will produce ISO C++11 warning, which can be fixed with C++20 `__VA_OPTS__(,)`.
#ifdef PART_PROFILE
    #define EXECUTE(tag, f, ...) (cuttlefish::profile([&](){ f(__VA_ARGS__); }, (tag)));
#else
    #define EXECUTE(tag, f, ...) (f(__VA_ARGS__));
#endif


namespace cuttlefish
{
    // `perf`-profiles the execution of the function `f` at file `record`.
    void profile(const std::function<void()>& f, const std::string& record);
};



#endif


#ifndef HYPERLOGLOG_HPP
#define HYPERLOGLOG_HPP



#include "utility.hpp"
#include "parlay/parallel.h"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <cmath>


namespace cuttlefish
{

// =============================================================================
// A class to estimate the cardinality of a data stream in parallel at accuracy
// i.e. SD 4.6%, with the HyperLogLog algorithm. Data must be provided as a
// stream of its hashes, which must be uniform for the accuracy bound to hold.
class HyperLogLog
{

private:

    // Number of substreams in the estimation process; SD of estimation is 1.03896 / sqrt(m).
    static constexpr uint32_t m = 512;
    static constexpr uint32_t log_m = 9;

    std::vector<Padded<uint8_t[m]>> M_w;    // `M_w[i]` contains the 'log'-registers for worker `i`.


public:

    // Constructs a Hyperloglog cardinality-estimator.
    HyperLogLog();

    // Adds the 32-bit hash `h` of a data item to the estimator.
    void add(uint32_t h);

    // Returns the cardinality estimation of the added stream of hashes.
    uint64_t estimate() const;
};


inline void HyperLogLog::add(const uint32_t h)
{
    constexpr auto substream_mask = m - 1;

    const auto stream = h & substream_mask;
    const auto h_proxy = h >> log_m;
    const auto tz = (h_proxy != 0 ? __builtin_ctz(h_proxy) : 32);
    auto& M = M_w[parlay::worker_id()].unwrap();
    M[stream] = std::max(M[stream], static_cast<uint8_t>(1 + tz));
}

}



#endif

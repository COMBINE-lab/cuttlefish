
#include "HyperLogLog.hpp"

#include <cstring>


namespace cuttlefish
{

HyperLogLog::HyperLogLog():
      M_w(parlay::num_workers())
{
    static_assert(m >= 128);
    static_assert((1lu << log_m) == m);

    for(auto& M : M_w)
        std::memset(M.unwrap(), 0, m);
}


uint64_t HyperLogLog::estimate() const
{
    uint8_t M[m];
    std::memset(M, 0, m);

    for(std::size_t i = 0; i < m; ++i)
        for(std::size_t w_id = 0; w_id < parlay::num_workers(); ++w_id)
            M[i] = std::max(M[i], M_w[w_id].unwrap()[i]);


    double est = 0;
    for(std::size_t j = 0; j < m; ++j)
        est += 1.0 / (1llu << M[j]);

    constexpr double alpha = 0.7213 / (1 + 1.079 / m);   // Factor to correct systemic multiplicative bias in estimation.

    est = alpha * m * m * (1 / est);
    if(est <= 2.5 * m)  // Small range.
    {
        uint32_t v = 0;
        for(std::size_t i = 0; i < m; ++i)
            v += (M[i] == 0);

        if(v > 0)
            est = m * std::log(m / static_cast<double>(v));
    }
    else if(est > ((1llu << 32) / 30.0)) // Long range.
        est = -1.0 * (1llu << 32) * std::log(1 - (est / (1llu << 32)));

    return est;
}

}

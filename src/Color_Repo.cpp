
#include "Color_Repo.hpp"
#include "parlay/parallel.h"

#include <cstdint>
#include <algorithm>


namespace cuttlefish
{

void Color_Repo::init(const std::string& path)
{
    B.reserve(parlay::num_workers());
    for(uint32_t w_id = 0; w_id < parlay::num_workers(); ++w_id)
        B.emplace_back(bucket_t(path + "." + std::to_string(w_id), 32 * 1024));
}


typename Color_Repo::bucket_t& Color_Repo::bucket()
{
    assert(B.size() == parlay::num_workers());
    return B[parlay::worker_id()].unwrap();
}


std::size_t Color_Repo::bytes() const
{
    std::size_t sz = 0;
    std::for_each(B.cbegin(), B.cend(), [&](const auto& b)
    {
        sz += b.unwrap().size();
    });

    return sz * sizeof(source_id_t);
}


void Color_Repo::close()
{
    parlay::parallel_for(0, parlay::num_workers(),
        [&](const std::size_t w)
        {
            B[w].unwrap().serialize();
        }, 1);
}

}

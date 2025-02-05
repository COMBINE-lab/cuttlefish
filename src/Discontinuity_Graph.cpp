
#include "Discontinuity_Graph.hpp"
#include "Minimizer_Iterator.hpp"
#include "Data_Logistics.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"


namespace cuttlefish
{


template <uint16_t k, bool Colored_>
Discontinuity_Graph<k, Colored_>::Discontinuity_Graph(const Build_Params& params, const Data_Logistics& logistics):
      min_len(params.min_len())
    , E_(params.vertex_part_count(), logistics.edge_matrix_path())
    , lmtigs(logistics.lmtig_buckets_path(), params.lmtig_bucket_count(), parlay::num_workers(), Colored_)
    , phantom_edge_count_(0)
    , max_source_id_(logistics.input_paths_collection().size())
{
    if constexpr(Colored_)
    {
        vertex_color_map_.reserve(lmtigs.bucket_count());
        vertex_color_map_.emplace_back(std::string());
        for(std::size_t b = 1; b < lmtigs.bucket_count(); ++b)
            vertex_color_map_.emplace_back(logistics.lmtig_buckets_path() + "/" + std::to_string(b) + ".col");
    }
}


template <uint16_t k, bool Colored_>
Discontinuity_Graph<k, Colored_>::Discontinuity_Graph(cereal::BinaryInputArchive& archive):
      min_len()
    , E_(archive)
    , lmtigs(archive)
    , max_source_id_()
{
    archive(*this);
}


template <uint16_t k, bool Colored_>
uint64_t Discontinuity_Graph<k, Colored_>::phantom_edge_upper_bound() const
{
    return phantom_edge_count_;
}


template <uint16_t k, bool Colored_>
void Discontinuity_Graph<k, Colored_>::close()
{
    E_.close();
    lmtigs.close();
}


template <uint16_t k, bool Colored_>
std::size_t Discontinuity_Graph<k, Colored_>::vertex_part_size_upper_bound() const
{
    std::size_t bound = 0;
    for(std::size_t j = 1; j <= E_.vertex_part_count(); ++j)
        // Each *original* edge from the non-diagonal blocks of column `j` corresponds to a
        // unique vertex of partition `j`, and in the worst-case, each original edge from the
        // diagonal block corresponds to two unique vertices.
        bound = std::max(bound, E_.col_size(j) + E_.block_size(j, j));

    return bound;
}


template <uint16_t k, bool Colored_>
bool Discontinuity_Graph<k, Colored_>::is_discontinuity(const char* const seq) const
{
    minimizer_t min_l, min_r;
    uint64_t h_l, h_r;
    std::size_t idx_l, idx_r;

    typedef Minimizer_Iterator<const char*, k - 1, true> min_it_t;
    min_it_t::minimizer(seq, min_len, min_seed, min_l, h_l, idx_l);
    min_it_t::minimizer(seq + 1, min_len, min_seed, min_r, h_r, idx_r);

    // const auto graph_id = [&](const auto h){ return h & (params.subgraph_count() - 1); };
    // TODO: update.
    const auto graph_id = [&](const auto h){ return h & (128 * 128 - 1); };
    return graph_id(h_l) != graph_id(h_r);
}


template <uint16_t k, bool Colored_>
std::size_t Discontinuity_Graph<k, Colored_>::RSS() const
{
    std::size_t v_c_map_bytes = 0;
    std::for_each(vertex_color_map_.cbegin(), vertex_color_map_.cend(), [&](auto& v){ v_c_map_bytes += v.unwrap().RSS(); });

    return E_.RSS() + lmtigs.RSS() + v_c_map_bytes;
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BOOL, cuttlefish::Discontinuity_Graph)

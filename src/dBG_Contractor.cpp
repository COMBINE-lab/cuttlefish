
#include "dBG_Contractor.hpp"
#include "Discontinuity_Graph.hpp"
#include "Graph_Partitioner.hpp"
#include "State_Config.hpp"
#include "Subgraphs_Manager.hpp"
#include "Discontinuity_Graph_Contractor.hpp"
#include "Contracted_Graph_Expander.hpp"
#include "Unitig_Collator.hpp"
#include "globals.hpp"
#include "profile.hpp"
#include "parlay/parallel.h"


namespace cuttlefish
{

template <uint16_t k>
dBG_Contractor<k>::dBG_Contractor(const Build_Params& params):
      params(params)
    , logistics(params)
    , op_buf(parlay::num_workers(), op_buf_t(output_sink.sink()))
{
    Edge_Frequency::set_edge_threshold(params.cutoff());
    std::cerr << "Edge frequency cutoff: " << params.cutoff() << ".\n";
}


template <uint16_t k>
template <bool Colored_>
void dBG_Contractor<k>::construct()
{
    // Clear the output file and initialize the output sink.
    const auto op_file_path = logistics.output_file_path();
    clear_file(op_file_path);
    output_sink.init_sink(op_file_path);

    Discontinuity_Graph<k, Colored_> gamma(params, logistics);  // The discontinuity graph.

    const auto t_0 = timer::now();
    decltype(timer::now()) t_part;

{
    Subgraphs_Manager<k, Colored_> G(logistics, params.min_len(), gamma, op_buf);

    if(params.is_read_graph())
    {
        EXECUTE("partition", (Graph_Partitioner<k, true, Colored_>(G, logistics, params.min_len())).partition)
    }
    else
    {
        EXECUTE("partition", (Graph_Partitioner<k, false, Colored_>(G, logistics, params.min_len())).partition)
    }

    G.finalize();

    t_part = timer::now();
    std::cerr << "Sequence splitting into subgraphs completed. Time taken: " << timer::duration(t_part - t_0) << " seconds.\n";

    {
        EXECUTE("subgraphs", G.process)

        std::cerr << "Trivial maximal unitig count: " << G.trivial_mtig_count() << ".\n";
        std::cerr << "Trivial ICC count: " << G.icc_count() << ".\n";
    }
}

    const auto t_subg = timer::now();
    std::cerr << "Subgraphs construction and contraction completed. Time taken: " << timer::duration(t_subg - t_part) << " seconds.\n";

    std::cerr << "Edge-matrix size: " << gamma.E().size() << "\n";
    std::cerr << "Phantom edge upper-bound: " << gamma.phantom_edge_upper_bound() << "\n";
    std::cerr << "Expecting at most " << ((gamma.E().row_size(0) + gamma.phantom_edge_upper_bound()) / 2) << " more non-DCC maximal unitigs\n";

    open_p_v();

    {
        Discontinuity_Graph_Contractor<k, Colored_> contractor(gamma, P_v, logistics);
        EXECUTE("contract", contractor.contract)

        gamma.close();
    }

    const auto t_c = timer::now();
    std::cerr << "Discontinuity-graph contraction completed. Time taken: " << timer::duration(t_c - t_subg) << " seconds.\n";

    open_p_e();

    {
        Contracted_Graph_Expander<k, Colored_> expander(gamma, P_v, P_e, logistics);
        EXECUTE("expand", expander.expand)
    }

    force_free(P_v);

    const auto t_e = timer::now();
    std::cerr << "Expansion of contracted graph completed. Time taken: " << timer::duration(t_e - t_c) << " seconds.\n";

    {
        Unitig_Collator<k, Colored_> collator(gamma, P_e, logistics, op_buf, params.gmtig_bucket_count());
        EXECUTE("collate", collator.collate)
    }

    // Flush data and close the output sink.
    parlay::parallel_for(0, parlay::num_workers(),
                        [&](const std::size_t idx){ op_buf[idx].unwrap().close(); }, 1);
    output_sink.close_sink();

    const auto t_uc = timer::now();
    std::cerr << "Unitigs-collation completed. Time taken: " << timer::duration(t_uc - t_e) << " seconds.\n";

    force_free(P_e);
}


template <uint16_t k>
void dBG_Contractor<k>::construct()
{
    params.color() ? construct<true>() : construct<false>();
}


template <uint16_t k>
void dBG_Contractor<k>::open_p_v()
{
    const auto p_v_path_pref = logistics.vertex_path_info_buckets_path();
    P_v.reserve(params.vertex_part_count() + 1);
    P_v.emplace_back(); // No vertex other than the ϕ-vertex belongs to partition 0.
    for(std::size_t j = 1; j <= params.vertex_part_count(); ++j)
        P_v.emplace_back(p_v_bucket_t(p_v_path_pref + "_" + std::to_string(j), 128 * 1024));    // 128 KB for `P_v` buffers.
}


template <uint16_t k>
void dBG_Contractor<k>::open_p_e()
{
    const auto p_e_path_pref = logistics.edge_path_info_buckets_path();
    P_e.reserve(params.lmtig_bucket_count() + 1);
    P_e.emplace_back(); // Using edge-partition 0 with edges that do not have any associated lm-tig (i.e. has weight > 1).
    for(std::size_t b = 1; b <= params.lmtig_bucket_count(); ++b)
        P_e.emplace_back(p_e_bucket_t(p_e_path_pref + "_" + std::to_string(b), 128 * 1024));    // 128 KB for `P_e` buffers.
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, cuttlefish::dBG_Contractor)

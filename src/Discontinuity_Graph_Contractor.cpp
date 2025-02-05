
#include "Discontinuity_Graph_Contractor.hpp"
#include "Data_Logistics.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <filesystem>
#include <cassert>


namespace cuttlefish
{

template <uint16_t k, bool Colored_>
Discontinuity_Graph_Contractor<k, Colored_>::Discontinuity_Graph_Contractor(Discontinuity_Graph<k, Colored_>& G, P_v_t& P_v, const Data_Logistics& logistics):
      G(G)
    , P_v(P_v)
    , compressed_diagonal_path(logistics.compressed_diagonal_path())
    , M(G.vertex_part_size_upper_bound())
    , D_c(parlay::num_workers())
    , phantom_count_(0)
    , icc_count(0)
{
    std::cerr << "Hash table capacity during contraction: " << M.capacity() << ".\n";
}


template <uint16_t k, bool Colored_>
void Discontinuity_Graph_Contractor<k, Colored_>::contract()
{
    // Debug
    double map_clr_time = 0;    // Time taken to clear the hash map.
    double edge_proc_time = 0;  // Time taken to process the non-diagonal edges.
    double diag_comp_time = 0;  // Time taken to compute the diagonal chains.
    double diag_cont_time = 0;  // Time taken to contract the diagonal chains.
    double phantom_filt_time = 0;   // Time taken to filter in the false-phantom edges.

    // TODO: document the phases.

    std::vector<Padded<Buffer<Discontinuity_Edge<k>>>> B(parlay::num_workers());    // Worker-local edge-read buffers.
    constexpr std::size_t buf_cap = 1 * 1024 * 1024 / sizeof(Discontinuity_Edge<k>);    // 1 MB read-capacity.
    parlay::parallel_for(0, B.size(), [&](const auto w){ B[w].unwrap().resize_uninit(buf_cap); });

    for(auto j = G.E().vertex_part_count(); j >= 1; --j)
    {
        std::cerr << "\rPart: " << j;

        auto t_s = now();
        M.clear();
        auto t_e = now();
        map_clr_time += duration(t_e - t_s);

        parlay::par_do(
            [&]()
            {
                const auto t_s = now();
                contract_diagonal_block(j, B[parlay::worker_id()].unwrap());
                const auto t_e = now();
                diag_comp_time += duration(t_e - t_s);
            }
            ,
            [&]()
            {
                const auto process_non_diagonal_edges = [&](auto)
                {
                    auto& buf = B[parlay::worker_id()].unwrap();
                    std::size_t read;
                    for(std::size_t i = 0; i < j; ++i)
                        while((read = G.E().read_block_buffered(i, j, buf, buf_cap)) > 0)
                            for(std::size_t idx = 0; idx < read; ++idx)
                            {
                                const auto& e = buf[idx];
                                Other_End* p_z;

                                assert(G.E().partition(e.y()) == j);
                                assert(e.x_is_phi() || G.E().partition(e.x()) < j);

                                if(M.insert(e.y(), Other_End(e.x(), e.s_x(), e.s_y(), e.x_is_phi(), e.w(), false), p_z))
                                {
                                    assert(M.find(e.y()));
                                    continue;
                                }

                                assert(M.find(e.y()));
                                auto& z = *p_z;
                                if(z.is_phi() && e.x_is_phi())
                                    form_meta_vertex(e.y(), j, e.s_y(), e.w(), z.w(), false);
                                else if(z.in_same_part())   // Corresponds to a compressed diagonal chain.
                                {
                                    assert(!z.is_phi());
                                    assert(M.find(z.v()));
                                    assert(G.E().partition(z.v()) == j);

                                    z = Other_End(e.x(), e.s_x(), e.s_y(), e.x_is_phi(), e.w(), false);
                                }
                                else
                                    G.add_edge(e.x(), e.s_x(), z.v(), z.s_v(), e.w() + z.w(), e.x_is_phi(), z.is_phi());

                                z.process();
                            }
                };

                const auto t_s = now();
                parlay::parallel_for(0, parlay::num_workers(), process_non_diagonal_edges, 1);
                const auto t_e = now();
                edge_proc_time += duration(t_e - t_s);
            }
        );


        t_s = now();
        std::vector<std::size_t> C(parlay::num_workers(), 0);
        for(std::size_t w = 1; w < D_c.size(); ++w)
            C[w] = C[w - 1] + D_c[w - 1].unwrap().size();

        const auto d_c_sz = C[D_c.size() - 1] + D_c.back().unwrap().size();
        D_c_flat.reserve_uninit(d_c_sz);
        parlay::parallel_for(0, D_c.size(),
            [&](const std::size_t w)
            {
                const auto& d = D_c[w].unwrap();
                if(!d.empty())
                    std::memcpy(reinterpret_cast<char*>(D_c_flat.data() + C[w]), reinterpret_cast<const char*>(d.data()), d.size() * sizeof(Discontinuity_Edge<k>));
            }
        , 1);

        const auto contract_diagonal_chain = [&](const std::size_t i)
        {
            const auto& e = D_c_flat[i];
            const auto w_xy = e.w();
            assert(!e.x_is_phi() && !e.y_is_phi());

            if(!M.find(e.x()))  // `e.x()` has a false-phantom edge.
            {
                assert(D.contains(e.x()) && D[e.x()].v() == e.y());
                phantom_count_++;
                G.add_edge(e.x(), inv_side(e.s_x()));
                M.insert(e.x(), Other_End(Discontinuity_Graph<k, Colored_>::phi(), side_t::back, inv_side(e.s_x()), true, 1, false));
            }

            if(!M.find(e.y()))  // `e.y()` has a false-phantom edge.
            {
                assert(D.contains(e.y()) && D[e.y()].v() == e.x());
                phantom_count_++;
                G.add_edge(e.y(), inv_side(e.s_y()));
                M.insert(e.y(), Other_End(Discontinuity_Graph<k, Colored_>::phi(), side_t::back, inv_side(e.s_y()), true, 1, false));
            }

            auto& m_x = *M.find(e.x());
            auto& m_y = *M.find(e.y());

            if(m_x.is_phi() && m_y.is_phi())
                form_meta_vertex(e.x(), j, inv_side(e.s_x()), m_x.w(), w_xy + m_y.w(), false);
            else
                G.add_edge(m_x.v(), m_x.s_v(), m_y.v(), m_y.s_v(), m_x.w() + w_xy + m_y.w(), m_x.is_phi(), m_y.is_phi());

            m_x.process(), m_y.process();
        };

        parlay::parallel_for(0, d_c_sz, contract_diagonal_chain);

        t_e = now();
        diag_cont_time += duration(t_e - t_s);


        t_s = now();
        const auto add_false_phantom_edges =
            [&](const std::size_t w_id)
            {
                auto it = M.iterator(parlay::num_workers(), w_id);
                Kmer<k> u;  // Vertex in the hash table.
                Other_End end;  // Other endpoint of `u` in the table.
                while(it.next(u, end))
                    if(!end.processed())    // `u` has a false-phantom edge.
                    {
                        phantom_count_++;
                        G.add_edge(u, inv_side(end.s_u()));

                        if(end.is_phi())
                            form_meta_vertex(u, j, end.s_u(), end.w(), 1, false);
                        else
                        {
                            assert(G.E().partition(end.v()) < j);
                            G.add_edge(Discontinuity_Graph<k, Colored_>::phi(), side_t::back, end.v(), end.s_v(), 1 + end.w(), true, false);
                        }
                    }
            };

        parlay::parallel_for(0, parlay::num_workers(), add_false_phantom_edges, 1);
        t_e = now();
        phantom_filt_time += duration(t_e - t_s);
    }

    std::cerr << "\n";


    uint64_t meta_v_c = 0;
    std::for_each(P_v.cbegin(), P_v.cend(), [&](const auto& b){ meta_v_c += b.unwrap().size(); });
    std::cerr << "Formed " << meta_v_c << " meta-vertices.\n";
    std::cerr << "Found " << icc_count << " ICCs.\n";
    std::cerr << "Found " << phantom_count_ << " phantoms.\n";
    std::cerr << "Map clearing time: " << map_clr_time << ".\n";
    std::cerr << "Non-diagonal edges contraction time: " << edge_proc_time << ".\n";
    std::cerr << "Diagonal-chain computation time: " << diag_comp_time << ".\n";
    std::cerr << "Diagonal-chain contraction time: " << diag_cont_time << ".\n";
    std::cerr << "Filtering in false-phantom edges time: " << phantom_filt_time << ".\n";
}


template <uint16_t k, bool Colored_>
void Discontinuity_Graph_Contractor<k, Colored_>::contract_diagonal_block(const std::size_t j, Buffer<Discontinuity_Edge<k>>& buf)
{
    D_j.clear();
    D.clear();

    constexpr std::size_t buf_cap = 1 * 1024 * 1024 / sizeof(Discontinuity_Edge<k>);
    buf.reserve_uninit(buf_cap);
    std::size_t read;
    while((read = G.E().read_block_buffered(j, j, buf, buf_cap)) > 0)
        for(std::size_t i = 0; i < read; ++i)
        {
            const auto& e = buf[i];
            const auto end_x = D.find(e.x());
            const auto end_y = D.find(e.y());

            const auto& u  = (end_x != D.end() ? end_x->second.v() : e.x());
            const auto s_u = (end_x != D.end() ? end_x->second.s_v() : e.s_x());
            const auto w_u = (end_x != D.end() ? end_x->second.w() : 0);
            const auto& v  = (end_y != D.end() ? end_y->second.v() : e.y());
            const auto s_v = (end_y != D.end() ? end_y->second.s_v() : e.s_y());
            const auto w_v = (end_y != D.end() ? end_y->second.w() : 0);
            const auto w   = w_u + e.w() + w_v;

            assert(G.E().partition(u) == j);
            assert(G.E().partition(v) == j);

            // When adding an `{u, v}` edge, if an `{u, v}` edge already exists,
            // then it's an Isolated Cordless Cycle (ICC). In this case, if `u ≠ v`,
            // then the current contracted form of the ICC resides entirely in this
            // diagonal block. Otherwise, the ICC has already been contracted down
            // to a single vertex that has a self-loop right now.
            // In both cases, form a meta-vertex with `u` having rank `1`.
            // In both cases, information-propagation through the side of `u`
            // "facing outward" of the linearized cycle needs to be discarded, but
            // the single lm-tig there has to be included in the maximal unitig
            // label. This forms a tricky and cumbersome special case.

            if(CF_UNLIKELY(u == v)) // The ICC has already been contracted to a single vertex.
            {
                assert(e.x() == e.y() && e.x() == u); assert(e.w() > 1);
                D.insert_or_assign(u, Other_End(Discontinuity_Graph<k, Colored_>::phi(), side_t::back, side_t::front, true, 1, false, true));
                form_meta_vertex(u, j, side_t::front, 1, w, true);
                icc_count++;
            }
            else if(CF_UNLIKELY(u == e.y() && v == e.x()))  // The currently contracted form of the actual ICC resides in this diagonal block.
            {
                assert(e.x() != e.y());
                assert(inv_side(s_u) == e.s_y()); assert(inv_side(s_v) == e.s_x());

                D.insert_or_assign(u, Other_End(Discontinuity_Graph<k, Colored_>::phi(), side_t::back, inv_side(s_u), true, 1, false, true));
                D.insert_or_assign(v, Other_End(Discontinuity_Graph<k, Colored_>::phi(), side_t::back, inv_side(s_v), true, 1, false, true));
                form_meta_vertex(u, j, inv_side(s_u), 1, true); // Propagation needs to go through `s_u`, since `u` has already been connected to the other cycle members through that side.
                icc_count++;
            }
            else
            {
                D.insert_or_assign(u, Other_End(v, s_v, s_u, false, w, true, true));
                D.insert_or_assign(v, Other_End(u, s_u, s_v, false, w, true, true));
                D_j.emplace_back(u, s_u, v, s_v, w, w == 1 ? e.b() : 0, w == 1 ? e.b_idx() : 0, false, false, side_t::unspecified);
            }
    }


    std::filesystem::create_directories(compressed_diagonal_path);
    const auto d_j_path = compressed_diagonal_path + "/" + std::to_string(j);
    std::ofstream output(d_j_path);
    output.write(reinterpret_cast<const char*>(D_j.data()), D_j.size() * sizeof(Discontinuity_Edge<k>));
    if(!output)
    {
        std::cerr << "Error writing compressed diagonal edge block at " << d_j_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    output.close();


    std::for_each(D_c.begin(), D_c.end(), [](auto& v){ v.unwrap().clear(); });
    const auto& kv = D.values();
    const auto collect_compressed_diagonal_chains =
        [&](const std::size_t i)
        {
            const auto& u = kv[i].first;
            const auto& end = kv[i].second;
            if(!end.is_phi())   // Not an ICC.
            {
                assert(D.find(end.v()) != D.end());
                if(u < end.v() && D[end.v()].v() == u)
                    D_c[parlay::worker_id()].unwrap().emplace_back(u, end.s_u(), end.v(), end.s_v(), end.w(), 0, 0, false, false, side_t::unspecified);
            }
        };

    parlay::parallel_for(0, D.size(), collect_compressed_diagonal_chains);
}


}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BOOL, cuttlefish::Discontinuity_Graph_Contractor)

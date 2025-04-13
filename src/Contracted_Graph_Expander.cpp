
#include "Contracted_Graph_Expander.hpp"
#include "Data_Logistics.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"

#include <fstream>
#include <cstdlib>
#include <cassert>


namespace cuttlefish
{

template <uint16_t k, bool Colored_>
Contracted_Graph_Expander<k, Colored_>::Contracted_Graph_Expander(Discontinuity_Graph<k, Colored_>& G, P_v_t& P_v, P_e_t& P_e, const Data_Logistics& logistics):
      G(G)
    , P_v(P_v)
    , P_e(P_e)
    , compressed_diagonal_path(logistics.compressed_diagonal_path())
    , M(G.vertex_part_size_upper_bound())
{
    std::cerr << "Hash table capacity during expansion: " << M.capacity() << ".\n";
}


template <uint16_t k, bool Colored_>
void Contracted_Graph_Expander<k, Colored_>::expand()
{
    // Debug
    double map_clr_time = 0;    // Time taken to clear the hash map.
    double p_v_load_time = 0;   // Time to load vertices' path-info.
    double diag_exp_time = 0;   // Time taken to expand the compressed diagonal chains.
    double edge_proc_time = 0;  // Time taken to process the non-diagonal edges.
    double spec_case_time = 0;  // Time taken to process the special edge-blocks.

    const auto v_part_c = G.E().vertex_part_count();    // Number of vertex-partitions.

    std::vector<Padded<Buffer<Obj_Path_Info_Pair<Kmer<k>, k>>>> B_p_v(parlay::num_workers());
    std::vector<Padded<Buffer<Discontinuity_Edge<k>>>> B_e(parlay::num_workers());
    constexpr std::size_t buf_cap_p_v = 1 * 1024 * 1024 / sizeof(Obj_Path_Info_Pair<Kmer<k>, k>);   // 1 MB read-capacity.
    constexpr std::size_t buf_cap_e = 1 * 1024 * 1024 / sizeof(Discontinuity_Edge<k>);  // 1 MB read-capacity.
    parlay::parallel_for(0, parlay::num_workers(),
        [&](const auto i) { B_p_v[i].unwrap().resize_uninit(buf_cap_p_v), B_e[i].unwrap().resize_uninit(buf_cap_e); });

    G.E().reset_read();

    for(std::size_t i = 1; i <= v_part_c; ++i)
    {
        std::cerr << "\rPart: " << i;

        auto t_s = now();
        M.clear();
        auto t_e = now();
        map_clr_time += duration(t_e - t_s);


        const auto load_path_info = [&](auto)
        {
            const auto& P_v_i = P_v[i].unwrap();
            auto& buf = B_p_v[parlay::worker_id()].unwrap();
            std::size_t read;
            while((read = P_v_i.read_buffered(buf, buf_cap_p_v)) > 0)
                for(std::size_t idx = 0; idx < read; ++idx)
                {
                    const auto& p_v = buf[idx];
                    if(!M.insert(p_v.obj(), p_v.path_info()))
                        assert(*M.find(p_v.obj()) == p_v.path_info());
                }
        };

        t_s = now();
        parlay::parallel_for(0, parlay::num_workers(), load_path_info, 1);
        t_e = now();
        p_v_load_time += duration(t_e - t_s);


        t_s = now();
        expand_diagonal_block(i);
        t_e = now();
        diag_exp_time += duration(t_e - t_s);


        const auto process_non_diagonal_edge = [&](auto)
        {
            auto& buf = B_e[parlay::worker_id()].unwrap();
            std::size_t read;
            for(std::size_t j = i + 1; j <= v_part_c; ++j)
                while((read = G.E().read_block_buffered(i, j, buf, buf_cap_e)) > 0)
                    for(std::size_t idx = 0; idx < read; ++idx)
                    {
                        const auto& e = buf[idx];

                        assert(M.find(e.x()));
                        const auto x_inf = *M.find(e.x());  // Path-information of the endpoint of the edge that is in the current partition, `i`.
                        const auto y_inf = infer(x_inf, e.s_x(), e.s_y(), e.w());
                        if(y_inf.r() > 0)   // `e` is not the back-edge to the rank-1 vertex in an ICC.
                        {
                            const auto j = G.E().partition(e.y());  // TODO: consider obtaining this info during the edge-reading process.
                            assert(j > i);
                            P_v[j].unwrap().emplace(e.y(), y_inf.p(), y_inf.r(), y_inf.o(), y_inf.is_cycle());
                        }

                        if(e.w() == 1)
                            add_edge_path_info(e, x_inf, y_inf);
                    }
        };

        t_s = now();
        parlay::parallel_for(0, parlay::num_workers(), process_non_diagonal_edge, 1);
        t_e = now();
        edge_proc_time += duration(t_e - t_s);


        t_s = now();
        parlay::parallel_for(0, parlay::num_workers(),
        [&](auto)
        {
            auto& buf = B_e[parlay::worker_id()].unwrap();
            std::size_t read;
            while((read = G.E().read_block_buffered(i, i, buf, buf_cap_e)) > 0)
                for(std::size_t idx = 0; idx < read; ++idx)
                {
                    const auto& e = buf[idx];
                    if(e.w() == 1)
                    {
                        assert(M.find(e.x()) && M.find(e.y()));
                        add_diagonal_edge_path_info(e, *M.find(e.x()));
                    }
                }
        }, 1);

        parlay::parallel_for(0, parlay::num_workers(),
        [&](auto)
        {
            auto& buf = B_e[parlay::worker_id()].unwrap();
            std::size_t read;
            while((read = G.E().read_block_buffered(0, i, buf, buf_cap_e)) > 0)
                for(std::size_t idx = 0; idx < read; ++idx)
                {
                    const auto& e = buf[idx];
                    if(e.w() == 1)
                    {
                        assert(M.find(e.y()));
                        add_edge_path_info(e, *M.find(e.y()));
                    }
                }
        }, 1);
        t_e = now();

        spec_case_time += duration(t_e - t_s);
    }

    std::cerr << "\n";

    const auto t_s = now();

    // Remove the vertices' path-information buckets.
    parlay::parallel_for(1, P_v.size(), [&](const auto i){ P_v[i].unwrap().remove(); });

    // Remove the edge-matrix.
    parlay::parallel_for(1, v_part_c + 1,
    [&](const auto j)
    {
        parlay::parallel_for(0, j + 1, [&](const auto i){ G.E().remove_block(i, j); }, 1);
    }, 1);

    const auto t_e = now();
    const auto rm_time = duration(t_e - t_s);


    // std::cerr << "Encountered " << og_edge_c << " original edges\n";
    std::cerr << "Map clearing time: " << map_clr_time << ".\n";
    std::cerr << "Path-info load time: " << p_v_load_time << ".\n";
    std::cerr << "Non-diagonal blocks expansion time: " << edge_proc_time << ".\n";
    std::cerr << "Diagonal block expansion time: " << diag_exp_time << ".\n";
    std::cerr << "Special case time: " << spec_case_time << ".\n";
    std::cerr << "Deletion time: " << rm_time << ".\n";
}


template <uint16_t k, bool Colored_>
void Contracted_Graph_Expander<k, Colored_>::expand_diagonal_block(const std::size_t i)
{
    const std::string d_i_path(compressed_diagonal_path + "/" + std::to_string(i));
    const auto file_sz = file_size(d_i_path);
    const auto edge_c = file_sz / sizeof(Discontinuity_Edge<k>);
    D_i.reserve_uninit(edge_c);

    std::ifstream input(d_i_path);
    input.read(reinterpret_cast<char*>(D_i.data()), file_sz);
    assert(static_cast<std::size_t>(input.gcount()) == file_sz);
    input.close();


    parlay::par_do(
        [&]()
        {
            if(!input || !remove_file(d_i_path))
            {
                std::cerr << "Error reading / removing of contracted edge block from " << d_i_path << ". Aborting.\n";
                std::exit(EXIT_FAILURE);
            }
        },
        [&]()
        {
            Path_Info<k> x_inf, y_inf;
            for(int64_t idx = static_cast<int64_t>(edge_c) - 1; idx >= 0; --idx)    // In reverse order of the newly introduced diagonal-edges to always ensure one endpoint having path-info ready.
            {
                const auto& e = D_i[idx];
                assert(M.find(e.x()) || M.find(e.y()));

                if(!M.find(e.y(), y_inf))
                {
                    assert(M.find(e.x()));
                    x_inf = *M.find(e.x()); // M.find(e.x(), x_inf);
                    y_inf = infer(x_inf, e.s_x(), e.s_y(), e.w());
                    if(y_inf.r() > 0)
                        M.insert(e.y(), y_inf);
                }
                else
                {
                    assert(M.find(e.y()));
                    x_inf = infer(y_inf, e.s_y(), e.s_x(), e.w());
                    if(x_inf.r() > 0)
                        M.insert(e.x(), x_inf);
                }
            }
        }
    );
}


/*
template <uint16_t k>
template <typename T_s_, typename T_d_>
uint64_t Contracted_Graph_Expander<k>::collate_w_local_bufs(T_s_& source, const std::size_t beg, const size_t end, T_d_& dest)
{
    std::vector<Padded_Data<uint64_t>> C(parlay::num_workers(), 0);

    parlay::parallel_for(beg, end,
        [&](const std::size_t idx)
        {
            for(auto& w_local_buf_list : source)
            {
                auto& buf = w_local_buf_list.data()[idx];
                dest[idx].add(buf.data(), buf.size());
                C[parlay::worker_id()].data() += buf.size();

                buf.clear();
            }
        }, 1);

    uint64_t c = 0;
    std::for_each(C.cbegin(), C.cend(), [&](auto& v){ c += v.data(); });
    return c;
}
*/

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BOOL, cuttlefish::Contracted_Graph_Expander)

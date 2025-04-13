
#include "Unitig_Collator.hpp"
#include "Unitig_File.hpp"
#include "FASTA_Record.hpp"
#include "Data_Logistics.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"
#include "xxHash/xxh3.h"

#include <vector>
#include <string>
#include <string_view>
#include <cstdlib>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <cassert>


// TODO: use more efficient data structures throughout the collator.


namespace cuttlefish
{

template <uint16_t k, bool Colored_>
Unitig_Collator<k, Colored_>::Unitig_Collator(Discontinuity_Graph<k, Colored_>& G, P_e_t& P_e, const Data_Logistics& logistics, op_buf_list_t& op_buf, const std::size_t gmtig_bucket_count):
      G(G)
    , P_e(P_e)
    , lmtig_buckets_path(logistics.lmtig_buckets_path())
    , unitig_coord_buckets_path(logistics.unitig_coord_buckets_path())
    , max_unitig_bucket_count(gmtig_bucket_count)
    , op_buf(op_buf)
    , phantom_c_(0)
{
     // TODO: fix better policy?
    if((max_unitig_bucket_count & (max_unitig_bucket_count - 1)) != 0)
    {
        std::cerr << "Maximal unitig bucket count needs to be a power of 2. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    max_bucket_sz = 0;
    std::size_t sum_bucket_sz = 0;
    std::for_each(P_e.cbegin(), P_e.cend(),
        [&](const auto& bucket){ max_bucket_sz = std::max(max_bucket_sz, bucket.unwrap().size()); sum_bucket_sz += bucket.unwrap().size(); });

    std::cerr << "Sum edge-bucket size: " << sum_bucket_sz << "\n";
    std::cerr << "Maximum edge-bucket size: " << max_bucket_sz << "\n";

    std::filesystem::create_directories(unitig_coord_buckets_path);
}


template <uint16_t k, bool Colored_>
void Unitig_Collator<k, Colored_>::collate()
{
    const auto t_0 = timer::now();
    map();

    const auto t_1 = timer::now();
    std::cerr << "Time taken in mapping: " << timer::duration(t_1 - t_0) << "s.\n";

    reduce();
    const auto t_2 = timer::now();
    std::cerr << "Time taken in reduction: " << timer::duration(t_2 - t_1) << "s.\n";

    if constexpr(Colored_)
    {
        emit_trivial_mtigs();
        const auto t_3 = timer::now();
        std::cerr << "Time taken in trivial-mtigs emission: " << timer::duration(t_3 - t_2) << "s.\n";

        std::cerr << "Found " << phantom_c_ << " phantom unitigs.\n";
    }

    // TODO: print meta-information over the maximal unitigs'.
}


template <uint16_t k, bool Colored_>
void Unitig_Collator<k, Colored_>::map()
{
    typedef Buffer<Path_Info<k>> map_t;
    typedef Buffer<unitig_path_info_t> buf_t;
    typedef Buffer<Vertex_Color_Mapping> v_c_map_t;

    std::vector<Padded<map_t>> M_vec(parlay::num_workers());   // Worker-local DATs of edge-index to edge path-info.
    std::vector<Padded<buf_t>> buf_vec(parlay::num_workers()); // Worker-local buffers to read in edge path-info.
    std::vector<Padded<v_c_map_t>> v_c_map_vec(parlay::num_workers());  // Worker-local buffers to read in vertex-color mappings.

    parlay::parallel_for(0, parlay::num_workers(),
        [&](const std::size_t w_id)
        {
            // TODO: no need to resize to the maximum bucket-size, as we use more suited buffer now instead of `vector`.
            M_vec[w_id].unwrap().resize_uninit(max_bucket_sz);   // TODO: thread-local allocation suits best here.
        }, 1);


    max_unitig_bucket.reserve(max_unitig_bucket_count);
    for(std::size_t i = 0; i < max_unitig_bucket_count; ++i)
        max_unitig_bucket.emplace_back(unitig_coord_buckets_path + "/" + std::to_string(i));

    std::atomic_uint64_t edge_c = 0;    // Number of edges (i.e. unitigs) found.
#ifndef NDEBUG
    std::atomic_uint64_t h_p_e = 0; // Hash of the edges' path-information.
#endif

    // Maps the edges (i.e. unitigs) from bucket `b` to maximal-unitig buckets.
    const auto map_to_max_unitig_bucket =
    [&](const std::size_t b)
    {
        const auto w_id = parlay::worker_id();
        auto& M = M_vec[w_id].unwrap();
        auto& buf = buf_vec[w_id].unwrap();
        auto& v_c_map = v_c_map_vec[w_id].unwrap();

        const auto b_sz = load_path_info(b, M.data(), buf);
        edge_c += b_sz;
        P_e[b].unwrap().remove();

#ifndef NDEBUG
        uint64_t h = 0;
        for(std::size_t i = 0; i < b_sz; ++i)
            h ^= M[i].hash();

        h_p_e ^= h;
#endif

        std::size_t v_c_map_sz = 0;
        if constexpr(Colored_)
        {
            v_c_map_sz = load_vertex_color_mapping(b, v_c_map);
            std::sort(v_c_map.data(), v_c_map.data() + v_c_map_sz); // TODO: replace.
        }

        const auto bucket_path = lmtig_buckets_path + "/" + std::to_string(b);
        assert(file_exists(bucket_path));
        Unitig_File_Reader unitig_reader(bucket_path);
        Buffer<char> unitig;    // Read-off unitig.
        uni_idx_t idx = 0;  // The unitig's sequential ID in the bucket.
        std::size_t uni_len;    // The unitig's length in bases.
        std::size_t color_idx = 0;  // The unitig's associated color-mappings' index into the sorted mappings.
        std::vector<Unitig_Color> color;    // Color-encodings of the unitig.
        for(; (uni_len = unitig_reader.read_next_unitig(unitig)) > 0; idx++)
        {
            assert(idx < b_sz);
            const auto p = M[idx].p();  // Path-ID of this unitig.
            const auto mapped_b_id = XXH3_64bits(&p, sizeof(p)) & (max_unitig_bucket_count - 1); // Hashed maximal unitig bucket.
            // TODO: use `wyhash`. Result: buckets mapped with wyhash (seed = 0) is unbalanced af. Unreliable.
            // TODO: interesting af issue: if we use the same hash function here and in `kmer.to_u64()`, it results in disastrous balancing in this case.

            if constexpr(!Colored_)
                max_unitig_bucket[mapped_b_id].unwrap().add(M[idx], unitig.data(), uni_len);
            else
            {
                color.clear();
                if(CF_UNLIKELY(color_idx == v_c_map_sz))    // Unitigs due to phantom k-mers do not have colors in this bucket.
                {
                    assert(uni_len == k);
                    phantom_c_++;
                }
                else
                {
                    for(; color_idx < v_c_map_sz; ++color_idx)
                    {
                        const auto& v_c = v_c_map[color_idx];
                        if(v_c.idx() != idx)
                            break;
                        else
                        {
                            if(!color.empty())
                                assert(v_c.off() > v_c_map[color_idx - 1].off());

                            color.emplace_back(v_c.off(), v_c.c());
                        }
                    }

                    assert(!color.empty());
                    assert(color.front().off() == 0 && color.back().off() < uni_len - k + 1);
                }

                max_unitig_bucket[mapped_b_id].unwrap().add(M[idx], unitig.data(), uni_len, color);
            }
        }

        assert(idx == b_sz);
        assert(color_idx == v_c_map_sz);

        unitig_reader.remove_files();
        if constexpr(Colored_)
            G.vertex_color_map(b).remove();
    };

    parlay::parallel_for(1, P_e.size(), map_to_max_unitig_bucket, 1);

    std::cerr << "Found " << edge_c << " edges.\n";
#ifndef NDEBUG
    std::cerr << "Edges' path-information signature: " << h_p_e << "\n";
#endif
}


template <uint16_t k, bool Colored_>
void Unitig_Collator<k, Colored_>::reduce()
{
    std::size_t max_max_uni_b_sz = 0;   // Maximum unitig-count in some maximal unitig bucket.
    std::size_t max_max_uni_b_label_len = 0;    // Maximum dump-string length in some maximal unitig bucket.
    std::size_t max_max_uni_b_color_c = 0;  // Maximum color-count in some maximal unitig bucket.
    std::for_each(max_unitig_bucket.cbegin(), max_unitig_bucket.cend(),
        [&](const auto& b)
        {
            max_max_uni_b_sz = std::max(max_max_uni_b_sz, b.unwrap().size()),
            max_max_uni_b_label_len = std::max(max_max_uni_b_label_len, b.unwrap().label_len());
            if constexpr(Colored_)
                max_max_uni_b_color_c = std::max(max_max_uni_b_color_c, b.unwrap().color_count());
        });

    std::cerr << "Maximum maximal unitig bucket size:   " << max_max_uni_b_sz << "\n";
    std::cerr << "Maximum label length in mtig-buckets: " << max_max_uni_b_label_len << "\n";
    if constexpr(Colored_)
        std::cerr << "Maximum colors in mtig-buckets:  " << max_max_uni_b_color_c << "\n",
        assert(max_max_uni_b_color_c >= max_max_uni_b_sz);


    typedef Buffer<Unitig_Coord<k, Colored_>> coord_buf_t;
    typedef Buffer<char> label_buf_t;
    typedef Buffer<Unitig_Color> color_buf_t;
    std::vector<Padded<coord_buf_t>> U_vec(parlay::num_workers()); // Worker-local buffers for unitig coordinate information.
    std::vector<Padded<label_buf_t>> L_vec(parlay::num_workers()); // Worker-local buffers for dump-strings in buckets.
    std::vector<Padded<color_buf_t>> C_vec(parlay::num_workers()); // Worker-local buffers for colors in buckets.

    parlay::parallel_for(0, parlay::num_workers(),
        [&](const std::size_t w_id)
        {
            // TODO: thread-local allocations here suit best.
            U_vec[w_id].unwrap().resize_uninit(max_max_uni_b_sz);
            L_vec[w_id].unwrap().resize_uninit(max_max_uni_b_label_len);
            if constexpr(Colored_)
                C_vec[w_id].unwrap().resize_uninit(max_max_uni_b_color_c);
        }, 1);

    // TODO: add per-worker progress tracker.

    const auto collate_max_unitig_bucket =
    [&](const std::size_t b)
    {
        const auto w_id = parlay::worker_id();
        auto const U = U_vec[w_id].unwrap().data(); // Coordinate info of the unitigs.
        auto const L = L_vec[w_id].unwrap().data();   // Dump-strings of the unitig labels.
        auto const C = C_vec[w_id].unwrap().data(); // Colors of the unitigs.
        auto& output = op_buf[w_id].unwrap();   // Output buffer for the maximal unitigs.

        const auto b_sz = max_unitig_bucket[b].unwrap().load_coords(U);
        const auto len = max_unitig_bucket[b].unwrap().load_labels(L);
        const auto color_c = (Colored_ ? max_unitig_bucket[b].unwrap().load_colors(C) : 0);
        max_unitig_bucket[b].unwrap().remove();
        (void)len; (void)color_c;

        std::sort(U, U + b_sz);

        Maximal_Unitig m_tig(*this);
        std::size_t i, j;
        for(i = 0; i < b_sz; i = j)
        {
            m_tig.clear();

            const bool is_cycle = U[i].is_cycle();
            const std::size_t s = i;
            std::size_t e = s;

            // Find the current maximal unitig's stretch in the bucket.
            while(e < b_sz && U[e].p() == U[s].p())
            {
                assert(U[e].o() != side_t::unspecified); assert(U[e].is_cycle() == is_cycle);

                assert(U[e].label_idx() + U[e].label_len() <= len);
                if constexpr(Colored_)
                    assert(U[e].color_idx() + U[e].color_c() <= color_c);

                e++;
            }

            if(e - s == 2 && !is_cycle) // Special case for handling orientation of path-traversals in the discontinuity graph.
            {
                assert(U[s].r() == 0); assert(U[s + 1].r() == 0);

                const std::string_view u0(L + U[s].label_idx(), U[s].label_len());
                const std::string_view u1(L + U[s + 1].label_idx(), U[s + 1].label_len());

                const bool rc_0 = (U[s].o() == side_t::front);
                bool rc_1 = (U[s + 1].o() == side_t::front);
                rc_1 = !rc_1;

                if constexpr(!Colored_)
                {
                    m_tig.init(u0, rc_0);
                    m_tig.append(u1, rc_1);
                }
                else
                {
                    m_tig.init(u0, rc_0, C + U[s].color_idx(), U[s].color_c());
                    m_tig.append(u1, rc_1, C + U[s + 1].color_idx(), U[s + 1].color_c());
                }
            }
            else
                for(j = s; j < e; ++j)
                {
                    const std::string_view u(L + U[j].label_idx(), U[j].label_len());
                    const bool rc = (U[j].o() == side_t::front);

                    if constexpr(!Colored_)
                        m_tig.empty() ?
                            m_tig.init(u, rc) :
                            m_tig.append(u, rc);
                    else
                        m_tig.empty() ?
                            m_tig.init(u, rc, C + U[j].color_idx(), U[j].color_c()) :
                            m_tig.append(u, rc, C + U[j].color_idx(), U[j].color_c());
                }


            // Cyclic maximal unitig traversals start and end at the same vertex, so one copy needs to be removed.
            if(is_cycle)
                m_tig.pop_back();


            // is_cycle ? m_tig.canonicalize_cycle(): m_tig.canonicalize();

            // TODO: decide record-ID choice.
            if constexpr(!Colored_)
            {
                is_cycle ? m_tig.canonicalize_cycle(): m_tig.canonicalize();
                output += FASTA_Record(0, std::string_view(m_tig.data(), m_tig.size()));
                // TODO: this forms a parallel-scaling bottleneck due to the shared sink of the buffers.
            }
            else
            {
                if(!is_cycle)
                    m_tig.canonicalize();

                output.template operator+=<true>(FASTA_Record(0, std::string_view(m_tig.data(), m_tig.size()), m_tig.color()));
                // TODO: this forms a parallel-scaling bottleneck due to the shared sink of the buffers.
            }

            j = e;
        }
    };

    std::cerr << "Peak-RAM before collation: " << process_peak_memory() / (1024.0 * 1024.0 * 1024.0) << "\n";
    parlay::parallel_for(0, max_unitig_bucket_count, collate_max_unitig_bucket, 1);
    std::cerr << "Peak-RAM after collation:  " << process_peak_memory() / (1024.0 * 1024.0 * 1024.0) << "\n";
}


template <uint16_t k, bool Colored_>
std::size_t Unitig_Collator<k, Colored_>::load_path_info(const std::size_t b, Path_Info<k>* const M, Buffer<unitig_path_info_t>& buf)
{
    buf.reserve_uninit(P_e[b].unwrap().size());  // TODO: perform one fixed resize beforehand, as the `P_e` buckets will not grow anymore.
    const std::size_t b_sz = P_e[b].unwrap().load(buf.data());
    assert(b_sz <= max_bucket_sz);

    for(std::size_t idx = 0; idx < b_sz; ++idx)
    {
        const auto& p_e = buf[idx];
        assert(p_e.obj() < max_bucket_sz);

        M[p_e.obj()] = p_e.path_info();
    }

    return b_sz;
}


template <uint16_t k, bool Colored_>
std::size_t Unitig_Collator<k, Colored_>::load_vertex_color_mapping(const std::size_t b, Buffer<Vertex_Color_Mapping>& buf)
{
    const auto& B = G.vertex_color_map(b);
    buf.reserve_uninit(B.size());
    return B.load(buf.data());
}


template <uint16_t k, bool Colored_>
void Unitig_Collator<k, Colored_>::emit_trivial_mtigs()
{
    assert(Colored_);

    parlay::parallel_for(0, parlay::num_workers(), [&](const std::size_t w)
    {
        Buffer<Vertex_Color_Mapping> v_c_map;
        const auto v_c_map_sz = load_vertex_color_mapping(P_e.size() + w, v_c_map);
        std::sort(v_c_map.data(), v_c_map.data() + v_c_map_sz);

        auto& output = op_buf[w].unwrap();  // Output buffer for the maximal unitigs.
        Unitig_File_Reader unitig_reader(lmtig_buckets_path + "/" + std::to_string(P_e.size() + w));
        Buffer<char> unitig;    // Read-off unitig.
        uni_idx_t idx = 0;  // The unitig's sequential ID in the bucket.
        std::size_t uni_len;    // The unitig's length in bases.
        std::size_t color_idx = 0;  // The unitig's associated color-mappings' index into the sorted mappings.
        std::vector<Unitig_Color> color;    // Color-encodings of the unitig.
        for(; (uni_len = unitig_reader.read_next_unitig(unitig)) > 0; idx++)
        {
            color.clear();
            for(; color_idx < v_c_map_sz; ++color_idx)
            {
                const auto& v_c = v_c_map[color_idx];
                if(v_c.idx() != idx)
                    break;
                else
                {
                    assert(v_c.off() <= uni_len - k);
                    if(color.empty())
                        assert(v_c.off() == 0);
                    else
                        assert(v_c.off() > v_c_map[color_idx - 1].off());

                    color.emplace_back(v_c.off(), v_c.c());
                }
            }

            assert(!color.empty());

            output.template operator+=<true>(FASTA_Record(0, std::string_view(unitig.data(), uni_len), color));
        }

        unitig_reader.remove_files();
        G.vertex_color_map(P_e.size() + w).remove();
    }, 1);
}

}



// Template-instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BOOL, cuttlefish::Unitig_Collator)

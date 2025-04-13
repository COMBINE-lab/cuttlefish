
#include "Subgraph.hpp"
#include "Super_Kmer_Bucket.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"
#include "color_sets/hybrid.hpp"

#include <vector>
#include <algorithm>
#include <filesystem>
#include <cassert>


namespace cuttlefish
{


template <uint16_t k, bool Colored_>
Subgraph<k, Colored_>::Subgraph(const Super_Kmer_Bucket<Colored_>& B, Discontinuity_Graph<k, Colored_>& G, op_buf_t& op_buf, Subgraphs_Scratch_Space<k, Colored_>& space):
      B(B)
    , work_space(space)
    , M(space.map())
    , kmer_count_(0)
    , edge_c(0)
    , label_sz(0)
    , disc_edge_c(0)
    , isolated(0)
    , G(G)
    , mtig_c(0)
    , trivial_mtig_c(0)
    , icc_count_(0)
    , color_shift_c(0)
    , v_new_col_c(0)
    , v_old_col_c(0)
    , color_rel_c(0)
    , op_buf(op_buf)
{
    M.clear();
}


template <uint16_t k, bool Colored_>
void Subgraph<k, Colored_>::construct()
{
    typename Super_Kmer_Bucket<Colored_>::Iterator super_kmer_it(B);    // Iterator over the weak super k-mers inducing this graph.

    typedef typename decltype(super_kmer_it)::label_unit_t label_unit_t;
    const auto word_count = super_kmer_it.super_kmer_word_count();  // Fixed number of words in a super k-mer label.

    Directed_Vertex<k> v;   // Current vertex in a scan over some super k-mer.
    Kmer<k> pred_v; // Previous vertex in a scan.

    Super_Kmer_Attributes<Colored_> att;
    const label_unit_t* label;
    source_id_t source = 0; // Source-ID of the current super k-mer.
    while(super_kmer_it.next(att, label))
    {
        const auto len = att.len();
        assert(len >= k);
        assert(len < 2 * (k - 1));
        kmer_count_ += len - (k - 1);

        if constexpr(Colored_)
        {
            // assert(att.source() >= source);
            source = att.source();
        }

        v.from_super_kmer(label, word_count);
        std::size_t kmer_idx = 0;
        while(true)
        {
            assert(kmer_idx + k - 1 < len);

            const auto is_canonical = v.in_canonical_form();
            const auto pred_base = (kmer_idx == 0 ? base_t::E : get_base(label, word_count, kmer_idx - 1));
            const auto succ_base = (kmer_idx + k == len ? base_t::E : get_base(label, word_count, kmer_idx + k));
            auto front = (is_canonical ? pred_base : DNA_Utility::complement(succ_base));
            auto back  = (is_canonical ? succ_base : DNA_Utility::complement(pred_base));

            if(CF_UNLIKELY(kmer_idx > 0 && v.canonical() == pred_v))    // Counter overcounting of self-loops.
                (is_canonical ? front : back) = base_t::E;

            edge_c += (succ_base != base_t::E);

            // Update hash table with the neighborhood info.
            ht_router::update(M, v.canonical(),
                                 front, back,
                                 kmer_idx == 0 && att.left_discontinuous() ? v.entrance_side() : side_t::unspecified,
                                 kmer_idx + k == len && att.right_discontinuous() ? v.exit_side() : side_t::unspecified,
                                 source);

            if(kmer_idx + k == len)
                break;
/*
            if((kmer_idx > 0 || !att.left_discontinuous()) && (kmer_idx + k < len || !att.right_discontinuous()))
                assert(!st.is_discontinuity());
*/

            pred_v = v.canonical();
            v.roll_forward(succ_base);
            kmer_idx++;
        }
    }

    ht_router::flush_updates(M);
}


/*
template <uint16_t k, bool Colored_>
void Subgraph<k, Colored_>::construct_loop_filtered()
{
    auto super_kmer_it = B.iterator();  // Iterator over the weak super k-mers inducing this graph.

    const auto word_count = super_kmer_it.super_kmer_word_count();  // Fixed number of words in a super k-mer label.

    Directed_Vertex<k> v_pre, v_suf;

    Super_Kmer_Attributes<Colored_> att;
    label_unit_t* label;
    while(super_kmer_it.next(att, label))
    {
        const auto len = att.len();
        assert(len >= k);
        assert(len < 2 * (k - 1));

        if(len == k)
            continue;

        std::size_t edge_idx = 0;
        v_pre.from_super_kmer(label, word_count);
        auto succ_base_pre = get_base(label, word_count, edge_idx + k);
        v_suf = std::as_const(v_pre).roll_forward(succ_base_pre);
        auto pred_base_suf = get_base(label, word_count, edge_idx);

        while(true)
        {
            assert(edge_idx + k < len);

            if(M.find(v_pre.canonical()) == M.end())
                M.emplace(v_pre.canonical(), State_Config());

            if(M.find(v_suf.canonical()) == M.end())
                M.emplace(v_suf.canonical(), State_Config());

            auto& st_pre = M.find(v_pre.canonical())->second;
            auto& st_suf = M.find(v_suf.canonical())->second;

            v_pre.in_canonical_form() ?
                st_pre.update_edge(side_t::back, succ_base_pre) :
                st_pre.update_edge(side_t::front, DNA_Utility::complement(succ_base_pre));
            if(!v_suf.is_same_vertex(v_pre))    // Avoid double-counting of self-loops.
                v_suf.in_canonical_form() ?
                    st_suf.update_edge(side_t::front, pred_base_suf) :
                    st_suf.update_edge(side_t::back, DNA_Utility::complement(pred_base_suf));
            edge_c++;

            if(edge_idx == 0 && att.left_discontinuous())
                st_pre.mark_discontinuous(v_pre.entrance_side());

            if(edge_idx + 1 + k == len)
            {
                if(att.right_discontinuous())
                    st_suf.mark_discontinuous(v_suf.exit_side());

                break;
            }


            edge_idx++;
            v_pre = v_suf;
            succ_base_pre = get_base(label, word_count, edge_idx + k);
            pred_base_suf = get_base(label, word_count, edge_idx);
            v_suf.roll_forward(succ_base_pre);
        }
    }
}
*/


template <uint16_t k, bool Colored_>
void Subgraph<k, Colored_>::contract()
{
    Maximal_Unitig_Scratch<k> maximal_unitig;   // Scratch space to be used to construct maximal unitigs.
    uint64_t vertex_count = 0;  // Count of vertices processed.

    std::vector<Kmer<k>> V; // Vertices in an lm-tig.
    std::vector<uint64_t> H;    // Color-set hashes of the vertices in an lm-tig.

    auto& C = work_space.color_map();   // Color-set map.
    auto& in_process = work_space.in_process_arr();
    auto& S = work_space.set();
    if constexpr(Colored_)
        in_process.clear(),
        S.clear();

    for(auto p = M.begin(); p != M.end(); ++p)
    {
        const auto& v = ht_router::get_key(p);
        const auto& v_st = ht_router::get_val(p);
        assert(!v_st.is_discontinuous(side_t::front) || !v_st.is_discontinuous(side_t::back));

        if(v_st.is_visited())   // The containing maximal unitig has already been outputted.
            continue;

        if(v_st.is_isolated())
        {
            if(v_st.is_discontinuity()) // A potential phantom-edge for the discontinuity graph is incident to `v`.
                G.inc_potential_phantom_edge();

            isolated++;
            continue;
        }


        std::size_t b;  // Bucket-ID of the produced unitig-label.
        std::size_t b_idx;  // Index of the unitig-label in the corresponding bucket.
        if(extract_maximal_unitig(v, maximal_unitig, b, b_idx))
        {
            vertex_count += maximal_unitig.size();
            label_sz += maximal_unitig.size() + k - 1;
            mtig_c++;

            if constexpr(Colored_)
            {
                maximal_unitig.get_vertices_and_hashes(V, H);
                assert(V.size() == H.size() && V.size() == maximal_unitig.size());

                uint64_t h_last = ~H[0];
                Color_Coordinate c;
                const auto w_id = parlay::worker_id();
                for(std::size_t i = 0; i < V.size(); ++i)
                    if(H[i] != h_last)  // This is either a color-shifting vertex, or the first vertex in the lm-tig.
                    {
                        const LMTig_Coord lmtig_coord(b, b_idx, i);
                        const auto color_status = C.mark_in_process(H[i], w_id, c);
                        switch(color_status)
                        {
                        case Color_Status::undiscovered:    // Mark this vertex's color as of interest to extract later, and keep the vertex as pending.
                        {
                            S.insert(V[i]);
                            in_process.emplace_back(lmtig_coord, H[i]);
                            v_new_col_c++;
                            break;
                        }

                        case Color_Status::in_process:  // Keep this vertex pending and revisit it once its color is available.
                            if(c.processing_worker() != w_id)   // The color is not new globally, but new to this worker.
                                S.insert(V[i]),
                                v_new_col_c++;

                            in_process.emplace_back(lmtig_coord, H[i]);
                            break;

                        case Color_Status::discovered:  // This vertex's color is available.
                            G.add_color(b, b_idx, i, c);
                            v_old_col_c++;
                            break;
                        }

                        h_last = H[i];
                        color_shift_c++;
                    }
            }
        }
    }

    assert(vertex_count + isolated == M.size());
}


template <uint16_t k, bool Colored_>
void Subgraph<k, Colored_>::extract_new_colors()
{
    const auto t_0 = timer::now();

    collect_color_rels();
    const auto t_1 = timer::now();
    t_collect_rels += timer::duration(t_1 - t_0);

    collect_color_sets();

    const auto t_3 = timer::now();
    attach_colors_to_vertices();
    const auto t_4 = timer::now();
    t_attach += timer::duration(t_4 - t_3);
}


template <uint16_t k, bool Colored_>
void Subgraph<k, Colored_>::collect_color_rels()
{
if constexpr(Colored_)
{
    constexpr auto color_rel_bucket_c = Subgraphs_Scratch_Space<k, Colored_>::color_rel_bucket_c();
    auto& color_rel_bucket_arr = work_space.color_rel_bucket_arr();
    std::for_each(color_rel_bucket_arr.begin(), color_rel_bucket_arr.end(), [](auto& b){ b.clear(); });

    typename Super_Kmer_Bucket<Colored_>::Iterator super_kmer_it(B);    // Iterator over the weak super k-mers inducing this graph.
    typedef typename decltype(super_kmer_it)::label_unit_t label_unit_t;
    const auto word_count = super_kmer_it.super_kmer_word_count();  // Fixed number of words in a super k-mer label.

    Directed_Vertex<k> v;   // Current vertex in a scan over some super k-mer.
    Super_Kmer_Attributes<Colored_> att;
    const label_unit_t* label;

    auto& S = work_space.set();

    while(super_kmer_it.next(att, label))
    {
        const auto len = att.len();
        assert(len >= k);
        assert(len < 2 * (k - 1));
        const auto source = att.source();

        v.from_super_kmer(label, word_count);
        std::size_t kmer_idx = 0;

        while(true)
        {
            assert(kmer_idx + k - 1 < len);

            if(S.contains(v.canonical()))
            {
                static_assert(is_pow_2(color_rel_bucket_c));
                color_rel_bucket_arr[v.canonical().to_u64() & (color_rel_bucket_c - 1)].emplace(v.canonical(), source);
                color_rel_c++;
            }

            if(kmer_idx + k == len)
                break;

            const auto succ_base = get_base(label, word_count, kmer_idx + k);
            v.roll_forward(succ_base);
            kmer_idx++;
        }
    }
}
}


template <uint16_t k, bool Colored_>
void Subgraph<k, Colored_>::semi_sort_color_rels(const color_rel_t* const x, color_rel_t* const y, const std::size_t sz)
{
    auto& count_map = work_space.count_map();
    count_map.clear();

    // TODO: skip `x`; add and use ext-mem-bucket iterator.
    for(std::size_t i = 0; i < sz; ++i)
        count_map[x[i].first]++;

    uint32_t pref_sum = 0;
    for(auto& kmer_count : count_map)
    {
        const auto temp = kmer_count.second;
        kmer_count.second = pref_sum;
        pref_sum += temp;
    }


    for(std::size_t i = 0; i < sz; ++i)
    {
        auto& off = count_map[x[i].first];
        assert(off < sz);
        y[off] = x[i];
        off++;
    }
}


template <uint16_t k, bool Colored_>
void Subgraph<k, Colored_>::sort_color_set(std::vector<source_id_t>& color)
{
    auto const wv = work_space.bv().data();
    const auto word_c = work_space.bv().capacity();

    for(std::size_t i = 0; i < color.size(); ++i)
        assert(color[i] / 64 < word_c),
        wv[color[i] / 64] |= (uint64_t(1) << (color[i] & 63));


    color.clear();
    for(std::size_t word_idx = 0; word_idx < word_c; ++word_idx)
        while(wv[word_idx] != 0)
        {
            const auto bit_idx = __builtin_ctzll(wv[word_idx]);
            color.push_back(word_idx * 64 + bit_idx);
            wv[word_idx] &= (wv[word_idx] - 1);
        }
}


template <uint16_t k, bool Colored_>
void Subgraph<k, Colored_>::collect_color_sets()
{
if constexpr(Colored_)
{
    // TODO: the following two need not be subgraph-local, rather worker-local.
    fulgor::color_set_builder builder(G.max_source_id() + 1);
    fulgor::bit_vector_builder bvb;

    auto& color_rel_bucket_arr = work_space.color_rel_bucket_arr();
    auto& color_rel = work_space.color_rel_arr();
    auto& color_rel_collated = work_space.color_rel_collate_arr();
    auto& C = work_space.color_map();
    auto& color_bucket = work_space.color_repo().bucket();

    std::vector<source_id_t> src;   // Sources of the current vertex.

    for(std::size_t b = 0; b < color_rel_bucket_arr.size(); ++b)
    {
        const auto t_0 = timer::now();
        auto& color_rel_bucket = color_rel_bucket_arr[b];
        const auto color_rel_c = color_rel_bucket.size();
        color_rel.reserve_uninit(color_rel_c);
        color_rel_bucket.load(color_rel.data());

        color_rel_collated.reserve_uninit(color_rel_c);
        semi_sort_color_rels(color_rel.data(), color_rel_collated.data(), color_rel_c);

        const auto t_1 = timer::now();
        t_sort += timer::duration(t_1 - t_0);

        for(std::size_t i = 0, j; i < color_rel_c; i = j)
        {
            const auto& v = color_rel_collated[i].first;
            src.clear();
            src.push_back(color_rel_collated[i].second);

            for(j = i + 1; j < color_rel_c; ++j)
            {
                if(color_rel_collated[j].first != v)
                    break;

                // assert(color_rel_collated[j].second >= color_rel_collated[j - 1].second);   // Ensure sortedness of source-IDs, for compression.

                // This is required to deduplicate relations, if the latter color-set sort does not deduplicate them.
                // This works because: super k-mers from some source going into the same subgraph bucket cluster
                // together due to the partitioning policy (and its buffering scheme). Hence color-relations from
                // those super k-mers cluster together in the relationship list. This list is semi-sorted with
                // counting sort, which is stable, and hence the clusters are retained in the sorted output.
                // if(color_rel_collated[j].second != color_rel_collated[j - 1].second)
                src.push_back(color_rel_collated[j].second);
            }

            const auto color_idx = color_bucket.size();
            if(C.update_if_in_process(M[v].color_hash(), Color_Coordinate(parlay::worker_id(), color_idx)))
            {
                sort_color_set(src);
                assert(std::is_sorted(src.cbegin(), src.cend()));

                bvb.clear();
                builder.process(src.data(), src.size(), bvb);
                auto& bit_vec_words = bvb.bits();
                color_bucket.add(bit_vec_words.data(), bit_vec_words.size());
            }
        }

        const auto t_2 = timer::now();
        t_collect_sets += timer::duration(t_2 - t_1);
    }
}
}


template <uint16_t k, bool Colored_>
void Subgraph<k, Colored_>::attach_colors_to_vertices()
{
    auto& in_process = work_space.in_process_arr();
    auto& C = work_space.color_map();
    for(const auto& p : in_process)
    {
        const auto& lmtig_coord = p.first;
        const auto& h = p.second;

        const auto c = C.get(h);
        G.add_color(lmtig_coord.b(), lmtig_coord.idx(), lmtig_coord.off(), c);
    }

    in_process.clear();
}


template <uint16_t k, bool Colored_>
Subgraphs_Scratch_Space<k, Colored_>::Subgraphs_Scratch_Space(const std::size_t max_sz, const std::string& color_rel_bucket_pref):
      in_process_arr_(parlay::num_workers())
{
    map_.reserve(parlay::num_workers());
    for(std::size_t i = 0; i < parlay::num_workers(); ++i)
        HT_Router<k, Colored_>::add_HT(map_, max_sz);

    if constexpr(Colored_)
    {
        color_rel_bucket_arr_.resize(parlay::num_workers());
        color_rel_arr_.resize(parlay::num_workers());
        color_rel_collate_arr_.resize(parlay::num_workers());
        count_map_.resize(parlay::num_workers());

        static_assert(is_pow_2(color_rel_bucket_c_));
        for(std::size_t w = 0; w < parlay::num_workers(); ++w)
        {
            const auto color_rel_dir = color_rel_bucket_pref + "/" + std::to_string(w);
            std::filesystem::create_directories(color_rel_dir);
            for(std::size_t b = 0; b < color_rel_bucket_c_; ++b)
                color_rel_bucket_arr_[w].unwrap().
                    emplace_back(color_rel_bucket_t(color_rel_dir + "/b_" + std::to_string(b),
                                    color_rel_buf_sz / sizeof(color_rel_t)));
        }

        bv_.resize(parlay::num_workers());
    }

    set_.resize(parlay::num_workers());
}


template <uint16_t k, bool Colored_>
auto Subgraphs_Scratch_Space<k, Colored_>::map() -> map_t&
{
    assert(map_.size() == parlay::num_workers());
    return map_[parlay::worker_id()].unwrap();
}


template <uint16_t k, bool Colored_>
auto Subgraphs_Scratch_Space<k, Colored_>::color_rel_bucket_arr() -> color_rel_bucket_arr_t&
{
    // assert(Colored_);
    assert(color_rel_bucket_arr_.size() == parlay::num_workers());
    return color_rel_bucket_arr_[parlay::worker_id()].unwrap();
}


template <uint16_t k, bool Colored_>
Color_Table& Subgraphs_Scratch_Space<k, Colored_>::color_map()
{
    // assert(Colored_);
    return M_c;
}


template <uint16_t k, bool Colored_>
auto Subgraphs_Scratch_Space<k, Colored_>::in_process_arr() -> in_process_arr_t&
{
    // assert(Colored_);
    assert(in_process_arr_.size() == parlay::num_workers());
    return in_process_arr_[parlay::worker_id()].unwrap();
}


template <uint16_t k, bool Colored_>
auto Subgraphs_Scratch_Space<k, Colored_>::color_rel_arr() -> color_rel_arr_t&
{
    // assert(Colored_);
    assert(color_rel_arr_.size() == parlay::num_workers());
    return color_rel_arr_[parlay::worker_id()].unwrap();
}


template <uint16_t k, bool Colored_>
auto Subgraphs_Scratch_Space<k, Colored_>::color_rel_collate_arr() -> color_rel_arr_t&
{
    // assert(Colored_);
    assert(color_rel_collate_arr_.size() == parlay::num_workers());
    return color_rel_collate_arr_[parlay::worker_id()].unwrap();
}


template <uint16_t k, bool Colored_>
auto Subgraphs_Scratch_Space<k, Colored_>::count_map() -> count_map_t&
{
    // assert(Colored_);
    assert(count_map_.size() == parlay::num_workers());
    return count_map_[parlay::worker_id()].unwrap();
}


template <uint16_t k, bool Colored_>
auto Subgraphs_Scratch_Space<k, Colored_>::bv() -> bit_vector_t&
{
    assert(Colored_);
    assert(bv_.size() == parlay::num_workers());
    return bv_[parlay::worker_id()].unwrap();
}


template <uint16_t k, bool Colored_>
auto Subgraphs_Scratch_Space<k, Colored_>::set() -> set_t&
{
    // assert(Colored_);
    assert(set_.size() == parlay::num_workers());
    return set_[parlay::worker_id()].unwrap();
}


template <uint16_t k, bool Colored_>
Color_Repo& Subgraphs_Scratch_Space<k, Colored_>::color_repo()
{
    assert(Colored_);
    return color_repo_;
}

}

/*
```
for b in B:    // B: collection of buckets
    for s in b:    // s: super k-mer
        i = 0
        for x in s:    // x: k-mer
            is_can = x.is_canonical()
            pn = pred_nuc(x) // or, pn = s[i - 1]
            sn = succ_nuc(x) // or, sn = s[i + k]
            front = (is_can ? pn : bar(sn))
            back  = (is_can ? sn : bar(pn))
            update HT for x.canonical() with front and back
            i = i + 1
```
*/



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BOOL, cuttlefish::Subgraph)
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_PER_BOOL, cuttlefish::Subgraphs_Scratch_Space)

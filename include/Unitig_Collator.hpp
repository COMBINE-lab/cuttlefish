
#ifndef UNITIG_COLLATOR_HPP
#define UNITIG_COLLATOR_HPP



#include "dBG_Contractor.hpp"
#include "Path_Info.hpp"
#include "Unitig_Coord_Bucket.hpp"
#include "DNA_Utility.hpp"
#include "Directed_Vertex.hpp"
#include "Discontinuity_Graph.hpp"
#include "Color_Encoding.hpp"
#include "globals.hpp"
#include "utility.hpp"

#include <cstddef>
#include <cstdint>
#include <atomic>
#include <cstring>
#include <string_view>
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>


class Data_Logistics;


namespace cuttlefish
{

// =============================================================================
// Collates locally-maximal unitigs from different unitig-buckets as per their
// path-information in a discontinuity graph of `k`-mers.
template <uint16_t k, bool Colored_>
class Unitig_Collator
{
private:

    Discontinuity_Graph<k, Colored_>& G;    // The discontinuity-graph.

    typedef typename dBG_Contractor<k>::unitig_path_info_t unitig_path_info_t;

    typedef typename dBG_Contractor<k>::P_e_t P_e_t;
    P_e_t& P_e; // `P_e[b]` contains path-info for edges in bucket `b`.

    const std::string lmtig_buckets_path;   // Path-prefix to the lm-tig buckets.
    const std::string unitig_coord_buckets_path;    // Path-prefix to the unitig-coordinate buckets produced in map-reduce.

    std::size_t max_bucket_sz;  // Maximum size of the edge-buckets.

    const std::size_t max_unitig_bucket_count;  // Number of buckets storing literal globally-maximal unitigs.
    // std::vector<Padded<Unitig_Coord_Bucket<k>>> max_unitig_bucket; // Key-value collation buckets for lm-unitigs.
    std::vector<Padded<Unitig_Coord_Bucket_Concurrent<k, Colored_>>> max_unitig_bucket; // Key-value collation buckets for lm-unitigs.

    typedef typename dBG_Contractor<k>::op_buf_list_t op_buf_list_t;
    op_buf_list_t& op_buf;  // Worker-specific output buffers.

    std::atomic_uint64_t phantom_c_;    // Number of phantom unitigs observed.

    class Maximal_Unitig;


    // Maps each locally-maximal unitig to its maximal unitig's corresponding
    // bucket.
    void map();

    // Reduces each maximal unitig bucket to its contained maximal unitigs.
    void reduce();

    // Loads the path-info of edges from bucket `b` into the table `M`, and
    // returns the size of the bucket. Uses the buffer `buf` to transfer the
    // information from the bucket to the table.
    std::size_t load_path_info(std::size_t b, Path_Info<k>* M, Buffer<unitig_path_info_t>& buf);

    // Loads the vertex-color mappings from bucket `b` into `buf`, and returns
    // the size of the bucket.
    std::size_t load_vertex_color_mapping(std::size_t b, Buffer<Vertex_Color_Mapping>& buf);

    // Emits the trivially maximal unitigs to the output stream. Only
    // applicable in the colored case.
    void emit_trivial_mtigs();


public:

    // Constructs a unitig-collator for unitigs with their associated path-info
    // at `P_e`, i.e. `P_e[b]` contains path-information of the unitigs'
    // corresponding edges at bucket `b`. `logistics` is the data logistics
    // manager for the algorithm execution. Worker-specific maximal unitigs are
    // written to the buffers in `op_buf`. `gmtig_bucket_count` many buckets
    // are used to partition the lm-tigs to their maximal unitigs. `G` is the
    // associated discontinuity graph.
    Unitig_Collator(Discontinuity_Graph<k, Colored_>& G, P_e_t& P_e, const Data_Logistics& logistics, op_buf_list_t& op_buf, std::size_t gmtig_bucket_count);

    // Collates the locally-maximal unitigs into global ones.
    void collate();
};


// Sequences associated to a maximal unitig.
template  <uint16_t k, bool Colored_>
class Unitig_Collator<k, Colored_>::Maximal_Unitig
{
private:

    const Unitig_Collator<k, Colored_>& collator;   // The unitig-collator using this maximal unitig.

    Buffer<char> label_;    // Label-sequence.
    std::size_t sz; // Size of the label.

    std::vector<Unitig_Color> color_;   // Color-sequence.

    Buffer<char> cycle_buf; // Working-space to process cyclic maximal unitigs.


    // Appends the `len_s`-length sequence `s` to the label. `RC_` specifies
    // whether `s` needs to be reverse-complemented or not.
    template <bool RC_> void append(const char* s, std::size_t len_s);

public:

    Maximal_Unitig(const Unitig_Collator& collator): collator(collator) {}

    // Returns the label sequence.
    auto data() const { return label_.data(); }

    // Returns the size of the label.
    auto size() const { return  sz; }

    // Clears the label.
    auto clear() { sz = 0; if constexpr(Colored_) color_.clear(); }

    // Returns `true` iff the label is empty.
    auto empty() const { return sz == 0; }

    // Returns the color-sequence.
    const auto& color() const { return color_; }

    // Initializes the label with the sequence `unitig`. `rc` specifies whether
    // `unitig` needs to be put in its reverse-complemented form.
    void init(const std::string_view& unitig, bool rc);

    // Initializes with the sequence `unitig` and its color sequence `color`
    // (of size `color_c`). `rc` specifies whether the sequences need to be put
    // in reverse-complemented forms.
    void init(const std::string_view& unitig, bool rc, const Unitig_Color* color, std::size_t color_c);

    // Appends the `k`-overlapping sequence `unitig` to the label. `rc`
    // specifies whether `unitig` needs to be added in its reverse-
    // complemented form.
    void append(const std::string_view& unitig, bool rc);

    // Appends the `k`-overlapping sequence `unitig` and the `color_c`-sized
    // color-sequence `color`. `rc` specifies whether the sequences need to be
    // added in reverse-complemented forms.
    void append(const std::string_view& unitig, bool rc, const Unitig_Color* color, std::size_t color_c);

    // Removes the last k-mer.
    void pop_back();

    // Transforms the label to its canonical form.
    void canonicalize();

    // Transforms the label to its canonical form given that the maximal unitig
    // is a cycle.
    void canonicalize_cycle();
};


template <uint16_t k, bool Colored_>
template <bool RC_>
inline void Unitig_Collator<k, Colored_>::Maximal_Unitig::append(const char* const s, const std::size_t len_s)
{
    if constexpr(!RC_)
        std::memcpy(label_.data() + sz, s, len_s);
    else
        for(std::size_t i = 0; i < len_s; ++i)
            label_[sz + i] = DNA_Utility::complement(s[len_s - 1 - i]);

    sz += len_s;
}


template <uint16_t k, bool Colored_>
inline void Unitig_Collator<k, Colored_>::Maximal_Unitig::init(const std::string_view& unitig, const bool rc)
{
    assert(!Colored_);
    assert(unitig.length() >= k);

    clear();

    label_.reserve_uninit(unitig.length());
    rc ?    append<true>(unitig.data(), unitig.length()) :
            append<false>(unitig.data(), unitig.length());
}


template <uint16_t k, bool Colored_>
inline void Unitig_Collator<k, Colored_>::Maximal_Unitig::init(const std::string_view& unitig, const bool rc, const Unitig_Color* const color, const std::size_t color_c)
{
    assert(Colored_);
    assert(unitig.length() >= k);
    if(color_c == 0)    // Unitig induced by a phantom k-mer—no color is observed.
        assert(unitig.length() == k && collator.G.is_discontinuity(unitig.data())); // Some necessary conditions.

    clear();

    label_.reserve_uninit(unitig.length());
    color_.insert(color_.end(), color, color + color_c);

    if(!rc)
        append<false>(unitig.data(), unitig.length());
    else
    {
        append<true>(unitig.data(), unitig.length());

        const auto vertex_c = unitig.length() - k + 1;
        color_.emplace_back(vertex_c, Color_Coordinate());
        for(std::size_t i = 0; i < color_.size() - 1; i++)
        {
            assert(color_[i].off() < color_[i + 1].off());
            const auto off_r_to_l = color_[i + 1].off() - 1;    // Offset if the colors' RL-encoding is done right-to-left.
            const auto off_rev_cmp = (vertex_c - off_r_to_l) - 1;   // Offset in the reverse-complemented form of the sequence.
            color_[i].set_off(off_rev_cmp);
        }
        color_.pop_back();

        std::reverse(color_.begin(), color_.end());
    }
}


template <uint16_t k, bool Colored_>
inline void Unitig_Collator<k, Colored_>::Maximal_Unitig::append(const std::string_view& unitig, const bool rc)
{
    assert(sz >= k); assert(unitig.length() >= k);

    label_.reserve(sz + unitig.length() - k);
    if(!rc)
        append<false>(unitig.data() + k, unitig.length() - k);
    else
        append<true>(unitig.data(), unitig.length() - k);
}


template <uint16_t k, bool Colored_>
inline void Unitig_Collator<k, Colored_>::Maximal_Unitig::append(const std::string_view& unitig, const bool rc, const Unitig_Color* color, std::size_t color_c)
{
    assert(Colored_);
    assert(sz >= k); assert(unitig.length() >= k);

    if(CF_UNLIKELY(color_c == 0))   // Unitig induced by a phantom k-mer—no color is observed.
    {
        assert(unitig.length() == k && collator.G.is_discontinuity(unitig.data()));
        append(unitig, rc);
        return;
    }

    const auto prev_v_c = sz - k + 1;
    const auto prev_col_c = color_.size();
    const auto unitig_v_c = unitig.length() - k + 1;
    label_.reserve(sz + unitig.length() - k);
    if(!rc)
    {
        append<false>(unitig.data() + k, unitig.length() - k);

        if(!color_.empty() && color_.back().c() == color[0].c())
            color++, color_c--;

        color_.insert(color_.end(), color, color + color_c);
        std::for_each(color_.begin() + prev_col_c, color_.end(),
        [&](auto& c)
        {
            assert(c.off() < unitig_v_c);
            c.set_off((prev_v_c + c.off()) - 1);    // Adjust based on the existing colors. Offsets from the new unitig gets left-shifted due to vertex-overlap.
        });
    }
    else
    {
        append<true>(unitig.data(), unitig.length() - k);

        color_.insert(color_.end(), color, color + color_c);
        color_.emplace_back(unitig_v_c, Color_Coordinate());
        for(std::size_t i = prev_col_c; i < color_.size() - 1; ++i)
        {
            assert(color_[i].off() < color_[i + 1].off());
            const auto off_r_to_l = color_[i + 1].off() - 1;    // Offset if the colors' RL-encoding is done right-to-left.
            const auto off_rev_cmp = (unitig_v_c - off_r_to_l) - 1; // Offset in the reverse-complemented form of the sequence.
            const auto off_adjusted = (prev_v_c + off_rev_cmp) - 1; // Offset adjusted based on the existing colors.
                                                                    // Offsets from the new unitig gets left-shifted due to vertex-overlap.

            color_[i].set_off(off_adjusted);
        }
        color_.pop_back();

        if(prev_col_c > 0 && color_[prev_col_c - 1].c() == color_.back().c())
            color_.pop_back();

        std::reverse(color_.begin() + prev_col_c, color_.end());
    }
}


template <uint16_t k, bool Colored_>
inline void Unitig_Collator<k, Colored_>::Maximal_Unitig::canonicalize()
{
    for(std::size_t i = 0; i < k; ++i)
    {
        const auto b_fw = label_[i];
        const auto b_bw = DNA_Utility::complement(label_[sz - 1 - i]);

        if(b_fw < b_bw) // Already in canonical form.
            return;

        if(b_fw > b_bw) // Reverse-complement is the canonical form.
        {
            for(std::size_t j = 0; j < sz / 2; ++j)
            {
                const auto b_l = label_[j];
                const auto b_r = label_[sz - 1 - j];
                label_[j] = DNA_Utility::complement(b_r),
                label_[sz - 1 - j] = DNA_Utility::complement(b_l);
            }

            if(sz & 1)
                label_[sz / 2] = DNA_Utility::complement(label_[sz / 2]);

            if constexpr(Colored_)
            {
                const int64_t vertex_c = sz - k + 1;
                color_.emplace_back(vertex_c, Color_Coordinate());
                for(std::size_t i = 0; i < color_.size() - 1; ++i)
                {
                    assert(color_[i].off() < color_[i + 1].off());
                    const auto off_r_to_l = color_[i + 1].off() - 1;
                    const auto off_rev_cmp = (vertex_c - off_r_to_l) - 1;
                    color_[i].set_off(off_rev_cmp);
                }
                color_.pop_back();

                std::reverse(color_.begin(), color_.end());
            }

            return;
        }
    }
}


template <uint16_t k, bool Colored_>
inline void Unitig_Collator<k, Colored_>::Maximal_Unitig::pop_back()
{
    assert(sz > k);

    const auto v_c = sz - k;
    sz--;
    if constexpr(Colored_)
        if(color_.back().off() == v_c - 1)
            color_.pop_back();
}


template <uint16_t k, bool Colored_>
inline void Unitig_Collator<k, Colored_>::Maximal_Unitig::canonicalize_cycle()
{
    Directed_Vertex<k> v(label_.data());
    Kmer<k> min_fw(v.kmer());   // Minimum k-mer in the forward strand.
    Kmer<k> min_bw(v.kmer_bar());   // Minimum k-mer in the backward strand.
    std::size_t min_idx_fw = 0;
    std::size_t min_idx_bw = k - 1;

    for(std::size_t i = 1; i + k <= sz; ++i)
    {
        v.roll_forward(DNA_Utility::map_base(label_[i + k - 1]));

        if(min_fw > v.kmer())
            min_fw = v.kmer(), min_idx_fw = i;

        if(min_bw > v.kmer_bar())
            min_bw = v.kmer_bar(), min_idx_bw = i + k - 1;
    }


    cycle_buf.reserve_uninit(sz);
    if(min_fw < min_bw)
    {
        const auto len_r = sz - min_idx_fw;
        const auto len_l = sz - len_r;
        std::memcpy(cycle_buf.data(), label_.data() + len_l, len_r);
        std::memcpy(cycle_buf.data() + len_r, label_.data() + k - 1, len_l);
    }
    else
    {
        const int64_t len_l = min_idx_bw + 1;
        const int64_t len_r = sz - len_l;
        std::size_t p = 0;

        for(int64_t i = len_l - 1; i >= 0; --i)
            cycle_buf[p++] = DNA_Utility::complement(label_[i]);

        for(int64_t i = len_r - 1; i >= 0; --i)
            cycle_buf[p++] = DNA_Utility::complement(label_[len_l + i - (k - 1)]);
    }

    std::memcpy(label_.data(), cycle_buf.data(), sz);

    // TODO: no fixing-policy for colors yet.
}

}



#endif

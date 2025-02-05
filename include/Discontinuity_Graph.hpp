
#ifndef DISCONTINUITY_GRAPH_HPP
#define DISCONTINUITY_GRAPH_HPP



#include "Kmer.hpp"
#include "Edge_Matrix.hpp"
#include "Unitig_File.hpp"
#include "Maximal_Unitig_Scratch.hpp"
#include "Ext_Mem_Bucket.hpp"
#include "Color_Encoding.hpp"
#include "Build_Params.hpp"
#include "globals.hpp"
#include "parlay/parallel.h"
#include "utility.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/archives/binary.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <utility>
#include <atomic>


class Data_Logistics;


namespace cuttlefish
{

class Vertex_Color_Mapping;


// =============================================================================
// A representation of a discontinuity graph of `k`-mers. `Colored_` denotes
// whether the edges have colors or not.
template <uint16_t k, bool Colored_>
class Discontinuity_Graph
{
private:

    // k-mer (super-)label of the ϕ-vertex in the discontinuity graph.  // TODO: revisit; almost sure we don't need this.
    static constexpr const char phi_label[] =   "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    static const Kmer<k> phi_;  // ϕ k-mer connected to each chain-end in the discontinuity graph.

    const uint16_t min_len; // Size of the l-minimizers.

    Edge_Matrix<k> E_;  // Edge-matrix of the discontinuity graph.

    Unitig_Write_Distributor lmtigs;    // Distribution-manager for the writes of locally maximal unitigs' labels.

    // TODO: check logs to see if this is a bottleneck.
    std::atomic_uint64_t phantom_edge_count_;   // Number of potential phantom edges identified.

    const source_id_t max_source_id_;   // Maximum source-ID, used for coloring.
    std::vector<Padded<Ext_Mem_Bucket<Vertex_Color_Mapping>>> vertex_color_map_;    // Buckets of vertex-color mappings.


public:

    // Constructs a discontinuity graph object that operates with the required
    // parameters in `params`. `logistics` is the data logistics manager for
    // the algorithm execution.
    Discontinuity_Graph(const Build_Params& params, const Data_Logistics& logistics);

    // Deserializes the discontinuity graph from the `cereal` archive `archive`.
    Discontinuity_Graph(cereal::BinaryInputArchive& archive);

    // Returns the ϕ k-mer connected to each chain-end in the discontinuity
    // graph.
    static const Kmer<k>& phi() { return phi_; }

    // Returns the edge-matrix of the graph.
    const Edge_Matrix<k>& E() const { return E_; }

    // Returns the edge-matrix of the graph.
    Edge_Matrix<k>& E() { return E_; }

    // Returns the maximum source-ID, used for coloring.
    auto max_source_id() const { return max_source_id_; }

    // Returns the number of potential phantom edges identified.
    uint64_t phantom_edge_upper_bound() const;

    // Adds the edge `({(u, s_u), (v, s_v)}, 1)` to the graph. `u_is_phi` and
    // `v_is_phi` denote whether `u` and `v` are `ϕ` respectively. The locally-
    // maximal unitig corresponding to the edge is `mtig`. The edge should be
    // an original edge of the graph. Returns `(b, b_idx)`, where the deposited
    // lm-tig is put into bucket `b` at index `b_idx`.
    std::pair<std::size_t, std::size_t> add_edge(const Kmer<k>& u, side_t s_u, const Kmer<k>& v, side_t s_v, bool u_is_phi, bool v_is_phi, const Maximal_Unitig_Scratch<k>& mtig);

    // Adds the edge `({ϕ, (v, s_v)}, 1)` to the graph. The edge should be an
    // original edge of the graph.
    void add_edge(const Kmer<k>& v, side_t s_v);

    // Adds the edge `({(u, s_u), (v, s_v)}, w)` to the graph. `u_is_phi` and
    // `v_is_phi` denote whether `u` and `v` are `ϕ` respectively. The edge
    // should be a contracted edge and not an original one.
    void add_edge(const Kmer<k>& u, side_t s_u, const Kmer<k>& v, side_t s_v, weight_t w, bool u_is_phi, bool v_is_phi);

    // Adds the trivial maximal unitig `mtig` to the graph. Nothing is added to
    // the graph per se, just the unitig label is stored. Returns `(b, b_idx)`,
    // where the deposited `mtig` is put into bucket `b` at index `b_idx`.
    std::pair<std::size_t, std::size_t> add_trivial_mtig(const Maximal_Unitig_Scratch<k>& mtig);

    // Adds the color-coordinate `c` to the `b`'th unitig bucket, where the
    // `b_idx`'th unitig has the corresponding color at offset `off`.
    void add_color(uint16_t b, uint32_t b_idx, uint16_t off, const Color_Coordinate& c);

    // Returns the `b`'th vertex-color mapping bucket.
    auto& vertex_color_map(const std::size_t b) { assert(b < vertex_color_map_.size()); return vertex_color_map_[b].unwrap(); }

    // Increments the potential phantom edge count.
    void inc_potential_phantom_edge() { phantom_edge_count_++; }

    // Closes and releases the streams depositing to the graph.
    void close();

    // Returns a tight upper bound of the maximum number of vertices in a
    // partition.
    std::size_t vertex_part_size_upper_bound() const;

    // Returns `true` iff the k-mer at `seq` is a discontinuity vertex.
    bool is_discontinuity(const char* seq) const;

    // Returns the resident set size of the space-dominant components of the
    // graph.
    std::size_t RSS() const;

    // (De)serializes the graph from / to the `cereal` archive `archive`.
    template <typename T_archive_> void serialize(T_archive_& archive);
};


template <uint16_t k, bool Colored_> const Kmer<k> Discontinuity_Graph<k, Colored_>::phi_(phi_label);


template <uint16_t k, bool Colored_>
inline std::pair<std::size_t, std::size_t> Discontinuity_Graph<k, Colored_>::add_edge(const Kmer<k>& u, const side_t s_u, const Kmer<k>& v, const side_t s_v, const bool u_is_phi, const bool v_is_phi, const Maximal_Unitig_Scratch<k>& mtig)
{
    const auto w_id = parlay::worker_id();
    const auto coord = lmtigs.add(w_id, mtig);
    E_.add(u, s_u, v, s_v, 1, coord.first, coord.second, u_is_phi, v_is_phi);

    return coord;
}


template <uint16_t k, bool Colored_>
inline void Discontinuity_Graph<k, Colored_>::add_edge(const Kmer<k>& v, const side_t s_v)
{
    const auto w_id = parlay::worker_id();
    const auto coord = lmtigs.add(w_id, s_v == side_t::front ? v : v.reverse_complement());
    E_.add(phi_, side_t::back, v, s_v, 1, coord.first, coord.second, true, false);
}


template <uint16_t k, bool Colored_>
inline void Discontinuity_Graph<k, Colored_>::add_edge(const Kmer<k>& u, const side_t s_u, const Kmer<k>& v, const side_t s_v, const weight_t w, const bool u_is_phi, const bool v_is_phi)
{
    // Edge-partition 0 associates to edges that do not have any corresponding lm-tig (i.e. has weight > 1).
    E_.add(u, s_u, v, s_v, w, 0, 0, u_is_phi, v_is_phi);
}


template <uint16_t k, bool Colored_>
inline std::pair<std::size_t, std::size_t> Discontinuity_Graph<k, Colored_>::add_trivial_mtig(const Maximal_Unitig_Scratch<k>& mtig)
{
    return lmtigs.add_trivial_mtig(parlay::worker_id(), mtig);
}


template <uint16_t k, bool Colored_>
inline void Discontinuity_Graph<k, Colored_>::add_color(const uint16_t b, const uint32_t b_idx, const uint16_t off, const Color_Coordinate& c)
{
    vertex_color_map_[b].unwrap().emplace(b_idx, off, c);
}


template <uint16_t k, bool Colored_>
template <typename T_archive_>
inline void Discontinuity_Graph<k, Colored_>::serialize(T_archive_& archive)
{
    uint64_t phantom_edge_c = phantom_edge_count_;
    archive(type::mut_ref(min_len), E_, lmtigs, phantom_edge_c, type::mut_ref(max_source_id_), vertex_color_map_);
    phantom_edge_count_ = phantom_edge_c;
}

}



#endif

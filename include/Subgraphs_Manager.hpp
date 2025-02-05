
#ifndef SUBGRAPHS_MANAGER_HPP
#define SUBGRAPHS_MANAGER_HPP



#include "Atlas.hpp"
#include "HyperLogLog.hpp"
#include "Directed_Vertex.hpp"
#include "DNA_Utility.hpp"
#include "Discontinuity_Graph.hpp"
#include "Async_Logger_Wrapper.hpp"
#include "Character_Buffer.hpp"
#include "utility.hpp"
#include "globals.hpp"

#include <cstdint>
#include <cstddef>
#include <deque>
#include <limits>
#include <atomic>
#include <string>
#include <vector>
#include <type_traits>


class Data_Logistics;


namespace cuttlefish
{

// =============================================================================
// Manager for the subgraphs of the de Bruijn graph—manages their super k-mer
// based sequence representations in buckets, constructs them from these
// representations, and contracts them into their compacted form. `Colored_`
// denotes whether the vertices in the graph have associated colors.
template <uint16_t k, bool Colored_>
class Subgraphs_Manager
{
private:

    const std::string path_pref;    // Path-prefix to the subgraph atlases.
    const std::string color_rel_path_pref;  // Path-prefix to color-relationship buckets.
    const uint16_t l;   // `l`-minimizer size to partition the graph.

    typedef Atlas<Colored_> atlas_t;
    static constexpr std::size_t chunk_bytes = 1024 * 1024; // 1 MB chunk capacity for each atlas.
    static constexpr std::size_t w_chunk_bytes = 64 * 1024; // 64 KB worker-local chunk capacity in each atlas.
    std::vector<Padded<atlas_t>> atlas; // Super k-mer buckets for the subgraph atlases.

    std::vector<Padded<HyperLogLog>> HLL;   // `HLL[g]` is the cardinality-estimator for subgraph `g`.

    Discontinuity_Graph<k, Colored_>& G_;   // The discontinuity graph.

    std::atomic_uint64_t trivial_mtig_count_;   // Number of trivial maximal unitigs in the subgraphs (i.e. also maximal unitigs in the supergraph).
    std::atomic_uint64_t icc_count_;    // Number of trivial maximal unitigs in the subgraphs that are ICCs.

    // TODO: move out the following to some CF3-centralized location.

    typedef Async_Logger_Wrapper sink_t;
    typedef Character_Buffer<sink_t> op_buf_t;
    typedef std::vector<Padded<op_buf_t>> op_buf_list_t;
    op_buf_list_t& op_buf; // Worker-specific output buffers.

    const std::string color_path_pref;  // Path-prefix to the output color buckets.

    // Adds the label `seq` and length `len` to the HLL estimate of the
    // subgraph `g` of the de Bruijn graph.
    void add_to_HLL(std::size_t g, const char* seq, std::size_t len);


public:

    // Constructs a manager for the subgraphs of a de Bruijn graph which is
    // partitioned according to `l`-minimizers. `logistics` is the data-
    // logistics manager for the algorithm execution. The discontinuity-graph
    // is produced at `G` without false-phantom edges. Worker-specific
    // trivially maximal unitigs are written to the buffers in `op_buf`.
    Subgraphs_Manager(const Data_Logistics& logistics, uint16_t l, Discontinuity_Graph<k, Colored_>& G, op_buf_list_t& op_buf);

    // Returns the number of subgraphs.
    auto graph_count() const { return Atlas<Colored_>::graph_count(); }

    // Returns the discontinuity graph.
    const auto& G() const { return G_; }

    // Adds a (weak) super k-mer to the subgraph `g` of the de Bruijn graph
    // with label `seq` and length `len`. The markers `l_disc` and `r_disc`
    // denote whether the left and the right ends of the (weak) super k-mer are
    // discontinuous or not.
    template <bool C_ = Colored_, std::enable_if_t<!C_, int> = 0>
    void add_super_kmer(std::size_t g, const char* seq, std::size_t len, bool l_disc, bool r_disc);

    // Adds a (weak) super k-mer to the subgraph `g` of the de Bruijn graph
    // with label `seq` and length `len` from source-ID `source`. The markers
    // `l_disc` and `r_disc` denote whether the left and the right ends of the
    // (weak) super k-mer are discontinuous or not.
    template <bool C_ = Colored_, std::enable_if_t<C_, int> = 0>
    void add_super_kmer(std::size_t g, const char* seq, std::size_t len, source_id_t source, bool l_disc, bool r_disc);

    // Collates the current super k-mer buffers in each subgraph per their
    // source-IDs into external-memory buckets. The source-IDs are supposed to
    // be in the range `[src_min, src_max]`.
    void collate_super_kmer_buffers(source_id_t src_min, source_id_t src_max);

    // Flushes the buffer of the `w`'th worker to the subgraphs in the atlas if
    // it is overflown.
    void flush_worker_if_req(std::size_t w);

    // Finalizes the subgraphs for iteration—no more content should be added
    // after this.
    void finalize();

    // Returns the largest estimated size of any subgraph.
    uint64_t estimate_size_max() const;

    // Constructs and contracts each subgraph.
    void process();

    // Returns the number of trivial maximal unitigs in the subgraphs (i.e. also
    // maximal unitigs in the supergraph).
    uint64_t trivial_mtig_count() const;

    // Returns the number of trivial maximal unitigs in the subgraphs that are
    // ICCs.
    uint64_t icc_count() const;

    // Returns the subgraph ID for a minimizer with 64-bit hash value `h`.
    uint64_t graph_ID(uint64_t h) const { return h & (Atlas<Colored_>::graph_count() - 1); }

    // Returns the resident set size of the space-dominant components of the
    // subgraphs-manager.
    std::size_t RSS() const;
};


template <uint16_t k, bool Colored_>
template <bool C_, std::enable_if_t<!C_, int>>
inline void Subgraphs_Manager<k, Colored_>::add_super_kmer(const std::size_t g, const char* const seq, const std::size_t len, const bool l_disc, const bool r_disc)
{
    assert(len >= k);

    const auto a = Atlas<Colored_>::atlas_ID(g);
    auto& bucket = atlas[a].unwrap();
    bucket.add(seq, len, l_disc, r_disc, g);

    // add_to_HLL(g, seq, len);
}


template <uint16_t k, bool Colored_>
template <bool C_, std::enable_if_t<C_, int>>
inline void Subgraphs_Manager<k, Colored_>::add_super_kmer(const std::size_t g, const char* const seq, const std::size_t len, const source_id_t source, const bool l_disc, const bool r_disc)
{
    assert(len >= k);

    const auto a = Atlas<Colored_>::atlas_ID(g);
    auto& bucket = atlas[a].unwrap();
    bucket.add(seq, len, source, l_disc, r_disc, g);

    // add_to_HLL(g, seq, len);
}


template <uint16_t k, bool Colored_>
inline void Subgraphs_Manager<k, Colored_>::add_to_HLL(const std::size_t g, const char* const seq, const std::size_t len)
{
    auto& hll = HLL[g].unwrap();
    Directed_Vertex<k> v{Kmer<k>(seq)};
    constexpr auto u32_mask = std::numeric_limits<uint32_t>::max();
    std::size_t next_idx = k;
    while(true)
    {
        hll.add(v.canonical().to_u64() & u32_mask);

        if(next_idx == len)
            break;

        v.roll_forward(DNA_Utility::map_base(seq[next_idx]));
        next_idx++;
    }
}

}



#endif

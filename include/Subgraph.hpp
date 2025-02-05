
#ifndef SUBGRAPH_HPP
#define SUBGRAPH_HPP



#include "dBG_Contractor.hpp"
#include "State_Config.hpp"
#include "Directed_Vertex.hpp"
#include "Kmer.hpp"
#include "Kmer_Hasher.hpp"
#include "DNA_Utility.hpp"
#include "Kmer_Hashtable.hpp"
#include "Unitig_Scratch.hpp"
#include "Maximal_Unitig_Scratch.hpp"
#include "Discontinuity_Graph.hpp"
#include "Ext_Mem_Bucket.hpp"
#include "Color_Table.hpp"
#include "Color_Repo.hpp"
#include "dBG_Utilities.hpp"
#include "utility.hpp"
#include "globals.hpp"
#include "unordered_dense/unordered_dense.h"

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <tuple>
#include <utility>
#include <unordered_map>


namespace cuttlefish
{

template <bool Colored_> class Super_Kmer_Bucket;

template <uint16_t k, bool Colored_> class HT_Router;

enum class Walk_Termination;    // Type of scenarios how a unitig-walk terminates in the subgraph.

class LMTig_Coord;  // lm-tig coordinate of a vertex (k-mer).


// Working space for workers processing different subgraphs.
template <uint16_t k, bool Colored_>
class Subgraphs_Scratch_Space
{
public:

    // typedef std::unordered_map<Kmer<k>, State_Config, Kmer_Hasher<k>> map_t;
    typedef ankerl::unordered_dense::map<Kmer<k>, State_Config<Colored_>, Kmer_Hasher<k>> map_t;
    // typedef Kmer_Hashtable<k, Colored_> map_t;

    typedef std::pair<LMTig_Coord, uint64_t> in_process_t;  // Vertex's lm-tig coordinate and color-hash.
    typedef std::vector<in_process_t> in_process_arr_t;

    typedef std::pair<Kmer<k>, source_id_t> color_rel_t;
    typedef Ext_Mem_Bucket<color_rel_t> color_rel_bucket_t;
    typedef std::vector<color_rel_bucket_t> color_rel_bucket_arr_t;
    typedef Buffer<color_rel_t> color_rel_arr_t;
    typedef ankerl::unordered_dense::map<Kmer<k>, uint32_t, Kmer_Hasher<k>> count_map_t;

    typedef Buffer<uint64_t> bit_vector_t;

    typedef ankerl::unordered_dense::set<Kmer<k>, Kmer_Hasher<k>> set_t;


    // Constructs working space for workers, supporting capacity of at least
    // `max_sz` vertices. For colored graphs, temporary color-relationship
    // buckets are stored at path-prefix `color_rel_bucket_pref`.
    Subgraphs_Scratch_Space(std::size_t max_sz, const std::string& color_rel_bucket_pref);

    // Returns the appropriate map for a worker.
    map_t& map();

    // Returns the appropriate container of in-process vertices, their lm-tig
    // coordinates and color-hashes, for a worker.
    in_process_arr_t& in_process_arr();

    // Returns the count of color-relationship buckets per worker.
    static constexpr auto color_rel_bucket_c() { return color_rel_bucket_c_; }

    // Returns the appropriate array of buckets for (vertex, source-ID)
    // relationships for a worker.
    color_rel_bucket_arr_t& color_rel_bucket_arr();

    // Returns the hashtable for color-sets.
    Color_Table& color_map();

    // Returns the appropriate container for (vertex, source-ID) relationships
    // for a worker.
    color_rel_arr_t& color_rel_arr();

    // Returns the appropriate container for collated (vertex, source-ID)
    // relationships for a worker.
    color_rel_arr_t& color_rel_collate_arr();

    // Returns the appropriate count map of (vertex, source-ID) relationships
    // for a worker.
    count_map_t& count_map();

    // Returns the appropriate color bit-vector of a worker.
    bit_vector_t& bv();

    // Returns the appropriate hashset for a worker.
    set_t& set();

    // Returns the external-memory color repository.
    Color_Repo& color_repo();

private:

    std::vector<Padded<map_t>> map_;   // Map collection for different workers.
    // TODO: try thread-local allocation for map-space, e.g. from parlay.

    Color_Table M_c;    // Hashtable for color-sets.

    // Collection of containers for in-process vertices: their lm-tig
    // coordinates and color-hashes, for different workers.
    std::vector<Padded<in_process_arr_t>> in_process_arr_;

    static constexpr std::size_t color_rel_bucket_c_ = 32;  // Count of color-relationship buckets per worker.
    static constexpr std::size_t color_rel_buf_sz = 1024 * 1024;    // 1 MB.

    // Collection of array of buckets for (vertex, source-ID) relationships, for
    // different workers.
    std::vector<Padded<color_rel_bucket_arr_t>> color_rel_bucket_arr_;

    // Collection of containers for (vertex, source-ID) relationships, for
    // different workers.
    std::vector<Padded<color_rel_arr_t>> color_rel_arr_;

    // Collection of containers for collated (vertex, source-ID) relationships,
    // for different workers.
    std::vector<Padded<color_rel_arr_t>> color_rel_collate_arr_;

    // Collection of count map of (vertex, source-ID) relationships, for
    // different workers.
    std::vector<Padded<count_map_t>> count_map_;

    // Collection of color bit-vectors of different workers.
    std::vector<Padded<bit_vector_t>> bv_;

    // Hashset collection for different workers.
    std::vector<Padded<set_t>> set_;    // Set collection for different workers.

    // External-memory color repository.
    Color_Repo color_repo_;
};


// =============================================================================
// A subgraph of a de Bruijn graph of `k`-mers. `Colored_` denotes whether the
// vertices have colors.
template <uint16_t k, bool Colored_>
class Subgraph
{
    typedef Walk_Termination termination_t;
    typedef typename Subgraphs_Scratch_Space<k, Colored_>::in_process_arr_t in_process_arr_t;
    typedef typename Subgraphs_Scratch_Space<k, Colored_>::color_rel_t color_rel_t;
    typedef typename Subgraphs_Scratch_Space<k, Colored_>::color_rel_arr_t color_rel_arr_t;

private:

    typedef HT_Router<k, Colored_> ht_router;

    const Super_Kmer_Bucket<Colored_>& B;   // The weak super k-mer bucket inducing this subgraph.

    Subgraphs_Scratch_Space<k, Colored_>& work_space;   // Collection of working space for various data structures, per worker.

    typename Subgraphs_Scratch_Space<k, Colored_>::map_t& M;    // Map to be used for this subgraph.

    uint64_t kmer_count_;   // Number of k-mer instances (copies) in the graph.

    uint64_t edge_c;    // Number of edges in the graph.
    uint64_t label_sz;  // Total number of characters in the literal representations of all the maximal unitigs.
    uint64_t disc_edge_c;   // Number of edges of the discontinuity graph induced from this subgraph.
    uint64_t isolated;  // Count of isolated vertices—not part of any edge.

    Discontinuity_Graph<k, Colored_>& G;    // The discontinuity graph.

    uint64_t mtig_c;    // Number of maximal unitigs in the graph.
    uint64_t trivial_mtig_c;    // Number of trivial maximal unitigs in the graph (i.e. also maximal unitigs in the supergraph).
    uint64_t icc_count_;    // Number of trivial maximal unitigs in the graph that are ICCs.

    uint64_t color_shift_c; // Number of vertices in the graph that either shift color or is the first vertex in an lm-tig.
    uint64_t v_new_col_c;   // Number of vertices in the graph attempting introduction of new colors to the global color-table.
    uint64_t v_old_col_c;   // Number of vertices in the graph with existing colors from the global color-table.
    uint64_t color_rel_c;   // Number of color-relationships (i.e. (k-mer, source) pairs) sorted in color-extraction.

    double t_collect_rels = 0;  // Time taken to collect color-relationships.
    double t_sort = 0;  // Time taken to semi-sort color-relationships.
    double t_collect_sets = 0;  // Time taken to collect color-sets.
    double t_attach = 0;    // Time taken to attach the color-sets to vertices appropriately.


    // TODO: move the following out to a central location.

    typedef uint64_t label_unit_t;

    typedef typename dBG_Contractor<k>::op_buf_t op_buf_t;
    op_buf_t& op_buf;   // Output buffer for trivially maximal unitigs of the underlying dBG.

    // Returns the `idx`'th base of the super k-mer label encoding `super_kmer`
    // that has `word_count` words.
    static base_t get_base(const label_unit_t* super_kmer, std::size_t word_count, std::size_t idx);

    // Extracts the maximal unitig containing the vertex `v_hat`, and
    // `maximal_unitig` is used as the working scratch for the extraction, i.e.
    // to build and store two unitigs connecting to the two sides of `v_hat`.
    // Returns `true` iff the containing maximal unitig has not been outputted
    // earlier. If succeeds, puts the produced lm-tig's bucket-ID to `b` and
    // the lm-tig's index in the bucket to `b_idx`.
    bool extract_maximal_unitig(const Kmer<k>& v_hat, Maximal_Unitig_Scratch<k>& maximal_unitig, std::size_t& b, std::size_t& b_idx);

    // Traverses a unitig starting from the vertex `v_hat`, exiting it through
    // the side `s_v_hat`. `unitig` is used as the scratch space to build the
    // unitig. Returns `true` iff the walk tried to exit the subgraph through a
    // discontinuous side; in which case that vertex is stored in `exit_v`.
    termination_t walk_unitig(const Kmer<k>& v_hat, side_t s_v_hat, Unitig_Scratch<k>& unitig, Directed_Vertex<k>& exit_v);

    // Collects color-relationships of vertices with potentially new colors.
    void collect_color_rels();

    // Semi-sorts the color-relationship array `x` of size `sz` to the array
    // `y`.
    void semi_sort_color_rels(const color_rel_t* x, color_rel_t* y, std::size_t sz);

    // Sorts the color-set (list) `color`.
    void sort_color_set(std::vector<source_id_t>& color);

    // Collates the color-sets of vertices from the collected color-relationship
    // array.
    void collect_color_sets();

    // Attaches extracted colors to vertices.
    void attach_colors_to_vertices();



public:

    // Constructs a subgraph object where the subgraph is induced by the weak
    // super k-mers in the bucket `B`. Updates the discontinuity graph `G` with
    // its edges observed from this subgraph and writes the trivially maximal
    // unitigs to `op_buf`. Uses scratch space for internal data structures
    // from `space`.
    Subgraph(const Super_Kmer_Bucket<Colored_>& B, Discontinuity_Graph<k, Colored_>& d_graph, op_buf_t& op_buf, Subgraphs_Scratch_Space<k, Colored_>& space);

    Subgraph(const Subgraph&) = delete;
    Subgraph(Subgraph&&) = delete;

    // Constructs the subgraph from the provided weak super k-mer bucket into
    // an internal navigable and membership data structure.
    void construct();

    // Constructs the subgraph from the provided weak super k-mer bucket into
    // an internal navigable and membership data structure. Addresses "exact"
    // loop-filtering opposed to `construct`.
    // void construct_loop_filtered();

    // Builds the compacted graph from the original graph.
    void contract();

    // Extracts the new color-sets available from this subgraph.
    void extract_new_colors();

    // Returns the size of the graph.
    std::size_t size() const { return M.size(); }

    // Returns the count of isolated vertices—not part of any edge.
    uint64_t isolated_vertex_count() const { return isolated; }

    // Returns the number of k-mer instances (copies) in the graph.
    uint64_t kmer_count() const { return kmer_count_; }

    // Returns the number of (multi-)edges in the graph.
    uint64_t edge_count() const { return edge_c; }

    // Returns the number of edges of the discontinuity graph produced from this
    // subgraph.
    uint64_t discontinuity_edge_count() const { return disc_edge_c; }

    // Returns the number of maximal unitigs in the graph.
    uint64_t mtig_count() const { return mtig_c; }

    // Returns the number of trivial maximal unitigs in the graph (i.e. also
    // maximal unitigs in the supergraph).
    uint64_t trivial_mtig_count() const { return trivial_mtig_c; }

    // Returns the number of trivial maximal unitigs in the graph that are ICCs.
    uint64_t icc_count() const { return icc_count_; }

    // Returns the number of vertices in the graph that either shift color or
    // is the first vertex in an lm-tig.
    uint64_t color_shift_count() const { return color_shift_c; }

    // Returns the number of vertices in the graph for which color-sets were
    // extracted (this may be larger than the count of unique colors).
    std::size_t color_extraction_count() const { return work_space.set().size(); }

    // Returns the number of vertices in the graph attempting introduction of
    // new colors to the global color-table.
    auto new_colored_vertex() const { return v_new_col_c; }

    // Returns the number of vertices in the graph with existing colors from
    // the global color-table.
    auto old_colored_vertex() const { return v_old_col_c; }

    // Returns the number of color-relationships (i.e. (k-mer, source) pairs)
    // sorted in color-extraction.
    auto color_rel_sorted() const { return color_rel_c; }

    // Returns the time taken to collect color-relationships.
    auto collect_rels_time() const { return t_collect_rels; }

    // Returns the time taken to sort color-relationships.
    auto sort_time() const { return t_sort; }

    // Returns the time taken to collect color-sets.
    auto collect_sets_time() const { return t_collect_sets; }

    // Returns the time taken to attach the color-sets to vertices
    // appropriately.
    auto attach_time() const { return t_attach; }

    // Returns the total number of characters in the literal representations of
    // all the maximal unitigs.
    uint64_t label_size() const { return label_sz; }
};


// Router class wrapping some hashtable methods to help switching map types.
template <uint16_t k, bool Colored_>
class HT_Router
{
    template <uint16_t, bool> friend class Subgraph;
    template <uint16_t, bool> friend class Subgraphs_Scratch_Space;

private:

    template <typename T_ht_> static void flush_updates(T_ht_& HT) { (void)HT; }
    static void flush_updates(Kmer_Hashtable<k, Colored_>& HT) { HT.flush_updates(); }

    template <typename T_ht_> static void add_HT(std::vector<Padded<T_ht_>>& vec, std::size_t sz) { vec.emplace_back(); (void)sz; }
    static void add_HT(std::vector<Padded<Kmer_Hashtable<k, Colored_>>>& vec, std::size_t sz) { vec.emplace_back(sz); }

    template <typename T_ht_> static void update(T_ht_& HT, const Kmer<k>& kmer, base_t front, base_t back, side_t disc_0, side_t disc_1, source_id_t source);
    static void update(Kmer_Hashtable<k, Colored_>& HT, const Kmer<k>& kmer, base_t front, base_t back, side_t disc_0, side_t disc_1);

    template <typename T_iter_> static const Kmer<k>& get_key(const T_iter_& it) { return it->first; }
    static const Kmer<k>& get_key(const typename Kmer_Hashtable<k, Colored_>::Iterator& it) { return it->key; }

    template <typename T_iter_> static State_Config<Colored_>& get_val(const T_iter_& it) { return it->second; }
    static State_Config<Colored_>& get_val(const typename Kmer_Hashtable<k, Colored_>::Iterator& it) { return it->val; }
};


// Type of scenarios how a unitig-walk terminates in the subgraph.
enum class Walk_Termination
{
    null,       // non-existent walk
    branched,   // branched off
    crossed,    // crossed to a different unitig, or looped / cycled back to the same unitig
    dead_ended, // no extension existed
    exitted,    // exitted the subgraph
};


// lm-tig coordinate of a vertex (k-mer).
class LMTig_Coord
{
private:

    uint16_t b_;    // Bucket-ID of the containing lm-tig: the `x` coordinate.
    uint16_t off_;  // Offset of the corresponding k-mer within the containing lm-tig label: the `z` coordinate.
    uint32_t idx_;  // Index of the containing lm-tig within its bucket: the `y` coordinate.


public:

    LMTig_Coord(uint16_t b, uint32_t idx, uint16_t off):
          b_(b)
        , off_(off)
        , idx_(idx)
    {}

    // Returns the bucket-ID of the containing lm-tig: the `x` coordinate.
    uint16_t b() const { return b_; }

    // Returns the index of the containing lm-tig within its bucket: the `y`
    // coordinate.
    uint32_t idx() const { return idx_; }

    // Returns the offset of the corresponding k-mer within the containing lm-
    // tig label: the `z` coordinate.
    uint16_t off() const { return off_; }
};


template <uint16_t k, bool Colored_>
inline base_t Subgraph<k, Colored_>::get_base(const label_unit_t* const super_kmer, const std::size_t word_count, const std::size_t idx)
{
    assert(idx / 32 < word_count);

    const auto word_idx = idx >> 5;
    const auto bit_idx  = (idx & 31) << 1;
    return base_t((super_kmer[(word_count - 1) - word_idx] >> (62 - bit_idx)) & 0b11lu);
}


template <uint16_t k, bool Colored_>
inline bool Subgraph<k, Colored_>::extract_maximal_unitig(const Kmer<k>& v_hat, Maximal_Unitig_Scratch<k>& maximal_unitig, std::size_t& b, std::size_t& b_idx)
{
    constexpr auto back = side_t::back;
    constexpr auto front = side_t::front;
    constexpr auto exitted = termination_t::exitted;

    assert(M.find(v_hat) != M.end());


    maximal_unitig.mark_linear();

    Directed_Vertex<k> v_l, v_r;    // Possible discontinuity ends of the maximal unitig at the left and the right extensions.
    termination_t walk_end_l(termination_t::null), walk_end_r(termination_t::null); // Whether the maximal unitig tried to exit the subgraph through the left and the right extensions.

    walk_end_r = walk_unitig(v_hat, back, maximal_unitig.unitig(back), v_r);
    if(maximal_unitig.unitig(back).is_cycle())
    {
        assert(walk_end_r == termination_t::crossed);
        maximal_unitig.mark_cycle(back);
    }
    else
    {
        walk_end_l = walk_unitig(v_hat, front, maximal_unitig.unitig(front), v_l);
        assert(!maximal_unitig.unitig(front).is_cycle());
    }

    if(walk_end_l == exitted || walk_end_r == exitted)  // The maximal unitig containing `v_hat` spans multiple subgraphs.
    {
        maximal_unitig.finalize_weak();
        std::tie(b, b_idx) = G.add_edge( walk_end_l == exitted ? v_l.canonical() : Discontinuity_Graph<k, Colored_>::phi(), walk_end_l == exitted ? v_l.entrance_side() : side_t::back,
                    walk_end_r == exitted ? v_r.canonical() : Discontinuity_Graph<k, Colored_>::phi(), walk_end_r == exitted ? v_r.entrance_side() : side_t::back,
                    walk_end_l != exitted, walk_end_r != exitted, maximal_unitig);
        disc_edge_c++;
    }
    else    // Extracted a trivial maximal unitig.
    {
        trivial_mtig_c++;
        if(maximal_unitig.is_cycle())
            icc_count_++;

        maximal_unitig.finalize();
        if constexpr(!Colored_)
            maximal_unitig.add_fasta_rec_to_buffer(op_buf);
        else
            std::tie(b, b_idx) = G.add_trivial_mtig(maximal_unitig);

    }

    return true;
}


template <uint16_t k, bool Colored_>
inline typename Subgraph<k, Colored_>::termination_t Subgraph<k, Colored_>::walk_unitig(const Kmer<k>& v_hat, const side_t s_v_hat, Unitig_Scratch<k>& unitig, Directed_Vertex<k>& exit_v)
{
    const auto s_icc_return = inv_side(s_v_hat);    // The side through which to return to `v_hat` if it's contained in an ICC.
    Directed_Vertex<k> v(s_v_hat == side_t::back ? v_hat : v_hat.reverse_complement()); // Current vertex being added to the unitig.
    side_t s_v = s_v_hat;   // The side of the current vertex `v_hat` through which to extend the unitig, i.e. to exit `v`.
    base_t b_ext;   // The nucleobase encoding the edge(s) incident to the side `s_v` of `v`.

    auto it = M.find(v.canonical());
    // auto it = M.find_positive(v.canonical());
    assert(it != M.end());
    State_Config state = ht_router::get_val(it);    // State of `v`.

    if constexpr(!Colored_) unitig.init(v);
    else                    unitig.init(v, state.color_hash());

    while(true)
    {
        ht_router::get_val(it).mark_visited();

        b_ext = state.edge_at(s_v);
        assert(!state.is_discontinuous(s_v) || b_ext == base_t::E); // If a side is discontinuous, it must be empty.
        if(b_ext == base_t::N)  // Reached a branching endpoint.
            return termination_t::branched;

        if(b_ext == base_t::E)
        {
            if(CF_UNLIKELY(!state.is_discontinuous(s_v)))   // Reached a truly empty side.
                return termination_t::dead_ended;

            // Trying to exit the subgraph through a discontinuity vertex.
            exit_v = v;
            return termination_t::exitted;
        }


        if(s_v == side_t::front)
            b_ext = DNA_Utility::complement(b_ext);

        v.roll_forward(b_ext);  // Walk to the next vertex.
        it = M.find(v.canonical());
        // it = M.find_positive(v.canonical());
        assert(it != M.end());
        state = ht_router::get_val(it);

        s_v = v.entrance_side();
        assert(!state.is_empty_side(s_v));
        if(state.is_branching_side(s_v))    // Crossed an endpoint and reached a different unitig.
            return termination_t::crossed;

        if(CF_UNLIKELY(state.is_visited())) // Hit the same unitig.
        {
            // The unitig is an ICC; crossed back to the same unitig.
            if(v.canonical() == v_hat && s_v == s_icc_return)
                unitig.mark_cycle();
            else
                // Otherwise, hit a looping edge—visiting the immediate predecessor vertex.
                // A special case; crossed to the same unitig from a different orientation.
                assert(std::as_const(v)
                        .roll_backward(s_v == side_t::front ? state.edge_at(s_v) : DNA_Utility::complement(state.edge_at(s_v)))
                        .is_same_vertex(v));

            return termination_t::crossed;
        }


        // Still within the unitig.
        bool e;
        if constexpr(!Colored_) e = unitig.extend(v, DNA_Utility::map_char(b_ext));
        else                    e = unitig.extend(v, state.color_hash(), DNA_Utility::map_char(b_ext));
        assert(e); (void)e;

        s_v = opposite_side(s_v);
    }

    return termination_t::null;
}


template <uint16_t k, bool Colored_>
template <typename T_ht_>
inline void HT_Router<k, Colored_>::update(T_ht_& HT, const Kmer<k>& kmer, const base_t front, const base_t back, const side_t disc_0, const side_t disc_1, const source_id_t source)
{
    auto& st = HT[kmer];
    st.update_edges(front, back);

    if(disc_0 != side_t::unspecified)
        st.mark_discontinuous(disc_0);
    if(disc_1 != side_t::unspecified)
        st.mark_discontinuous(disc_1);

    if constexpr(Colored_)
        st.add_source(source);
}


template <uint16_t k, bool Colored_>
inline void HT_Router<k, Colored_>::update(Kmer_Hashtable<k, Colored_>& HT, const Kmer<k>& kmer, const base_t front, const base_t back, const side_t disc_0, const side_t disc_1)
{
    HT.update(kmer, front, back, disc_0, disc_1);
}

}



#endif

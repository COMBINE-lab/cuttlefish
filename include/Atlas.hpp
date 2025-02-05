
#ifndef ATLAS_HPP
#define ATLAS_HPP



#include "Super_Kmer_Chunk.hpp"
#include "Super_Kmer_Bucket.hpp"
#include "Spin_Lock.hpp"
#include "globals.hpp"
#include "utility.hpp"
#include "parlay/parallel.h"

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <memory>


namespace cuttlefish
{

template <bool Colored_>
class Atlas
{
private:

    typedef typename Super_Kmer_Chunk<Colored_>::attribute_t attribute_t;
    typedef typename Super_Kmer_Chunk<Colored_>::label_unit_t label_unit_t;


    static constexpr uint64_t atlas_count_ = 128;   // Number of subgraph atlases.
    static constexpr uint64_t graph_per_atlas_ = 128;    // Number of subgraphs per atlas.
    static constexpr uint64_t graph_count_ = atlas_count_ * graph_per_atlas_;   // Number of subgraphs.


    const std::string path_;    // Directory of the external-memory bucket.

    uint64_t size_; // Number of super k-mers in the atlas. It's not necessarily correct before closing it.

    typedef Super_Kmer_Chunk<Colored_> chunk_t;
    const std::size_t chunk_cap;    // Capacity of the super k-mer chunk of the bucket.
    const std::size_t w_local_chunk_cap;    // Capacity of the worker-local super k-mer chunks.
    std::unique_ptr<chunk_t> chunk; // Super k-mer chunk of the bucket for the worker-local chunks to dump to.
    std::unique_ptr<chunk_t> flush_buf; // Super k-mer chunk acting as buffer between the main chunk and the subgraphs.
    std::vector<Padded<chunk_t>> chunk_w;   // `chunk_w[i]` is the specific super k-mer chunk for worker `i`.

    const std::size_t rec_size; // Size of a super k-mer record in bytes.

    Spin_Lock chunk_lock;   // Lock to the chunk.
    Spin_Lock flush_lock;   // Lock to the flush buffer.

    // Byte-capacity of the chunk of each subgraph in the atlas: 64KB.
    static constexpr std::size_t subgraph_chunk_cap_bytes = 64 * 1024;
    std::vector<Super_Kmer_Bucket<Colored_>> subgraph;  // Subgraphs in the atlas.

    std::vector<uint32_t> src_hist; // Frequency histogram of super k-mer sources currently in the chunk.


    // Empties the local chunk of worker `w_id` to the chunk of the bucket in a
    // thread-safe manner.
    void empty_w_local_chunk(std::size_t w_id);

    // Flushes the super k-mers from the chunk `c` to the appropriate subgraphs.
    void flush_chunk(chunk_t& c);

public:

    // Returns the number of subgraph atlases.
    static constexpr auto atlas_count() { return atlas_count_; }

    // Returns the number of subgraphs per atlas.
    static constexpr auto graph_per_atlas() { return graph_per_atlas_; }

    // Returns the number of subgraphs.
    static constexpr auto graph_count() { return graph_count_; }

    // Returns the atlas-ID of the `g`'th subgraph.
    static auto atlas_ID(const uint64_t g) { return g >> log_2(graph_per_atlas()); }

    // Returns the graph-ID of the `g`'th subgraph within its atlas.
    static auto graph_ID(const uint64_t g) { return g & (graph_per_atlas() - 1); }

    // Constructs a super k-mer atlas for `k`-mers and `l`-minimizers, at
    // external-memory path-prefix `path`. The super chunk buffer of the atlas
    // will have a soft capacity of `chunk_cap` and each worker-local buffer
    // will have a hard capacity of `chunk_cap_per_w`.
    Atlas(uint16_t k, uint16_t l, const std::string& path, std::size_t chunk_cap, std::size_t chunk_cap_per_w);

    Atlas(Atlas&&);

    Atlas(const Atlas&) = delete;
    Atlas& operator=(const Atlas&) = delete;
    Atlas& operator=(Atlas&&) = delete;

    // Returns the number of super k-mers in the atlas. It's not necessarily
    // correct before closing it.
    auto size() const { return size_; }

    // Returns the size of the atlas in bytes. It's not necessarily correct
    // before closing the bucket.
    auto bytes() const { return size() * rec_size; }

    // Adds a super k-mer to the atlas with label `seq` and length `len`. The
    // markers `l_disc` and `r_disc` denote whether the left and the right ends
    // of the (weak) super k-mer are discontinuous or not. The associated super
    // k-mer is to reside in the `g_id`'th subgraph.
    void add(const char* seq, std::size_t len, bool l_disc, bool r_disc, uint16_t g_id);

    // Adds a super k-mer to the atlas with label `seq` and length `len` from
    // source-ID `source`. The markers `l_disc` and `r_disc` denote whether the
    // left and the right ends of the (weak) super k-mer are discontinuous or
    // not. The associated super k-mer is to reside in the `g_id`'th subgraph.
    void add(const char* seq, std::size_t len, source_id_t source, bool l_disc, bool r_disc, uint16_t g_id);

    // Collates the worker-local super k-mers in the bucket per their source-ID
    // and flushes them to the subgraphs in the atlas. The source-IDs are
    // supposed to be in the range `[src_min, src_max]`.
    void flush_collated(source_id_t src_min, source_id_t src_max);

    // Flushes the buffer of the `w`'th worker to the subgraphs in the atlas if
    // it is overflown.
    void flush_worker_if_req(std::size_t w);

    // Closes the atlas—no more content should be added afterwards.
    void close();

    // Returns the super k-mer bucket of the `g`'th subgraph in the atlas.
    auto& bucket(const std::size_t g) { return subgraph[g]; }

    // Returns the resident set size of the space-dominant components of the
    // atlas.
    std::size_t RSS() const;
};


template <>
inline void Atlas<false>::add(const char* const seq, const std::size_t len, const bool l_disc, const bool r_disc, const uint16_t g_id)
{
    const auto w_id = parlay::worker_id();
    auto& c_w = chunk_w[w_id].unwrap(); // Worker-specific chunk.

    c_w.add(seq, len, l_disc, r_disc, g_id);
    if(c_w.full())
        empty_w_local_chunk(w_id);
}


template <>
inline void Atlas<true>::add(const char* const seq, const std::size_t len, const source_id_t source, const bool l_disc, const bool r_disc, const uint16_t g_id)
{
    const auto w_id = parlay::worker_id();
    auto& c_w = chunk_w[w_id].unwrap(); // Worker-specific chunk.

    c_w.add(seq, len, source, l_disc, r_disc, g_id);
    // No flush until collation / flush is invoked explicitly from outside.
}

}



#endif


#ifndef SUPER_KMER_BUCKET_HPP
#define SUPER_KMER_BUCKET_HPP



#include "Super_Kmer_Chunk.hpp"
#include "globals.hpp"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <cassert>


namespace cuttlefish
{

// =============================================================================
// A bucket of super k-mers corresponding to a subgraph of the underlying de
// Bruijn graph. `Colored_` denotes whether the super k-mers in the bucket each
// has an associated source ID.
template <bool Colored_>
class Super_Kmer_Bucket
{
    typedef typename Super_Kmer_Chunk<Colored_>::attribute_t attribute_t;
    typedef typename Super_Kmer_Chunk<Colored_>::label_unit_t label_unit_t;

private:

    const std::string path_;    // Path to the external-memory bucket.
    std::ofstream output;   // Output stream to the external-memory bucket.

    uint64_t size_; // Number of super k-mers in the bucket. It's not necessarily correct before closing the bucket.

    typedef Super_Kmer_Chunk<Colored_> chunk_t;
    std::size_t chunk_cap;  // Capacity (in number of super k-mers) of the chunk of the bucket.
    mutable chunk_t chunk;  // Super k-mer chunk for the bucket.

    std::vector<uint32_t> chunk_sz; // Sizes of the flushed chunks.
    std::vector<std::pair<int32_t, int32_t>> cmp_bytes; // Sizes (in bytes) of the compressed chunks' attributes and labels.

    std::size_t bytes_; // Total number of bytes in the bucket. Not necessarily exact before closing.
    std::size_t compressed_bytes_;  // Total number of bytes in the compressed bucket. Not necessarily exact before closing.


    // Flushes the super k-mer chunk to the external-memory bucket.
    void flush_chunk();

public:

    class Iterator;
    friend class Iterator;

    // Constructs a super k-mer bucket for `k`-mers and `l`-minimizers, at
    // external-memory path `path`. The super k-mer chunk of the bucket has
    // soft capacity of `chunk_cap`.
    Super_Kmer_Bucket(uint16_t k, uint16_t l, const std::string& path, std::size_t chunk_cap);

    Super_Kmer_Bucket(Super_Kmer_Bucket&& rhs) = default;

    Super_Kmer_Bucket(const Super_Kmer_Bucket&) = delete;
    Super_Kmer_Bucket& operator=(const Super_Kmer_Bucket&) = delete;
    Super_Kmer_Bucket&& operator=(Super_Kmer_Bucket&&) = delete;

    // Returns the number of super k-mers in the bucket. It's not necessarily
    // correct before closing the bucket.
    auto size() const { return size_; }

    // Returns the total number of bytes in the bucket. Not necessarily exact
    // before closing.
    auto bytes() const { return bytes_; }

    // Returns the total number of bytes in the compressed bucket. Not
    // necessarily exact before closing.
    auto compressed_bytes() const { return compressed_bytes_; }

    // Issues prefetch request for the end of the chunk.
    void fetch_end() { chunk.fetch_end(); }

    // Adds a super k-mer directly to the chunk with encoding `seq` and
    // attributes `att`.
    void add(const label_unit_t* seq, const attribute_t& att);

    // Closes the bucket—no more content should be added afterwards.
    void close();

    // Removes the bucket.
    void remove();

    // Returns the resident set size of the space-dominant components of the
    // bucket.
    std::size_t RSS() const;
};


// Iterator over super k-mer buckets.
template <bool Colored_>
class Super_Kmer_Bucket<Colored_>::Iterator
{
private:

    const Super_Kmer_Bucket& B; // Bucket to iterate over.
    std::ifstream input;    // Input stream from the external-memory bucket.

    std::size_t idx;    // Current slot-index the iterator is in, i.e. next super k-mer to access.
    std::size_t chunk_start_idx;    // Index into the bucket where the current in-memory chunk starts.
    std::size_t chunk_end_idx;  // Non-inclusive index into the bucket where the current in-memory chunk ends.
    std::size_t chunk_id;   // Sequential-ID of the chunk being processed right now.


    // Reads in the next super k-mer chunk from the bucket and returns the
    // number of super k-mers read.
    std::size_t read_chunk();

public:

    typedef typename Super_Kmer_Chunk<Colored_>::attribute_t attribute_t;
    typedef typename Super_Kmer_Chunk<Colored_>::label_unit_t label_unit_t;

    // Constructs an iterator for the super k-mer bucket `B`.
    Iterator(const Super_Kmer_Bucket& B);

    Iterator(const Iterator&) = delete;
    Iterator& operator=(const Iterator&) = delete;
    Iterator(Iterator&&) = delete;
    Iterator& operator=(Iterator&&) = delete;

    // Return the number of 64-bit words in super k-mer encodings.
    auto super_kmer_word_count() const { return B.chunk.super_kmer_word_count(); }

    // Moves the iterator to the next super k-mer in the bucket. Iff the bucket
    // is not depleted, the associated super k-mer's attribute and label-
    // encoding are put in `att` and `label` respectively, and returns `true`.
    bool next(attribute_t& att, const label_unit_t*& label);
};


template <bool Colored_>
inline void Super_Kmer_Bucket<Colored_>::add(const label_unit_t* const seq, const attribute_t& att)
{
    chunk.add(seq, att);
    size_++;
    if(chunk.full())
        flush_chunk();
}


template <bool Colored_>
inline bool Super_Kmer_Bucket<Colored_>::Iterator::next(attribute_t& att, const label_unit_t*& label)
{
    assert(idx <= B.size());

    if(CF_UNLIKELY(idx == B.size()))
    {
        B.chunk.clear();
        return false;
    }

    if(idx == chunk_end_idx)
    {
        chunk_start_idx = chunk_end_idx;
        chunk_end_idx += read_chunk();
    }

    assert(idx >= chunk_start_idx && idx < chunk_end_idx);
    B.chunk.get_super_kmer(idx - chunk_start_idx, att, label);
    idx++;

    return true;
}

}



#endif

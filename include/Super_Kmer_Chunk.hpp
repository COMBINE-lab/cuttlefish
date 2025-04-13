
#ifndef SUPER_KMER_CHUNK_HPP
#define SUPER_KMER_CHUNK_HPP



#include "Super_Kmer_Attributes.hpp"
#include "Kmer_Utility.hpp"
#include "globals.hpp"
#include "utility.hpp"

#include "lz4.h"
#include <cstdint>
#include <cstddef>
#include <cstring>
#include <utility>
#include <fstream>
#include <iostream>
#include <cassert>


namespace cuttlefish
{

// =============================================================================
// A chunk of weak super k-mers: their attributes and labels. `Colored_` denotes
// whether each super k-mer has an associated source ID.
template <bool Colored_>
class Super_Kmer_Chunk
{
public:

    typedef Super_Kmer_Attributes<Colored_> attribute_t;
    typedef uint64_t label_unit_t;

    class Iterator;

private:

    std::size_t max_sup_kmer_len;   // Maximum length of the (weak) super k-mers.
    std::size_t sup_kmer_word_c;    // Number of 64-bit words in super k-mer encodings.

    std::size_t cap_ = 0;   // Maximum capacity of the chunk in number of super k-mers.
    std::size_t size_ = 0;  // Size of the chunk in number of super k-mers.

    Buffer<attribute_t> att_buf;    // Buffer of attributes of the super k-mers.
    Buffer<label_unit_t> label_buf; // Buffer of concatenated labels of the super k-mers.

    mutable Buffer<uint8_t> cmp_buf;    // Buffer to (de)compress data.


    // Adds the 2-bit encoded form of the label `seq` with length `len` to the
    // chunk.
    void add_encoded_label(const char* seq, std::size_t len);


public:

    // Constructs a placeholder chunk.
    Super_Kmer_Chunk(){}

    // Constructs a super k-mer chunk for `k`-mers and `l`-minimizers, with
    // maximum capacity `cap` in number of super k-mers.
    Super_Kmer_Chunk(uint16_t k, uint16_t l, std::size_t cap);

    Super_Kmer_Chunk(Super_Kmer_Chunk&&);

    Super_Kmer_Chunk& operator=(Super_Kmer_Chunk&&);

    Super_Kmer_Chunk(const Super_Kmer_Chunk&) = delete;
    Super_Kmer_Chunk& operator=(const Super_Kmer_Chunk&) = delete;

    // Returns the number of 64-bit words in super k-mer encodings.
    auto super_kmer_word_count() const { return sup_kmer_word_c; }

    // Returns the number of super k-mers in the chunk.
    auto size() const { return size_; }

    // Returns whether the chunk is empty.
    bool empty() const { return size_ == 0; }

    // Returns the maximum capacity of the chunk in number of super k-mers.
    auto capacity() const {return cap_;  }

    // Returns the free capacity of the chunk in number of super k-mers.
    auto free_capacity() const { return capacity() - size(); }

    // Returns whether the chunk is full or not.
    auto full() const { return size() == capacity(); }

    // Returns the number of units, i.e. 64-bit words, in the label buffer.
    auto label_units() const { return size() * sup_kmer_word_c; }

    // Returns the total number of bytes in the chunk.
    auto bytes() const { return size() * sizeof(attribute_t) + label_units() * sizeof(label_unit_t); }

    // Reserves sufficient space for at least `cap` many super k-mers.
    void reserve(std::size_t cap);

    // Reserves sufficient space for at least `cap` many super k-mers. No
    // guarantees are made for the existing elements.
    void reserve_uninit(std::size_t cap);

    // Resizes the chunk to `n` many super k-mers.
    void resize(std::size_t n);

    // Resizes the chunk to `n` many super k-mers. No guarantees are made for
    // the existing elements.
    void resize_uninit(std::size_t n);

    // Clears the chunk.
    void clear() { size_ = 0; }

    // Frees up the memory used by this chunk.
    void free();

    // Returns the size of a super k-mer record in bytes, that is over `k`-mers
    // and `l`-minimizers.
    static std::size_t record_size(uint16_t k, uint16_t l);

    // Returns the size of a super k-mer record in bytes.
    std::size_t record_size() const { return sizeof(attribute_t) + sup_kmer_word_c * sizeof(label_unit_t); }

    // Serializes the chunk to the stream `os`.
    template <typename T_os_>
    void serialize(T_os_& os) const;

    // Serializes the chunk in a compressed format to the stream `os` and
    // returns the compressed sizes of the attributes and the labels.
    auto serialize_compressed(std::ofstream& os) const -> std::pair<int32_t, int32_t>;

    // Deserializes a chunk from the stream `is` with `sz` super k-mers.
    template <typename T_is_>
    void deserialize(T_is_& is, std::size_t sz);

    // Deserializes a compressed chunk with `sz` super k-mers that has size
    // `cmp_bytes` in the compressed form, from the stream `is`.
    void deserialize_decompressed(std::ifstream& is, std::size_t sz, std::pair<int32_t, int32_t> cmp_bytes);

    // Issues prefetch request for the end of the chunk.
    void fetch_end() const;

    // Adds a super k-mer to the chunk with label `seq` and length `len`. The
    // markers `l_disc` and `r_disc` denote whether the left and the right ends
    // of the (weak) super k-mer are discontinuous or not. The associated super
    // k-mer is to reside in the `g_id`'th subgraph.
    void add(const char* seq, std::size_t len, bool l_disc, bool r_disc, uint16_t g_id);

    // Adds a super k-mer to the chunk with label `seq` and length `len` from
    // source-ID `source`. The markers `l_disc` and `r_disc` denote whether the
    // left and the right ends of the (weak) super k-mer are discontinuous or
    // not. The associated super k-mer is to reside in the `g_id`'th subgraph.
    void add(const char* seq, std::size_t len, source_id_t source, bool l_disc, bool r_disc, const uint16_t g_id);

    // Adds a super k-mer to the chunk with encoding `seq` and attributes
    // `att`.
    void add(const label_unit_t* seq, const attribute_t& att);

    // Appends the chunk `c`'s contents in the indices `[l, r)` to this chunk.
    void append(const Super_Kmer_Chunk& c, std::size_t l, std::size_t r);

    // Appends the chunk `c` to the end of this chunk.
    void append(const Super_Kmer_Chunk& c);

    // Copies `n` super k-mers from the chunk `c`'s index `src_idx` to the
    // index `dest_idx` of this chunk. The indices `[dest_idx, dest_idx + n)`
    // are overwritten.
    void copy(std::size_t dest_idx, const Super_Kmer_Chunk& c, std::size_t src_idx, std::size_t n);

    // Moves `n` super k-mers from index `src_idx` to index `dest_idx`. The
    // indices `[dest_idx, dest_idx + n)` are overwritten.
    void move(std::size_t dest_idx, std::size_t src_idx, std::size_t n);

    // Gets the `idx`'th super k-mer's (in the chunk) attributes to `att` and
    // label to `label`.
    void get_super_kmer(std::size_t idx, attribute_t& att, const label_unit_t*& label);

    // Returns the attribute of the super k-mer at index `i`.
    const attribute_t& att_at(const std::size_t i) const { assert(i < size()); return att_buf[i]; }

    // Returns the attribute of the super k-mer at the front of the chunk.
    const attribute_t& front_att() const { return att_at(0); }

    // Returns the attribute of the super k-mer at the back of the chunk.
    const attribute_t& back_att() const { return att_at(size() - 1); }

    // Returns the location of the label of the super k-mer at index `i`.
    const label_unit_t* label_at(const std::size_t i) const { assert(i < size()); return label_buf.data() + i * sup_kmer_word_c; }

    // Returns an iterator over the super k-mers in the chunk.
    Iterator iterator() const { return Iterator(*this); }

    // Returns the resident set size of the space-dominant components of the
    // chunk.
    std::size_t RSS() const;
};


template <bool Colored_>
class Super_Kmer_Chunk<Colored_>::Iterator
{
private:

    const Super_Kmer_Chunk& chunk;  // Chunk to iterate over.
    std::size_t idx;    // Current slot-index the iterator is in, i.e. next super k-mer to access.
    std::size_t label_off;  // Current label-offset the iterator is in, i.e. next super k-mer's label's starting offset.

public:

    // Constructs an iterator for the super k-mer chunk `chunk`.
    Iterator(const Super_Kmer_Chunk& chunk);

    // Moves the iterator to the next super k-mer in the chunk. Iff the chunk
    // is not depleted, the associated super k-mer's attribute and label-
    // encoding are put in `att` and `label` respectively, and returns `true`.
    bool next(attribute_t& att, const label_unit_t*& label);
};


template <bool Colored_>
inline void Super_Kmer_Chunk<Colored_>::reserve(const std::size_t cap)
{
    if(capacity() >= cap)
        return;

    att_buf.reserve(cap);
    label_buf.reserve(cap * sup_kmer_word_c);
    cap_ = att_buf.capacity();
}


template <bool Colored_>
inline void Super_Kmer_Chunk<Colored_>::reserve_uninit(const std::size_t cap)
{
    if(capacity() >= cap)
        return;

    att_buf.reserve_uninit(cap);
    label_buf.reserve_uninit(cap * sup_kmer_word_c);
    cap_ = att_buf.capacity();
}


template <bool Colored_>
inline void Super_Kmer_Chunk<Colored_>::resize(const std::size_t n)
{
    if(n > capacity())
        reserve(n);

    size_ = n;
}



template <bool Colored_>
inline void Super_Kmer_Chunk<Colored_>::resize_uninit(const std::size_t n)
{
    if(n > capacity())
        reserve_uninit(n);

    size_ = n;
}


template <bool Colored_>
template <typename T_os_>
inline void Super_Kmer_Chunk<Colored_>::serialize(T_os_& os) const
{
    os.write(reinterpret_cast<const char*>(att_buf.data()), size() * sizeof(attribute_t));
    os.write(reinterpret_cast<const char*>(label_buf.data()), label_units() * sizeof(label_unit_t));

    if(!os)
    {
        std::cerr << "Serialization of super k-mer chunk of size " << size() << " failed. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <bool Colored_>
inline auto Super_Kmer_Chunk<Colored_>::serialize_compressed(std::ofstream& os) const -> std::pair<int32_t, int32_t>
{
    const auto max_att_bytes = LZ4_compressBound(size() * sizeof(attribute_t));
    const auto max_label_bytes = LZ4_compressBound(label_units() * sizeof(label_unit_t));
    assert(max_att_bytes > 0 && max_label_bytes > 0);

    cmp_buf.reserve_uninit(max_att_bytes + max_label_bytes);
    auto* const sink = reinterpret_cast<char*>(cmp_buf.data());

    auto* const sink_att = sink;
    const auto att_bytes = LZ4_compress_default(reinterpret_cast<const char*>(att_buf.data()), sink_att, size() * sizeof(attribute_t), cmp_buf.capacity());
    assert(att_bytes > 0);

    auto* const sink_label = sink + att_bytes;
    const auto label_bytes = LZ4_compress_default(reinterpret_cast<const char*>(label_buf.data()), sink_label, label_units() * sizeof(label_unit_t), cmp_buf.capacity() - att_bytes);
    assert(label_bytes > 0);

    os.write(sink, att_bytes + label_bytes);
    if(!os)
    {
        std::cerr << "Serialization of compressed super k-mer chunk of size " << size() << " failed. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    return {att_bytes, label_bytes};
}


template <bool Colored_>
template <typename T_is_>
inline void Super_Kmer_Chunk<Colored_>::deserialize(T_is_& is, const std::size_t sz)
{
    assert(sz <= cap_);

    size_ = sz;

    is.read(reinterpret_cast<char*>(att_buf.data()), size() * sizeof(attribute_t));
    assert(static_cast<std::size_t>(is.gcount()) == size() * sizeof(attribute_t));

    is.read(reinterpret_cast<char*>(label_buf.data()), label_units() * sizeof(label_unit_t));
    assert(static_cast<std::size_t>(is.gcount()) == label_units() * sizeof(label_unit_t));

    if(!is)
    {
        std::cerr << "Deserialization of super k-mer chunk of size " << size() << " failed. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <bool Colored_>
inline void Super_Kmer_Chunk<Colored_>::deserialize_decompressed(std::ifstream& is, const std::size_t sz, const std::pair<int32_t, int32_t> cmp_bytes)
{
    assert(sz <= cap_);
    size_ = sz;

    cmp_buf.reserve_uninit(cmp_bytes.first + cmp_bytes.second);
    auto* const src = reinterpret_cast<char*>(cmp_buf.data());

    is.read(src, cmp_bytes.first + cmp_bytes.second);
    assert(is.gcount() == static_cast<std::streamsize>(cmp_bytes.first + cmp_bytes.second));

    const auto src_att = src;
    const auto att_bytes = LZ4_decompress_safe(src_att, reinterpret_cast<char*>(att_buf.data()), cmp_bytes.first, att_buf.capacity() * sizeof(attribute_t));
    assert(att_bytes >= 0); (void)att_bytes;
    assert(static_cast<std::size_t>(att_bytes) == size() * sizeof(attribute_t));

    const auto src_label = src + cmp_bytes.first;
    const auto label_bytes = LZ4_decompress_safe(src_label, reinterpret_cast<char*>(label_buf.data()), cmp_bytes.second, label_buf.capacity() * sizeof(label_unit_t));
    assert(label_bytes >= 0);   (void)label_bytes;
    assert(static_cast<std::size_t>(label_bytes) == label_units() * sizeof(label_unit_t));
}


template <>
inline void Super_Kmer_Chunk<false>::add(const char* const seq, const std::size_t len, const bool l_disc, const bool r_disc, const uint16_t g_id)
{
    assert(len <= max_sup_kmer_len);
    assert(size() < cap_);

    att_buf[size()] = Super_Kmer_Attributes<false>(len, l_disc, r_disc, g_id);
    add_encoded_label(seq, len);
    size_++;
}


template <>
inline void Super_Kmer_Chunk<true>::add(const char* const seq, const std::size_t len, const source_id_t source, const bool l_disc, const bool r_disc, const uint16_t g_id)
{
    assert(len <= max_sup_kmer_len);

    reserve(size() + 1);

    att_buf[size()] = Super_Kmer_Attributes<true>(len, source, l_disc, r_disc, g_id);
    add_encoded_label(seq, len);
    size_++;
}


template <bool Colored_>
inline void Super_Kmer_Chunk<Colored_>::add(const label_unit_t* const seq, const attribute_t& att)
{
    if constexpr(!Colored_)
        assert(size() < cap_);
    else
        reserve(size() + 1);

    att_buf[size()] = att;
    std::memcpy(label_buf.data() + label_units(), seq, sup_kmer_word_c * sizeof(label_unit_t));
    size_++;
}


template <bool Colored_>
inline void Super_Kmer_Chunk<Colored_>::add_encoded_label(const char* const seq, const std::size_t len)
{
    const auto label_off = label_units();   // Offset into the packed-encoding concatenation where to put the label.
    int64_t word_idx = sup_kmer_word_c - 1; // Index of the current word being encoded from the label; the encoding is MSB-boundary aligned for now.
    // TODO: use vectorized encoding.
    for(std::size_t b_idx = 0; b_idx < len; b_idx += 32)
    {
        label_buf[label_off + word_idx] = Kmer_Utility::encode_checked<32>(seq + b_idx);
        word_idx--;
    }
}


template <bool Colored_>
inline void Super_Kmer_Chunk<Colored_>::append(const Super_Kmer_Chunk& c, const std::size_t l, const std::size_t r)
{
    assert(r >= l);

    const auto n = r - l;
    reserve(size() + n);
    assert(size() + n <= cap_);

    std::memcpy(att_buf.data() + size(), c.att_buf.data() + l, n * sizeof(attribute_t));
    std::memcpy(label_buf.data() + label_units(), c.label_buf.data() + l * sup_kmer_word_c, n * sup_kmer_word_c * sizeof(label_unit_t));

    size_ += n;
}


template <bool Colored_>
inline void Super_Kmer_Chunk<Colored_>::append(const Super_Kmer_Chunk& c)
{
    append(c, 0, c.size());
}


template <bool Colored_>
inline void Super_Kmer_Chunk<Colored_>::copy(const std::size_t dest_idx, const Super_Kmer_Chunk& c, const std::size_t src_idx, const std::size_t n)
{
    assert(dest_idx + n <= size());

    std::memcpy(att_buf.data() + dest_idx, c.att_buf.data() + src_idx, n * sizeof(attribute_t));
    std::memcpy(label_buf.data() + dest_idx * sup_kmer_word_c, c.label_buf.data() + src_idx * sup_kmer_word_c, n * sup_kmer_word_c * sizeof(label_unit_t));
}


template <bool Colored_>
inline void Super_Kmer_Chunk<Colored_>::move(const std::size_t dest_idx, const std::size_t src_idx, const std::size_t n)
{
    assert(dest_idx + n <= size());

    std::memmove(att_buf.data() + dest_idx, att_buf.data() + src_idx, n * sizeof(attribute_t));
    std::memmove(label_buf.data() + dest_idx * sup_kmer_word_c, label_buf.data() + src_idx * sup_kmer_word_c, n * sup_kmer_word_c * sizeof(label_unit_t));
}


template <bool Colored_>
inline void Super_Kmer_Chunk<Colored_>::get_super_kmer(std::size_t idx, attribute_t& att, const label_unit_t*& label)
{
    assert(idx < size());

    att = att_buf[idx];
    label = label_buf.data() + idx * sup_kmer_word_c;
}


template <bool Colored_>
inline Super_Kmer_Chunk<Colored_>::Iterator::Iterator(const Super_Kmer_Chunk& chunk):
      chunk(chunk)
    , idx(0)
    , label_off(0)
{}


template <bool Colored_>
inline bool Super_Kmer_Chunk<Colored_>::Iterator::next(attribute_t& att, const label_unit_t*& label)
{
    assert(idx <= chunk.size());

    if(CF_UNLIKELY(idx == chunk.size()))
        return false;

    att = chunk.att_buf[idx];
    label = chunk.label_buf.data() + label_off;

    idx++;
    label_off += chunk.super_kmer_word_count();

    return true;
}

}



#endif

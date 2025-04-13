
#ifndef COLOR_ENCODING_HPP
#define COLOR_ENCODING_HPP



#include <cstdint>
#include <cassert>


namespace cuttlefish
{

// =============================================================================
// Coordinate of a color in the actual color-collection.
class Color_Coordinate
{
    typedef uint64_t pack_t;

private:

    // Flag to denote whether the corresponding color is in the process of
    // extraction or not.
    static constexpr pack_t in_process = (pack_t(1) << (sizeof(pack_t) * 8 - 1));

    static constexpr uint32_t idx_pos = 8;  // Position of the index (in worker-local bucket) of a color-set.

    static constexpr pack_t w_limit = (pack_t(1) << idx_pos);   // Maximum worker count: 2^8.
    static constexpr pack_t idx_limit  = (pack_t(1) << 32); // Maximum worker-local bucket size: 2^32.

    pack_t bit_pack;    // Packed representation of the color-coordinate.

public:

    // Constructs an empty coordinate.
    Color_Coordinate(): bit_pack(0)
    {}

    // Constructs an empty coordinate marked as being processed by worker-ID
    // `w`.
    explicit Color_Coordinate(const pack_t w_id): bit_pack(w_id | in_process) { assert(w_id < w_limit); }

    // Constructs the coordinate `(w_id, idx)`.
    Color_Coordinate(const pack_t w_id, const pack_t idx): bit_pack(w_id | (pack_t(idx) << idx_pos)) { assert(w_id < w_limit); assert(idx < idx_limit); }

    // Returns whether the corresponding color is in the process of extraction
    // or not.
    bool is_in_process() const { return bit_pack & in_process; }

    // Returns the worker-ID that marked this coordinate as processing.
    pack_t processing_worker() const { assert(is_in_process()); return bit_pack & (~in_process); }

    // Returns 40-bit packing of the color-coordinate.
    auto as_u40() const { assert(bit_pack < (pack_t(1) << 40)); return bit_pack; }

    // (De)serializes the coordinate from / to the `cereal` archive `archive`.
    template <typename T_archive_> void serialize(T_archive_& archive) { archive(bit_pack); }
};


// =============================================================================
// Mapping between a vertex (in a given unitig bucket) and its color.
class Vertex_Color_Mapping
{
private:

    // TODO: consider packing `idx` and `off`?
    uint32_t idx_;  // Index of the vertex's containing unitig in its bucket.
    uint16_t off_;  // Offset of the vertex in the unitig label.
    Color_Coordinate c_;    // Coordinate of the vertex's color in the color-repository.

public:

    // For some given unitig bucket, constructs a vertex-color mapping between
    // the vertex at offset `off` in the unitig at index `idx` in the bucket
    // and the color-coordinate `c`.
    Vertex_Color_Mapping(const uint32_t idx, const uint16_t off, const Color_Coordinate& c):
        idx_(idx), off_(off), c_(c)
    {}

    // Returns the index of the vertex's containing unitig in its bucket.
    auto idx() const { return idx_; }

    // Returns the offset of the vertex in the unitig label.
    auto off() const { return off_; }

    // Returns the coordinate of the vertex's color in the color-repository.
    const auto& c() const { return c_; }

    // Returns `true` iff this vertex precedes the associated vertex of `rhs`
    // in their bucket.
    bool operator<(const Vertex_Color_Mapping& rhs) const { return idx_ != rhs.idx_ ? (idx_ < rhs.idx_) : (off_ < rhs.off_); }

    // (De)serializes the mapping from / to the `cereal` archive `archive`.
    template <typename T_archive_> void serialize(T_archive_& archive) { archive(idx_, off_, c_); }
};


// ============================================================================
// Encoding of a color in a unitig: the offset in the unitig where the color
// is, and the color's coordinate in the global color-repository.
class Unitig_Color
{
private:

    uint64_t bit_pack;  // Encoding of the offset and the color.

public:

    // Constructs a color-encoding for a unitig at its offset `off` and color-
    // coordinate `c`.
    Unitig_Color(const std::size_t off, const Color_Coordinate c):
          bit_pack((c.as_u40() << 24) | off)
    {
        assert(off <= 0xFF'FF'FF);
    }

    // Returns the offset of the color in the unitig.
    uint32_t off() const { return bit_pack & 0xFF'FF'FF; }

    // Returns the coordinate of the color in the global color-repository.
    uint64_t c() const { return bit_pack >> 24; }

    // Sets the offset of thee color in the unitig to `s`.
    void set_off(const uint32_t o) { assert(o <= 0xFF'FF'FF); bit_pack = (bit_pack & ~0xFF'FF'FFlu) | o; }

    // Returns the 64-bit representation of the unitig-color.
    uint64_t to_u64() const { return bit_pack; }
};

}



#endif

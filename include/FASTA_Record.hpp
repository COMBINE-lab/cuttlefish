
#ifndef FASTA_RECORD_HPP
#define FASTA_RECORD_HPP



#include "Color_Encoding.hpp"
#include "fmt/format.h"

#include <cstdint>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>
#include <cassert>


// =============================================================================
// A class wrapping a basic FASTA record. The class is specifically designed for
// writing purposes of output maximal unitigs in the FASTA format.
class FASTA_Record
{
private:

    const fmt::format_int id_;  // Identifier for the FASTA sequence.
    const std::string_view seq_;    // The FASTA sequence.
    const std::string_view seq_add_;    // Additional FASTA sequence (in case the original sequence `*seq` is broken into two parts).
    const std::size_t offset_;  // Offset position into the sequence `seq_`—data before this index will be skipped in the record.
    const std::size_t offset_add_;  // Offset position into the additional sequence `seq_add`—data before this index will be skipped in the record.

    typedef std::vector<cuttlefish::Unitig_Color> color_list_t;
    const color_list_t* const color_;   // Color-encodings associated to the FASTA sequence.


    // Constructs a FASTA header with identifier `id`, along with the sequences
    // `seq` and `seq_add` (onward their indices `offset` and `offset_add`,
    // respectively) and the color-list `color`.
    FASTA_Record(uint64_t id, const std::string_view& seq, const std::string_view& seq_add, std::size_t offset = 0, std::size_t offset_add = 0, const color_list_t* color = nullptr);


public:

    // Constructs a FASTA header with identifier `id` and the sequence `seq`.
    // Only a constant reference to the sequence is captured, so the record's
    // correctness holds as long as the referred sequence itself remains
    // unaltered.
    FASTA_Record(uint64_t id, const std::string& str);

    // Constructs a FASTA header with identifier `id` and the sequence `seq`.
    // Only a constant reference to the sequence is captured, so the record's
    // correctness holds as long as the referred sequence itself remains
    // unaltered.
    FASTA_Record(uint64_t id, const std::string_view& str);

    // Constructs a FASTA header with identifier `id`, the sequence `seq`, and
    // color-list `color`. Only a constant reference to the sequence is
    // captured, so the record's correctness holds as long as the referred
    // sequence itself remains unaltered.
    FASTA_Record(uint64_t id, const std::string_view& str, const color_list_t& color);

    // Constructs a FASTA header with identifier `id`, along with the sequences
    // `seq` and `seq_add` (onward their indices `offset` and `offset_add`,
    // respectively). Only constant references to the sequences are captured,
    // so the record's correctness holds as long as the referred sequences themselves
    // remains unaltered.
    FASTA_Record(uint64_t id, const std::string& seq, const std::string& seq_add, std::size_t offset, std::size_t offset_add);

    // Returns the length of the header line of the record.
    std::size_t header_size() const { return id_.size() + 1; /* One additional byte for `>`. */ }

    // Returns the length of the sequence of the record.
    std::size_t seq_size() const { return (seq_.size() - offset_) + (!seq_add_.empty() ? (seq_add_.size() - offset_add_) : 0); }

    // Returns the size of the color-list.
    std::size_t color_list_size() const { assert(color_ != nullptr); return color_->size(); }

    // Appends the header line to the `buffer`.
    void append_header(std::string& buffer) const;

    // Appends the FASTA sequence to the `buffer`.
    void append_seq(std::string& buffer) const;

    // Appends the FASTA sequence to `buffer` in a rotated form — the
    // added sequence is supposed to be a cycle in a de Bruijn graph `G(·, k)`,
    // and it is right rotated so that the character at index `pivot` is at
    // index 0 finally.
    template <uint16_t k>
    void append_rotated_cycle(std::string& buffer, std::size_t pivot) const;

    // Appends the color-list to the `buf`.
    void append_color_list(std::string& buf) const;
};


inline FASTA_Record::FASTA_Record(const uint64_t id, const std::string& seq):
    FASTA_Record(id, std::string_view(seq), std::string_view())
{}


inline FASTA_Record::FASTA_Record(const uint64_t id, const std::string_view& seq):
    FASTA_Record(id, seq, std::string_view())
{}


inline FASTA_Record::FASTA_Record(uint64_t id, const std::string_view& str, const color_list_t& color):
    FASTA_Record(id, str, std::string_view(), 0, 0, &color)
{}


inline FASTA_Record::FASTA_Record(const uint64_t id, const std::string& seq, const std::string& seq_add, const std::size_t offset, const std::size_t offset_add):
    FASTA_Record(id, std::string_view(seq), std::string_view(seq_add), offset, offset_add)
{}


inline FASTA_Record::FASTA_Record(const uint64_t id, const std::string_view& seq, const std::string_view& seq_add, const std::size_t offset, const std::size_t offset_add, const color_list_t* const color):
    id_(id),
    seq_(seq),
    seq_add_(seq_add),
    offset_(offset),
    offset_add_(offset_add),
    color_(color)
{}


inline void FASTA_Record::append_header(std::string& buffer) const
{
    buffer.push_back('>');

    buffer.insert(buffer.end(), id_.data(), id_.data() + id_.size());
}


inline void FASTA_Record::append_seq(std::string& buffer) const
{
    buffer.append(seq_.cbegin() + offset_, seq_.cend());
    if(!seq_add_.empty())
        buffer.append(seq_add_.cbegin() + offset_add_, seq_add_.cend());
}


template <uint16_t k>
inline void FASTA_Record::append_rotated_cycle(std::string& buffer, const std::size_t pivot) const
{
    buffer.append(seq_.cbegin() + pivot, seq_.cend());
    buffer.append(seq_.cbegin() + k - 1, seq_.cbegin() + k - 1 + pivot);
}


inline void FASTA_Record::append_color_list(std::string& buf) const
{
    assert(color_ != nullptr);
    for(const auto& c : *color_)
    {
        const auto v = fmt::format_int(c.to_u64());
        buf.push_back(' ');
        buf.append(v.data(), v.size());
    }
}



#endif

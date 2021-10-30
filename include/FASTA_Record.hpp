
#ifndef FASTA_RECORD_HPP
#define FASTA_RECORD_HPP



#include "fmt/format.h"

#include <cstdint>


// =============================================================================
// A class wrapping a basic FASTA record: the sequence of type `T_seq_` and its
// header/identifier of type `T_id`. The class is specifically designed for
// writing purposed of output maximal unitigs in the FASTA format.
template <typename T_seq_, typename T_id_ = fmt::format_int>
class FASTA_Record
{
private:

    const T_id_ id_; // Identifier for the FASTA sequence.
    const T_seq_* const seq_;   // Pointer to the FASTA sequence.
    const T_seq_* const seq_add_;   // Additional FASTA sequence (in case the original sequence `*seq` is broken into two parts).
    const std::size_t offset_;  // Offset position into the sequence `seq_`—data before this index will be skipped in the record.
    const std::size_t offset_add_;  // Offset position into the additional sequence `seq_add`—data before this index will be skipped in the record.


    // Constructs a FASTA header with identifier `id`, along with the sequences
    // `seq` and `seq_add` (onward their indices `offset` and `offset_add`,
    // respectively). Only constant references to the sequences are captured,
    // so the record's correctness holds as long as the referred sequences themselves
    // remains unaltered.
    FASTA_Record(uint64_t id, const T_seq_* seq, const T_seq_* seq_add, std::size_t offset = 0, std::size_t offset_add = 0);


public:

    // Constructs a FASTA header with identifier `id` and the sequence `seq`
    // (onward its index `offset`). Only a constant reference to the sequence
    // is captured, so the record's correctness holds as long as the referred
    // sequence itself remains unaltered.
    FASTA_Record(uint64_t id, const T_seq_& str, std::size_t offset = 0);

    // Constructs a FASTA header with identifier `id`, along with the sequences
    // `seq` and `seq_add` (onward their indices `offset` and `offset_add`,
    // respectively). Only constant references to the sequences are captured,
    // so the record's correctness holds as long as the referred sequences themselves
    // remains unaltered.
    FASTA_Record(uint64_t id, const T_seq_& seq, const T_seq_& seq_add, std::size_t offset = 0, std::size_t offset_add = 0);

    // Returns the length of the header line of the record.
    std::size_t header_size() const;

    // Returns the length of the sequence of the record.
    std::size_t seq_size() const;

    // Appends the header line to the vector `buffer`.
    void append_header(std::vector<char>& buffer) const;

    // Appends the FASTA sequence to the vector `buffer`.
    void append_seq(std::vector<char>& buffer) const;

    // Appends the FASTA sequence to the vector `buffer` in a rotated form — the
    // added sequence is supposed to be a cycle in a de Bruijn graph `G(·, k)`,
    // and it is right rotated so that the character at index `pivot` is at
    // index 0 finally.
    template <uint16_t k>
    void append_rotated_cycle(std::vector<char>& buffer, std::size_t pivot) const;
};


template <typename T_seq_, typename T_id_>
inline FASTA_Record<T_seq_, T_id_>::FASTA_Record(const uint64_t id, const T_seq_& seq, const std::size_t offset):
    FASTA_Record(id, &seq, nullptr, offset)
{}


template <typename T_seq_, typename T_id_>
inline FASTA_Record<T_seq_, T_id_>::FASTA_Record(const uint64_t id, const T_seq_& seq, const T_seq_& seq_add, const std::size_t offset, const std::size_t offset_add):
    FASTA_Record(id, &seq, &seq_add, offset, offset_add)
{}


template <typename T_seq_, typename T_id_>
inline FASTA_Record<T_seq_, T_id_>::FASTA_Record(const uint64_t id, const T_seq_* const seq, const T_seq_* const seq_add, const std::size_t offset, const std::size_t offset_add):
    id_(id),
    seq_(seq),
    seq_add_(seq_add),
    offset_(offset),
    offset_add_(offset_add)
{}


template <typename T_seq_, typename T_id_>
inline std::size_t FASTA_Record<T_seq_, T_id_>::header_size() const
{
    return  id_.size() + static_cast<std::size_t>(1U);  // One additional byte for `>`.
}


template <typename T_seq_, typename T_id_>
inline std::size_t FASTA_Record<T_seq_, T_id_>::seq_size() const
{
    return (seq_->size() - offset_) + (seq_add_ != nullptr ? (seq_add_->size() - offset_add_) : 0);
}


template <typename T_seq_, typename T_id_>
inline void FASTA_Record<T_seq_, T_id_>::append_header(std::vector<char>& buffer) const
{
    buffer.emplace_back('>');

    buffer.insert(buffer.end(), id_.data(), id_.data() + id_.size());
}


template <typename T_seq_, typename T_id_>
inline void FASTA_Record<T_seq_, T_id_>::append_seq(std::vector<char>& buffer) const
{
    // `std::memcpy` at the end of `buffer` does not update the size of the vector `buffer`.
    buffer.insert(buffer.end(), seq_->begin() + offset_, seq_->end());
    if(seq_add_ != nullptr)
        buffer.insert(buffer.end(), seq_add_->begin() + offset_add_, seq_add_->end());
}


template <typename T_seq_, typename T_id_>
template <uint16_t k>
inline void FASTA_Record<T_seq_, T_id_>::append_rotated_cycle(std::vector<char>& buffer, const std::size_t pivot) const
{
    buffer.insert(buffer.end(), seq_->begin() + pivot, seq_->end());
    buffer.insert(buffer.end(), seq_->begin() + k - 1, seq_->begin() + k - 1 + pivot);
}



#endif

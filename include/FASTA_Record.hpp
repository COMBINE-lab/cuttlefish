
#ifndef FASTA_RECORD_HPP
#define FASTA_RECORD_HPP



#include "fmt/format.h"

// =============================================================================
// A class wrapping a basic FASTA record: the sequence of type `T_seq_` and its
// header/identifier of type `T_id`. The class is specifically designed for
// writing purposed of output maximal unitigs in the FASTA format.
template <typename T_seq_, typename T_id_ = fmt::format_int>
class FASTA_Record
{
private:

    const T_id_ id_; // Identifier for the FASTA sequence.
    const T_seq_& seq_; // The FASTA sequence.


public:

    // Constructs a FASTA header with identifier `id` and the sequence `seq`.
    // Only a constant reference to the sequence is captured, so the record's
    // correctness holds as long as the referred sequence itself remains unaltered.
    FASTA_Record(uint64_t id, const T_seq_& str);

    // Returns the length of the header line of the record.
    std::size_t header_size() const;

    // Returns the length of the sequence of the record.
    std::size_t seq_size() const;

    // Appends the header line to the vector `buffer`.
    void append_header(std::vector<char>& buffer) const;

    // Appends the FASTA sequence to the vector `buffer`.
    void append_seq(std::vector<char>& buffer) const;
};


template <typename T_seq_, typename T_id_>
inline FASTA_Record<T_seq_, T_id_>::FASTA_Record(const uint64_t id, const T_seq_& seq):
    id_(id),
    seq_(seq)
{}


template <typename T_seq_, typename T_id_>
inline std::size_t FASTA_Record<T_seq_, T_id_>::header_size() const
{
    return  id_.size() + static_cast<std::size_t>(1U);  // One additional byte for `>`.
}


template <typename T_seq_, typename T_id_>
inline std::size_t FASTA_Record<T_seq_, T_id_>::seq_size() const
{
    return seq_.size();
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
    buffer.insert(buffer.end(), seq_.begin(), seq_.end());
}



#endif

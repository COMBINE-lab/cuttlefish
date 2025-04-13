
// TODO: Remove

#ifndef MINIMIZER_INSTANCE_HPP
#define MINIMIZER_INSTANCE_HPP



#include "Min_Heap.hpp"
#include "globals.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <utility>


// =============================================================================


// A `(minimizer, offset)` pair, where `minimizer` is the l-minimizer (for a
// pre-defined `l`) of some k-mer `x` in some sequence `seq`, where `minimizer`
// is at the index `offset` in `seq`.
class Minimizer_Instance
{
private:

    typedef cuttlefish::minimizer_t minimizer_t;

    minimizer_t minimizer_; // The minimizer.
    std::size_t offset_;    // Offset of the minimizer in the underlying sequence.

public:

    // Constructs an empty minimizer instance.
    Minimizer_Instance()
    {}

    // Constructs an instance of the minimizer `minimizer`, situated at the
    // offset `offset` some sequence.
    Minimizer_Instance(const minimizer_t minimizer, const std::size_t offset):
        minimizer_(minimizer),
        offset_(offset)
    {}

    // Returns the minimizer.
    cuttlefish::minimizer_t minimizer() const { return minimizer_; }

    // Returns the offset.
    std::size_t offset() const { return offset_; }

    // Shifts (to the right) the offset of the minimizer by `offset_shit`.
    void shift(const std::size_t offset_shift)  { offset_ += offset_shift; }

    // Returns `true` iff this instance has a lesser minimizer than `rhs`; and
    // in the case of same minimizers, returns `true` iff it has a lower offset.
    bool operator<(const Minimizer_Instance& rhs) const
    {
        return minimizer_ != rhs.minimizer_ ? (minimizer_ < rhs.minimizer_) : (offset_ < rhs.offset_);
    }
};



#endif

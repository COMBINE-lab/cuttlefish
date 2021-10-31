
#ifndef UNITIG_SCRATCH_HPP
#define UNITIG_SCRATCH_HPP



#include "Directed_Vertex.hpp"
#include "dBG_Utilities.hpp"

#include <cstdint>
#include <vector>


// =============================================================================
// A class to keep scratch data (i.e. working space) for unitigs.
template <uint16_t k>
class Unitig_Scratch
{
private:

    // 100K (soft limit) unitig vertices can be retained in memory, at most, before reallocations.
    static constexpr std::size_t BUFF_SZ = 100 * 1024UL;

    Directed_Vertex<k> endpoint_;   // The current end of the unitig through which farther extensions can be done.
                                    // (The side for the extension is to be handled by the client code, although can
                                    // also be inferred from the "directed" vertex.)
    std::vector<char> label_;       // Literal label of the unitig.
    std::vector<uint64_t> hash_;    // Hashes of the constituent vertices of the unitig.


    // Clears the scratch data.
    void clear();

public:

    // Constructs an empty unitig scratch.
    Unitig_Scratch();

    // Initializes the unitig scratch with the vertex `v`.
    void init(const Directed_Vertex<k>& v);

    // Extends the unitig scratch with the vertex `v`, and its literal form
    // with the symbol `b`.
    void extend(const Directed_Vertex<k>& v, char b);

    // Reverse complements the label sequence (literal form) of the unitig.
    void rev_compl_label();

    // Returns the literal label of the unitig.
    const std::vector<char>& label() const;

    // Returns the hash collection of the unitig vertices.
    const std::vector<uint64_t>& hash() const;

    // Returns the current extension-end vertex of the unitig.
    const Directed_Vertex<k>& endpoint() const;

    // Returns the count of vertices in this unitig.
    std::size_t size() const;
};


template <uint16_t k>
inline void Unitig_Scratch<k>::clear()
{
    label_.clear();
    hash_.clear();
}


template <uint16_t k>
inline void Unitig_Scratch<k>::init(const Directed_Vertex<k>& v)
{
    clear();

    endpoint_ = v;
    endpoint_.kmer().get_label(label_);
    hash_.emplace_back(endpoint_.hash());
}


template <uint16_t k>
inline void Unitig_Scratch<k>::extend(const Directed_Vertex<k>& v, const char b)
{
    endpoint_ = v;

    label_.emplace_back(b);
    hash_.emplace_back(endpoint_.hash());
}


template <uint16_t k>
inline void Unitig_Scratch<k>::rev_compl_label()
{
    cuttlefish::reverse_complement(label_);
}


template <uint16_t k>
inline const std::vector<char>& Unitig_Scratch<k>::label() const
{
    return label_;
}


template <uint16_t k>
inline const std::vector<uint64_t>& Unitig_Scratch<k>::hash() const
{
    return hash_;
}


template <uint16_t k>
inline const Directed_Vertex<k>& Unitig_Scratch<k>::endpoint() const
{
    return endpoint_;
}


template <uint16_t k>
inline std::size_t Unitig_Scratch<k>::size() const
{
    return hash_.size();
}



#endif

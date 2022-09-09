
#ifndef MAXIMAL_UNITIG_SCRATCH_HPP
#define MAXIMAL_UNITIG_SCRATCH_HPP



#include "Unitig_Scratch.hpp"
#include "FASTA_Record.hpp"
#include "Character_Buffer.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>


// =============================================================================
// A class to keep scratch data for building maximal unitigs from two of its
// constituent unitigs that cover it and overlap at a meeting-point vertex.
// That is, the maximal unitig is split into two unitigs `u_b` and `u_f`, at
// some vertex `v`—`u_b` and `u_f` are connected to the front and to the back
// of `v`, respectively. The unitigs are built such that the paths start from
// `v`. Thus, the maximal unitig in literal form is `\bar(u_f) \glue_k u_b`
// (or its reverse complement).
template <uint16_t k>
class Maximal_Unitig_Scratch
{
private:

    Unitig_Scratch<k> unitig_back;  // The unitig `u_b` (see note above class body).
    Unitig_Scratch<k> unitig_front; // The unitig `u_f` (see note above class body).

    uint64_t id_;   // The unique ID of the maximal unitig.

    Unitig_Scratch<k>* cycle;   // Pointer to either `u_b` or `u_f`, whichever contains the maximal unitig in case of it's a cycle.


    // Returns whether the maximal unitig is in canonical form.
    bool is_canonical() const;


public:

    // Constructs an empty scratch space for the unitig.
    Maximal_Unitig_Scratch();

    // Returns the unitig scratch `u_b` or `u_f`, based on `s` (see note above
    // class body).
    Unitig_Scratch<k>& unitig(const cuttlefish::side_t s);

    // Returns the unique ID of the maximal unitig.
    uint64_t id() const;

    // Returns whether the maximal unitig is linear, i.e. it is a linear path
    // and not a Detached Chordless Cycle (DCC) in the underlying graph.
    bool is_linear() const;

    // Returns the hashes of the vertices of the unitig at side `s`.
    const std::vector<uint64_t>& unitig_hash(cuttlefish::side_t s) const;

    // Returns the hashes of the vertices in the maximal unitig in case it's
    // a DCC.
    const std::vector<uint64_t>& cycle_hash() const;

    // Returns the count of vertices in the maximal unitig.
    std::size_t size() const;

    // Returns the signature vertex of the maximal unitig, which is the first
    // vertex in the canonical form of the unitig.
    const Directed_Vertex<k>& sign_vertex() const;

    // Marks the maximal unitig as linear, i.e not a DCC.
    void mark_linear();

    // Marks the maximal unitig as a DCC, and signals that the cycle has been
    // extracted in the unitig scratch at side `s`.
    void mark_cycle(cuttlefish::side_t s);

    // Signals the scratch that the unitig pieces `u_b` and `u_f` are in their
    // final forms and will not be modified anymore. So it restructures the
    // maximal unitig so as to put its label in canonical form and sets its
    // unique ID.
    void finalize();

    // Returns `true` iff the maximal unitig has been marked as a cycle.
    bool is_cycle() const;

    // Returns a FASTA record of the maximal unitig (in canonical form).
    // Applicable when the maximal unitig is linear.
    const FASTA_Record<std::vector<char>> fasta_rec() const;

    // Adds a corresponding FASTA record for the maximal unitig into `buffer`.
    template <std::size_t CAPACITY, typename T_sink_> void add_fasta_rec_to_buffer(Character_Buffer<CAPACITY, T_sink_>& buffer) const;
};


template <uint16_t k>
inline Unitig_Scratch<k>& Maximal_Unitig_Scratch<k>::unitig(const cuttlefish::side_t s)
{
    return s == cuttlefish::side_t::back ? unitig_back : unitig_front;
}


template <uint16_t k>
inline bool Maximal_Unitig_Scratch<k>::is_canonical() const
{
    return unitig_front.endpoint().kmer_bar() < unitig_back.endpoint().kmer_bar();
}


template <uint16_t k>
inline bool Maximal_Unitig_Scratch<k>::is_linear() const
{
    return cycle == nullptr;
}


template <uint16_t k>
inline uint64_t Maximal_Unitig_Scratch<k>::id() const
{
    return id_;
}


template <uint16_t k>
inline void Maximal_Unitig_Scratch<k>::mark_linear()
{
    cycle = nullptr;
}


template <uint16_t k>
inline const std::vector<uint64_t>& Maximal_Unitig_Scratch<k>::unitig_hash(const cuttlefish::side_t s) const
{
    return (s == cuttlefish::side_t::back ? unitig_back.hash() : unitig_front.hash());
}


template <uint16_t k>
inline const std::vector<uint64_t>& Maximal_Unitig_Scratch<k>::cycle_hash() const
{
    return cycle->hash();
}


template <uint16_t k>
inline std::size_t Maximal_Unitig_Scratch<k>::size() const
{
    return is_linear() ?    (unitig_back.size() + unitig_front.size() - 1) :
                            cycle->size();
}


template <uint16_t k>
inline const Directed_Vertex<k>& Maximal_Unitig_Scratch<k>::sign_vertex() const
{
    return is_linear() ?   (is_canonical() ? unitig_front.endpoint() : unitig_back.endpoint()) :
                            cycle->min_vertex();
}


template <uint16_t k>
inline void Maximal_Unitig_Scratch<k>::mark_cycle(const cuttlefish::side_t s)
{
    cycle = &(s == cuttlefish::side_t::back ? unitig_back : unitig_front);
}


template <uint16_t k>
inline void Maximal_Unitig_Scratch<k>::finalize()
{
    if(is_linear())
    {
        if(is_canonical())
            id_ = unitig_front.endpoint().hash(),
            unitig_front.reverse_complement();
        else
            id_ = unitig_back.endpoint().hash(),
            unitig_back.reverse_complement();
    }
    else
    {
        id_ = cycle->min_vertex().hash();
        if(!cycle->min_vertex().in_canonical_form())
            cycle->reverse_complement();
    }
}


template <uint16_t k>
inline bool Maximal_Unitig_Scratch<k>::is_cycle() const
{
    return !is_linear();
}


template <uint16_t k>
inline const FASTA_Record<std::vector<char>> Maximal_Unitig_Scratch<k>::fasta_rec() const
{
    return is_canonical() ?
            FASTA_Record<std::vector<char>>(id(), unitig_front.label(), unitig_back.label(), 0, k) :
            FASTA_Record<std::vector<char>>(id(), unitig_back.label(), unitig_front.label(), 0, k);
}


template <uint16_t k>
template <std::size_t CAPACITY, typename T_sink_>
inline void Maximal_Unitig_Scratch<k>::add_fasta_rec_to_buffer(Character_Buffer<CAPACITY, T_sink_>& buffer) const
{
    if(is_linear())
        buffer += fasta_rec();
    else
        buffer.template rotate_append_cycle<k>(FASTA_Record<std::vector<char>>(id(), cycle->label()), cycle->min_vertex_idx());
}



#endif

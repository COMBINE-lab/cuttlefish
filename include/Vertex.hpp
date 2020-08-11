
#ifndef VERTEX_HPP
#define VERTEX_HPP



#include "globals.hpp"

#include <cstdint>
#include <cstdlib>
#include <iostream>


template <uint16_t k> class CdBG;


class Vertex
{
    template <uint16_t k>
    friend class CdBG;

private:

    cuttlefish::Vertex_Class vertex_class_;
    cuttlefish::nucleotide_t enter_;
    cuttlefish::nucleotide_t exit_;
    bool visited_;
    bool outputted_;


public:

    // Constructs an unvisited vertex.
    Vertex(): visited_(false)
    {}

    // Constructs a vertex of class `single_in_single_out`.
    Vertex(cuttlefish::Vertex_Class vertex_class, cuttlefish::nucleotide_t enter, cuttlefish::nucleotide_t exit);

    // Constructs a vertex of class either `multi_in_single_out` or `single_in_multi_out`.
    Vertex(cuttlefish::Vertex_Class vertex_class, cuttlefish::nucleotide_t nucl);

    // Constructs a vertex of class `vertex_class`, with the provided outputted status.
    Vertex(cuttlefish::Vertex_Class vertex_class, bool outputted = false);

    // Returns the vertex class.
    cuttlefish::Vertex_Class vertex_class() const;

    // Returns the entering nucleotide.
    cuttlefish::nucleotide_t enter() const;

    // Returns the exitting nucleotide.
    cuttlefish::nucleotide_t exit() const;

    // Returns the outputted status.
    bool outputted() const;

    // For debugging purposes: prints the vertex information.
    friend std::ostream& operator <<(std::ostream& out, const Vertex& vertex);
};



inline cuttlefish::Vertex_Class Vertex::vertex_class() const
{
    return vertex_class_;
}


inline cuttlefish::nucleotide_t Vertex::enter() const
{
    return enter_;
}


inline cuttlefish::nucleotide_t Vertex::exit() const
{
    return exit_;
}


inline bool Vertex::outputted() const
{
    return outputted_;
}



#endif

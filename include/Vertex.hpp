
#ifndef VERTEX_HPP
#define VERTEX_HPP


#include "globals.hpp"

#include <cstdint>
#include <cstdlib>
#include <iostream>


class Vertex
{
public:
    cuttlefish::Vertex_Class vertex_class;
    cuttlefish::nucleotide_t enter, exit;
    bool visited;
    bool outputted;


    // Constructs an unvisited vertex.
    Vertex(): visited(false)
    {}

    // Constructs a vertex of class `single_in_single_out`.
    Vertex(const cuttlefish::Vertex_Class vertex_class, const cuttlefish::nucleotide_t enter, const cuttlefish::nucleotide_t exit);

    // Constructs a vertex of class either `multi_in_single_out` or `single_in_multi_out`.
    Vertex(const cuttlefish::Vertex_Class vertex_class, const cuttlefish::nucleotide_t nucl);

    // Constructs a vertex of class `vertex_class`, with the provided outputted status.
    Vertex(const cuttlefish::Vertex_Class vertex_class, const bool outputted = false);

    // For debugging purposes: prints the vertex information.
    friend std::ostream& operator <<(std::ostream& out, const Vertex& vertex);
};


#endif

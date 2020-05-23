
#ifndef VERTEX_HPP
#define VERTEX_HPP


#include "globals.hpp"

#include <cstdint>
#include <cstdlib>
#include <iostream>


class Vertex
{
public:
    cuttlefish::state_t state;
    cuttlefish::nucleotide_t enter, exit;
    bool visited;
    bool outputted;


    // Constructs an unvisited vertex.
    Vertex(): visited(false)
    {}

    // Constructs a vertex of type SINGLE_IN_SINGLE_OUT.
    Vertex(const cuttlefish::state_t state, const cuttlefish::nucleotide_t enter, const cuttlefish::nucleotide_t exit);

    // Constructs a vertex of type MULTI_IN_SINGLE_OUT or SINGLE_IN_MULTI_OUT.
    Vertex(const cuttlefish::state_t state, const cuttlefish::nucleotide_t nucl);

    // Constructs a vertex of type `state`, with the provided outputted status.
    Vertex(const cuttlefish::state_t state, const bool outputted = false);

    // For debugging purposes: prints the vertex information.
    friend std::ostream& operator <<(std::ostream& out, const Vertex& vertex);
};


#endif

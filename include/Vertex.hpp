
#ifndef VERTEX_HPP
#define VERTEX_HPP



#include "globals.hpp"

#include <cstdint>
#include <iostream>


template <uint16_t k> class CdBG;


class Vertex
{
    template <uint16_t k>
    friend class CdBG;

private:

    cuttlefish::State_Class state_class_;
    cuttlefish::base_t front_;
    cuttlefish::base_t back_;
    bool visited_;
    bool outputted_;


public:

    // Constructs an unvisited vertex.
    Vertex(): visited_(false)
    {}

    // Constructs a vertex of state-class `single_in_single_out`.
    Vertex(cuttlefish::State_Class state_class, cuttlefish::base_t front, cuttlefish::base_t back);

    // Constructs a vertex of state-class either `multi_in_single_out` or `single_in_multi_out`.
    Vertex(cuttlefish::State_Class state_class, cuttlefish::base_t base);

    // Constructs a vertex of state-class `state_class`, with the provided outputted status.
    Vertex(cuttlefish::State_Class state_class, bool outputted = false);

    // Returns the state-class.
    cuttlefish::State_Class state_class() const;

    // Returns the base at the front of the vertex.
    cuttlefish::base_t front() const;

    // Returns the base at the back of the vertex.
    cuttlefish::base_t back() const;

    // Returns the outputted status.
    bool outputted() const;

    // For debugging purposes: prints the vertex information.
    friend std::ostream& operator <<(std::ostream& out, const Vertex& vertex);
};



inline cuttlefish::State_Class Vertex::state_class() const
{
    return state_class_;
}


inline cuttlefish::base_t Vertex::front() const
{
    return front_;
}


inline cuttlefish::base_t Vertex::back() const
{
    return back_;
}


inline bool Vertex::outputted() const
{
    return outputted_;
}



#endif

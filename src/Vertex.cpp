
#include "Vertex.hpp"
#include <cassert>


Vertex::Vertex(const cuttlefish::state_t state, const cuttlefish::nucleotide_t enter, const cuttlefish::nucleotide_t exit):
    state(state), enter(enter), exit(exit), visited(true), outputted(false)
{
    assert(state == cuttlefish::SINGLE_IN_SINGLE_OUT);
}


Vertex::Vertex(const cuttlefish::state_t state, const cuttlefish::nucleotide_t nucl):
    state(state), visited(true), outputted(false)
{
    assert(state == cuttlefish::MULTI_IN_SINGLE_OUT || state == cuttlefish::SINGLE_IN_MULTI_OUT);


    if(state == cuttlefish::MULTI_IN_SINGLE_OUT)
    {
        enter = cuttlefish::PLACEHOLDER_NUCLEOTIDE;
        exit = nucl;
    }
    else
    {
        enter = nucl;
        exit = cuttlefish::PLACEHOLDER_NUCLEOTIDE;
    }
}


Vertex::Vertex(const cuttlefish::state_t state, const bool outputted):
    state(state), visited(true), outputted(outputted)
{
    enter = cuttlefish::PLACEHOLDER_NUCLEOTIDE;
    exit = cuttlefish::PLACEHOLDER_NUCLEOTIDE;
}


std::ostream& operator <<(std::ostream& out, const Vertex& vertex)
{
    std::string label;

    switch (vertex.state)
    {
    case cuttlefish::SINGLE_IN_SINGLE_OUT:
        label = "Single_In_Single_Out";
        break;
    
    case cuttlefish::MULTI_IN_SINGLE_OUT:
        label = "Multi_In_Single_Out";
        break;

    case cuttlefish::SINGLE_IN_MULTI_OUT:
        label = "Single_In_Multi_Out";
        break;
    
    case cuttlefish::MULTI_IN_MULTI_OUT:
        label = "Multi_In_Multi_Out";
        break;
    
    default:
        std::cerr << "Invalid state encountered for vertices.\n";
        break;
    }


    out << label;
    return out;
}

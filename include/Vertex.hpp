
#ifndef VERTEX_HPP
#define VERTEX_HPP


#include <cstdint>
#include <cstdlib>
#include <iostream>

#include "globals.hpp"


class Vertex
{
public:
    uint8_t state;
    cuttlefish::nucleotide_t enter, exit;
    bool outputted;


    Vertex()
    {}


    Vertex(const cuttlefish::state_t state, const cuttlefish::nucleotide_t enter, const cuttlefish::nucleotide_t exit):
            state(state), enter(enter), exit(exit), outputted(false)
    {}


    Vertex(const cuttlefish::state_t state, const cuttlefish::nucleotide_t nucl):
            state(state), outputted(false)
    {
        if(state == cuttlefish::MULTI_IN_SINGLE_OUT)
        {
            enter = cuttlefish::PLACEHOLDER_NUCLEOTIDE;
            exit = nucl;
        }
        else if(state == cuttlefish::SINGLE_IN_MULTI_OUT)
        {
            enter = nucl;
            exit = cuttlefish::PLACEHOLDER_NUCLEOTIDE;
        }
        else
        {
            std::cerr << "Invalid construction of a vertex. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
    }


    Vertex(const cuttlefish::state_t state):
            state(state), outputted(false)
    {
        if(state == cuttlefish::MULTI_IN_MULTI_OUT)
        {
            enter = cuttlefish::PLACEHOLDER_NUCLEOTIDE;
            exit = cuttlefish::PLACEHOLDER_NUCLEOTIDE;
        }
    }


    friend std::ostream& operator <<(std::ostream& out, const Vertex& vertex)
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
};


#endif

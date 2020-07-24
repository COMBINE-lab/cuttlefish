
#include "Vertex.hpp"
#include <cassert>


Vertex::Vertex(const cuttlefish::Vertex_Class vertex_class, const cuttlefish::nucleotide_t enter, const cuttlefish::nucleotide_t exit):
    vertex_class_(vertex_class), enter_(enter), exit_(exit), visited_(true), outputted_(false)
{
    assert(vertex_class == cuttlefish::Vertex_Class::single_in_single_out);
}


Vertex::Vertex(const cuttlefish::Vertex_Class vertex_class, const cuttlefish::nucleotide_t nucl):
    vertex_class_(vertex_class), visited_(true), outputted_(false)
{
    assert(vertex_class == cuttlefish::Vertex_Class::multi_in_single_out || vertex_class == cuttlefish::Vertex_Class::single_in_multi_out);


    if(vertex_class == cuttlefish::Vertex_Class::multi_in_single_out)
    {
        enter_ = cuttlefish::PLACEHOLDER_NUCLEOTIDE;
        exit_ = nucl;
    }
    else
    {
        enter_ = nucl;
        exit_ = cuttlefish::PLACEHOLDER_NUCLEOTIDE;
    }
}


Vertex::Vertex(const cuttlefish::Vertex_Class vertex_class, const bool outputted):
    vertex_class_(vertex_class), visited_(true), outputted_(outputted)
{
    enter_ = cuttlefish::PLACEHOLDER_NUCLEOTIDE;
    exit_ = cuttlefish::PLACEHOLDER_NUCLEOTIDE;
}


std::ostream& operator <<(std::ostream& out, const Vertex& vertex)
{
    std::string label;

    switch (vertex.vertex_class_)
    {
    case cuttlefish::Vertex_Class::single_in_single_out:
        label = "Single_In_Single_Out";
        break;
    
    case cuttlefish::Vertex_Class::multi_in_single_out:
        label = "Multi_In_Single_Out";
        break;

    case cuttlefish::Vertex_Class::single_in_multi_out:
        label = "Single_In_Multi_Out";
        break;
    
    case cuttlefish::Vertex_Class::multi_in_multi_out:
        label = "Multi_In_Multi_Out";
        break;
    
    default:
        std::cerr << "Invalid state encountered for vertices.\n";
        break;
    }


    out << label;
    return out;
}

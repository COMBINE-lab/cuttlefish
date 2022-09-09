
#include "Vertex.hpp"
#include "DNA.hpp"

#include <cassert>


Vertex::Vertex(const cuttlefish::State_Class state_class, const cuttlefish::base_t front, const cuttlefish::base_t back):
    state_class_(state_class), front_(front), back_(back), visited_(true), outputted_(false)
{
    assert(state_class == cuttlefish::State_Class::single_in_single_out);
}


Vertex::Vertex(const cuttlefish::State_Class state_class, const cuttlefish::base_t base):
    state_class_(state_class), visited_(true), outputted_(false)
{
    assert(state_class == cuttlefish::State_Class::multi_in_single_out || state_class == cuttlefish::State_Class::single_in_multi_out);


    if(state_class == cuttlefish::State_Class::multi_in_single_out)
    {
        front_ = DNA::N;
        back_ = base;
    }
    else
    {
        front_ = base;
        back_ = DNA::N;
    }
}


Vertex::Vertex(const cuttlefish::State_Class state_class, const bool outputted):
    state_class_(state_class), visited_(true), outputted_(outputted)
{
    front_ = DNA::N;
    back_ = DNA::N;
}


std::ostream& operator <<(std::ostream& out, const Vertex& vertex)
{
    std::string label;

    switch (vertex.state_class_)
    {
    case cuttlefish::State_Class::single_in_single_out:
        label = "Single_In_Single_Out";
        break;
    
    case cuttlefish::State_Class::multi_in_single_out:
        label = "Multi_In_Single_Out";
        break;

    case cuttlefish::State_Class::single_in_multi_out:
        label = "Single_In_Multi_Out";
        break;
    
    case cuttlefish::State_Class::multi_in_multi_out:
        label = "Multi_In_Multi_Out";
        break;
    
    default:
        std::cerr << "Invalid state encountered for vertices.\n";
        break;
    }


    out << label;
    return out;
}

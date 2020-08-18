
#include "State.hpp"
#include "globals.hpp"

#include <cstdlib>
#include <iostream>


State::State(): code(0b00000)
{}


State::State(const Vertex& vertex)
{
    if(!vertex.outputted())
        switch(vertex.vertex_class())
        {
        case cuttlefish::Vertex_Class::single_in_single_out:
            code = 0b10000;
            set_nibble_upper_half(vertex.enter());
            set_nibble_lower_half(vertex.exit());

            break;

        case cuttlefish::Vertex_Class::multi_in_single_out:
            code = 0b00100;
            set_nibble_lower_half(vertex.exit());

            break;

        case cuttlefish::Vertex_Class::single_in_multi_out:
            code = 0b01000;
            set_nibble_lower_half(vertex.enter());

            break;

        case cuttlefish::Vertex_Class::multi_in_multi_out:
            code = 0b00011;

            break;

        default:
            std::cerr << "Invalid vertex class " << (uint16_t)vertex.vertex_class() << " encountered during state construction. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
    else
        switch(vertex.vertex_class())
        {
        case cuttlefish::Vertex_Class::single_in_single_out:
            code = 0b01100;

            break;

        case cuttlefish::Vertex_Class::multi_in_single_out:
            code = 0b01101;

            break;

        case cuttlefish::Vertex_Class::single_in_multi_out:
            code = 0b01110;

            break;

        case cuttlefish::Vertex_Class::multi_in_multi_out:
            code = 0b01111;

            break;

        default:
            std::cerr << "Invalid vertex class " << (uint16_t)vertex.vertex_class() << " encountered during state construction. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
}


inline void State::set_nibble_lower_half(const cuttlefish::base_t base)
{
    switch(base)
    {
    case DNA_Base::A:   // A --> 00
        break;
    
    case DNA_Base::C:   // C --> 01
        code |= uint8_t(0b0001);
        break;

    case DNA_Base::G:   // G --> 10
        code |= uint8_t(0b0010);
        break;

    case DNA_Base::T:   // T --> 11
        code |= uint8_t(0b0011);
        break;

    default:
        std::cerr << "Invalid base " << base << " encountered during state construction. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


inline void State::set_nibble_upper_half(const cuttlefish::base_t base)
{
    switch(base)
    {
    case DNA_Base::A:   // A --> 00
        break;
    
    case DNA_Base::C:   // C --> 01
        code |= uint8_t(0b0100);
        break;

    case DNA_Base::G:   // G --> 10
        code |= uint8_t(0b1000);
        break;

    case DNA_Base::T:   // T --> 11
        code |= uint8_t(0b1100);
        break;

    default:
        std::cerr << "Invalid base " << base << " encountered during state construction. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


Vertex State::decode() const
{
    // TODO: Define a class member const static array and return elements from those,
    // instead of constructing values everytime.

    switch(code)
    {
    case 0b00000:   // 0
        return Vertex();

    case 0b00001:   // 1
    case 0b00010:   // 2
        std::cerr << "Invalid state " << (uint16_t)code << " encountered during state decoding. Aborting.\n";
        std::exit(EXIT_FAILURE);

    case 0b00011:   // 3
        return Vertex(cuttlefish::Vertex_Class::multi_in_multi_out);

    case 0b00100:   // 4
        return Vertex(cuttlefish::Vertex_Class::multi_in_single_out, A);

    case 0b00101:   // 5
        return Vertex(cuttlefish::Vertex_Class::multi_in_single_out, C);

    case 0b00110:   // 6
        return Vertex(cuttlefish::Vertex_Class::multi_in_single_out, G);

    case 0b00111:   // 7
        return Vertex(cuttlefish::Vertex_Class::multi_in_single_out, T);

    case 0b01000:   // 8
        return Vertex(cuttlefish::Vertex_Class::single_in_multi_out, A);

    case 0b01001:   // 9
        return Vertex(cuttlefish::Vertex_Class::single_in_multi_out, C);

    case 0b01010:   // 10
        return Vertex(cuttlefish::Vertex_Class::single_in_multi_out, G);

    case 0b01011:   // 11
        return Vertex(cuttlefish::Vertex_Class::single_in_multi_out, T);

    case 0b01100:   // 12
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, true);

    case 0b01101:   // 13
        return Vertex(cuttlefish::Vertex_Class::multi_in_single_out, true);

    case 0b01110:   // 14
        return Vertex(cuttlefish::Vertex_Class::single_in_multi_out, true);

    case 0b01111:   // 15
        return Vertex(cuttlefish::Vertex_Class::multi_in_multi_out, true);

    case 0b10000:   // 16
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, A, A);

    case 0b10001:   // 17
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, A, C);

    case 0b10010:   // 18
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, A, G);

    case 0b10011:   // 19
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, A, T);

    case 0b10100:   // 20
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, C, A);

    case 0b10101:   // 21
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, C, C);

    case 0b10110:   // 22
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, C, G);

    case 0b10111:   // 23
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, C, T);

    case 0b11000:   // 24
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, G, A);

    case 0b11001:   // 25
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, G, C);

    case 0b11010:   // 26
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, G, G);

    case 0b11011:   // 27
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, G, T);

    case 0b11100:   // 28
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, T, A);

    case 0b11101:   // 29
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, T, C);

    case 0b11110:   // 30
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, T, G);

    case 0b11111:   // 31
        return Vertex(cuttlefish::Vertex_Class::single_in_single_out, T, T);

    default:
        std::cerr << "Invalid state " << (uint16_t)code << " encountered during state decoding. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


bool State::is_visited() const
{
    switch(code)
    {
    case 0b00000:   // 0
        return false;

    case 0b00001:   // 1
    case 0b00010:   // 2
        std::cerr << "Invalid state " << (uint16_t)code << " encountered during checking visited status for vertices. Aborting.\n";
        std::exit(EXIT_FAILURE);

    case 0b00011:   // 3
    case 0b00100:   // 4
    case 0b00101:   // 5
    case 0b00110:   // 6
    case 0b00111:   // 7
    case 0b01000:   // 8
    case 0b01001:   // 9
    case 0b01010:   // 10
    case 0b01011:   // 11
    case 0b01100:   // 12
    case 0b01101:   // 13
    case 0b01110:   // 14
    case 0b01111:   // 15
    case 0b10000:   // 16
    case 0b10001:   // 17
    case 0b10010:   // 18
    case 0b10011:   // 19
    case 0b10100:   // 20
    case 0b10101:   // 21
    case 0b10110:   // 22
    case 0b10111:   // 23
    case 0b11000:   // 24
    case 0b11001:   // 25
    case 0b11010:   // 26
    case 0b11011:   // 27
    case 0b11100:   // 28
    case 0b11101:   // 29
    case 0b11110:   // 30
    case 0b11111:   // 31
        return true;

    default:
        std::cerr << "Invalid state " << (uint16_t)code << " encountered during checking visited status for vertices. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


State State::outputted() const
{
    switch(code)
    {
    case 0b00000:   // 0
        std::cerr << "Invalid output attempt for an unvisited vertex. Aborting.\n";
        std::exit(EXIT_FAILURE);

    case 0b00001:   // 1
    case 0b00010:   // 2
        std::cerr << "Invalid state " << (uint16_t)code << " encountered during output attempt. Aborting.\n";
        std::exit(EXIT_FAILURE);

    case 0b00011:   // 3
        return State(0b01111);

    case 0b00100:   // 4
    case 0b00101:   // 5
    case 0b00110:   // 6
    case 0b00111:   // 7
        return State(0b01101);

    case 0b01000:   // 8
    case 0b01001:   // 9
    case 0b01010:   // 10
    case 0b01011:   // 11
        return State(0b01110);

    case 0b01100:   // 12
        return State(0b01100);

    case 0b01101:   // 13
        return State(0b01101);

    case 0b01110:   // 14
        return State(0b01110);

    case 0b01111:   // 15
        return State(0b01111);

    case 0b10000:   // 16
    case 0b10001:   // 17
    case 0b10010:   // 18
    case 0b10011:   // 19
    case 0b10100:   // 20
    case 0b10101:   // 21
    case 0b10110:   // 22
    case 0b10111:   // 23
    case 0b11000:   // 24
    case 0b11001:   // 25
    case 0b11010:   // 26
    case 0b11011:   // 27
    case 0b11100:   // 28
    case 0b11101:   // 29
    case 0b11110:   // 30
    case 0b11111:   // 31
        return State(0b01100);

    default:
        std::cerr << "Invalid state " << (uint16_t)code << " encountered during output attempt. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


cuttlefish::Vertex_Class State::vertex_class() const
{
    switch(code)
    {
    case 0b00000:   // 0
        std::cerr << "No class for an unvisited vertex. Aborting.\n";
        std::exit(EXIT_FAILURE);

    case 0b00001:   // 1
    case 0b00010:   // 2
        std::cerr << "Invalid state " << (uint16_t)code << " encountered during vertex class decoding. Aborting.\n";
        std::exit(EXIT_FAILURE);

    case 0b00011:   // 3
        return cuttlefish::Vertex_Class::multi_in_multi_out;

    case 0b00100:   // 4
    case 0b00101:   // 5
    case 0b00110:   // 6
    case 0b00111:   // 7
        return cuttlefish::Vertex_Class::multi_in_single_out;

    case 0b01000:   // 8
    case 0b01001:   // 9
    case 0b01010:   // 10
    case 0b01011:   // 11
        return cuttlefish::Vertex_Class::single_in_multi_out;

    case 0b01100:   // 12
        return cuttlefish::Vertex_Class::single_in_single_out;

    case 0b01101:   // 13
        return cuttlefish::Vertex_Class::multi_in_single_out;

    case 0b01110:   // 14
        return cuttlefish::Vertex_Class::single_in_multi_out;

    case 0b01111:   // 15
        return cuttlefish::Vertex_Class::multi_in_multi_out;

    case 0b10000:   // 16
    case 0b10001:   // 17
    case 0b10010:   // 18
    case 0b10011:   // 19
    case 0b10100:   // 20
    case 0b10101:   // 21
    case 0b10110:   // 22
    case 0b10111:   // 23
    case 0b11000:   // 24
    case 0b11001:   // 25
    case 0b11010:   // 26
    case 0b11011:   // 27
    case 0b11100:   // 28
    case 0b11101:   // 29
    case 0b11110:   // 30
    case 0b11111:   // 31
        return cuttlefish::Vertex_Class::single_in_single_out;

    default:
        std::cerr << "Invalid state " << (uint16_t)code << " encountered during vertex class decoding. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


bool State::is_outputted() const
{
    switch(code)
    {
    case 0b00000:   // 0
        return false;

    case 0b00001:   // 1
    case 0b00010:   // 2
        std::cerr << "Invalid state " << (uint16_t)code << " encountered during checking ouuputted status. Aborting.\n";
        std::exit(EXIT_FAILURE);

    case 0b00011:   // 3
    case 0b00100:   // 4
    case 0b00101:   // 5
    case 0b00110:   // 6
    case 0b00111:   // 7
    case 0b01000:   // 8
    case 0b01001:   // 9
    case 0b01010:   // 10
    case 0b01011:   // 11
        return false;

    case 0b01100:   // 12
    case 0b01101:   // 13
    case 0b01110:   // 14
    case 0b01111:   // 15
        return true;

    case 0b10000:   // 16
    case 0b10001:   // 17
    case 0b10010:   // 18
    case 0b10011:   // 19
    case 0b10100:   // 20
    case 0b10101:   // 21
    case 0b10110:   // 22
    case 0b10111:   // 23
    case 0b11000:   // 24
    case 0b11001:   // 25
    case 0b11010:   // 26
    case 0b11011:   // 27
    case 0b11100:   // 28
    case 0b11101:   // 29
    case 0b11110:   // 30
    case 0b11111:   // 31
        return false;

    default:
        std::cerr << "Invalid state " << (uint16_t)code << " encountered during checking outputted status. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


std::ostream& operator<<(std::ostream& out, const State& state)
{
    out << (uint16_t)state.code;

    return out;
}
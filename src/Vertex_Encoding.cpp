
#include "Vertex_Encoding.hpp"
#include "globals.hpp"

#include <cstdlib>
#include <iostream>


Vertex_Encoding::Vertex_Encoding():
    vertex_code(0b00000)
{}


Vertex_Encoding::Vertex_Encoding(const Vertex& vertex)
{
    if(!vertex.outputted)
        switch(vertex.state)
        {
        case cuttlefish::SINGLE_IN_SINGLE_OUT:
            vertex_code = 0b10000;
            set_nibble_upper_half(vertex.enter);
            set_nibble_lower_half(vertex.exit);

            break;

        case cuttlefish::MULTI_IN_SINGLE_OUT:
            vertex_code = 0b00100;
            set_nibble_lower_half(vertex.exit);

            break;

        case cuttlefish::SINGLE_IN_MULTI_OUT:
            vertex_code = 0b01000;
            set_nibble_lower_half(vertex.enter);

            break;

        case cuttlefish::MULTI_IN_MULTI_OUT:
            vertex_code = 0b00011;

            break;

        default:
            std::cerr << "Invalid vertex state " << (uint16_t)vertex.state << " encountered. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
    else
        switch(vertex.state)
        {
        case cuttlefish::SINGLE_IN_SINGLE_OUT:
            vertex_code = 0b01100;

            break;

        case cuttlefish::MULTI_IN_SINGLE_OUT:
            vertex_code = 0b01101;

            break;

        case cuttlefish::SINGLE_IN_MULTI_OUT:
            vertex_code = 0b01110;

            break;

        case cuttlefish::MULTI_IN_MULTI_OUT:
            vertex_code = 0b01111;

            break;

        default:
            std::cerr << "Invalid vertex state " << (uint16_t)vertex.state << " encountered. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
}


inline void Vertex_Encoding::set_nibble_lower_half(const cuttlefish::nucleotide_t nucl)
{
    switch(nucl)
    {
    case 'A':   // A --> 00
        break;
    
    case 'C':   // C --> 01
        vertex_code |= uint8_t(0b0001);
        break;

    case 'G':   // G --> 10
        vertex_code |= uint8_t(0b0010);
        break;

    case 'T':   // T --> 11
        vertex_code |= uint8_t(0b0011);
        break;

    default:
        // Placeholder rule to handle `N` nucleotides.
        // TODO: Need to make an informed rule for this.
        break;
        
        // std::cerr << "Invalid nucelotide " << nucl << " encountered. Aborting.\n";
        // std::exit(EXIT_FAILURE);
    }
}


inline void Vertex_Encoding::set_nibble_upper_half(const cuttlefish::nucleotide_t nucl)
{
    switch(nucl)
    {
    case 'A':   // A --> 00
        break;
    
    case 'C':   // C --> 01
        vertex_code |= uint8_t(0b0100);
        break;

    case 'G':   // G --> 10
        vertex_code |= uint8_t(0b1000);
        break;

    case 'T':   // T --> 11
        vertex_code |= uint8_t(0b1100);
        break;

    default:
        // Placeholder rule to handle `N` nucleotides.
        // TODO: Need to make an informed rule for this.
        break;

        // std::cerr << "Invalid nucelotide " << nucl << " encountered. Aborting.\n";
        // std::exit(EXIT_FAILURE);
    }
}


Vertex Vertex_Encoding::decode() const
{
    switch(vertex_code)
    {
    case 0b00000:   // 0
        return Vertex();

    case 0b00001:   // 1
    case 0b00010:   // 2
        std::cerr << "Invalid vertex encoding " << (uint16_t)vertex_code << " encountered. Aborting.\n";
        std::exit(EXIT_FAILURE);

    case 0b00011:   // 3
        return Vertex(cuttlefish::MULTI_IN_MULTI_OUT);

    case 0b00100:   // 4
        return Vertex(cuttlefish::MULTI_IN_SINGLE_OUT, 'A');

    case 0b00101:   // 5
        return Vertex(cuttlefish::MULTI_IN_SINGLE_OUT, 'C');

    case 0b00110:   // 6
        return Vertex(cuttlefish::MULTI_IN_SINGLE_OUT, 'G');

    case 0b00111:   // 7
        return Vertex(cuttlefish::MULTI_IN_SINGLE_OUT, 'T');

    case 0b01000:   // 8
        return Vertex(cuttlefish::SINGLE_IN_MULTI_OUT, 'A');

    case 0b01001:   // 9
        return Vertex(cuttlefish::SINGLE_IN_MULTI_OUT, 'C');

    case 0b01010:   // 10
        return Vertex(cuttlefish::SINGLE_IN_MULTI_OUT, 'G');

    case 0b01011:   // 11
        return Vertex(cuttlefish::SINGLE_IN_MULTI_OUT, 'T');

    case 0b01100:   // 12
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, true);

    case 0b01101:   // 13
        return Vertex(cuttlefish::MULTI_IN_SINGLE_OUT, true);

    case 0b01110:   // 14
        return Vertex(cuttlefish::SINGLE_IN_MULTI_OUT, true);

    case 0b01111:   // 15
        return Vertex(cuttlefish::MULTI_IN_MULTI_OUT, true);

    case 0b10000:   // 16
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'A', 'A');

    case 0b10001:   // 17
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'A', 'C');

    case 0b10010:   // 18
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'A', 'G');

    case 0b10011:   // 19
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'A', 'T');

    case 0b10100:   // 20
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'C', 'A');

    case 0b10101:   // 21
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'C', 'C');

    case 0b10110:   // 22
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'C', 'G');

    case 0b10111:   // 23
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'C', 'T');

    case 0b11000:   // 24
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'G', 'A');

    case 0b11001:   // 25
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'G', 'C');

    case 0b11010:   // 26
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'G', 'G');

    case 0b11011:   // 27
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'G', 'T');

    case 0b11100:   // 28
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'T', 'A');

    case 0b11101:   // 29
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'T', 'C');

    case 0b11110:   // 30
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'T', 'G');

    case 0b11111:   // 31
        return Vertex(cuttlefish::SINGLE_IN_SINGLE_OUT, 'T', 'T');

    default:
        std::cerr << "Invalid vertex encoding " << vertex_code << " encountered. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


Vertex_Encoding Vertex_Encoding::outputted() const
{
    switch(vertex_code)
    {
    case 0b00000:   // 0
        std::cerr << "Invalid output attempt for an unvisited vertex. Aborting.\n";
        std::exit(EXIT_FAILURE);

    case 0b00001:   // 1
    case 0b00010:   // 2
        std::cerr << "Invalid vertex encoding " << vertex_code << " encountered. Aborting.\n";
        std::exit(EXIT_FAILURE);

    case 0b00011:   // 3
        return Vertex_Encoding(0b01111);

    case 0b00100:   // 4
    case 0b00101:   // 5
    case 0b00110:   // 6
    case 0b00111:   // 7
        return Vertex_Encoding(0b01101);

    case 0b01000:   // 8
    case 0b01001:   // 9
    case 0b01010:   // 10
    case 0b01011:   // 11
        return Vertex_Encoding(0b01110);

    case 0b01100:   // 12
        return Vertex_Encoding(0b01100);

    case 0b01101:   // 13
        return Vertex_Encoding(0b01101);

    case 0b01110:   // 14
        return Vertex_Encoding(0b01110);

    case 0b01111:   // 15
        return Vertex_Encoding(0b01111);

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
        return Vertex_Encoding(0b01100);

    default:
        std::cerr << "Invalid vertex encoding " << vertex_code << " encountered. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


cuttlefish::state_t Vertex_Encoding::state() const
{
    switch(vertex_code)
    {
    case 0b00000:   // 0
        std::cerr << "No state for an unvisited vertex. Aborting.\n";
        std::exit(EXIT_FAILURE);

    case 0b00001:   // 1
    case 0b00010:   // 2
        std::cerr << "Invalid vertex encoding " << vertex_code << " encountered. Aborting.\n";
        std::exit(EXIT_FAILURE);

    case 0b00011:   // 3
        return cuttlefish::MULTI_IN_MULTI_OUT;

    case 0b00100:   // 4
    case 0b00101:   // 5
    case 0b00110:   // 6
    case 0b00111:   // 7
        return cuttlefish::MULTI_IN_SINGLE_OUT;

    case 0b01000:   // 8
    case 0b01001:   // 9
    case 0b01010:   // 10
    case 0b01011:   // 11
        return cuttlefish::SINGLE_IN_MULTI_OUT;

    case 0b01100:   // 12
        return cuttlefish::SINGLE_IN_SINGLE_OUT;

    case 0b01101:   // 13
        return cuttlefish::MULTI_IN_SINGLE_OUT;

    case 0b01110:   // 14
        return cuttlefish::SINGLE_IN_MULTI_OUT;

    case 0b01111:   // 15
        return cuttlefish::MULTI_IN_MULTI_OUT;

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
        return cuttlefish::SINGLE_IN_SINGLE_OUT;

    default:
        std::cerr << "Invalid vertex encoding " << vertex_code << " encountered. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


bool Vertex_Encoding::is_outputted() const
{
    switch(vertex_code)
    {
    case 0b00000:   // 0
        return false;

    case 0b00001:   // 1
    case 0b00010:   // 2
        std::cerr << "Invalid vertex encoding " << vertex_code << " encountered. Aborting.\n";
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
        std::cerr << "Invalid vertex encoding " << vertex_code << " encountered. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


std::ostream& operator <<(std::ostream& out, const Vertex_Encoding& vertex_encoding)
{
    out << (uint16_t)vertex_encoding.vertex_code;

    return out;
}

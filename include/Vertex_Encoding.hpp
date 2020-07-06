
#ifndef VERTEX_ENCODING_HPP
#define VERTEX_ENCODING_HPP



#include "Vertex.hpp"
#include "compact_vector/compact_iterator.hpp"

#include <cstdint>


class Kmer_Hash_Table;
class Kmer_Hash_Entry_API;


class Vertex_Encoding
{
    friend class Kmer_Hash_Table;
    friend class Kmer_Hash_Entry_API;

private:

    // The encoded information of a vertex.
    cuttlefish::vertex_code_t vertex_code;


    // Constructs an encoding with the provided vertex code.
    Vertex_Encoding(const cuttlefish::vertex_code_t vertex_code);

    // Constructs an encoding from the state stored at bitvector entry `bv_entry`.
    Vertex_Encoding(const cuttlefish::bitvector_entry_t& bv_entry);

    // Sets the nucleotide 2-bit encoding at the bits b1 and b0 of `vertex_code`.
    // Requirement: the two bits must be zero before the call, for consistent
    // behavior.
    void set_nibble_lower_half(const cuttlefish::nucleotide_t nucl);

    // Sets the nucleotide 2-bit encoding at the bits b3 and b2 of `vertex_code`.
    // Requirement: the two bits must be zero before the call, for consistent
    // behavior.
    void set_nibble_upper_half(const cuttlefish::nucleotide_t nucl);

    // Returns the wrapped vertex code value.
    cuttlefish::vertex_code_t get_vertex_code() const
    {
        return vertex_code;
    }


public:

    // Constructs an encoding of an univisisted vertex.
    Vertex_Encoding();

    // Constructs an encoding of the provided `vertex`.
    Vertex_Encoding(const Vertex& vertex);

    // Returns a vertex by decoding the encoded information.
    Vertex decode() const;

    // Returns a Boolean denoting whether the underlying vertex has been
    // visited.
    bool is_visited() const;

    // Returns an encoding of the underlying vertex if it had been outputted.
    Vertex_Encoding outputted() const;

    // Returns the state of the encoded vertex.
    cuttlefish::Vertex_Class vertex_class() const;

    // Returns the output status of the encoded vertex.
    bool is_outputted() const;

    // TODO
    bool is_dead_end() const;

    // Returns true iff the underlying codes are different.
    bool operator!=(const Vertex_Encoding& rhs) const
    {
        return vertex_code != rhs.vertex_code;
    }

    // For debugging.
    friend std::ostream& operator <<(std::ostream& out, const Vertex_Encoding& vertex_encoding);



private:
    // TODO:
    // Are lookup-tables faster than switch-cases? If yes, define the following arrays,
    // and modify the public function definitions.

    const static Vertex decoded_vertex[32];

    const static Vertex_Encoding encoded_outputted_vertex[32];

    const static cuttlefish::Vertex_Class decoded_vertex_class[32];

    const static bool decoded_output_status[32];
};


inline Vertex_Encoding::Vertex_Encoding(const cuttlefish::vertex_code_t vertex_code):
    vertex_code(vertex_code)
{
    if(vertex_code == 0b00001 || vertex_code == 0b00010)
    {
        std::cerr << "Invalid vertex encoding " << (uint16_t)vertex_code << " encountered during construction from code. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


inline Vertex_Encoding::Vertex_Encoding(const cuttlefish::bitvector_entry_t& bv_entry)
{
    // CAS vector `fetch` does not work yet.
    // bv_entry.fetch_val(vertex_code);

    vertex_code = bv_entry;
}



#endif
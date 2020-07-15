
#ifndef KMER_HASH_ENTRY_API_HPP
#define KMER_HASH_ENTRY_API_HPP



#include "globals.hpp"
#include "Vertex_Encoding.hpp"


class Kmer_Hash_Table;


// Wrapper class acting as an API to the entries of the bitvector used as hash table for k-mers.
class Kmer_Hash_Entry_API
{
    friend class Kmer_Hash_Table;


private:

    // Position information (base pointer and offset) for the bitvector entry.
    cuttlefish::bitvector_entry_t bv_entry;

    // Value read from the bitvector entry when the object is constructed; is immutable.
    const Vertex_Encoding vertex_encoding_read;

    // Value read from the bitvector entry when the object is constrcuted; is mutable.
    Vertex_Encoding vertex_encoding;


    // Constructs an API to the bitvector entry `bv_entry`.
    Kmer_Hash_Entry_API(const cuttlefish::bitvector_entry_t& bv_entry):
        bv_entry(bv_entry), vertex_encoding_read(bv_entry)
    {
        vertex_encoding = vertex_encoding_read;
    }


    // Returns the vertex encoding value read when the object was constructed.
    uint8_t get_read_encoding() const;

    
    // Returns the value of the mutable vertex encoding wrapped inside the API,
    // i.e. the encoding value that had been read at the object creation, and then
    // possibly have been modified.
    uint8_t get_current_encoding() const;


public:

    // Returns a reference to the mutable copy of the wrapped vertex encoding.
    Vertex_Encoding& get_vertex_encoding();
};



inline uint8_t Kmer_Hash_Entry_API::get_read_encoding() const
{
    return vertex_encoding_read.get_vertex_code();
}


inline uint8_t Kmer_Hash_Entry_API::get_current_encoding() const
{
    return vertex_encoding.get_vertex_code();
}


inline Vertex_Encoding& Kmer_Hash_Entry_API::get_vertex_encoding()
{
    return vertex_encoding;
}



#endif

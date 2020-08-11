
#ifndef KMER_HASH_ENTRY_API_HPP
#define KMER_HASH_ENTRY_API_HPP



#include "globals.hpp"
#include "State.hpp"


template <uint16_t k> class Kmer_Hash_Table;


// Wrapper class acting as an API to the entries of the bitvector used as hash table for k-mers.
class Kmer_Hash_Entry_API
{
    template <uint16_t k>
    friend class Kmer_Hash_Table;


private:

    // Position information (base pointer and offset) for the bitvector entry.
    cuttlefish::bitvector_entry_t bv_entry;

    // Value read from the bitvector entry when the object is constructed; is immutable.
    const State state_read;

    // Value read from the bitvector entry when the object is constrcuted; is mutable.
    State state;


    // Constructs an API to the bitvector entry `bv_entry`.
    Kmer_Hash_Entry_API(const cuttlefish::bitvector_entry_t& bv_entry);

    // Returns the state value read when the object was constructed.
    cuttlefish::state_code_t get_read_state() const;
    
    // Returns the value of the mutable state value wrapped inside the API,
    // i.e. the state value that had been read at the object creation, and then
    // possibly have been modified.
    cuttlefish::state_code_t get_current_state() const;


public:

    // Returns a reference to the mutable copy of the wrapped state value.
    State& get_state();
};


inline Kmer_Hash_Entry_API::Kmer_Hash_Entry_API(const cuttlefish::bitvector_entry_t& bv_entry):
    bv_entry(bv_entry), state_read(bv_entry)
{
    state = state_read;
}


inline cuttlefish::state_code_t Kmer_Hash_Entry_API::get_read_state() const
{
    return state_read.get_state();
}


inline cuttlefish::state_code_t Kmer_Hash_Entry_API::get_current_state() const
{
    return state.get_state();
}


inline State& Kmer_Hash_Entry_API::get_state()
{
    return state;
}



#endif

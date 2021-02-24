
#ifndef KMER_HASH_ENTRY_API_HPP
#define KMER_HASH_ENTRY_API_HPP



#include "globals.hpp"
#include "State.hpp"
#include "State_Read_Space.hpp"


template <uint16_t k, uint8_t BITS_PER_KEY> class Kmer_Hash_Table;


// Wrapper class acting as an API to the entries of the bitvector used as hash table for k-mers.
template <uint8_t BITS_PER_KEY>
class Kmer_Hash_Entry_API
{};


// Instantiation of the API class used in reference-dBG compaction.
template <>
class Kmer_Hash_Entry_API<cuttlefish::BITS_PER_REF_KMER>
{
    template <uint16_t k, uint8_t BITS_PER_KEY>
    friend class Kmer_Hash_Table;

    typedef compact::iterator_imp::lhs_setter<cuttlefish::state_code_t, cuttlefish::BITS_PER_REF_KMER, uint64_t, true, 64U> bitvector_entry_t;


private:

    // Position information (base pointer and offset) for the bitvector entry.
    bitvector_entry_t bv_entry;

    // Value read from the bitvector entry when the object is constructed; is immutable.
    const State state_read;

    // Value read from the bitvector entry when the object is constructed; is mutable.
    State state;


    // Constructs an API to the bitvector entry `bv_entry`.
    Kmer_Hash_Entry_API(const bitvector_entry_t& bv_entry):
        bv_entry(bv_entry), state_read(bv_entry)
    {
        state = state_read;
    }

    // Returns the state value read when the object was constructed.
    cuttlefish::state_code_t get_read_state() const
    {
        return state_read.get_state();
    }
    
    // Returns the value of the mutable state value wrapped inside the API,
    // i.e. the state value that had been read at the object creation, and then
    // possibly have been modified.
    cuttlefish::state_code_t get_current_state() const
    {
        return state.get_state();
    }


public:

    // Returns a reference to the mutable copy of the wrapped state value.
    State& get_state()
    {
        return state;
    }
};


// Instantiation of the API class used in read-dBG compaction.
template <>
class Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER>
{
    template <uint16_t k, uint8_t BITS_PER_KMER>
    friend class Kmer_Hash_Table;

    typedef compact::iterator_imp::lhs_setter<cuttlefish::state_code_t, cuttlefish::BITS_PER_READ_KMER, uint64_t, true, 64U> bitvector_entry_t;


private:

    // Position information (base pointer and offset) for the bitvector entry.
    bitvector_entry_t bv_entry;

    // Value read from the bitvector entry when the object is constructed; is immutable.
    const State_Read_Space state_read;

    // Value read from the bitvector entry when the object is constructed; is mutable.
    State_Read_Space state;


    // Constructs an API to the bitvector entry `bv_entry`.
    Kmer_Hash_Entry_API(const bitvector_entry_t& bv_entry):
        bv_entry(bv_entry), state_read(bv_entry)
    {
        state = state_read;
    }

    // Returns the state value read when the object was constructed.
    cuttlefish::state_code_t get_read_state() const
    {
        return state_read.get_state();
    }

    // Returns the value of the mutable state value wrapped inside the API,
    // i.e. the state value that had been read at the object creation, and then
    // possibly have been modified.
    cuttlefish::state_code_t get_current_state() const
    {
        return state.get_state();
    }


public:

    // Returns a reference to the mutable copy of the wrapped state value.
    State_Read_Space& get_state()
    {
        return state;
    }
};



#endif


#ifndef STATE_HPP
#define STATE_HPP



#include "globals.hpp"
#include "Vertex.hpp"
#include "compact_vector/compact_iterator.hpp"

#include <cstdint>
#include <cstdlib>
#include <iostream>


template <uint16_t k, uint8_t BITS_PER_KEY> class Kmer_Hash_Table;
template <uint8_t BITS_PER_KEY> class Kmer_Hash_Entry_API;


class State
{
    template <uint16_t k, uint8_t BITS_PER_KEY>
    friend class Kmer_Hash_Table;

    friend class Kmer_Hash_Entry_API<cuttlefish::BITS_PER_REF_KMER>;

    typedef compact::iterator_imp::lhs_setter<cuttlefish::state_code_t, cuttlefish::BITS_PER_REF_KMER, uint64_t, true, 64U> ref_bitvector_entry_t;

private:

    // The code of the state.
    cuttlefish::state_code_t code;


    // Constructs a `State` with the provided code `state`.
    State(cuttlefish::state_code_t code);

    // Constructs a `State` from the state stored at the bitvector entry `bv_entry`.
    State(const ref_bitvector_entry_t& bv_entry);

    // Sets the DNA base 2-bit encoding at the bits b1 and b0 of `code`.
    // Requirement: the two bits must be zero before the call, for consistent behavior.
    void set_nibble_lower_half(cuttlefish::base_t base);

    // Sets the DNA base 2-bit encoding at the bits b3 and b2 of `code`.
    // Requirement: the two bits must be zero before the call, for consistent behavior.
    void set_nibble_upper_half(cuttlefish::base_t base);

    // Returns the wrapped state code value.
    cuttlefish::state_code_t get_state() const;


public:

    // Constructs the state of an univisisted vertex.
    State();

    // Constructs the state of the provided `vertex`.
    State(const Vertex& vertex);

    // Returns a vertex by decoding the state code value.
    Vertex decode() const;

    // Returns true iff if this state can correspond to a vertex that has already
    // been visited by the compaction algorithm.
    bool is_visited() const;

    // Returns the state of the underlying vertex simulating if it had been outputted.
    State outputted() const;

    // Returns the class of this state.
    cuttlefish::State_Class state_class() const;

    // Returns the output status of a vertex that has this state.
    bool is_outputted() const;

    // Returns `true` iff this state qualifies as a dead-end in the state-transition model
    // during the states classification step of the algorithm.
    bool is_dead_end() const;

    // TODO?
    // const State& operator=();

    // Returns true iff the underlying codes are the same.
    bool operator==(const State& rhs) const;

    // For debugging.
    friend std::ostream& operator <<(std::ostream& out, const State& state);



private:
    // TODO:
    // Are lookup-tables faster than switch-cases? If yes, define the following arrays,
    // and modify the public function definitions.

    const static Vertex decoded_vertex[32];

    const static State outputted_vertex_state[32];

    const static cuttlefish::State_Class decoded_vertex_class[32];

    const static bool decoded_output_status[32];
};



inline State::State(const cuttlefish::state_code_t code):
    code(code)
{
    if(code == 0b00001 || code == 0b00010)
    {
        std::cerr << "Invalid state " << (uint16_t)code << " encountered during construction from code. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


inline State::State(const ref_bitvector_entry_t& bv_entry)
{
    // CAS vector `fetch` does not work.
    // bv_entry.fetch_val(vertex_code);

    code = bv_entry;

    if(code == 0b00001 || code == 0b00010)
    {
        std::cerr << "Invalid state " << (uint16_t)code << " encountered during construction from bitvector entry. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


inline cuttlefish::state_code_t State::get_state() const
{
    return code;
}


inline bool State::is_dead_end() const
{
    return is_visited() && state_class() == cuttlefish::State_Class::multi_in_multi_out;
}


inline bool State::operator==(const State& rhs) const
{
    return code == rhs.code;
}



#endif
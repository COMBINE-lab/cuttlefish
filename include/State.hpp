
#ifndef STATE_HPP
#define STATE_HPP



#include "Vertex.hpp"
#include "compact_vector/compact_iterator.hpp"

#include <cstdint>


class Kmer_Hash_Table;
class Kmer_Hash_Entry_API;


class State
{
    friend class Kmer_Hash_Table;
    friend class Kmer_Hash_Entry_API;

private:

    // The code of the state.
    cuttlefish::state_code_t code;


    // Constructs a `State` with the provided code `state`.
    State(const cuttlefish::state_code_t code);

    // Constructs a `State` from the state stored at the bitvector entry `bv_entry`.
    State(const cuttlefish::bitvector_entry_t& bv_entry);

    // Sets the nucleotide 2-bit encoding at the bits b1 and b0 of `code`.
    // Requirement: the two bits must be zero before the call, for consistent behavior.
    void set_nibble_lower_half(const cuttlefish::nucleotide_t nucl);

    // Sets the nucleotide 2-bit encoding at the bits b3 and b2 of `code`.
    // Requirement: the two bits must be zero before the call, for consistent behavior.
    void set_nibble_upper_half(const cuttlefish::nucleotide_t nucl);

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

    // Returns the class (or, type) that any vertex with this state has.
    cuttlefish::Vertex_Class vertex_class() const;

    // Returns the output status of a vertex that has this state.
    bool is_outputted() const;

    // Returns `true` iff this state qualifies as a dead-end in the state-transition model
    // during the states classification step of the algorithm.
    bool is_dead_end() const;

    // Returns true iff the underlying codes are different.
    bool operator!=(const State& rhs) const
    {
        return code != rhs.code;
    }

    // For debugging.
    friend std::ostream& operator <<(std::ostream& out, const State& state);



private:
    // TODO:
    // Are lookup-tables faster than switch-cases? If yes, define the following arrays,
    // and modify the public function definitions.

    const static Vertex decoded_vertex[32];

    const static State outputted_vertex_state[32];

    const static cuttlefish::Vertex_Class decoded_vertex_class[32];

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


inline State::State(const cuttlefish::bitvector_entry_t& bv_entry)
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
    return is_visited() && vertex_class() == cuttlefish::Vertex_Class::multi_in_multi_out;
}



#endif
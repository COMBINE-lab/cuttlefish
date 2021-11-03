
#ifndef STATE_READ_SPACE_HPP
#define STATE_READ_SPACE_HPP



#include "globals.hpp"


template <uint8_t BITS_PER_KEY> class Kmer_Hash_Entry_API;


// Class for a state in the state-space of the automata in read de Bruijn graphs.
class State_Read_Space
{
    friend class Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER>;

private:

    cuttlefish::state_code_t code;  // Numeric code of the state.

    static constexpr uint8_t BITS_PER_SIDE = 3; // Number of bits required to `Extended_Base`-encode edges incident to a side.
    static constexpr uint8_t FRONT_IDX = BITS_PER_SIDE; // Starting index of the three bits encoding the front-incident edge.
    static constexpr uint8_t BACK_IDX = 0;  // Starting index of the three bits encoding the back-incident edge.

    // Bitmask to extract the edge-encoding from some side. Has to be shifted to appropriate index before extraction.
    static constexpr uint8_t SIDE_MASK = (1 << BITS_PER_SIDE) - 1;

    // Bitmask used to extract the 'Extended_Base`-encoding of the edge(s) incident to the front side of a vertex.
    static constexpr cuttlefish::state_code_t FRONT_MASK = SIDE_MASK << FRONT_IDX;

    // Bitmask used to extract the 'Extended_Base`-encoding of the edge(s) incident to the back side of a vertex.
    static constexpr cuttlefish::state_code_t BACK_MASK = SIDE_MASK << BACK_IDX;


    // Constructs a state that wraps the provided numeric value `code`.
    State_Read_Space(cuttlefish::state_code_t code);

    // Sets the back-encoding of the state to the `Extended_Base`-encoding `edge`.
    void set_back_encoding(cuttlefish::edge_encoding_t edge);

    // Sets the front-encoding of the state to the `Extended_Base`-encoding `edge`.
    void set_front_encoding(cuttlefish::edge_encoding_t edge);


public:

    // Constructs the state of a vertex having both its sides unvisited.
    constexpr State_Read_Space();

    // Returns the wrapped state-code value.
    cuttlefish::state_code_t get_state() const;

    // Returns `true` iff some vertex having this state has been outputted.
    bool is_outputted() const;

    // Returns the `Extended_Base`-encoding of the edge(s) incident to the side
    // `side` of a vertex having this state.
    cuttlefish::edge_encoding_t edge_at(cuttlefish::side_t side) const;

    // Returns `true` iff some vertex having this state is branching (i.e. has
    // multiple incident edges) at its side `side`, and hasn't been outputted yet.
    bool is_branching_side(cuttlefish::side_t side) const;

    // Returns `true` iff some vertex having this state is branching (i.e. has
    // multiple incident edges) at its side `side`, and has already been outputted.
    bool was_branching_side(cuttlefish::side_t side) const;

    // Updates the `Extended_Base` encoding of the side `side` of this state, with
    // `edge`.
    void update_edge_at(cuttlefish::side_t side, cuttlefish::edge_encoding_t edge);

    // Marks the state as already been outputted.
    void mark_outputted();

    // Returns `true` iff the underlying code is the same as that one of `rhs`.
    bool operator==(const State_Read_Space& rhs) const;

    // Returns the state for the vertices that have been marked as outputted.
    static const State_Read_Space& get_outputted_state();

    // For the given code `code` of some state `s`, returns the code of the
    // state `s_op` which is the corresponding state where the vertices having
    // the DFA state `s` in the underlying graph transition to when outputted. 
    static cuttlefish::state_code_t mark_outputted(cuttlefish::state_code_t code);
};


inline constexpr State_Read_Space::State_Read_Space():
    code{(static_cast<cuttlefish::state_code_t>(cuttlefish::edge_encoding_t::E) << FRONT_IDX) | static_cast<cuttlefish::state_code_t>(cuttlefish::edge_encoding_t::E)}
{}


inline State_Read_Space::State_Read_Space(const cuttlefish::state_code_t code):
    code(code)
{}


inline void State_Read_Space::set_back_encoding(cuttlefish::edge_encoding_t edge)
{
    code = (code & FRONT_MASK) | (static_cast<cuttlefish::state_code_t>(edge) << BACK_IDX);
}


inline void State_Read_Space::set_front_encoding(cuttlefish::edge_encoding_t edge)
{
    code = (code & BACK_MASK) | (static_cast<cuttlefish::state_code_t>(edge) << FRONT_IDX);
}


inline cuttlefish::state_code_t State_Read_Space::get_state() const
{
    return code;
}


inline cuttlefish::edge_encoding_t State_Read_Space::edge_at(const cuttlefish::side_t side) const
{
    return static_cast<cuttlefish::edge_encoding_t>(side == cuttlefish::side_t::front ? (code & FRONT_MASK) >> FRONT_IDX : (code & BACK_MASK) >> BACK_IDX);
}


inline bool State_Read_Space::is_branching_side(const cuttlefish::side_t side) const
{
    return edge_at(side) == cuttlefish::edge_encoding_t::N;
}


inline bool State_Read_Space::was_branching_side(const cuttlefish::side_t side) const
{
    return edge_at(side) == cuttlefish::edge_encoding_t::OP_branching;
}


inline void State_Read_Space::update_edge_at(const cuttlefish::side_t side, const cuttlefish::edge_encoding_t edge)
{
    side == cuttlefish::side_t::front ? set_front_encoding(edge) : set_back_encoding(edge);
}


inline void State_Read_Space::mark_outputted()
{
    static constexpr cuttlefish::edge_encoding_t OP_non_branch = cuttlefish::edge_encoding_t::OP_non_branch;
    static constexpr cuttlefish::edge_encoding_t OP_branching = cuttlefish::edge_encoding_t::OP_branching;
    
    if(!is_outputted())
    {
        set_back_encoding(is_branching_side(cuttlefish::side_t::back) ? OP_branching : OP_non_branch);
        set_front_encoding(is_branching_side(cuttlefish::side_t::front) ? OP_branching : OP_non_branch);
    }
}


inline bool State_Read_Space::is_outputted() const
{
    static constexpr uint8_t OP_non_branch = static_cast<uint8_t>(cuttlefish::edge_encoding_t::OP_non_branch);
    static constexpr uint8_t OP_branching = static_cast<uint8_t>(cuttlefish::edge_encoding_t::OP_branching);

    return  code == ((OP_non_branch << FRONT_IDX)   |   (OP_non_branch << BACK_IDX))    ||
            code == ((OP_non_branch << FRONT_IDX)   |   (OP_branching << BACK_IDX))     ||
            code == ((OP_branching << FRONT_IDX)    |   (OP_non_branch << BACK_IDX))    ||
            code == ((OP_branching << FRONT_IDX)    |   (OP_branching << BACK_IDX)); 
}


inline bool State_Read_Space::operator==(const State_Read_Space& rhs) const
{
    return code == rhs.code;
}


inline cuttlefish::state_code_t State_Read_Space::mark_outputted(const cuttlefish::state_code_t code)
{
    State_Read_Space state(code);
    state.mark_outputted();

    return state.get_state();
}



#endif

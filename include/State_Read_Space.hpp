
#ifndef STATE_READ_SPACE_HPP
#define STATE_READ_SPACE_HPP



#include "globals.hpp"


template <uint8_t BITS_PER_KEY> class Kmer_Hash_Entry_API;


// Class for a state in the state-space of the automata in read de Bruijn graphs.
class State_Read_Space
{
    friend class Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER>;
    
    typedef DNA::Extended_Base edge_encoding_t; // TODO: replace this with `cuttlefish::edge_encoding_t`.

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

    // State code for vertices that has been outputted.
    // TODO: Use a well-thought-out value as the marker.
    static constexpr cuttlefish::state_code_t OUTPUTTED = static_cast<cuttlefish::state_code_t>((0b101 << FRONT_IDX) | 0b101 << BACK_IDX);


    // Constructs a state that wraps the provided numeric value `code`.
    State_Read_Space(cuttlefish::state_code_t code);

    // Sets the back-encoding of the state to the `Extended_Base`-encoding `edge`.
    // Requirement: except while for setting `Extended_Base::N`, the bits must be zero beforehand.
    void set_back_encoding(edge_encoding_t edge);

    // Sets the front-encoding of the state to the `Extended_Base`-encoding `edge`.
    // Requirement: except while for setting `Extended_Base::N`, the bits must be zero beforehand.
    void set_front_encoding(edge_encoding_t edge);

    // Returns the wrapped state-code value.
    cuttlefish::state_code_t get_state() const;


public:

    // Constructs the state of a vertex having both its sides unvisited.
    constexpr State_Read_Space();

    // Returns `true` iff some vertex having this state has been outputted.
    bool is_outputted() const;

    // Returns the `Extended_Base`-encoding of the edge(s) incident to the side
    // `side` of a vertex having this state.
    edge_encoding_t edge_at(cuttlefish::side_t side) const;

    // Updates the `Extended_Base` encoding of the side `side` of this state, with
    // `edge`. For optimization purposes, only certain edge-updates have defined
    // behavior: empty-to-rest and unique-to-multi.
    void update_edge_at(cuttlefish::side_t side, edge_encoding_t edge);

    // Returns `true` iff the underlying code is the same as that one of `rhs`.
    bool operator==(const State_Read_Space& rhs) const;
};


inline constexpr State_Read_Space::State_Read_Space():
    code{(static_cast<cuttlefish::state_code_t>(edge_encoding_t::E) << FRONT_IDX) | static_cast<cuttlefish::state_code_t>(edge_encoding_t::E)}
{}


inline State_Read_Space::State_Read_Space(const cuttlefish::state_code_t code):
    code(code)
{}


inline void State_Read_Space::set_back_encoding(edge_encoding_t edge)
{
    code |= (static_cast<cuttlefish::state_code_t>(edge) << BACK_IDX);
}


inline void State_Read_Space::set_front_encoding(edge_encoding_t edge)
{
    code |= (static_cast<cuttlefish::state_code_t>(edge) << FRONT_IDX);
}


inline cuttlefish::state_code_t State_Read_Space::get_state() const
{
    return code;
}


inline bool State_Read_Space::is_outputted() const
{
    return code == OUTPUTTED;
}


inline State_Read_Space::edge_encoding_t State_Read_Space::edge_at(const cuttlefish::side_t side) const
{
    return static_cast<edge_encoding_t>(side == cuttlefish::side_t::front ? (code & FRONT_MASK) >> FRONT_IDX : (code & BACK_MASK) >> BACK_IDX);
}


inline void State_Read_Space::update_edge_at(const cuttlefish::side_t side, const edge_encoding_t edge)
{
    side == cuttlefish::side_t::front ? set_front_encoding(edge) : set_back_encoding(edge);
}


inline bool State_Read_Space::operator==(const State_Read_Space& rhs) const
{
    return code == rhs.code;
}



#endif

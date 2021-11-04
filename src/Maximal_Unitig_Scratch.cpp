
#include "Maximal_Unitig_Scratch.hpp"


template <uint16_t k>
Maximal_Unitig_Scratch<k>::Maximal_Unitig_Scratch()
{}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Maximal_Unitig_Scratch)

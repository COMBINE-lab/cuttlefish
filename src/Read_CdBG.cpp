
#include "Read_CdBG.hpp"


template <uint16_t k>
Read_CdBG<k>::Read_CdBG(const Build_Params& params):
    params(params)
{}


template <uint16_t k>
void Read_CdBG<k>::construct()
{}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG)

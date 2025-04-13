
#include "State_Config.hpp"


namespace cuttlefish
{

uint8_t Edge_Frequency::f_th;

void Edge_Frequency::set_edge_threshold(const uint8_t f_th)
{
    assert(f_th <= max_f);
    Edge_Frequency::f_th = f_th;
}

}

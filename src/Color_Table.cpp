
#include "Color_Table.hpp"


namespace cuttlefish
{

Color_Table::Color_Table():
#ifndef USE_PARLAY_HASH
      M(map_sz_init)
#else
      M(map_sz_init, false)
#endif
{}

}

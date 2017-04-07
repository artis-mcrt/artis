#ifndef UPDATE_PACKETS_H
#define UPDATE_PACKETS_H

#include "update_grid.h"

void update_packets(int nts);

inline void update_cell(const int mgi)
///=calculate_levelpops for non isothermal homogeneous grids
///
{
  updatecellcounter++;

  cellhistory_reset(mgi, false);
}

#endif //UPDATE_PACKETS_H

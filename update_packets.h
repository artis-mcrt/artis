#ifndef UPDATE_PACKETS_H
#define UPDATE_PACKETS_H

#include "update_grid.h"

void update_packets(int nts);

inline void update_cell(const int cellnumber)
///=calculate_levelpops for non isothermal homogeneous grids
///
{
  updatecellcounter++;

  cellhistory_reset(cellnumber, true);
}

#endif //UPDATE_PACKETS_H

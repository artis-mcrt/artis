#ifndef UPDATE_GRID_H
#define UPDATE_GRID_H

#include <cstdio>
#include <ctime>

void update_grid(FILE *estimators_file, int nts, int nts_prev, int my_rank, int nstart, int ndo, int titer,
                 time_t real_time_start);
void cellhistory_reset(int modelgridindex);

#endif  // UPDATE_GRID_H

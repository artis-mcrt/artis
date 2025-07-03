#ifndef UPDATE_GRID_H
#define UPDATE_GRID_H

#include <ctime>
#include <fstream>

void update_grid(std::fstream &estimators_file, int nts, int nts_prev, int titer, time_t real_time_start);
void cellcache_change_cell(int nonemptymgi);

#endif  // UPDATE_GRID_H

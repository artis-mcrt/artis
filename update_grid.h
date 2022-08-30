#ifndef UPDATE_GRID_H
#define UPDATE_GRID_H

#include <cstdio>
#include <ctime>

#include "cuda.h"

void update_grid(FILE *estimators_file, int nts, int nts_prev, int my_rank, int nstart, int ndo, int titer,
                 const time_t real_time_start);
void precalculate_partfuncts(int modelgridindex);
__host__ __device__ void cellhistory_reset(int cellnumber, bool set_population);
double calculate_populations(int modelgridindex);
double calculate_electron_densities(int modelgridindex);

#endif  // UPDATE_GRID_H

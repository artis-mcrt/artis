#ifndef UPDATE_GRID_H
#define UPDATE_GRID_H

#include <stdio.h>
#include "sn3d.h"

void update_grid(FILE *estimators_file, int nts, int nts_prev, int my_rank, int nstart, int ndo, int titer);
void precalculate_partfuncts(int modelgridindex);
void cellhistory_validate_or_reset(int cellnumber, int timestep);
double calculate_populations(int modelgridindex);
double calculate_electron_densities(int modelgridindex);
void write_grid_restart_data(const int timestep);

inline double get_abundance(const int modelgridindex, const int element)
{
  return modelgrid[modelgridindex].composition[element].abundance;
}

#endif //UPDATE_GRID_H

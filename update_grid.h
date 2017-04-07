#ifndef UPDATE_GRID_H
#define UPDATE_GRID_H

#include "sn3d.h"

void update_grid(const int nts, const int my_rank, const int nstart, const int ndo, const int titer);
void precalculate_partfuncts(int modelgridindex);
void cellhistory_reset(const int cellnumber, const bool set_population);
double calculate_populations(int modelgridindex);
double calculate_electron_densities(int modelgridindex);
void write_grid_restart_data(void);

inline double get_abundance(int modelgridindex, int element)
{
  return modelgrid[modelgridindex].composition[element].abundance;
}

#endif //UPDATE_GRID_H

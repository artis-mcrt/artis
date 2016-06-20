#ifndef UPDATE_GRID_H
#define UPDATE_GRID_H

#include "sn3d.h"

void update_grid(int m, int my_rank, int nstart, int nblock, int titer);
void precalculate_partfuncts(int modelgridindex);
double calculate_populations(int modelgridindex, int first_nonempty_cell);
double calculate_electron_densities(int modelgridindex);
void write_grid_restart_data(void);

inline double get_abundance(int modelgridindex, int element)
{
  return modelgrid[modelgridindex].composition[element].abundance;
}

#endif //UPDATE_GRID_H

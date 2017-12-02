#ifndef UPDATE_GRID_H
#define UPDATE_GRID_H

void update_grid(int nts, int my_rank, int nstart, int ndo, int titer);
void precalculate_partfuncts(int modelgridindex);
void cellhistory_reset(int cellnumber, bool set_population);
double calculate_populations(int modelgridindex);
double calculate_electron_densities(int modelgridindex);
void write_grid_restart_data(void);

inline double get_abundance(int modelgridindex, int element)
{
  return modelgrid[modelgridindex].composition[element].abundance;
}

#endif //UPDATE_GRID_H

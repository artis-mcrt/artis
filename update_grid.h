#ifndef UPDATE_GRID_H
#define UPDATE_GRID_H

//void update_grid(FILE *estimators_file, int nts, int nts_prev, int my_rank, int nstart, int ndo, int titer);
void update_grid(FILE *estimators_file, int nts, int nts_prev, int my_rank, int *indices, int titer);
void precalculate_partfuncts(int modelgridindex);
void cellhistory_reset(int cellnumber, bool set_population);
double calculate_populations(int modelgridindex);
double calculate_electron_densities(int modelgridindex);
void write_grid_restart_data(const int timestep);

inline double get_abundance(const int modelgridindex, const int element)
{
  return modelgrid[modelgridindex].composition[element].abundance;
}

#endif //UPDATE_GRID_H

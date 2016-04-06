#ifndef UPDATE_GRID_H
#define UPDATE_GRID_H

int update_grid(int m, int my_rank, int nstart, int nblock, int titer);
int get_cell(double x, double y, double z, double t);
double get_abundance(int modelgridindex, int element);
void update_abundances(int modelgridindex, double t_current);
double calculate_populations(int modelgridindex, int first_nonempty_cell);
double calculate_electron_densities(int modelgridindex);
void precalculate_partfuncts(int modelgridindex);
void get_radfield_params(double J, double nuJ, int modelgridindex, double *T_J, double *T_R, double *W);
void write_grid_restart_data(void);

#endif //UPDATE_GRID_H

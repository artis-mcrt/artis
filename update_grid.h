#ifndef UPDATE_GRID_H
#define UPDATE_GRID_H

int update_grid(int m, int my_rank, int nstart, int nblock, int titer);
double get_abundance(int modelgridindex, int element);
double calculate_populations(int modelgridindex, int first_nonempty_cell);
double calculate_electron_densities(int modelgridindex);
void write_grid_restart_data(void);

#endif //UPDATE_GRID_H

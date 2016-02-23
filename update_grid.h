#ifndef UPDATE_GRID_H
  #define UPDATE_GRID_H

  int get_cell(double x, double y, double z, double t);
  double get_abundance(int modelgridindex, int element);
  double calculate_populations(int modelgridindex, int first_nonempty_cell);
  void precalculate_partfuncts(int modelgridindex);

#endif //UPDATE_GRID_H

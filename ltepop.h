#ifndef LTEPOP_H
  #define LTEPOP_H

  double nne_solution_f(double x, void *paras);
  double ionfract(int element, int ion, int modelgridindex, double nne);
  double phi(int element, int ion, int modelgridindex);
  double calculate_partfunct(int element, int ion, int modelgridindex);
  double get_groundlevelpop(int modelgridindex, int element, int ion);
  double calculate_exclevelpop(int modelgridindex, int element, int ion, int level);
  double ionstagepop(int modelgridindex, int element, int ion);
  void calculate_levelpops(int modelgridindex);
  double get_levelpop(int element, int ion, int level);
  double calculate_sahafact(int element, int ion, int level, int phixstargetindex, double T, double E_threshold);
  double get_sahafact(int element, int ion, int level, int phixstargetindex, double T, double E_threshold);
  void initialise_photoionestimators();

#endif //LTEPOP_H

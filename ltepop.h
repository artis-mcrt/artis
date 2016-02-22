#ifndef LTEPOP_H
  #define LTEPOP_H

  double get_groundlevelpop(int modelgridindex, int element, int ion);
  double calculate_exclevelpop(int modelgridindex, int element, int ion, int level);
  double ionstagepop(int modelgridindex, int element, int ion);
  void calculate_levelpops(int modelgridindex);
  double get_levelpop(int element, int ion, int level);
  double calculate_sahafact(int element, int ion, int level, int phixstargetindex, double T, double E_threshold);
  double get_sahafact(int element, int ion, int level, int phixstargetindex, double T, double E_threshold);
  double ionstagepop(int modelgridindex, int element, int ion);

#endif //LTEPOP_H

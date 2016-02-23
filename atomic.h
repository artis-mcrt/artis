#ifndef ATOMIC_H
  #define ATOMIC_H

  int get_element(int element);
  int get_elementindex(int Z);
  int get_nions(int element);
  int get_ionstage(int element, int ion);
  int get_nlevels(int element, int ion);
  int get_nlevels_nlte(int element, int ion);
  int get_ionisinglevels(int element, int ion);
  int get_bfcontinua(int element, int ion);
  double epsilon(int element, int ion, int level);
  double stat_weight(int element, int ion, int level);
  short is_nlte(int element, int ion, int level);
  int get_continuumindex(int element, int ion, int level);
  int get_nphixstargets(int element, int ion, int level);
  int get_phixsupperlevel(int element, int ion, int level, int phixstargetindex);
  float get_phixsprobability(int element, int ion, int level, int phixstargetindex);
  int transitioncheck(int upper, int lower);
  double einstein_spontaneous_emission(int lineindex);
  double osc_strength(int lineindex);
  double coll_str(int lineindex);
  double statw_up(int lineindex);
  double statw_down(int lineindex);
  double photoionization_crosssection(double nu_edge, double nu);
#endif //ATOMIC_H

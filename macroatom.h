#ifndef MACROATOM_H
  #define MACROATOM_H

  double rad_deexcitation(PKT *pkt_ptr, int lower, double epsilon_trans, double statweight_target, int lineindex, double t_current);
  double rad_recombination(int modelgridindex, int lower, double epsilon_trans);
  double rad_excitation(PKT *pkt_ptr, int upper, double epsilon_trans, double statweight_target, int lineindex, double t_current);//, double T_R, double W);
  double photoionization(int modelgridindex, int phixstargetindex, double epsilon_trans);

  double col_excitation(int modelgridindex, int upper, int lineindex, double epsilon_trans);
  double col_ionization(int modelgridindex, int phixstargetindex, double epsilon_trans);
  double col_deexcitation(int modelgridindex, int lower, double epsilon_trans, double statweight_target, int lineindex);
  double col_recombination(int modelgridindex, int lower, double epsilon_trans);

  double radfield(double nu, int modelgridindex);
  double radfield2(double nu, double T, double W);

  double get_individ_rad_deexc(int i);
  double get_individ_internal_down_same(int i);
  double get_individ_internal_up_same(int i);

#endif //MACROATOM_H

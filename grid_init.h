#ifndef GRIDINIT_H
  #define GRIDINIT_H

  int grid_init();
  double vol_init(CELL *grid_ptr);
  float get_rhoinit(int modelgridindex);
  float get_rho(int modelgridindex);
  float get_nne(int modelgridindex);
  float get_nnetot(int modelgridindex);
  float get_fni(int modelgridindex);
  float get_fco(int modelgridindex);
  float get_f52fe(int modelgridindex);
  float get_f48cr(int modelgridindex);
  float get_ffe(int modelgridindex);
  float get_fnistable(int modelgridindex);
  float get_fcostable(int modelgridindex);
  float get_ffestable(int modelgridindex);
  float get_fmnstable(int modelgridindex);
  float get_fcrstable(int modelgridindex);
  float get_fvstable(int modelgridindex);
  float get_ftirstable(int modelgridindex);
  float get_kappagrey(int modelgridindex);
  float get_Te(int modelgridindex);
  float get_TR(int modelgridindex);
  float get_TJ(int modelgridindex);
  float get_W(int modelgridindex);
  void set_rhoinit(int modelgridindex, float x);
  void set_rho(int modelgridindex, float x);
  void set_nne(int modelgridindex, float x);
  void set_nnetot(int modelgridindex, float x);
  void set_fni(int modelgridindex, float x);
  void set_fco(int modelgridindex, float x);
  void set_f52fe(int modelgridindex, float x);
  void set_f48cr(int modelgridindex, float x);
  void set_ffe(int modelgridindex, float x);
  void set_fnistable(int modelgridindex, float x);
  void set_fcostable(int modelgridindex, float x);
  void set_ffestable(int modelgridindex, float x);
  void set_fmnstable(int modelgridindex, float x);
  void set_fcrstable(int modelgridindex, float x);
  void set_fvstable(int modelgridindex, float x);
  void set_ftistable(int modelgridindex, float x);
  void set_kappagrey(int modelgridindex, float x);
  void set_Te(int modelgridindex, float x);
  void set_TR(int modelgridindex, float x);
  void set_TJ(int modelgridindex, float x);
  void set_W(int modelgridindex, float x);

#endif //GRIDINIT_H

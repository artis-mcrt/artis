#ifndef LIGHT_CURVE_H
  #define LIGHT_CURVE_H

  void init_light_curve(void);
  int write_light_curve(FILE *lc_file, int current_abin);
  int gather_light_curve(void);
  int gather_light_curve_res(int current_abin);

#endif //LIGHT_CURVE_H

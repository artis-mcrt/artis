#ifndef LIGHT_CURVE_H
#define LIGHT_CURVE_H

#include "exspec.h"

void init_light_curve(void);
int write_light_curve(FILE *lc_file, int current_abin);
int gather_light_curve(void);
int add_to_lc(const EPKT *pkt_ptr);
int gather_light_curve_res(int current_abin);
int add_to_lc_res(const EPKT *pkt_ptr, int current_abin);

#endif //LIGHT_CURVE_H

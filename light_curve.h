#ifndef LIGHT_CURVE_H
#define LIGHT_CURVE_H

#include <stdio.h>

void init_light_curve(void);
void write_light_curve(FILE *lc_file, int current_abin);
void gather_light_curve(void);
void gather_light_curve_res(int current_abin);

#endif //LIGHT_CURVE_H

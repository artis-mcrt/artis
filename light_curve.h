#ifndef LIGHT_CURVE_H
#define LIGHT_CURVE_H

#include <stdio.h>

void init_light_curve(void);
void write_light_curve(char lc_filename[], int current_abin);
void gather_light_curve(EPKT *epkts, int nepkts);
void gather_light_curve_res(EPKT *epkts, int nepkts, int current_abin);

#endif //LIGHT_CURVE_H

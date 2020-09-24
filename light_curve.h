#ifndef LIGHT_CURVE_H
#define LIGHT_CURVE_H

#include <stdio.h>
#include "exspec.h"

void gather_light_curve(EPKT *epkts, int nepkts, double *light_curve_lum, double *light_curve_lumcmf);

void gather_light_curve_res(EPKT *epkts, int nepkts, int current_abin, double *light_curve_lum);

void write_light_curve(char lc_filename[], int current_abin, const double *light_curve_lum,
const double *light_curve_lumcmf);

#endif //LIGHT_CURVE_H

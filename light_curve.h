#ifndef LIGHT_CURVE_H
#define LIGHT_CURVE_H

#include <stdio.h>
#include "exspec.h"

void add_to_lc_res(const PKT *pkt_ptr, int current_abin, double *light_curve_lum, double *light_curve_lumcmf);

void write_light_curve(char lc_filename[], int current_abin, const double *light_curve_lum,
const double *light_curve_lumcmf);

#endif //LIGHT_CURVE_H

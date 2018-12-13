#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <stdio.h>

void write_spectrum(char spec_filename[], bool do_emission_res, char emission_filename[], char trueemission_filename[], char absorption_filename[]);
void gather_spectrum(EPKT *epkts, int nepkts, int depth, bool do_emission_res);
void gather_spectrum_res(EPKT *epkts, int nepkts, int current_abin);
void init_spectrum(void);

#endif //SPECTRUM_H

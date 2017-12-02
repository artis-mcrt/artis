#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <stdio.h>

void write_spectrum(FILE *spec_file, FILE *emission_file, FILE *trueemission_file, FILE *absorption_file);
void gather_spectrum(int depth);
void gather_spectrum_res(int current_abin);

#endif //SPECTRUM_H

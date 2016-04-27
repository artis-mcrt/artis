#ifndef SPECTRUM_H
#define SPECTRUM_H

#include "exspec.h"
#include <stdio.h>

int write_spectrum(FILE *spec_file, FILE *emission_file, FILE *absorption_file);
int gather_spectrum(int depth);
int gather_spectrum_res(int current_abin);

#endif //SPECTRUM_H

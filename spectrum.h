#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <stdio.h>
#include "exspec.h"
#include "types.h"

void write_spectrum(char spec_filename[], bool do_emission_res, char emission_filename[], char trueemission_filename[], char absorption_filename[], struct spec *spectra);
void write_specpol(char spec_filename[], bool do_emission_res, char emission_filename[], char absorption_filename[]);
void add_to_spec_res(const PKT *const pkt_ptr, int current_abin, const bool do_emission_res, struct spec *spectra, const double nu_min, const double nu_max);
void init_spectrum(struct spec *spectra, const double nu_min, const double nu_max, const bool do_emission_res);

#endif //SPECTRUM_H

#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <stdio.h>
#include "exspec.h"

void write_spectrum(char spec_filename[], bool do_emission_res, char emission_filename[], char trueemission_filename[], char absorption_filename[]);
void add_to_spec(const PKT *const pkt_ptr, const bool do_emission_res);
void add_to_spec_res(const PKT *const pkt_ptr, int current_abin);
void init_spectrum(void);

#endif //SPECTRUM_H

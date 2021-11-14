#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <cstdio>
#include "exspec.h"
#include "types.h"

void write_spectrum(
  char spec_filename[], char emission_filename[], char trueemission_filename[],
  char absorption_filename[], struct spec *spectra);

void write_specpol(
  char spec_filename[], char emission_filename[], char absorption_filename[],
  struct spec *stokes_i, struct spec *stokes_q, struct spec *stokes_u);

void add_to_spec_res(
  const PKT *const pkt_ptr, int current_abin, struct spec *spectra,
  struct spec *stokes_i, struct spec *stokes_q, struct spec *stokes_u);

struct spec *alloc_spectra(const bool do_emission_res);
void init_spectra(struct spec *spectra, const double nu_min, const double nu_max, const bool do_emission_res);
void init_spectrum_trace(void);
void free_spectra(struct spec *spectra);
void write_partial_lightcurve_spectra(int my_rank, int nts, PKT *pkts);

#endif //SPECTRUM_H

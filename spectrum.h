#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <string>

struct timestepspec {
  double *flux = nullptr;
  double *absorption = nullptr;
  double *emission = nullptr;
  double *trueemission = nullptr;
};

struct spec {
  double nu_min = -1.;
  double nu_max = -1.;
  float *lower_freq = nullptr;
  float *delta_freq = nullptr;
  double *fluxalltimesteps = nullptr;
  double *absorptionalltimesteps = nullptr;
  double *emissionalltimesteps = nullptr;
  double *trueemissionalltimesteps = nullptr;
  struct timestepspec *timesteps = nullptr;
  bool do_emission_res = true;
};

void write_spectrum(const std::string &spec_filename, const char *emission_filename, const char *trueemission_filename,
                    const char *absorption_filename, struct spec *spectra, int num_timesteps);

void write_specpol(const std::string &specpol_filename, const std::string &emission_filename,
                   const std::string &absorption_filename, struct spec *stokes_i, struct spec *stokes_q,
                   struct spec *stokes_u);

void add_to_spec_res(const struct packet *const pkt_ptr, int current_abin, struct spec *spectra, struct spec *stokes_i,
                     struct spec *stokes_q, struct spec *stokes_u);

struct spec *alloc_spectra(const bool do_emission_res);
void init_spectra(struct spec *spectra, const double nu_min, const double nu_max, const bool do_emission_res);
void init_spectrum_trace();
void free_spectra(struct spec *spectra);
void write_partial_lightcurve_spectra(int my_rank, int nts, struct packet *pkts);

#endif  // SPECTRUM_H

#pragma once
#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <string>
#include <vector>

struct TimeStepstepspec {
  double *flux = nullptr;
  double *absorption = nullptr;
  double *emission = nullptr;
  double *trueemission = nullptr;
};

struct Spectra {
  double nu_min = -1.;
  double nu_max = -1.;
  std::vector<float> lower_freq;
  std::vector<float> delta_freq;
  std::vector<double> fluxalltimesteps;
  std::vector<double> absorptionalltimesteps;
  std::vector<double> emissionalltimesteps;
  std::vector<double> trueemissionalltimesteps;
  std::vector<struct TimeStepstepspec> timesteps;
  bool do_emission_res = false;
};

void write_spectrum(const std::string &spec_filename, const std::string &emission_filename,
                    const std::string &trueemission_filename, const std::string &absorption_filename,
                    const struct Spectra &spectra, int numtimesteps);

void write_specpol(const std::string &specpol_filename, const std::string &emission_filename,
                   const std::string &absorption_filename, const struct Spectra *stokes_i,
                   const struct Spectra *stokes_q, const struct Spectra *stokes_u);

void add_to_spec_res(const struct Packet &pkt_ptr, int current_abin, struct Spectra &spectra,
                     const struct Spectra *stokes_i, const struct Spectra *stokes_q, const struct Spectra *stokes_u);

void init_spectra(struct Spectra &spectra, double nu_min, double nu_max, bool do_emission_res);
void init_spectrum_trace();
void write_partial_lightcurve_spectra(int my_rank, int nts, struct Packet *pkts);

#endif  // SPECTRUM_H

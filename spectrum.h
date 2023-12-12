#pragma once
#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <string>
#include <vector>

struct timestepspec {
  double *flux = nullptr;
  double *absorption = nullptr;
  double *emission = nullptr;
  double *trueemission = nullptr;
};

struct spec {
  double nu_min = -1.;
  double nu_max = -1.;
  std::vector<float> lower_freq;
  std::vector<float> delta_freq;
  std::vector<double> fluxalltimesteps;
  std::vector<double> absorptionalltimesteps;
  std::vector<double> emissionalltimesteps;
  std::vector<double> trueemissionalltimesteps;
  std::vector<struct timestepspec> timesteps;
  bool do_emission_res = false;
};

void write_spectrum(const std::string &spec_filename, const std::string &emission_filename,
                    const std::string &trueemission_filename, const std::string &absorption_filename,
                    const struct spec &spectra, int numtimesteps);

void write_specpol(const std::string &specpol_filename, const std::string &emission_filename,
                   const std::string &absorption_filename, const struct spec *stokes_i, const struct spec *stokes_q,
                   const struct spec *stokes_u);

void add_to_spec_res(const struct packet *pkt_ptr, int current_abin, struct spec &spectra, const struct spec *stokes_i,
                     const struct spec *stokes_q, const struct spec *stokes_u);

void init_spectra(struct spec &spectra, double nu_min, double nu_max, bool do_emission_res);
void init_spectrum_trace();
void write_partial_lightcurve_spectra(int my_rank, int nts, struct packet *pkts);

#endif  // SPECTRUM_H

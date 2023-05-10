#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <memory>
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
  std::unique_ptr<float[]> lower_freq = nullptr;
  std::unique_ptr<float[]> delta_freq = nullptr;
  std::unique_ptr<double[]> fluxalltimesteps = nullptr;
  std::unique_ptr<double[]> absorptionalltimesteps = nullptr;
  std::unique_ptr<double[]> emissionalltimesteps = nullptr;
  std::unique_ptr<double[]> trueemissionalltimesteps = nullptr;
  std::unique_ptr<struct timestepspec[]> timesteps = nullptr;
  bool do_emission_res = true;
};

void write_spectrum(const std::string &spec_filename, const char *emission_filename, const char *trueemission_filename,
                    const char *absorption_filename, struct spec *spectra, int numtimesteps);

void write_specpol(const std::string &specpol_filename, const std::string &emission_filename,
                   const std::string &absorption_filename, struct spec *stokes_i, struct spec *stokes_q,
                   struct spec *stokes_u);

void add_to_spec_res(const struct packet *pkt_ptr, int current_abin, struct spec *spectra, struct spec *stokes_i,
                     struct spec *stokes_q, struct spec *stokes_u);

struct spec *alloc_spectra(bool do_emission_res);
void init_spectra(struct spec *spectra, double nu_min, double nu_max, bool do_emission_res);
void init_spectrum_trace();
void free_spectra(struct spec *spectra);
void write_partial_lightcurve_spectra(int my_rank, int nts, struct packet *pkts);

#endif  // SPECTRUM_H

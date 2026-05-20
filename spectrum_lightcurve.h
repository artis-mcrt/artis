#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <array>
#include <cstddef>
#include <span>
#include <string>

#include "exspec.h"
#include "mpi_logging.h"
#include "packet.h"

struct Spectra {
  double nu_min = -1.;
  double nu_max = -1.;
  std::array<float, MNUBINS> lower_freq{};
  std::array<float, MNUBINS> delta_freq{};

  MPI_shared_array<double> fluxalltimesteps;
  MPI_shared_array<double> absorptionalltimesteps;
  MPI_shared_array<double> emissionalltimesteps;
  MPI_shared_array<double> trueemissionalltimesteps;

  bool do_emission_absorption = false;

  [[nodiscard]] auto mem_usage_bytes() const -> size_t {
    auto mem_usage = sizeof(Spectra);
    mem_usage += sizeof(float) * (lower_freq.size() + delta_freq.size());
    mem_usage += sizeof(double) * (fluxalltimesteps.size() + absorptionalltimesteps.size() +
                                   emissionalltimesteps.size() + trueemissionalltimesteps.size());
    // Note: Allocator overhead is not included in this calculation.
    return mem_usage;
  }
};

void write_spectra(const std::string& spec_filename, const std::string& emission_filename,
                   const std::string& trueemission_filename, const std::string& absorption_filename,
                   const Spectra& spectra, int numtimesteps);

void write_specpol(const std::string& specpol_filename, const std::string& emission_filename,
                   const std::string& absorption_filename, const Spectra& spectra_I, const Spectra& spectra_Q,
                   const Spectra& spectra_U);

void add_to_spec_res(const Packet& pkt, int dirbin, Spectra& spectra_I, Spectra* spectra_Q, Spectra* spectra_U);

void init_spectra(Spectra& spectra, double nu_min, double nu_max, bool do_emission_absorption);
void init_spectrum_trace();
void write_partial_lightcurve_spectra(int nts, std::span<const Packet> pkts);

void add_to_lc_res(const Packet& pkt, int dirbin, std::span<double> light_curve_lum,
                   std::span<double> light_curve_lumcmf);

void write_light_curve(const std::string& lc_filename, std::span<const double> light_curve_lum,
                       std::span<const double> light_curve_lumcmf, int numtimesteps);

#endif  // SPECTRUM_H

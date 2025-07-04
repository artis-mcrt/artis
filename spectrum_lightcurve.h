#ifndef SPECTRUM_H
#define SPECTRUM_H

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include <array>
#include <cstddef>
#include <span>
#include <string>
#include <vector>

#include "exspec.h"
#include "globals.h"
#include "packet.h"

struct Spectra {
  double nu_min = -1.;
  double nu_max = -1.;
  std::array<float, MNUBINS> lower_freq{};
  std::array<float, MNUBINS> delta_freq{};

  std::span<double> fluxalltimesteps;
  MPI_Win win_fluxalltimesteps{MPI_WIN_NULL};

  std::span<double> absorptionalltimesteps;
  MPI_Win win_absorptionalltimesteps{MPI_WIN_NULL};

  std::span<double> emissionalltimesteps;
  MPI_Win win_emissionalltimesteps{MPI_WIN_NULL};

  std::span<double> trueemissionalltimesteps;
  MPI_Win win_trueemissionalltimesteps{MPI_WIN_NULL};

  bool do_emission_absorption = false;

  [[nodiscard]] auto mem_usage_bytes() const -> size_t {
    auto mem_usage = sizeof(Spectra);
    mem_usage += sizeof(float) * (lower_freq.size() + delta_freq.size());
    mem_usage += sizeof(double) * (fluxalltimesteps.size() + absorptionalltimesteps.size() +
                                   emissionalltimesteps.size() + trueemissionalltimesteps.size());
    // Note: Allocator overhead is not included in this calculation.
    return mem_usage;
  }

  auto operator=(Spectra &&) -> Spectra & = delete;
  auto operator=(const Spectra &) -> Spectra & = delete;
  Spectra(Spectra &&) = delete;
  Spectra(const Spectra &) = delete;

  // constructor to allocate MPI windows
  Spectra() = default;

  // destructor to free MPI windows
  ~Spectra() {
    if (globals::mpi_finalized) {
      return;  // do not free MPI windows after MPI_Finalize
    }
    if (win_fluxalltimesteps != MPI_WIN_NULL) {
      fluxalltimesteps = {};
      MPI_Win_free(&win_fluxalltimesteps);
      win_fluxalltimesteps = MPI_WIN_NULL;
    }
    if (win_absorptionalltimesteps != MPI_WIN_NULL) {
      absorptionalltimesteps = {};
      MPI_Win_free(&win_absorptionalltimesteps);
      win_absorptionalltimesteps = MPI_WIN_NULL;
    }
    if (win_emissionalltimesteps != MPI_WIN_NULL) {
      emissionalltimesteps = {};
      MPI_Win_free(&win_emissionalltimesteps);
      win_emissionalltimesteps = MPI_WIN_NULL;
    }
    if (win_trueemissionalltimesteps != MPI_WIN_NULL) {
      trueemissionalltimesteps = {};
      MPI_Win_free(&win_trueemissionalltimesteps);
      win_trueemissionalltimesteps = MPI_WIN_NULL;
    }
  }
};

void write_spectrum(const std::string &spec_filename, const std::string &emission_filename,
                    const std::string &trueemission_filename, const std::string &absorption_filename,
                    const Spectra &spectra, int numtimesteps);

void write_specpol(const std::string &specpol_filename, const std::string &emission_filename,
                   const std::string &absorption_filename, const Spectra *stokes_i, const Spectra *stokes_q,
                   const Spectra *stokes_u);

void add_to_spec_res(const Packet &pkt, int dirbin, Spectra &spectra, Spectra *stokes_i, Spectra *stokes_q,
                     Spectra *stokes_u);

void init_spectra(Spectra &spectra, double nu_min, double nu_max, bool do_emission_absorption);
void init_spectrum_trace();
void write_partial_lightcurve_spectra(int nts, std::span<const Packet> pkts);

void add_to_lc_res(const Packet &pkt, int dirbin, std::span<double> light_curve_lum,
                   std::span<double> light_curve_lumcmf);

void write_light_curve(const std::string &lc_filename, int dirbin, const std::vector<double> &light_curve_lum,
                       const std::vector<double> &light_curve_lumcmf, int numtimesteps);

#endif  // SPECTRUM_H

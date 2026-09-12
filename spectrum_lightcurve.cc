// Binning of escaped packets into spectra and light curves for each observer direction bin,
// optionally decomposed into the emission and absorption contributions of each atomic process.

#include "spectrum_lightcurve.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <format>
#include <ios>
#include <iterator>
#include <print>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "exspec.h"
#include "globals.h"
#include "grid.h"
#include "mpi_logging.h"
#include "packet.h"
#include "sn3d.h"
#include "vectors.h"

namespace {

struct Spectra {
  double dlognu = -1.;
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
    auto mem_usage = sizeof(Spectra);  // includes the inline lower_freq and delta_freq arrays
    mem_usage += sizeof(double) * (fluxalltimesteps.size() + absorptionalltimesteps.size() +
                                   emissionalltimesteps.size() + trueemissionalltimesteps.size());
    // Note: Allocator overhead is not included in this calculation.
    return mem_usage;
  }
};

Spectra rpkt_spectra_I;
Spectra rpkt_spectra_Q;
Spectra rpkt_spectra_U;
Spectra gamma_spectra;

// the other "atomicadd" function is atomic only for multithreaded modes (STDPAR or OpenMP), but here we need it to be
// atomic for node-shared memory between processes even in single-threaded mode
template <typename T, typename U>
constexpr void atomicadd_always(T& var, U&& val) {
  std::atomic_ref<T>(var).fetch_add(std::forward<U>(val), std::memory_order_relaxed);
}

// number of different emission processes (bf and bb for each ion, and free-free)
auto get_proccount() -> int { return (2 * get_nelements() * get_max_nions()) + 1; }

auto columnindex_from_emissiontype(const int et) -> int {
  if (et >= 0) {
    // bb-emission
    const int element = globals::linelist.elementindex[et];
    const int ion = globals::linelist.ionindex[et];
    return (element * get_max_nions()) + ion;
  }
  if (et == EMTYPE_FREEFREE) {
    // ff-emission

    const int contindex = get_bflistindex_from_emtype_continuum(et);
    assert_always(contindex >= globals::nbfcontinua);  // make sure the special value didn't collide with a real process

    return 2 * get_nelements() * get_max_nions();
  }
  if (et == EMTYPE_NOTSET) {
    return -1;
  }
  // bf-emission
  const int bfindex = get_bflistindex_from_emtype_continuum(et);
  if (globals::nbfcontinua == 0) {
    // no bf continua are in use, so a bf emission type should be impossible; count it in the free-free column
    return 2 * get_nelements() * get_max_nions();
  }
  assert_always(bfindex < globals::nbfcontinua);
  const int element = globals::bflist[bfindex].elementindex;
  const int ion = globals::bflist[bfindex].ionindex;
  const int level = globals::bflist[bfindex].levelindex;
  const int phixstargetindex = globals::bflist[bfindex].phixstargetindex;

  assert_always(get_emtype_continuum(element, ion, level, phixstargetindex) == et);

  return (get_nelements() * get_max_nions()) + (element * get_max_nions()) + ion;
}

[[nodiscard]] auto get_absorption_spectrum_index(const ptrdiff_t nts, const ptrdiff_t nnu_abs) -> ptrdiff_t {
  const ptrdiff_t nelements = get_nelements();
  const ptrdiff_t max_nions = get_max_nions();
  return (nnu_abs * globals::ntimesteps * nelements * max_nions) + (nts * nelements * max_nions);
}

[[nodiscard]] inline auto get_timestep(const double time) -> int {
  assert_always(time >= globals::tmin);
  assert_always(time < globals::tmax);
  for (int nts = 0; nts < globals::ntimesteps; nts++) {
    const double tsend = (nts < (globals::ntimesteps - 1)) ? globals::timesteps[nts + 1].start : globals::tmax;
    if (time >= globals::timesteps[nts].start && time < tsend) {
      return nts;
    }
  }
  assert_always(false);  // could not find matching timestep

  return -1;
}

// Only one rank writes each file. Different file numbers go to different ranks on different nodes, if available.
auto this_rank_writes_file(const int filenum) -> bool {
  return (filenum % globals::node_count == globals::node_id) &&
         (filenum % globals::node_nprocs == globals::rank_in_node);
}

void write_spectrum_file(const std::string& spec_filename, const Spectra& spectra, const int numtimesteps) {
  auto spec_file = fstream_required(spec_filename, std::ios::out | std::ios::trunc);
  std::print(spec_file, "0 ");
  for (int p = 0; p < numtimesteps; p++) {
    std::print(spec_file, "{:g} ", globals::timesteps[p].mid / DAY);
  }
  std::println(spec_file, "");

  const auto ntimesteps_all = static_cast<ptrdiff_t>(globals::ntimesteps);
  for (auto nubin = 0Z; nubin < MNUBINS; nubin++) {
    std::print(spec_file, "{:g} ", (spectra.lower_freq[nubin] + (spectra.delta_freq[nubin] / 2)));

    for (auto nts = 0Z; nts < numtimesteps; nts++) {
      std::print(spec_file, "{:g} ", spectra.fluxalltimesteps[(nubin * ntimesteps_all) + nts]);
    }
    std::println(spec_file, "");
  }
}

// Write an emission-type spectrum (emission or true emission) with a line for each frequency bin of
// each timestep, holding one column per emission process (see get_proccount).
void write_emission_spectrum_file(const std::string& emission_filename,
                                  const std::span<const double> emission_alltimesteps, const int numtimesteps) {
  assert_always(numtimesteps <= globals::ntimesteps);
  assert_always(!emission_filename.empty());
  auto emission_file = fstream_required(emission_filename, std::ios::out | std::ios::trunc);
  const auto ntimesteps_all = static_cast<ptrdiff_t>(globals::ntimesteps);
  const auto proccount = static_cast<ptrdiff_t>(get_proccount());
  for (auto nubin = 0Z; nubin < MNUBINS; nubin++) {
    for (auto nts = 0Z; nts < numtimesteps; nts++) {
      const auto emindex_nts_nubin = (nubin * ntimesteps_all * proccount) + (nts * proccount);
      for (int nproc = 0; nproc < proccount; nproc++) {
        std::print(emission_file, "{:g} ", emission_alltimesteps[emindex_nts_nubin + nproc]);
      }
      std::println(emission_file, "");
    }
  }
}

void write_absorption_spectrum_file(const std::string& absorption_filename, const Spectra& spectra,
                                    const int numtimesteps) {
  assert_always(numtimesteps <= globals::ntimesteps);
  assert_always(!absorption_filename.empty());
  auto absorption_file = fstream_required(absorption_filename, std::ios::out | std::ios::trunc);
  const int ioncount = get_nelements() * get_max_nions();  // may be higher than the true included ion count
  for (auto nubin = 0Z; nubin < MNUBINS; nubin++) {
    for (auto nts = 0Z; nts < numtimesteps; nts++) {
      for (int i = 0; i < ioncount; i++) {
        std::print(absorption_file, "{:g} ",
                   spectra.absorptionalltimesteps[get_absorption_spectrum_index(nts, nubin) + i]);
      }
      std::println(absorption_file, "");
    }
  }
}

void write_spectra(const std::string& spec_filename, const std::string& emission_filename,
                   const std::string& trueemission_filename, const std::string& absorption_filename,
                   const Spectra& spectra, const int numtimesteps) {
  assert_always(numtimesteps <= globals::ntimesteps);

  if (this_rank_writes_file(1)) {
    write_spectrum_file(spec_filename, spectra, numtimesteps);
  }

  if (spectra.do_emission_absorption) {
    if (this_rank_writes_file(2)) {
      write_emission_spectrum_file(emission_filename, spectra.emissionalltimesteps.span(), numtimesteps);
    }

    if (this_rank_writes_file(3)) {
      write_emission_spectrum_file(trueemission_filename, spectra.trueemissionalltimesteps.span(), numtimesteps);
    }

    if (this_rank_writes_file(4)) {
      write_absorption_spectrum_file(absorption_filename, spectra, numtimesteps);
    }
  }
}

void write_specpol(const std::string& specpol_filename, const std::string& emission_filename,
                   const std::string& absorption_filename, const Spectra& spectra_I, const Spectra& spectra_Q,
                   const Spectra& spectra_U, const int numtimesteps) {
  assert_always(numtimesteps <= globals::ntimesteps);
  assert_always(std::ssize(spectra_I.delta_freq) == MNUBINS);
  assert_always(std::ssize(spectra_I.lower_freq) == MNUBINS);
  const auto stokes_spectra = {&spectra_I, &spectra_Q, &spectra_U};
  const auto ntimesteps_all = static_cast<ptrdiff_t>(globals::ntimesteps);

  if (this_rank_writes_file(5)) {
    printlnlog("Writing {}", specpol_filename);
    auto specpol_file = fstream_required(specpol_filename, std::ios::out | std::ios::trunc);
    std::print(specpol_file, "{:g}", 0.0);
    for (int l = 0; l < 3; l++) {
      for (int p = 0; p < numtimesteps; p++) {
        std::print(specpol_file, " {:g}", globals::timesteps[p].mid / DAY);
      }
    }
    std::println(specpol_file, "");

    for (auto nnu = 0Z; nnu < MNUBINS; nnu++) {
      std::print(specpol_file, "{:g}", (spectra_I.lower_freq[nnu] + (spectra_I.delta_freq[nnu] / 2)));
      for (const auto* spec : stokes_spectra) {
        for (auto nts = 0Z; nts < numtimesteps; nts++) {
          std::print(specpol_file, " {:g}", spec->fluxalltimesteps[(nnu * ntimesteps_all) + nts]);
        }
      }
      std::println(specpol_file, "");
    }
  }

  if (!spectra_I.do_emission_absorption) {
    return;
  }

  if (this_rank_writes_file(6)) {
    printlnlog("Writing {}", emission_filename);
    auto emissionpol_file = fstream_required(emission_filename, std::ios::out | std::ios::trunc);
    const auto proccount = static_cast<ptrdiff_t>(get_proccount());
    for (auto nnu = 0Z; nnu < MNUBINS; nnu++) {
      for (const auto* spec : stokes_spectra) {
        for (auto nts = 0Z; nts < numtimesteps; nts++) {
          for (auto nproc = 0Z; nproc < proccount; nproc++) {
            const auto emindex = (nnu * ntimesteps_all * proccount) + (nts * proccount) + nproc;
            if (nproc > 0) {
              std::print(emissionpol_file, " ");
            }
            std::print(emissionpol_file, "{:g}", spec->emissionalltimesteps[emindex]);
          }
          std::println(emissionpol_file, "");
        }
      }
    }
  }

  if (this_rank_writes_file(7)) {
    printlnlog("Writing {}", absorption_filename);
    auto absorptionpol_file = fstream_required(absorption_filename, std::ios::out | std::ios::trunc);
    const int ioncount = get_nelements() * get_max_nions();  // may be higher than the true included ion count
    for (auto nnu = 0Z; nnu < MNUBINS; nnu++) {
      for (const auto* spec : stokes_spectra) {
        for (auto nts = 0Z; nts < numtimesteps; nts++) {
          for (int i = 0; i < ioncount; i++) {
            if (i > 0) {
              std::print(absorptionpol_file, " ");
            }
            std::print(absorptionpol_file, "{:g}",
                       spec->absorptionalltimesteps[get_absorption_spectrum_index(nts, nnu) + i]);
          }
          std::println(absorptionpol_file, "");
        }
      }
    }
  }
}

// resize and initialize the spectra object
void init_spectra(Spectra& spectra, const double nu_min, const double nu_max, const bool do_emission_absorption) {
  // setup the time and frequency bins using a logarithmic spacing in both t and nu

  assert_always(MNUBINS > 0);
  const double dlognu = (log(nu_max) - log(nu_min)) / MNUBINS;

  spectra.dlognu = dlognu;
  spectra.nu_min = nu_min;
  spectra.nu_max = nu_max;
  spectra.do_emission_absorption = do_emission_absorption;
  const bool print_memusage =
      (spectra.fluxalltimesteps.empty() || (do_emission_absorption && spectra.absorptionalltimesteps.empty()));

  for (auto nnu = 0Z; nnu < MNUBINS; nnu++) {
    spectra.lower_freq[nnu] = static_cast<float>(get_loggrid_edge(nu_min, dlognu, static_cast<double>(nnu)));
    spectra.delta_freq[nnu] =
        static_cast<float>(get_loggrid_edge(nu_min, dlognu, static_cast<double>(nnu + 1)) - spectra.lower_freq[nnu]);
  }

  if (spectra.fluxalltimesteps.empty()) {
    spectra.fluxalltimesteps = MPI_shared_array<double>(globals::ntimesteps * MNUBINS);
  }
  assert_always(std::ssize(spectra.fluxalltimesteps) == globals::ntimesteps * MNUBINS);
  std::ranges::fill(spectra.fluxalltimesteps, 0.0);
  MPI_Barrier_node();

  if (do_emission_absorption) {
    if (spectra.absorptionalltimesteps.empty()) {
      spectra.absorptionalltimesteps =
          MPI_shared_array<double>(globals::ntimesteps * MNUBINS * get_nelements() * get_max_nions());
    }
    assert_always(std::ssize(spectra.absorptionalltimesteps) ==
                  globals::ntimesteps * MNUBINS * get_nelements() * get_max_nions());
    std::ranges::fill(spectra.absorptionalltimesteps, 0.0);

    if (spectra.emissionalltimesteps.empty()) {
      spectra.emissionalltimesteps = MPI_shared_array<double>(globals::ntimesteps * MNUBINS * get_proccount());
    }
    assert_always(std::ssize(spectra.emissionalltimesteps) == globals::ntimesteps * MNUBINS * get_proccount());
    std::ranges::fill(spectra.emissionalltimesteps, 0.0);

    if (spectra.trueemissionalltimesteps.empty()) {
      spectra.trueemissionalltimesteps = MPI_shared_array<double>(globals::ntimesteps * MNUBINS * get_proccount());
    }
    assert_always(std::ssize(spectra.trueemissionalltimesteps) == globals::ntimesteps * MNUBINS * get_proccount());
    std::ranges::fill(spectra.trueemissionalltimesteps, 0.0);
  }
  MPI_Barrier_allranks();

  if (print_memusage) {
    printlnlog("[info] mem_usage: set of spectra{} occupy {:.3f} MB (node shared memory)",
               do_emission_absorption ? " (with emission/absorption tracing)" : "",
               spectra.mem_usage_bytes() / 1024. / 1024.);
  }
}

// Add a packet to the outgoing spectrum.
void add_to_spec_res(const Packet& pkt, const int dirbin, Spectra& spectra_I, Spectra* spectra_Q, Spectra* spectra_U) {
  if (dirbin != -1 && get_escapedirectionbin(pkt.dir) != dirbin) {
    return;  // do not add to the spectrum if the direction bin does not match
  }

  // Need to (1) decide which time bin to put it in and (2) which frequency bin.

  // specific angle bins contain fewer packets than the full sphere, so must be normalised to match
  const double nu_min = spectra_I.nu_min;
  const double nu_max = spectra_I.nu_max;
  const double dlognu = spectra_I.dlognu;
  const double t_arrive = pkt.escape_time - (dot(pkt.pos, pkt.dir) / CLIGHT_PROP);
  if (t_arrive > globals::tmin && t_arrive < globals::tmax && pkt.nu_rf > nu_min && pkt.nu_rf < nu_max) {
    const auto nts = get_timestep(t_arrive);

    // a binary search into freq_lower would probably be faster than this double logarithm
    const auto nnu = get_logbinindex(pkt.nu_rf, nu_min, dlognu, MNUBINS);

    const double solidanglefactor = (dirbin >= 0) ? MABINS : 1.;
    const double deltaE = pkt.e_rf / globals::timesteps[nts].width / spectra_I.delta_freq[nnu] / 4.e12 / PI / PARSEC /
                          PARSEC / globals::nprocs_exspec * solidanglefactor;

    const auto fluxindex = (nnu * static_cast<ptrdiff_t>(globals::ntimesteps)) + nts;
    atomicadd_always(spectra_I.fluxalltimesteps[fluxindex], deltaE);

    if (spectra_Q != nullptr) {
      atomicadd_always(spectra_Q->fluxalltimesteps[fluxindex], pkt.stokes_q * deltaE);
    }
    if (spectra_U != nullptr) {
      atomicadd_always(spectra_U->fluxalltimesteps[fluxindex], pkt.stokes_u * deltaE);
    }

    if (spectra_I.do_emission_absorption) {
      const auto proccount = get_proccount();

      const auto truenproc = columnindex_from_emissiontype(pkt.trueemissiontype);
      assert_always(truenproc < proccount);
      if (truenproc >= 0) {
        const auto emindex = (nnu * globals::ntimesteps * proccount) + (nts * proccount) + truenproc;
        atomicadd_always(spectra_I.trueemissionalltimesteps[emindex], deltaE);
      }

      const auto nproc = columnindex_from_emissiontype(pkt.emissiontype);
      assert_always(nproc < proccount);
      if (nproc >= 0) {  // -1 means EMTYPE_NOTSET
        const auto emindex = (nnu * globals::ntimesteps * proccount) + (nts * proccount) + nproc;
        atomicadd_always(spectra_I.emissionalltimesteps[emindex], deltaE);

        if (spectra_Q != nullptr && spectra_Q->do_emission_absorption) {
          atomicadd_always(spectra_Q->emissionalltimesteps[emindex], pkt.stokes_q * deltaE);
        }
        if (spectra_U != nullptr && spectra_U->do_emission_absorption) {
          atomicadd_always(spectra_U->emissionalltimesteps[emindex], pkt.stokes_u * deltaE);
        }
      }

      if (pkt.absorptionfreq > nu_min && pkt.absorptionfreq < nu_max) {
        const auto nnu_abs = get_logbinindex(pkt.absorptionfreq, nu_min, dlognu, MNUBINS);
        const double deltaE_absorption = pkt.e_rf / globals::timesteps[nts].width / spectra_I.delta_freq[nnu_abs] /
                                         4.e12 / PI / PARSEC / PARSEC / globals::nprocs_exspec * solidanglefactor;
        const int at = pkt.absorptiontype;
        if (at >= 0) {
          // bb-absorption
          const int element = globals::linelist.elementindex[at];
          const int ion = globals::linelist.ionindex[at];
          const auto absindex = get_absorption_spectrum_index(nts, nnu_abs) + (element * get_max_nions()) + ion;
          atomicadd_always(spectra_I.absorptionalltimesteps[absindex], deltaE_absorption);

          if (spectra_Q != nullptr && spectra_Q->do_emission_absorption) {
            atomicadd_always(spectra_Q->absorptionalltimesteps[absindex], pkt.stokes_q * deltaE_absorption);
          }
          if (spectra_U != nullptr && spectra_U->do_emission_absorption) {
            atomicadd_always(spectra_U->absorptionalltimesteps[absindex], pkt.stokes_u * deltaE_absorption);
          }
        }
      }
    }
  }
}

void write_light_curve(const std::string& lc_filename, const std::span<const double> light_curve_lum,
                       const std::span<const double> light_curve_lumcmf, const int numtimesteps) {
  if (globals::node_id != 0 || globals::rank_in_node != 0) {
    return;
  }
  assert_always(numtimesteps <= globals::ntimesteps);

  auto lc_file = fstream_required(lc_filename, std::ios::out | std::ios::trunc);

  // UVOIR bolometric light curve
  for (int nts = 0; nts < numtimesteps; nts++) {
    std::println(lc_file, "{:g} {:g} {:g}", globals::timesteps[nts].mid / DAY, light_curve_lum[nts] / LSUN,
                 light_curve_lumcmf[nts] / LSUN);
  }
}

// add a packet to the outgoing light-curve.
void add_to_lc_res(const Packet& pkt, const int dirbin, std::span<double> light_curve_lum,
                   std::span<double> light_curve_lumcmf) {
  if (dirbin >= 0 && get_escapedirectionbin(pkt.dir) != dirbin) {
    return;
  }
  const double solidanglefactor = (dirbin >= 0) ? MABINS : 1.;
  // dirbin -1 means all full 4π angle average (no angle filtering)

  const double t_arrive = pkt.escape_time - (dot(pkt.pos, pkt.dir) / CLIGHT_PROP);
  if (t_arrive > globals::tmin && t_arrive < globals::tmax) {
    const int nts = get_timestep(t_arrive);
    atomicadd_always(light_curve_lum[nts],
                     pkt.e_rf / globals::timesteps[nts].width * solidanglefactor / globals::nprocs_exspec);
  }

  const double inverse_gamma = std::sqrt(1. - (globals::vmax * globals::vmax / CLIGHTSQUARED));

  // Now do the cmf light curve. Unlike t_arrive above, this has no light-travel-time term:
  // it is the escape time transformed (time-dilated) into the comoving frame.
  const double t_escape_cmf = pkt.escape_time * inverse_gamma;

  if (t_escape_cmf > globals::tmin && t_escape_cmf < globals::tmax) {
    const int nts = get_timestep(t_escape_cmf);
    atomicadd_always(light_curve_lumcmf[nts],
                     pkt.e_cmf / globals::timesteps[nts].width * solidanglefactor / globals::nprocs_exspec);
  }
}

// The packets of the spectrum are a subset of the packets of the light curve, because the spectrum has a
// limited frequency range. The frequency-integrated spectrum must therefore stay below the light curve. This
// writes a warning to the log and changes no output file.
void check_spectrum_lightcurve_consistency(const Spectra& spectra_I, const std::span<const double> light_curve_lum,
                                           const int numtimesteps) {
  assert_always(numtimesteps <= globals::ntimesteps);
  if (globals::my_rank != 0) {
    return;
  }
  const auto ntimesteps_all = static_cast<ptrdiff_t>(globals::ntimesteps);
  for (int nts = 0; nts < numtimesteps; nts++) {
    double lum_from_spec = 0.;
    for (auto nnu = 0Z; nnu < MNUBINS; nnu++) {
      lum_from_spec += spectra_I.fluxalltimesteps[(nnu * ntimesteps_all) + nts] * spectra_I.delta_freq[nnu];
    }
    // undo the flux normalisation of add_to_spec_res() to get back to a luminosity
    lum_from_spec *= 4.e12 * PI * PARSEC * PARSEC;
    const double lum_lightcurve = light_curve_lum[nts];
    if (lum_lightcurve > 0. && lum_from_spec > (lum_lightcurve * 1.001)) {
      printlnlog(
          "[warning] consistency check failed for timestep {}: frequency-integrated spec.out luminosity {:g} "
          "[erg/s] exceeds the light_curve.out luminosity {:g} [erg/s], but the spectrum's packets should be a "
          "subset of the light curve's packets",
          nts, lum_from_spec, lum_lightcurve);
    }
  }
}

// Sum the spectra of the nodes. The arrays are node-shared memory, so only one rank of each node takes part.
void mpi_allreduce_spectra(Spectra& spectra) {
  assert_always(globals::rank_in_node == 0);
  MPI_Allreduce_safe(spectra.fluxalltimesteps, MPI_SUM, globals::mpi_comm_internode);
  if (spectra.do_emission_absorption) {
    MPI_Allreduce_safe(spectra.absorptionalltimesteps, MPI_SUM, globals::mpi_comm_internode);
    MPI_Allreduce_safe(spectra.emissionalltimesteps, MPI_SUM, globals::mpi_comm_internode);
    MPI_Allreduce_safe(spectra.trueemissionalltimesteps, MPI_SUM, globals::mpi_comm_internode);
  }
}

void write_partial_lightcurve_spectra_dirbin(const int nts, std::span<const Packet> packets,
                                             const bool do_emission_absorption, const int dirbin) {
  THREADLOCALONHOST std::vector<double> rpkt_light_curve_lum;
  THREADLOCALONHOST std::vector<double> rpkt_light_curve_lumcmf;
  THREADLOCALONHOST std::vector<double> gamma_light_curve_lum;
  THREADLOCALONHOST std::vector<double> gamma_light_curve_lumcmf;
  reserve_resize(rpkt_light_curve_lum, globals::ntimesteps);
  std::ranges::fill(rpkt_light_curve_lum, 0.);
  reserve_resize(rpkt_light_curve_lumcmf, globals::ntimesteps);
  std::ranges::fill(rpkt_light_curve_lumcmf, 0.);
  if constexpr (KEEP_ESCAPED_GAMMAS) {
    reserve_resize(gamma_light_curve_lum, globals::ntimesteps);
    std::ranges::fill(gamma_light_curve_lum, 0.);
    reserve_resize(gamma_light_curve_lumcmf, globals::ntimesteps);
    std::ranges::fill(gamma_light_curve_lumcmf, 0.);
  }

  init_spectra(rpkt_spectra_I, NU_MIN_R, NU_MAX_R, do_emission_absorption);
  if constexpr (POL_ON) {
    init_spectra(rpkt_spectra_Q, NU_MIN_R, NU_MAX_R, do_emission_absorption);
    init_spectra(rpkt_spectra_U, NU_MIN_R, NU_MAX_R, do_emission_absorption);
  }
  // the gamma packets go only into the angle-averaged spectrum and light curve
  const bool do_gamma_spectrum = KEEP_ESCAPED_GAMMAS && (dirbin == -1);
  if (do_gamma_spectrum) {
    init_spectra(gamma_spectra, NU_MIN_GAMMA, NU_MAX_GAMMA, false);
  }

  MPI_Barrier_node();
#if defined REPRODUCIBLE && REPRODUCIBLE
  for (int node_rank = 0; node_rank < globals::node_nprocs; node_rank++) {
    // do one rank at a time to keep the results reproducible (instead of simultaneous atomic adds to shared memory)
#else
  {
    // all ranks on the node simultaneously contribute to the light curves and spectra in shared memory using
    // atomic operations
    const int node_rank = globals::rank_in_node;
#endif
    if (node_rank == globals::rank_in_node) {
      for (const auto& pkt : packets) {
        if (pkt.type == TYPE_ESCAPE) {
          if (pkt.escape_type == TYPE_RPKT) {
            add_to_lc_res(pkt, dirbin, rpkt_light_curve_lum, rpkt_light_curve_lumcmf);
            add_to_spec_res(pkt, dirbin, rpkt_spectra_I, POL_ON ? &rpkt_spectra_Q : nullptr,
                            POL_ON ? &rpkt_spectra_U : nullptr);
          } else if (do_gamma_spectrum && pkt.escape_type == TYPE_GAMMA) {
            add_to_lc_res(pkt, dirbin, gamma_light_curve_lum, gamma_light_curve_lumcmf);
            add_to_spec_res(pkt, dirbin, gamma_spectra, nullptr, nullptr);
          }
        }
      }
    }
    MPI_Barrier_node();
  }

  const int numtimesteps = nts + 1;  // only produce spectra and light curves up to one past nts
  assert_always(numtimesteps <= globals::ntimesteps);

  MPI_Barrier_allranks();
  if (globals::rank_in_node == 0) {
    mpi_allreduce_spectra(rpkt_spectra_I);
    if constexpr (POL_ON) {
      mpi_allreduce_spectra(rpkt_spectra_Q);
      mpi_allreduce_spectra(rpkt_spectra_U);
    }
    if (do_gamma_spectrum) {
      mpi_allreduce_spectra(gamma_spectra);
    }
  }
  MPI_Allreduce_safe(rpkt_light_curve_lum, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce_safe(rpkt_light_curve_lumcmf, MPI_SUM, MPI_COMM_WORLD);
  if constexpr (KEEP_ESCAPED_GAMMAS) {
    MPI_Allreduce_safe(gamma_light_curve_lum, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce_safe(gamma_light_curve_lumcmf, MPI_SUM, MPI_COMM_WORLD);
  }
  MPI_Barrier_allranks();

  if (dirbin == -1) {
    write_light_curve("light_curve.out", rpkt_light_curve_lum, rpkt_light_curve_lumcmf, numtimesteps);
    if constexpr (KEEP_ESCAPED_GAMMAS) {
      write_light_curve("gamma_light_curve.out", gamma_light_curve_lum, gamma_light_curve_lumcmf, numtimesteps);
      write_spectra("gamma_spec.out", "", "", "", gamma_spectra, numtimesteps);
    }
    write_spectra("spec.out", "emission.out", "emissiontrue.out", "absorption.out", rpkt_spectra_I, numtimesteps);
    if constexpr (POL_ON) {
      write_specpol("specpol.out", "emissionpol.out", "absorptionpol.out", rpkt_spectra_I, rpkt_spectra_Q,
                    rpkt_spectra_U, numtimesteps);
    }
    if (do_emission_absorption) {
      check_spectrum_lightcurve_consistency(rpkt_spectra_I, rpkt_light_curve_lum, numtimesteps);
    }
  } else {
    if (globals::my_rank == 0 && !std::filesystem::exists(outdir_resfiles)) {
      std::filesystem::create_directory(outdir_resfiles);
    }
    MPI_Barrier_allranks();

    write_light_curve(std::format("{}light_curve_res_{:02d}.out", outdir_resfiles, dirbin), rpkt_light_curve_lum,
                      rpkt_light_curve_lumcmf, numtimesteps);
    write_spectra(std::format("{}spec_res_{:02d}.out", outdir_resfiles, dirbin),
                  std::format("{}emission_res_{:02d}.out", outdir_resfiles, dirbin),
                  std::format("{}emissiontrue_res_{:02d}.out", outdir_resfiles, dirbin),
                  std::format("{}absorption_res_{:02d}.out", outdir_resfiles, dirbin), rpkt_spectra_I, numtimesteps);

    if constexpr (POL_ON) {
      write_specpol(std::format("{}specpol_res_{:02d}.out", outdir_resfiles, dirbin),
                    std::format("{}emissionpol_res_{:02d}.out", outdir_resfiles, dirbin),
                    std::format("{}absorptionpol_res_{:02d}.out", outdir_resfiles, dirbin), rpkt_spectra_I,
                    rpkt_spectra_Q, rpkt_spectra_U, numtimesteps);
    }
  }
  MPI_Barrier_allranks();
}

}  // anonymous namespace

void write_partial_lightcurve_spectra(const int nts, std::span<const Packet> pkts) {
  // sn3d calls this with the packets of its rank, and exspec calls it with the packets of all ranks
  const bool simulation_complete = (nts >= globals::timestep_finish - 1);

  // the emission resolved spectra are slow to generate, and require a lot of memory. The code
  // therefore makes them only at the end of the simulation.
  const bool do_emission_absorption = simulation_complete;

  const bool multdimensional = grid::get_modelgridtype() != GridType::SPHERICAL1D;
  const int dirbinend = (multdimensional && simulation_complete) ? MABINS : 0;

  const auto time_func_start = std::chrono::steady_clock::now();

  for (int dirbin = -1; dirbin < dirbinend; dirbin++) {
    write_partial_lightcurve_spectra_dirbin(nts, pkts, do_emission_absorption, dirbin);
    if (dirbin >= 0 && globals::my_rank == 0) {
      printlnlog("timestep {}: wrote the files of direction bin {} (the last bin is {})", nts, dirbin, dirbinend - 1);
    }
  }

  const auto duration_write_spectra =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - time_func_start).count();
  printlnlog("timestep {}: Saving light curves and {}spectra took {:.1f}s", nts,
             do_emission_absorption ? "emission/absorption " : "", duration_write_spectra);
}

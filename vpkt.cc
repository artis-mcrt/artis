// Virtual packets (optional VPKT mode): at real-packet emission and scattering events, virtual
// packets are launched toward predefined observer directions and their optical-depth-weighted
// energies are accumulated into viewing-angle-dependent (and optionally polarised) spectra.

#include "vpkt.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <iterator>
#include <print>
#include <sstream>
#include <string>
#include <system_error>
#include <tuple>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "ltepop.h"
#include "mpi_logging.h"
#include "packet.h"
#include "rpkt.h"
#include "sn3d.h"
#include "vectors.h"

namespace vpkt {
namespace {

static_assert(!VPKT_ON || POL_ON,
              "POL_ON must be true if VPKT_ON is true because vpkt needs stokes parameters to be tracked");

struct StokesParams {
  double I = 0.;
  double Q = 0.;
  double U = 0.;
};

struct VSpecPol {
  std::array<StokesParams, VSPEC_NUBINS> flux;
  float lower_time{NAN};
  float delta_t{NAN};
};

std::vector<std::vector<VSpecPol>> vspecpol{};

std::array<float, VSPEC_NUBINS> lower_freq_vspec;
std::array<float, VSPEC_NUBINS> delta_freq_vspec;

int nobsdirections = 0;  // Number of observer directions
int nspectraperobsdir = 0;  // Number of virtual packet spectra per observer direction (total + elements switched off)
std::vector<double> obsdirs_costheta;
std::vector<double> obsdirs_phi;
double vspec_timemin_input;
double vspec_timemax_input;
int nwavelengthranges = 0;  // Number of wavelength ranges

std::vector<double> vspec_numin_input;
std::vector<double> vspec_numax_input;
double tau_max_vpkt;

std::vector<int> opacityexclusions;  // vector of opacity contribution setups:
                                     // 0: full opacity
                                     // -1: no line opacity; -2: no bf opacity; -3: no ff opacity; -4: no es opacity,
                                     // +ve: exclude element with atomic number's contribution to bound-bound opacity
// per-thread accumulated optical depths (one per opacity-choice spectrum), since vpkts are traced
// from within the OpenMP/stdpar-parallel packet loop
THREADLOCALONHOST std::vector<double> tau_vpkt;

std::ofstream vpkt_contrib_file;
PaddedMutex vpkt_contrib_file_mutex;

struct VGrid {
  std::vector<std::vector<StokesParams>> flux;
  double yvel{NAN};
  double zvel{NAN};
};

std::array<std::array<VGrid, VGRID_NZ>, VGRID_NY> vgrid;

int grid_nwavelengthranges{};
double tmin_grid;
double tmax_grid;
std::vector<double> nu_grid_min;
std::vector<double> nu_grid_max;
bool vgrid_on;

const double dlogt_vspec = (std::log(VSPEC_TIMEMAX) - std::log(VSPEC_TIMEMIN)) / VSPEC_TIMEBINS;
const double dlognu_vspec = (std::log(VSPEC_NUMAX) - std::log(VSPEC_NUMIN)) / VSPEC_NUBINS;

// Virtual packet is killed when tau reaches tau_max_vpkt for ALL the different setups
// E.g. imagine that a packet in the first setup (all elements included) reaches tau = tau_max_vpkt
// because of the element Zi. If we remove Zi, tau now could be lower than tau_max_vpkt and could
// thus contribute to the spectrum.
constexpr auto all_taus_past_taumax(std::vector<double>& tau, const double tau_max) -> bool {
  return std::ranges::all_of(tau, [tau_max](const double tau_i) { return tau_i > tau_max; });
}

// Add a packet to the outgoing spectrum
void add_to_vspecpol(const double nu_rf, const double e_rf, const double prob, const double q_rf, const double u_rf,
                     const int obsdirindex, const int opachoiceindex, const double t_arrive) {
  // Need to decide in which (1) time and (2) frequency bin the vpkt is escaping

  if (t_arrive <= VSPEC_TIMEMIN || t_arrive >= VSPEC_TIMEMAX || nu_rf <= VSPEC_NUMIN || nu_rf >= VSPEC_NUMAX) {
    return;
  }

  const auto nt = static_cast<int>(get_logbinindex(t_arrive, VSPEC_TIMEMIN, dlogt_vspec, VSPEC_TIMEBINS));
  const auto nnu = static_cast<int>(get_logbinindex(nu_rf, VSPEC_NUMIN, dlognu_vspec, VSPEC_NUBINS));

  const int ind_comb = (nspectraperobsdir * obsdirindex) + opachoiceindex;
  const double pktcontrib = e_rf / vspecpol[nt][ind_comb].delta_t / delta_freq_vspec[nnu] / 4.e12 / PI / PARSEC /
                            PARSEC / globals::nprocs * 4 * PI;

  atomicadd(vspecpol[nt][ind_comb].flux[nnu].I, prob * pktcontrib);
  atomicadd(vspecpol[nt][ind_comb].flux[nnu].Q, prob * q_rf * pktcontrib);
  atomicadd(vspecpol[nt][ind_comb].flux[nnu].U, prob * u_rf * pktcontrib);
}

// Add a packet to the outgoing spectrum
void add_to_vpkt_grid(const double nu_rf, const double e_rf, const double prob, const double stokes_q,
                      const double stokes_u, const Vec3d& vel, const int wlbin, const int obsdirindex,
                      const Vec3d& obsdir) {
  double vref1{NAN};
  double vref2{NAN};

  // obs is the observer orientation

  // Packet velocity

  if (obsdir[0] == 1) {
    // if obsdir = +x , vref1 = vy and vref2 = vz
    vref1 = vel[1];
    vref2 = vel[2];
  } else if (obsdir[0] == -1) {
    // if obsdir = -x , vref1 = -vy and vref2 = -vz
    vref1 = -vel[1];
    vref2 = -vel[2];
  } else {
    // Rotate velocity into projected area seen by the observer (see notes)
    // Rotate velocity from (x,y,z) to (obsdir,vref1,vref2) so that x corresponds to obsdir
    const double crossterm = obsdir[1] * obsdir[2] / (1 + obsdir[0]);
    vref1 = (-obsdir[1] * vel[0]) + ((obsdir[0] + (pow2(obsdir[2]) / (1 + obsdir[0]))) * vel[1]) - (crossterm * vel[2]);
    vref2 = (-obsdir[2] * vel[0]) - (crossterm * vel[1]) + ((obsdir[0] + (pow2(obsdir[1]) / (1 + obsdir[0]))) * vel[2]);
  }

  // Outside the grid
  if (fabs(vref1) >= globals::vmax || fabs(vref2) >= globals::vmax) {
    return;
  }

  // Bin size

  // vgrid cell (can be different to propagation cell size)
  const int ny = static_cast<int>((globals::vmax - vref1) / (2 * globals::vmax / VGRID_NY));
  const int nz = static_cast<int>((globals::vmax - vref2) / (2 * globals::vmax / VGRID_NZ));

  // Add contribution
  if (nu_rf > nu_grid_min[wlbin] && nu_rf < nu_grid_max[wlbin]) {
    atomicadd(vgrid[ny][nz].flux[wlbin][obsdirindex].I, prob * e_rf);
    atomicadd(vgrid[ny][nz].flux[wlbin][obsdirindex].Q, prob * stokes_q * e_rf);
    atomicadd(vgrid[ny][nz].flux[wlbin][obsdirindex].U, prob * stokes_u * e_rf);
  }
}

auto trace_vpkt_direction(const Packet& rpkt, const double t_arrive, const double nu_rf, const double e_rf,
                          const double rpkt_doppler, const int obsdirindex, const Vec3d& obsdir,
                          const enum packet_type type_before_rpkt, std::string& vpkt_contrib_row) -> bool {
  int mgi = 0;

  auto cellindex = rpkt.cellindex;
  auto next_trans = rpkt.next_trans;  // should be -1 since doppler factor changed it?
  auto e_cmf = rpkt.e_cmf;
  auto nu_cmf = rpkt.nu_cmf;
  auto vpktpos = rpkt.pos;

  bool end_packet = false;
  const double t_start = rpkt.prop_time;
  double t_future = t_start;

  // The cell's material state (level populations, rho, and the binned expansion opacities) was calculated
  // for the middle of the current timestep, so t_gridstate is the reference time when scaling those
  // quantities to the packet's own time as it propagates outwards through the expanding ejecta.
  const double t_gridstate = globals::timesteps[globals::timestep].mid;

  tau_vpkt.assign(nspectraperobsdir, 0.);  // sizes this thread's copy on first use and zeroes it

  atomicadd(nvpkt_created, 1);

  const auto vel_vec = get_velocity(rpkt.pos, t_start);

  // Scattering event: dipole function

  // MACROATOM and KPKT: isotropic emission
  double pn{1 / (4 * PI)};
  // normalised Stokes parameters in the rest frame (RF)
  double q_rf{0.};
  double u_rf{0.};
  if (type_before_rpkt == TYPE_RPKT) {
    // Transform Stokes Parameters from the RF to the CMF
    const auto [old_dir_cmf, q_i_cmf, u_i_cmf] = POL_ON
                                                     ? frame_transform(rpkt.dir, rpkt.stokes_q, rpkt.stokes_u, vel_vec)
                                                     : std::make_tuple(angle_ab(rpkt.dir, vel_vec), 0., 0.);

    // Need to rotate Stokes Parameters in the scattering plane

    const auto obs_cmf = angle_ab(obsdir, vel_vec);
    std::tie(std::ignore, q_rf, u_rf, pn) = scatter_polarisation_to_rf(old_dir_cmf, obs_cmf, q_i_cmf, u_i_cmf, vel_vec);

  } else {
    assert_testmodeonly(type_before_rpkt == TYPE_KPKT || type_before_rpkt == TYPE_MA);
  }

  // pn is the emission probability per unit solid angle in the CMF, but the vpkt escapes along a fixed
  // rest-frame direction. Convert with the solid-angle transformation dOmega_cmf / dOmega_rf = 1 / doppler^2
  // (relativistic beaming), consistent with real packets, whose CMF-sampled directions are aberrated to the RF
  pn /= pow2(rpkt_doppler);

  // compute the optical depth to boundary

  mgi = grid::get_propcell_modelgridindex(cellindex);
  THREADLOCALONHOST auto chi_vpkt_cont = ContinuumOpacity{false};

  while (!end_packet) {
    // distance to the next cell
    const auto [boundarydist, next_cellindex] = grid::boundary_distance(obsdir, vpktpos, t_future, cellindex);
    if (mgi < 0) {
      next_trans = -1;
    } else if (boundarydist > 0) {
      const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);
      const auto doppler = calculate_doppler_nucmf_on_nurf(vpktpos, obsdir, t_future);
      calculate_chi_rpkt_cont<false>(nu_cmf, chi_vpkt_cont, nonemptymgi);

      // ratio of the number densities when the packet crosses this cell to those of the grid state the
      // opacities were calculated from
      const double densityscalefactor = pow3(t_gridstate / t_future);

      // Scale each continuum process by its own power of the density rather than by a single common
      // factor. Electron scattering (sigma_T * nne) and bound-free (nnlevel * sigma_bf) are linear in a
      // number density, but free-free (nnion * nne) is quadratic. NB: the stimulated recombination that
      // is subtracted inside chi_boundfree is also quadratic in density, but it is not separable from the
      // sum that calculate_chi_bf_gammacontr() returns, so bound-free is scaled as purely linear here.
      const double chi_escatter = chi_vpkt_cont.chi_escatter * densityscalefactor;
      const double chi_bf = chi_vpkt_cont.chi_boundfree * densityscalefactor;
      const double chi_ff = chi_vpkt_cont.chi_freefree_heat * pow2(densityscalefactor);

      const double chi_cont = chi_escatter + chi_bf + chi_ff;

      for (int opacchoiceindex = 0; opacchoiceindex < nspectraperobsdir; opacchoiceindex++) {
        double chi_cont_thischoice = chi_cont;
        if (opacityexclusions[opacchoiceindex] == -2) {
          chi_cont_thischoice -= chi_bf;
        } else if (opacityexclusions[opacchoiceindex] == -3) {
          chi_cont_thischoice -= chi_ff;
        } else if (opacityexclusions[opacchoiceindex] == -4) {
          chi_cont_thischoice -= chi_escatter;
        }
        tau_vpkt[opacchoiceindex] += chi_cont_thischoice * boundarydist * doppler;
      }

      // kill vpkt with high optical depth
      if (all_taus_past_taumax(tau_vpkt, tau_max_vpkt)) {
        return false;
      }

      const Vec3d pos_boundary{
          vpktpos[0] + (obsdir[0] * boundarydist),
          vpktpos[1] + (obsdir[1] * boundarydist),
          vpktpos[2] + (obsdir[2] * boundarydist),
      };

      const double nu_cmf_boundary = std::min(
          nu_rf * calculate_doppler_nucmf_on_nurf(pos_boundary, obsdir, t_future + (boundarydist / CLIGHT_PROP)),
          nu_cmf);

      const double dnu_on_dl = (nu_cmf_boundary - nu_cmf) / boundarydist;
      // Trace individual lines from nu_cmf until dist_limit. Returns false when all vpkt opacity setups exceed tau_max.
      const auto trace_lines_to_dist = [&](const double dist_limit) -> bool {
        while (true) {
          const int lineindex = closest_transition(nu_cmf, next_trans, globals::linelist.nu);

          if (lineindex < 0) {
            next_trans = globals::nlines + 1;
            break;
          }
          const double nutrans = globals::linelist.nu[lineindex];

          next_trans = lineindex + 1;

          const auto ldist = get_linedistance(t_future, nu_cmf, nutrans, dnu_on_dl);

          if (ldist > dist_limit) {
            next_trans--;
            break;
          }

          const double t_line = t_future + (ldist / CLIGHT_PROP);

          const int element = globals::linelist.elementindex[lineindex];
          const int ion = globals::linelist.ionindex[lineindex];
          const auto ionuniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);
          const auto lower_uniquelevelindex = globals::linelist.uniquelevelindex_lower[lineindex];
          const auto upper_uniquelevelindex = globals::linelist.uniquelevelindex_upper[lineindex];
          const int lower = lower_uniquelevelindex - ionuniquelevelindexstart;
          const int upper = upper_uniquelevelindex - ionuniquelevelindexstart;

          const double B_ul = globals::linelist.B_ul[lineindex];
          const double B_lu = globals::linelist.B_lu[lineindex];

          const auto n_u = calculate_levelpop(nonemptymgi, element, ion, upper);
          const auto n_l = calculate_levelpop(nonemptymgi, element, ion, lower);
          // The level populations belong to the grid state at t_gridstate, but the packet only reaches this
          // line resonance at t_line, by which time the ejecta have expanded further. Scale them with
          // n ∝ t^-3, exactly as densityscalefactor does for the continuum opacity. Combined with the factor of
          // t_line in the Sobolev optical depth this gives the expected tau_sobolev ∝ t^-2.
          const double popscalefactor = pow3(t_gridstate / t_line);
          const double tau_line =
              std::max(0., ((B_lu * n_l) - (B_ul * n_u)) * popscalefactor * HCLIGHTOVERFOURPI * t_line);

          // Check on the element to exclude (or -1 for no line opacity)
          const int Z = get_atomicnumber(element);
          for (int opacchoiceindex = 0; opacchoiceindex < nspectraperobsdir; opacchoiceindex++) {
            if (opacityexclusions[opacchoiceindex] != -1 && (opacityexclusions[opacchoiceindex] != Z)) {
              tau_vpkt[opacchoiceindex] += tau_line;
            }
          }

          // kill vpkt with high optical depth
          if (all_taus_past_taumax(tau_vpkt, tau_max_vpkt)) {
            return false;
          }
        }

        return true;
      };
      if constexpr (VPKT_USE_EXPANSION_OPACITIES) {
        // -1 means below the wavelength grid, in which case there is continuum opacity but no expansion opacity
        const auto binindex_start =
            std::max(get_linearbinindex(1e8 * CLIGHT / nu_cmf, expopac_lambdamin, expopac_deltalambda), -1Z);

        if (binindex_start < expopac_nbins) {
          // trace line-by-line from nu_cmf to the next bin edge, because the expansion opacity bin at nu_cmf includes
          // lines with nu > nu_cmf
          const auto first_bin_edge_nu =
              (binindex_start < 0) ? get_expopac_bin_nu_upper(0) : get_expopac_bin_nu_lower(binindex_start);
          const auto first_bin_edge_dist = get_linedistance(t_future, nu_cmf, first_bin_edge_nu, dnu_on_dl);
          const double line_by_line_limit = std::min(first_bin_edge_dist, boundarydist);
          next_trans = -1;  // trigger binary search from nu_cmf
          if (!trace_lines_to_dist(line_by_line_limit)) {
            return false;
          }

          double dist = line_by_line_limit;

          if (dist < boundarydist) {
            // continue with expansion opacities for subsequent full bins
            for (auto binindex = binindex_start + 1; binindex < expopac_nbins; binindex++) {
              const auto next_bin_edge_nu = get_expopac_bin_nu_lower(binindex);
              const auto binedgedist = get_linedistance(t_future, nu_cmf, next_bin_edge_nu, dnu_on_dl);

              const auto kappa = expansionopacities[(nonemptymgi * expopac_nbins) + binindex];
              // kappa_exp * rho = (1 / (c t)) * sum_lines (lambda_line / delta_lambda) * (1 - exp(-tau_sobolev))
              // and was tabulated at t_gridstate, so scale it to the packet's time. In the optically thin limit
              // (1 - exp(-tau_sobolev)) -> tau_sobolev ∝ t^-2, which together with the explicit 1/(c t)
              // prefactor gives kappa_exp * rho ∝ t^-3, i.e. the same density scaling as the linear continuum
              // terms above. Saturated lines keep (1 - exp(-tau_sobolev)) ~ 1 and so fall off only as 1/t, but
              // their individual tau_sobolev cannot be recovered from the binned kappa, so the thin limit is
              // used for all bins.
              const double chi_bb_expansionopac = kappa * grid::get_rho(nonemptymgi) * densityscalefactor;

              const double tau_bin = chi_bb_expansionopac * (std::min(binedgedist, boundarydist) - dist);
              dist = std::min(binedgedist, boundarydist);

              for (int opacchoiceindex = 0; opacchoiceindex < nspectraperobsdir; opacchoiceindex++) {
                assert_testmodeonly(opacityexclusions[opacchoiceindex] <=
                                    0);  // expansion opacities include all elements, so cannot
                                         // be used with custom lists that exclude some elements
                if (opacityexclusions[opacchoiceindex] != -1) {
                  tau_vpkt[opacchoiceindex] += tau_bin;
                }
              }

              // kill vpkt with high optical depth
              if (all_taus_past_taumax(tau_vpkt, tau_max_vpkt)) {
                return false;
              }

              if (dist >= boundarydist) {
                break;
              }
            }
          }
        }  // if (binindex_start < expopac_nbins)
      } else {
        if (!trace_lines_to_dist(boundarydist)) {
          return false;
        }
      }
    }

    // virtual packet is still at the starting position
    // move it to cell boundary and go to next cell

    move_pkt_withtime(vpktpos, obsdir, t_future, nu_rf, nu_cmf, e_rf, e_cmf, boundarydist);

    if (next_cellindex >= 0) {
      if (next_cellindex != cellindex) {
        // snap the position onto the crossed cell boundary so that rounding errors cannot
        // accumulate over many crossings (same-cell returns from boundary_distance() due to
        // max_path_step leave the packet mid-cell and must not be snapped)
        grid::snap_pos_to_cell(vpktpos, t_future, next_cellindex);
      }
      cellindex = next_cellindex;
      mgi = grid::get_propcell_modelgridindex(cellindex);
      if (mgi >= 0) {
        const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);

        // kill vpkt with pass through a thick cell
        if (grid::thick_allcells[nonemptymgi] != 0) {
          return false;
        }
      }
    } else {
      end_packet = true;
    }
  }

  // increment the number of escaped virtual packet in the given timestep
  if (type_before_rpkt == TYPE_RPKT) {
    atomicadd(nvpkt_esc_from_rpkt, 1);
  } else if (type_before_rpkt == TYPE_KPKT) {
    atomicadd(nvpkt_esc_from_kpkt, 1);
  } else if (type_before_rpkt == TYPE_MA) {
    atomicadd(nvpkt_esc_from_macroatom, 1);
  }

  // Final Stokes vector

  if (VPKT_WRITE_CONTRIBS) {
    std::format_to(std::back_inserter(vpkt_contrib_row), " {:g} {:g}", t_arrive / DAY, nu_rf);
  }

  for (int opacchoiceindex = 0; opacchoiceindex < nspectraperobsdir; opacchoiceindex++) {
    const double prob = pn * std::exp(-tau_vpkt[opacchoiceindex]);

    assert_always(std::isfinite(prob));

    add_to_vspecpol(nu_rf, e_rf, prob, q_rf, u_rf, obsdirindex, opacchoiceindex, t_arrive);

    if constexpr (VPKT_WRITE_CONTRIBS) {
      std::format_to(std::back_inserter(vpkt_contrib_row), " {:g}", e_rf * prob);
    }
  }

  // vpkt grid

  if (vgrid_on) {
    const double prob = pn * exp(-tau_vpkt[0]);

    for (int wlbin = 0; wlbin < grid_nwavelengthranges; wlbin++) {
      if ((nu_rf > nu_grid_min[wlbin] && nu_rf < nu_grid_max[wlbin]) &&
          (t_arrive > tmin_grid && t_arrive < tmax_grid))  // Frequency selection
      {  // Time selection
        add_to_vpkt_grid(nu_rf, e_rf, prob, q_rf, u_rf, vel_vec, wlbin, obsdirindex, obsdir);
      }
    }
  }
  return true;  // true if we added columns to vpkt_contrib_row
}

void init_vspecpol() {
  vspecpol.resize(VSPEC_TIMEBINS, {});

  const int indexmax = nspectraperobsdir * nobsdirections;
  for (int p = 0; p < VSPEC_TIMEBINS; p++) {
    vspecpol[p].resize(indexmax, {});
  }

  for (int m = 0; m < VSPEC_NUBINS; m++) {
    lower_freq_vspec[m] = static_cast<float>(std::exp(log(VSPEC_NUMIN) + (m * dlognu_vspec)));
    delta_freq_vspec[m] =
        static_cast<float>(std::exp(log(VSPEC_NUMIN) + ((m + 1) * dlognu_vspec)) - lower_freq_vspec[m]);
  }

  // start by setting up the time and frequency bins.
  // it is all done interms of a logarithmic spacing in both t and nu - get the
  // step sizes first.
  for (int n = 0; n < VSPEC_TIMEBINS; n++) {
    for (int ind_comb = 0; ind_comb < indexmax; ind_comb++) {
      vspecpol[n][ind_comb].lower_time = static_cast<float>(std::exp(log(VSPEC_TIMEMIN) + (n * dlogt_vspec)));
      vspecpol[n][ind_comb].delta_t =
          static_cast<float>(std::exp(log(VSPEC_TIMEMIN) + ((n + 1) * dlogt_vspec)) - vspecpol[n][ind_comb].lower_time);

      for (auto& flux : vspecpol[n][ind_comb].flux) {
        flux.I = 0.;
        flux.Q = 0.;
        flux.U = 0.;
      }
    }
  }
}

void write_vspecpol(const std::string& filename) {
  auto vspecpol_file = fstream_required(filename, std::ios::out | std::ios::trunc);
  for (int ind_comb = 0; ind_comb < (nobsdirections * nspectraperobsdir); ind_comb++) {
    std::print(vspecpol_file, "{:g}", 0.);

    for (int l = 0; l < 3; l++) {
      for (int p = 0; p < VSPEC_TIMEBINS; p++) {
        std::print(vspecpol_file, " {:g}",
                   (vspecpol[p][ind_comb].lower_time + (vspecpol[p][ind_comb].delta_t / 2.)) / DAY);
      }
    }

    std::println(vspecpol_file, "");

    for (int m = 0; m < VSPEC_NUBINS; m++) {
      std::print(vspecpol_file, "{:g}", (lower_freq_vspec[m] + (delta_freq_vspec[m] / 2.)));

      // Stokes I
      for (int p = 0; p < VSPEC_TIMEBINS; p++) {
        std::print(vspecpol_file, " {:g}", vspecpol[p][ind_comb].flux[m].I);
      }

      // Stokes Q
      for (int p = 0; p < VSPEC_TIMEBINS; p++) {
        std::print(vspecpol_file, " {:g}", vspecpol[p][ind_comb].flux[m].Q);
      }

      // Stokes U
      for (int p = 0; p < VSPEC_TIMEBINS; p++) {
        std::print(vspecpol_file, " {:g}", vspecpol[p][ind_comb].flux[m].U);
      }

      std::println(vspecpol_file, "");
    }
  }
}

void read_vspecpol(const int my_rank, const int nts) {
  const auto filename = std::format("vspecpol_{:04d}_ts{}.tmp", my_rank, nts);
  printlnlog("Reading {}", filename);

  auto vspecpol_file = fstream_required(filename, std::ios::in);
  std::string line;

  for (int ind_comb = 0; ind_comb < (nobsdirections * nspectraperobsdir); ind_comb++) {
    get_noncommentline(vspecpol_file, line);  // first line is the header of times

    // Initialise I,Q,U fluxes from temporary files
    for (int j = 0; j < VSPEC_NUBINS; j++) {
      get_noncommentline(vspecpol_file, line);
      auto ssline = std::stringstream(line);

      float c{};
      assert_always(ssline >> c);  // frequency bin center

      // Stokes I
      for (int p = 0; p < VSPEC_TIMEBINS; p++) {
        assert_always(ssline >> vspecpol[p][ind_comb].flux[j].I);
      }

      // Stokes Q
      for (int p = 0; p < VSPEC_TIMEBINS; p++) {
        assert_always(ssline >> vspecpol[p][ind_comb].flux[j].Q);
      }

      // Stokes U
      for (int p = 0; p < VSPEC_TIMEBINS; p++) {
        assert_always(ssline >> vspecpol[p][ind_comb].flux[j].U);
      }
    }
  }
}

void init_vpkt_grid() {
  const double ybin = 2 * globals::vmax / VGRID_NY;
  const double zbin = 2 * globals::vmax / VGRID_NZ;

  for (int n = 0; n < VGRID_NY; n++) {
    for (int m = 0; m < VGRID_NZ; m++) {
      const double yvel = globals::vmax - ((n + 0.5) * ybin);
      const double zvel = globals::vmax - ((m + 0.5) * zbin);

      vgrid[n][m].yvel = yvel;
      vgrid[n][m].zvel = zvel;

      reserve_resize(vgrid[n][m].flux, grid_nwavelengthranges);

      for (int wlbin = 0; wlbin < grid_nwavelengthranges; wlbin++) {
        reserve_resize(vgrid[n][m].flux[wlbin], nobsdirections);
        std::ranges::fill(vgrid[n][m].flux[wlbin], StokesParams{.I = 0., .Q = 0., .U = 0.});
      }
    }
  }
}

void write_vpkt_grid(const std::string& filename) {
  auto vpkt_grid_file = fstream_required(filename, std::ios::out | std::ios::trunc);

  for (int obsdirindex = 0; obsdirindex < nobsdirections; obsdirindex++) {
    for (int wlbin = 0; wlbin < grid_nwavelengthranges; wlbin++) {
      for (int n = 0; n < VGRID_NY; n++) {
        for (int m = 0; m < VGRID_NZ; m++) {
          std::println(vpkt_grid_file, "{:g} {:g} {:g} {:g} {:g}", vgrid[n][m].yvel, vgrid[n][m].zvel,
                       vgrid[n][m].flux[wlbin][obsdirindex].I, vgrid[n][m].flux[wlbin][obsdirindex].Q,
                       vgrid[n][m].flux[wlbin][obsdirindex].U);
        }
      }
    }
  }
}

void read_vpkt_grid(const int my_rank, const int nts) {
  const std::string filename = std::format("vpkt_grid_{:04d}_ts{}.tmp", my_rank, nts);
  printlnlog("Reading {}", filename);
  auto vpkt_grid_file = fstream_required(filename, std::ios::in);

  for (int obsdirindex = 0; obsdirindex < nobsdirections; obsdirindex++) {
    for (int wlbin = 0; wlbin < grid_nwavelengthranges; wlbin++) {
      for (int n = 0; n < VGRID_NY; n++) {
        for (int m = 0; m < VGRID_NZ; m++) {
          assert_always(vpkt_grid_file >> vgrid[n][m].yvel >> vgrid[n][m].zvel >>
                        vgrid[n][m].flux[wlbin][obsdirindex].I >> vgrid[n][m].flux[wlbin][obsdirindex].Q >>
                        vgrid[n][m].flux[wlbin][obsdirindex].U);
        }
      }
    }
  }
}

}  // anonymous namespace

void remove_temp_vpkt_file(const int nts, const int my_rank) {
  const std::array filenames{
      std::format("vspecpol_{:04d}_ts{}.tmp", my_rank, nts),
      std::format("vpkt_grid_{:04d}_ts{}.tmp", my_rank, nts),
      std::format("vpackets_{:04d}_ts{}.tmp", my_rank, nts),
  };

  for (const auto& filename : filenames) {
    // these temp files are not needed again, so deletion failures are ignored (the error_code
    // overload of remove() reports the failure here instead of throwing)
    std::error_code remove_error;
    std::filesystem::remove(filename, remove_error);
  }
}

// read the virtual packet configuration from vpkt.txt: the observer directions, wavelength
// ranges, emission-time window, and the optical-depth and cell-selection limits
void read_vpktparameterfile() {
  FILE* input_file = fopen_required("vpkt.txt", "r");

  assert_always(fscanf(input_file, "%d", &nobsdirections) == 1);

  printlnlog("vpkt.txt: nobsdirections {}", nobsdirections);

  // nz_obs_vpkt. Cos(theta) to the observer. A list in the case of many observers
  obsdirs_costheta.resize(nobsdirections);
  for (int i = 0; i < nobsdirections; i++) {
    assert_always(fscanf(input_file, "%lg", &obsdirs_costheta[i]) == 1);

    if (fabs(obsdirs_costheta[i]) > 1) {
      printlnlog("[error] vpkt.txt observer direction {} has costheta {:g}, which is outside [-1, 1]. aborting", i,
                 obsdirs_costheta[i]);
      std::abort();
    }
    if (obsdirs_costheta[i] == 1) {
      obsdirs_costheta[i] = 0.9999;
    } else if (obsdirs_costheta[i] == -1) {
      obsdirs_costheta[i] = -0.9999;
    }
  }

  // phi to the observer (degrees). A list in the case of many observers
  obsdirs_phi.resize(nobsdirections);
  for (int i = 0; i < nobsdirections; i++) {
    double phi_degrees = 0.;
    assert_always(fscanf(input_file, "%lg", &phi_degrees) == 1);
    obsdirs_phi[i] = phi_degrees * PI / 180.;
    const double theta_degrees = std::acos(obsdirs_costheta[i]) / PI * 180.;

    printlnlog(
        "vpkt.txt:   direction {}: theta {:.1f} [degrees] (costheta {:g}), phi {:.1f} [degrees] ({:g} [radians])", i,
        theta_degrees, obsdirs_costheta[i], phi_degrees, obsdirs_phi[i]);
  }

  // Nspectra opacity choices (i.e. Nspectra spectra for each observer)
  int nspectra_customlist_flag = 0;
  assert_always(fscanf(input_file, "%d ", &nspectra_customlist_flag) == 1);

  if (nspectra_customlist_flag != 1) {
    nspectraperobsdir = 1;
    opacityexclusions.resize(nspectraperobsdir, 0);

    opacityexclusions[0] = 0;
  } else {
    assert_always(fscanf(input_file, "%d ", &nspectraperobsdir) == 1);
    opacityexclusions.resize(nspectraperobsdir, 0);

    for (int opacchoiceindex = 0; opacchoiceindex < nspectraperobsdir; opacchoiceindex++) {
      assert_always(fscanf(input_file, "%d ", &opacityexclusions[opacchoiceindex]) == 1);

      // The first number should be equal to zero!
      assert_always(opacityexclusions[0] == 0);  // The first spectrum should allow for all opacities (exclude[i]=0)
      assert_always(opacityexclusions[opacchoiceindex] <= 0 ||
                    !VPKT_USE_EXPANSION_OPACITIES);  // expansion opacities include all elements, so cannot be used
                                                     // with custom lists that exclude some elements
    }
  }

  printlnlog("vpkt.txt: Nspectra {} per observer", nspectraperobsdir);

  // time window. If dum4=1 it restrict vpkt to time windown (dum5,dum6)
  int override_tminmax = 0;
  double vspec_tmin_in_days = 0.;
  double vspec_tmax_in_days = 0.;
  assert_always(fscanf(input_file, "%d %lg %lg", &override_tminmax, &vspec_tmin_in_days, &vspec_tmax_in_days) == 3);

  printlnlog("vpkt: compiled with VSPEC_TIMEMIN {:.1f} [d] VSPEC_TIMEMAX {:.1f} [d] VSPEC_TIMEBINS {}",
             VSPEC_TIMEMIN / DAY, VSPEC_TIMEMAX / DAY, VSPEC_TIMEBINS);
  if (override_tminmax == 1) {
    vspec_timemin_input = vspec_tmin_in_days * DAY;
    vspec_timemax_input = vspec_tmax_in_days * DAY;
    printlnlog("vpkt.txt: VSPEC_TIMEMIN_input {:.1f} [d], VSPEC_TIMEMAX_input {:.1f} [d]", vspec_timemin_input / DAY,
               vspec_timemax_input / DAY);
  } else {
    vspec_timemin_input = VSPEC_TIMEMIN;
    vspec_timemax_input = VSPEC_TIMEMAX;
    printlnlog(
        "vpkt.txt: VSPEC_TIMEMIN_input {:.1f} [d], VSPEC_TIMEMAX_input {:.1f} [d] (inherited from VSPEC_TIMEMIN and "
        "VSPEC_TIMEMAX)",
        vspec_timemin_input / DAY, vspec_timemax_input / DAY);
  }

  assert_always(vspec_timemin_input >= VSPEC_TIMEMIN);
  assert_always(vspec_timemax_input <= VSPEC_TIMEMAX);
  assert_always(vspec_timemin_input >= globals::tmin);
  assert_always(vspec_timemax_input <= globals::tmax);

  // frequency window. dum4 restrict vpkt to a frequency range, dum5 indicates the number of ranges,
  // followed by a list of ranges (dum6,dum7)
  int flag_custom_freq_ranges = 0;
  assert_always(fscanf(input_file, "%d ", &flag_custom_freq_ranges) == 1);

  printlnlog("vpkt: compiled with VSPEC_NUBINS {}", VSPEC_NUBINS);
  assert_always(VSPEC_NUMAX > VSPEC_NUMIN);
  printlnlog("vpkt: compiled with VSPEC_NUMAX {:g} [Hz] lambda_min {:g} [Angstroms]", VSPEC_NUMAX,
             1e8 * CLIGHT / VSPEC_NUMAX);
  printlnlog("vpkt: compiled with VSPEC_NUMIN {:g} [Hz] lambda_max {:g} [Angstroms]", VSPEC_NUMIN,
             1e8 * CLIGHT / VSPEC_NUMIN);

  if (flag_custom_freq_ranges == 1) {
    assert_always(fscanf(input_file, "%d ", &nwavelengthranges) == 1);
    vspec_numin_input.resize(nwavelengthranges, 0.);
    vspec_numax_input.resize(nwavelengthranges, 0.);

    printlnlog("vpkt.txt: Nrange {} frequency intervals per spectrum per observer", nwavelengthranges);

    for (int i = 0; i < nwavelengthranges; i++) {
      double lmin_vspec_input = 0.;
      double lmax_vspec_input = 0.;
      assert_always(fscanf(input_file, "%lg %lg", &lmin_vspec_input, &lmax_vspec_input) == 2);

      vspec_numin_input[i] = CLIGHT / (lmax_vspec_input * 1e-8);
      vspec_numax_input[i] = CLIGHT / (lmin_vspec_input * 1e-8);

      assert_always(vspec_numin_input[i] >= VSPEC_NUMIN);
      assert_always(vspec_numax_input[i] <= VSPEC_NUMAX);
    }
  } else {
    nwavelengthranges = 1;

    vspec_numin_input.push_back(VSPEC_NUMIN);
    vspec_numax_input.push_back(VSPEC_NUMAX);

    printlnlog("vpkt.txt: Nrange 1 frequency interval (inherited from VSPEC_NUMIN and VSPEC_NUMAX)");
  }

  for (int i = 0; i < nwavelengthranges; i++) {
    printlnlog("vpkt.txt:   range {} lambda [{:g}, {:g}] [Angstroms]", i, 1e8 * CLIGHT / vspec_numax_input[i],
               1e8 * CLIGHT / vspec_numin_input[i]);
  }

  // if dum7=1, vpkt are not created when cell optical depth is larger than optical_depth_is_thick_vpkt
  int override_thickcell_tau = 0;
  assert_always(fscanf(input_file, "%d %lg", &override_thickcell_tau, &optical_depth_is_thick_vpkt) == 2);

  if (override_thickcell_tau == 1) {
    printlnlog("vpkt.txt: optical_depth_is_thick_vpkt {:g}", optical_depth_is_thick_vpkt);
  } else {
    optical_depth_is_thick_vpkt = globals::optical_depth_is_thick;
    printlnlog("vpkt.txt: optical_depth_is_thick_vpkt {:g} (inherited from optical_depth_is_thick)",
               optical_depth_is_thick_vpkt);
  }

  // Maximum optical depth. If a vpkt reaches dum7 is thrown away
  assert_always(fscanf(input_file, "%lg", &tau_max_vpkt) == 1);
  printlnlog("vpkt.txt: tau_max_vpkt {:g}", tau_max_vpkt);

  // Produce velocity grid map if =1
  int in_vgrid_on = 0;
  assert_always(fscanf(input_file, "%d", &in_vgrid_on) == 1);
  vgrid_on = in_vgrid_on != 0;
  printlnlog("vpkt.txt: velocity grid map {}", vgrid_on ? "ENABLED" : "DISABLED");

  if (vgrid_on) {
    double tmin_grid_in_days{NAN};
    double tmax_grid_in_days{NAN};
    // Specify time range for velocity grid map
    assert_always(fscanf(input_file, "%lg %lg", &tmin_grid_in_days, &tmax_grid_in_days) == 2);
    tmin_grid = tmin_grid_in_days * DAY;
    tmax_grid = tmax_grid_in_days * DAY;

    printlnlog("vpkt.txt: velocity grid time range tmin_grid {:g} [d] tmax_grid {:g} [d]", tmin_grid / DAY,
               tmax_grid / DAY);

    // Specify wavelength range: number of intervals (dum9) and limits (dum10,dum11)
    assert_always(fscanf(input_file, "%d ", &grid_nwavelengthranges) == 1);

    printlnlog("vpkt.txt: velocity grid frequency intervals {}", grid_nwavelengthranges);

    nu_grid_max.resize(grid_nwavelengthranges, 0.);
    nu_grid_min.resize(grid_nwavelengthranges, 0.);
    for (int i = 0; i < grid_nwavelengthranges; i++) {
      double range_lambda_min = 0.;
      double range_lambda_max = 0.;
      assert_always(fscanf(input_file, "%lg %lg", &range_lambda_min, &range_lambda_max) == 2);

      nu_grid_max[i] = CLIGHT / (range_lambda_min * 1e-8);
      nu_grid_min[i] = CLIGHT / (range_lambda_max * 1e-8);

      printlnlog("vpkt.txt:   velgrid range {} lambda [{:g}, {:g}] [Angstroms]", i, 1e8 * CLIGHT / nu_grid_max[i],
                 1e8 * CLIGHT / nu_grid_min[i]);
    }
  }

  fclose(input_file);
}

void write_timestep(const int nts, const bool is_final) {
  if constexpr (!VPKT_ON) {
    return;
  }
  const int my_rank = globals::my_rank;
  // write specpol of the virtual packets
  const auto filename_vspecpol =
      is_final ? std::format("vspecpol_{:04d}.out", my_rank) : std::format("vspecpol_{:04d}_ts{}.tmp", my_rank, nts);
  write_vspecpol(filename_vspecpol);

  if (vgrid_on) {
    const auto filename_vpktgrid = is_final ? std::format("vpkt_grid_{:04d}.out", my_rank)
                                            : std::format("vpkt_grid_{:04d}_ts{}.tmp", my_rank, nts);
    write_vpkt_grid(filename_vpktgrid);
  }

  if constexpr (VPKT_WRITE_CONTRIBS) {
    vpkt_contrib_file.close();
    const auto filename_source = std::format("vpackets_{:04d}_ts{}.tmp", my_rank, is_final ? nts + 1 : nts);
    const auto filename_dest = is_final ? std::format("vpackets_{:04d}.out", my_rank)
                                        : std::format("vpackets_{:04d}_ts{}.tmp", my_rank, nts + 1);

    std::filesystem::copy_file(filename_source, filename_dest, std::filesystem::copy_options::overwrite_existing);

    if (!is_final) {
      vpkt_contrib_file = std::ofstream(filename_dest, std::ios::app);
    }
  }

  printlnlog("timestep {}: wrote {}vpkt output files", nts, is_final ? "final " : "temporary ");
}

void init(const int nts, const bool continued_from_saved) {
  if constexpr (!VPKT_ON) {
    return;
  }

  init_vspecpol();
  if (vgrid_on) {
    init_vpkt_grid();
  }

  if constexpr (VPKT_WRITE_CONTRIBS) {
    const auto filename = std::format("vpackets_{:04d}_ts{}.tmp", globals::my_rank, nts + 1);

    if (continued_from_saved) {
      const auto filename_prev = std::format("vpackets_{:04d}_ts{}.tmp", globals::my_rank, nts);
      printlnlog("Copying {} to {}", filename_prev, filename);
      std::filesystem::copy_file(filename_prev, filename, std::filesystem::copy_options::overwrite_existing);
    } else {
      // Create new file with header line
      vpkt_contrib_file = std::ofstream(filename, std::ios::trunc);

      std::print(vpkt_contrib_file, "#emissiontype trueemissiontype absorption_type absorption_freq");

      for (int obsdirindex = 0; obsdirindex < nobsdirections; obsdirindex++) {
        std::print(vpkt_contrib_file, " dir{}_t_arrive_d dir{}_nu_rf", obsdirindex, obsdirindex);
        for (int opacchoiceindex = 0; opacchoiceindex < nspectraperobsdir; opacchoiceindex++) {
          std::print(vpkt_contrib_file, " dir{}_e_rf_{}", obsdirindex, opacchoiceindex);
        }
      }

      std::println(vpkt_contrib_file, "");
      vpkt_contrib_file.flush();
      vpkt_contrib_file.close();
    }

    vpkt_contrib_file = std::ofstream(filename, std::ios::app);
  }

  if (continued_from_saved) {
    // Continue simulation: read into temporary files

    read_vspecpol(globals::my_rank, nts);

    if (vgrid_on) {
      read_vpkt_grid(globals::my_rank, nts);
    }
  }
}

auto trace_vpkts(const Packet& pkt, const enum packet_type type_before_rpkt) -> void {
  assert_testmodeonly(VPKT_ON);

  // Cut on vpkts
  const auto nonemptymgi = grid::get_propcell_nonemptymgi(pkt.cellindex);

  if (grid::thick_allcells[nonemptymgi] != 0) {
    return;
  }

  THREADLOCALONHOST std::string vpkt_contrib_row;
  vpkt_contrib_row.clear();

  bool any_dir_escaped = false;
  for (int obsdirindex = 0; obsdirindex < nobsdirections; obsdirindex++) {
    // loop over different observer directions

    const auto obsdir = Vec3d{
        sqrt(1 - (obsdirs_costheta[obsdirindex] * obsdirs_costheta[obsdirindex])) * cos(obsdirs_phi[obsdirindex]),
        sqrt(1 - (obsdirs_costheta[obsdirindex] * obsdirs_costheta[obsdirindex])) * sin(obsdirs_phi[obsdirindex]),
        obsdirs_costheta[obsdirindex],
    };

    const double t_arrive = pkt.prop_time - (dot(pkt.pos, obsdir) / CLIGHT_PROP);

    bool dir_escaped = false;
    if (t_arrive >= vspec_timemin_input && t_arrive <= vspec_timemax_input) {
      // time selection

      const double doppler = calculate_doppler_nucmf_on_nurf(pkt.pos, obsdir, pkt.prop_time);
      const double nu_rf = pkt.nu_cmf / doppler;
      const double e_rf = pkt.e_cmf / doppler;

      // Loop over frequency intervals for this observer and check if the vpkt frequency falls in any of them. If it
      // does, trace the vpkt to see if it escapes in this direction.
      for (int i = 0; i < nwavelengthranges; i++) {
        if ((nu_rf > vspec_numin_input[i] && nu_rf < vspec_numax_input[i]) ||
            (pkt.absorptionfreq > vspec_numin_input[i] && pkt.absorptionfreq < vspec_numax_input[i])) {
          // frequency selection
          dir_escaped = dir_escaped || trace_vpkt_direction(pkt, t_arrive, nu_rf, e_rf, doppler, obsdirindex, obsdir,
                                                            type_before_rpkt, vpkt_contrib_row);
          break;  // we only need to match one frequency interval to trace the vpkt
        }
      }
    }

    if (dir_escaped) {
      any_dir_escaped = true;
    } else {
      vpkt_contrib_row += " -1. -1.";  // t_arrive_d nu_rf
      for (int opacchoiceindex = 0; opacchoiceindex < nspectraperobsdir; opacchoiceindex++) {
        vpkt_contrib_row += " 0.";  // e_rf_diri_j
      }
    }
  }
  if (VPKT_WRITE_CONTRIBS && any_dir_escaped) {
    // serialise writes so that rows from concurrent threads cannot interleave
    [[maybe_unused]] const ScopedMutex lock{vpkt_contrib_file_mutex};
    std::println(vpkt_contrib_file, "{} {} {} {:g}{}", pkt.emissiontype, pkt.trueemissiontype, pkt.absorptiontype,
                 pkt.absorptionfreq, vpkt_contrib_row);
  }
}
}  // namespace vpkt

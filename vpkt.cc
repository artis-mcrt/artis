#include "vpkt.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <ios>
#include <sstream>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "ltepop.h"
#include "packet.h"
#include "rpkt.h"
#include "sn3d.h"
#include "vectors.h"

namespace {

struct StokesParams {
  double i = 0.;
  double q = 0.;
  double u = 0.;
};

struct VSpecPol {
  std::array<StokesParams, VMNUBINS> flux;
  float lower_time{NAN};
  float delta_t{NAN};
};

std::vector<std::vector<VSpecPol>> vspecpol{};

std::array<float, VMNUBINS> lower_freq_vspec;
std::array<float, VMNUBINS> delta_freq_vspec;

// --------- INPUT PARAMETERS -----------

int Nobs = 0;      // Number of observer directions
int Nspectra = 0;  // Number of virtual packet spectra per observer direction (total + elements switched off)
std::vector<double> nz_obs_vpkt;
std::vector<double> phiobs;
double VSPEC_TIMEMIN_input;
double VSPEC_TIMEMAX_input;
int Nrange = 0;  // Number of wavelength ranges

std::vector<double> VSPEC_NUMIN_input;
std::vector<double> VSPEC_NUMAX_input;
double tau_max_vpkt;

std::vector<int> exclude;  // vector of opacity contribution setups:
                           // 0: full opacity
                           // -1: no line opacity; -2: no bf opacity; -3: no ff opacity; -4: no es opacity,
                           // +ve: exclude element with atomic number's contribution to bound-bound opacity
std::vector<double> tau_vpkt;

std::ofstream vpkt_contrib_file;

// --------- VPacket GRID -----------

struct VGrid {
  std::vector<std::vector<StokesParams>> flux;
  double yvel{NAN};
  double zvel{NAN};
};

std::array<std::array<VGrid, VGRID_NZ>, VGRID_NY> vgrid;

int Nrange_grid;
double tmin_grid;
double tmax_grid;
std::vector<double> nu_grid_min;
std::vector<double> nu_grid_max;
bool vgrid_on;

const double dlogt_vspec = (std::log(VSPEC_TIMEMAX) - std::log(VSPEC_TIMEMIN)) / VMTBINS;
const double dlognu_vspec = (std::log(VSPEC_NUMAX) - std::log(VSPEC_NUMIN)) / VMNUBINS;

// Virtual packet is killed when tau reaches tau_max_vpkt for ALL the different setups
// E.g. imagine that a packet in the first setup (all elements included) reaches tau = tau_max_vpkt
// because of the element Zi. If we remove Zi, tau now could be lower than tau_max_vpkt and could
// thus contribute to the spectrum.
constexpr auto all_taus_past_taumax(std::vector<double> &tau, const double tau_max) -> bool {
  return std::ranges::all_of(tau, [tau_max](const double tau_i) { return tau_i > tau_max; });
}

// Routine to add a packet to the outcoming spectrum.
void add_to_vspecpol(const Packet &vpkt, const int obsdirindex, const int opachoiceindex, const double t_arrive) {
  // Need to decide in which (1) time and (2) frequency bin the vpkt is escaping

  const int nt = static_cast<int>((log(t_arrive) - log(VSPEC_TIMEMIN)) / dlogt_vspec);
  const int nnu = static_cast<int>((log(vpkt.nu_rf) - log(VSPEC_NUMIN)) / dlognu_vspec);
  if (nt < 0 || nt >= VMTBINS || nnu < 0 || nnu >= VMNUBINS) {
    return;
  }

  const int ind_comb = (Nspectra * obsdirindex) + opachoiceindex;
  const double pktcontrib = vpkt.e_rf / vspecpol[nt][ind_comb].delta_t / delta_freq_vspec[nnu] / 4.e12 / PI / PARSEC /
                            PARSEC / globals::nprocs * 4 * PI;

  atomicadd(vspecpol[nt][ind_comb].flux[nnu].i, vpkt.stokes[0] * pktcontrib);
  atomicadd(vspecpol[nt][ind_comb].flux[nnu].q, vpkt.stokes[1] * pktcontrib);
  atomicadd(vspecpol[nt][ind_comb].flux[nnu].u, vpkt.stokes[2] * pktcontrib);
}

// Routine to add a packet to the outcoming spectrum.
void add_to_vpkt_grid(const Packet &vpkt, const std::array<double, 3> &vel, const int wlbin, const int obsdirindex,
                      const std::array<double, 3> &obs) {
  double vref1{NAN};
  double vref2{NAN};

  // obs is the observer orientation

  // Packet velocity

  // if nobs = x , vref1 = vy and vref2 = vz
  if (obs[0] == 1) {
    vref1 = vel[1];
    vref2 = vel[2];
  }
  // if nobs = -x , vref1 = -vy and vref2 = -vz
  else if (obs[0] == -1) {
    vref1 = -vel[1];
    vref2 = -vel[2];
  }

  // Rotate velocity into projected area seen by the observer (see notes)
  else {
    // Rotate velocity from (x,y,z) to (n_obs,ref1,ref2) so that x correspond to n_obs (see notes)
    vref1 = -obs[1] * vel[0] + (obs[0] + obs[2] * obs[2] / (1 + obs[0])) * vel[1] -
            obs[1] * obs[2] * (1 - obs[0]) / sqrt(1 - (obs[0] * obs[0])) * vel[2];
    vref2 = -obs[2] * vel[0] - obs[1] * obs[2] * (1 - obs[0]) / sqrt(1 - (obs[0] * obs[0])) * vel[1] +
            (obs[0] + obs[1] * obs[1] / (1 + obs[0])) * vel[2];
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
  if (vpkt.nu_rf > nu_grid_min[wlbin] && vpkt.nu_rf < nu_grid_max[wlbin]) {
    atomicadd(vgrid[ny][nz].flux[wlbin][obsdirindex].i, vpkt.stokes[0] * vpkt.e_rf);
    atomicadd(vgrid[ny][nz].flux[wlbin][obsdirindex].q, vpkt.stokes[1] * vpkt.e_rf);
    atomicadd(vgrid[ny][nz].flux[wlbin][obsdirindex].u, vpkt.stokes[2] * vpkt.e_rf);
  }
}

auto rlc_emiss_vpkt(const Packet &pkt, const double t_current, const double t_arrive, const double nu_rf,
                    const double e_rf, const int obsdirindex, const std::array<double, 3> &obsdir,
                    const enum packet_type type_before_rpkt, std::stringstream &vpkt_contrib_row) -> bool {
  int mgi = 0;

  Packet vpkt = pkt;
  vpkt.nu_rf = nu_rf;
  vpkt.e_rf = e_rf;
  vpkt.dir = obsdir;
  vpkt.last_cross = BOUNDARY_NONE;

  bool end_packet = false;
  double ldist = 0;
  double t_future = t_current;

  for (int opacindex = 0; opacindex < Nspectra; opacindex++) {
    tau_vpkt[opacindex] = 0;
  }

  atomicadd(nvpkt, 1);  // increment the number of virtual packet in the given timestep

  const auto vel_vec = get_velocity(pkt.pos, t_current);
  double Qi = vpkt.stokes[1];
  double Ui = vpkt.stokes[2];

  // ------------ SCATTERING EVENT: dipole function --------------------

  double pn{NAN};
  constexpr double I = 1.;
  double Q{NAN};
  double U{NAN};
  if (type_before_rpkt == TYPE_RPKT) {
    // Transform Stokes Parameters from the RF to the CMF

    const auto old_dir_cmf = frame_transform(pkt.dir, &Qi, &Ui, vel_vec);

    // Need to rotate Stokes Parameters in the scattering plane

    const auto obs_cmf = angle_ab(vpkt.dir, vel_vec);

    const auto [ref1_old, ref2_old] = meridian(old_dir_cmf);

    // This is the i1 angle of Bulla+2015, obtained by computing the angle between the
    // reference axes ref1 and ref2 in the meridian frame and the corresponding axes
    // ref1_sc and ref2_sc in the scattering plane.
    const double i1 = get_rot_angle(old_dir_cmf, obs_cmf, ref1_old, ref2_old);
    const double cos2i1 = cos(2 * i1);
    const double sin2i1 = sin(2 * i1);

    const double Qold = (Qi * cos2i1) - (Ui * sin2i1);
    const double Uold = (Qi * sin2i1) + (Ui * cos2i1);

    // Scattering

    const double mu = dot(old_dir_cmf, obs_cmf);

    pn = 3. / (16. * PI) * (1 + pow(mu, 2.) + (pow(mu, 2.) - 1) * Qold);

    const double Inew = 0.75 * ((mu * mu + 1.0) + Qold * (mu * mu - 1.0));
    const double Qnew = (0.75 * ((mu * mu - 1.0) + Qold * (mu * mu + 1.0))) / Inew;
    const double Unew = (1.5 * mu * Uold) / Inew;

    // Need to rotate Stokes Parameters out of the scattering plane to the meridian frame

    const auto [ref1, ref2] = meridian(obs_cmf);

    // This is the i2 angle of Bulla+2015, obtained from the angle THETA between the
    // reference axes ref1_sc and ref2_sc in the scattering plane and ref1 and ref2 in the
    // meridian frame. NB: we need to add PI to transform THETA to i2
    const double i2 = PI + get_rot_angle(obs_cmf, old_dir_cmf, ref1, ref2);
    const double cos2i2 = cos(2 * i2);
    const double sin2i2 = sin(2 * i2);

    Q = Qnew * cos2i2 + Unew * sin2i2;
    U = -Qnew * sin2i2 + Unew * cos2i2;

    // Transform Stokes Parameters from the CMF to the RF

    const auto vel_rev = std::array<double, 3>{-vel_vec[0], -vel_vec[1], -vel_vec[2]};

    frame_transform(obs_cmf, &Q, &U, vel_rev);

  } else if (type_before_rpkt == TYPE_KPKT || type_before_rpkt == TYPE_MA) {
    // MACROATOM and KPKT: isotropic emission
    Q = 0;
    U = 0;
    pn = 1 / (4 * PI);
  }

  // compute the optical depth to boundary

  mgi = grid::get_cell_modelgridindex(vpkt.where);
  Rpkt_continuum_absorptioncoeffs chi_vpkt_cont{};

  while (!end_packet) {
    // distance to the next cell
    const auto [sdist, snext] =
        grid::boundary_distance(vpkt.dir, vpkt.pos, vpkt.prop_time, vpkt.where, &vpkt.last_cross);
    const double s_cont = sdist * t_current * t_current * t_current / (t_future * t_future * t_future);

    if (mgi == grid::get_npts_model()) {
      vpkt.next_trans = -1;
    } else {
      calculate_chi_rpkt_cont(vpkt.nu_cmf, chi_vpkt_cont, mgi);

      const double chi_cont = chi_vpkt_cont.total;

      for (int ind = 0; ind < Nspectra; ind++) {
        if (exclude[ind] == -2) {
          const double chi_cont_nobf = chi_cont - chi_vpkt_cont.bf;
          tau_vpkt[ind] += chi_cont_nobf * s_cont;
        } else if (exclude[ind] == -3) {
          const double chi_cont_noff = chi_cont - chi_vpkt_cont.ffheat;
          tau_vpkt[ind] += chi_cont_noff * s_cont;
        } else if (exclude[ind] == -4) {
          const double chi_cont_noes = chi_cont - chi_vpkt_cont.ffescat;
          tau_vpkt[ind] += chi_cont_noes * s_cont;
        } else {
          tau_vpkt[ind] += chi_cont * s_cont;
        }
      }

      // kill vpkt with high optical depth
      if (all_taus_past_taumax(tau_vpkt, tau_max_vpkt)) {
        return false;
      }

      Packet dummypkt_abort = vpkt;
      move_pkt_withtime(dummypkt_abort, sdist);
      const double nu_cmf_abort = dummypkt_abort.nu_cmf;
      assert_testmodeonly(nu_cmf_abort <= vpkt.nu_cmf);
      const double d_nu_on_d_l = (nu_cmf_abort - vpkt.nu_cmf) / sdist;

      ldist = 0;
      while (ldist < sdist) {
        const int lineindex = closest_transition(vpkt.nu_cmf, vpkt.next_trans);

        if (lineindex < 0) {
          // no more lines below the current frequency
          vpkt.next_trans = globals::nlines + 1;
          break;
        }
        const double nutrans = globals::linelist[lineindex].nu;

        vpkt.next_trans = lineindex + 1;

        ldist = get_linedistance(vpkt.prop_time, vpkt.nu_cmf, nutrans, d_nu_on_d_l);

        if (ldist > sdist) {
          // exit the while loop if you reach the boundary; go back to the previous transition to start next cell with
          // the excluded line

          vpkt.next_trans -= 1;
          // printout("ldist > sdist : line in the next cell\n");
          break;
        }

        const double t_line = vpkt.prop_time + (ldist / CLIGHT);

        const int element = globals::linelist[lineindex].elementindex;
        const int ion = globals::linelist[lineindex].ionindex;
        const int upper = globals::linelist[lineindex].upperlevelindex;
        const int lower = globals::linelist[lineindex].lowerlevelindex;
        const auto A_ul = globals::linelist[lineindex].einstein_A;

        const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nutrans, 3) * A_ul;
        const double B_lu = stat_weight(element, ion, upper) / stat_weight(element, ion, lower) * B_ul;

        const auto n_u = calculate_levelpop(mgi, element, ion, upper);
        const auto n_l = calculate_levelpop(mgi, element, ion, lower);
        const double tau_line = std::max(0., (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_line);

        // Check on the element to exclude (or -1 for no line opacity)
        const int anumber = get_atomicnumber(element);
        for (int ind = 0; ind < Nspectra; ind++) {
          if (exclude[ind] != -1 && (exclude[ind] != anumber)) {
            tau_vpkt[ind] += tau_line;
          }
        }

        // kill vpkt with high optical depth
        if (all_taus_past_taumax(tau_vpkt, tau_max_vpkt)) {
          return false;
        }
      }
    }

    // virtual packet is still at the starting position
    // move it to cell boundary and go to next cell

    t_future += (sdist / CLIGHT_PROP);
    move_pkt_withtime(vpkt, sdist);
    vpkt.prop_time = t_future;

    grid::change_cell(vpkt, snext);
    end_packet = (vpkt.type == TYPE_ESCAPE);

    mgi = grid::get_cell_modelgridindex(vpkt.where);

    // kill vpkt with pass through a thick cell
    if (grid::modelgrid[mgi].thick != 0) {
      return false;
    }
  }

  // increment the number of escaped virtual packet in the given timestep
  if (type_before_rpkt == TYPE_RPKT) {
    atomicadd(nvpkt_esc1, 1);
  } else if (type_before_rpkt == TYPE_KPKT) {
    atomicadd(nvpkt_esc2, 1);
  } else if (type_before_rpkt == TYPE_MA) {
    atomicadd(nvpkt_esc3, 1);
  }

  // -------------- final stokes vector ---------------

  if (VPKT_WRITE_CONTRIBS) {
    vpkt_contrib_row << " " << t_arrive / DAY << " " << vpkt.nu_rf;
  }

  for (int ind = 0; ind < Nspectra; ind++) {
    const double prob = pn * std::exp(-tau_vpkt[ind]);

    assert_always(std::isfinite(prob));

    vpkt.stokes = {I * prob, Q * prob, U * prob};

    for (const auto stokeval : vpkt.stokes) {
      assert_always(std::isfinite(stokeval));
    }

    add_to_vspecpol(vpkt, obsdirindex, ind, t_arrive);

    if constexpr (VPKT_WRITE_CONTRIBS) {
      vpkt_contrib_row << " " << vpkt.e_rf * prob;
    }
  }

  // vpkt grid

  if (vgrid_on) {
    const double prob = pn * exp(-tau_vpkt[0]);

    vpkt.stokes = {I * prob, Q * prob, U * prob};

    for (int wlbin = 0; wlbin < Nrange_grid; wlbin++) {
      if (vpkt.nu_rf > nu_grid_min[wlbin] && vpkt.nu_rf < nu_grid_max[wlbin]) {  // Frequency selection
        if (t_arrive > tmin_grid && t_arrive < tmax_grid) {                      // Time selection
          add_to_vpkt_grid(vpkt, vel_vec, wlbin, obsdirindex, obsdir);
        }
      }
    }
  }
  return true;  // true if we added columns to vpkt_contrib_row
}

void init_vspecpol() {
  vspecpol.resize(VMTBINS);

  const int indexmax = Nspectra * Nobs;
  for (int p = 0; p < VMTBINS; p++) {
    vspecpol[p].resize(indexmax);
  }

  for (int m = 0; m < VMNUBINS; m++) {
    lower_freq_vspec[m] = exp(log(VSPEC_NUMIN) + (m * (dlognu_vspec)));
    delta_freq_vspec[m] = exp(log(VSPEC_NUMIN) + ((m + 1) * (dlognu_vspec))) - lower_freq_vspec[m];
  }

  // start by setting up the time and frequency bins.
  // it is all done interms of a logarithmic spacing in both t and nu - get the
  // step sizes first.
  for (int n = 0; n < VMTBINS; n++) {
    for (int ind_comb = 0; ind_comb < indexmax; ind_comb++) {
      vspecpol[n][ind_comb].lower_time = exp(log(VSPEC_TIMEMIN) + (n * (dlogt_vspec)));
      vspecpol[n][ind_comb].delta_t =
          exp(log(VSPEC_TIMEMIN) + ((n + 1) * (dlogt_vspec))) - vspecpol[n][ind_comb].lower_time;

      for (auto &flux : vspecpol[n][ind_comb].flux) {
        flux.i = 0.;
        flux.q = 0.;
        flux.u = 0.;
      }
    }
  }
}

void write_vspecpol(FILE *specpol_file) {
  for (int ind_comb = 0; ind_comb < (Nobs * Nspectra); ind_comb++) {
    fprintf(specpol_file, "%g ", 0.);

    for (int l = 0; l < 3; l++) {
      for (int p = 0; p < VMTBINS; p++) {
        fprintf(specpol_file, "%g ", (vspecpol[p][ind_comb].lower_time + (vspecpol[p][ind_comb].delta_t / 2.)) / DAY);
      }
    }

    fprintf(specpol_file, "\n");

    for (int m = 0; m < VMNUBINS; m++) {
      fprintf(specpol_file, "%g ", (lower_freq_vspec[m] + (delta_freq_vspec[m] / 2.)));

      // Stokes I
      for (int p = 0; p < VMTBINS; p++) {
        fprintf(specpol_file, "%g ", vspecpol[p][ind_comb].flux[m].i);
      }

      // Stokes Q
      for (int p = 0; p < VMTBINS; p++) {
        fprintf(specpol_file, "%g ", vspecpol[p][ind_comb].flux[m].q);
      }

      // Stokes U
      for (int p = 0; p < VMTBINS; p++) {
        fprintf(specpol_file, "%g ", vspecpol[p][ind_comb].flux[m].u);
      }

      fprintf(specpol_file, "\n");
    }
  }
}

void read_vspecpol(const int my_rank, const int nts) {
  char filename[MAXFILENAMELENGTH];

  snprintf(filename, MAXFILENAMELENGTH, "vspecpol_%.4d_ts%d.tmp", my_rank, nts);
  printout("Reading %s\n", filename);

  FILE *vspecpol_file = fopen_required(filename, "r");

  float a{NAN};
  float b{NAN};
  float c{NAN};

  for (int ind_comb = 0; ind_comb < (Nobs * Nspectra); ind_comb++) {
    // Initialise I,Q,U fluxes from temporary files
    assert_always(fscanf(vspecpol_file, "%g ", &a) == 1);

    for (int l = 0; l < 3; l++) {
      for (int p = 0; p < VMTBINS; p++) {
        assert_always(fscanf(vspecpol_file, "%g ", &b) == 1);
      }
    }

    assert_always(fscanf(vspecpol_file, "\n") == 0);

    for (int j = 0; j < VMNUBINS; j++) {
      assert_always(fscanf(vspecpol_file, "%g ", &c) == 1);

      // Stokes I
      for (int p = 0; p < VMTBINS; p++) {
        assert_always(fscanf(vspecpol_file, "%lg ", &vspecpol[p][ind_comb].flux[j].i) == 1);
      }

      // Stokes Q
      for (int p = 0; p < VMTBINS; p++) {
        assert_always(fscanf(vspecpol_file, "%lg ", &vspecpol[p][ind_comb].flux[j].q) == 1);
      }

      // Stokes U
      for (int p = 0; p < VMTBINS; p++) {
        assert_always(fscanf(vspecpol_file, "%lg ", &vspecpol[p][ind_comb].flux[j].u) == 1);
      }

      assert_always(fscanf(vspecpol_file, "\n") == 0);
    }
  }

  fclose(vspecpol_file);
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

      vgrid[n][m].flux.resize(Nrange_grid, {});

      for (int wlbin = 0; wlbin < Nrange_grid; wlbin++) {
        vgrid[n][m].flux[wlbin].resize(Nobs, {0., 0., 0.});
      }
    }
  }
}

void write_vpkt_grid(FILE *vpkt_grid_file) {
  for (int obsdirindex = 0; obsdirindex < Nobs; obsdirindex++) {
    for (int wlbin = 0; wlbin < Nrange_grid; wlbin++) {
      for (int n = 0; n < VGRID_NY; n++) {
        for (int m = 0; m < VGRID_NZ; m++) {
          fprintf(vpkt_grid_file, "%g %g %g %g %g \n", vgrid[n][m].yvel, vgrid[n][m].zvel,
                  vgrid[n][m].flux[wlbin][obsdirindex].i, vgrid[n][m].flux[wlbin][obsdirindex].q,
                  vgrid[n][m].flux[wlbin][obsdirindex].u);
        }
      }
    }
  }
}

void read_vpkt_grid(const int my_rank, const int nts) {
  if (!vgrid_on) {
    return;
  }

  char filename[MAXFILENAMELENGTH];
  snprintf(filename, MAXFILENAMELENGTH, "vpkt_grid_%.4d_ts%d.tmp", my_rank, nts);
  printout("Reading vpkt grid file %s\n", filename);
  FILE *vpkt_grid_file = fopen_required(filename, "r");

  for (int obsdirindex = 0; obsdirindex < Nobs; obsdirindex++) {
    for (int wlbin = 0; wlbin < Nrange_grid; wlbin++) {
      for (int n = 0; n < VGRID_NY; n++) {
        for (int m = 0; m < VGRID_NZ; m++) {
          assert_always(fscanf(vpkt_grid_file, "%lg %lg %lg %lg %lg \n", &vgrid[n][m].yvel, &vgrid[n][m].zvel,
                               &vgrid[n][m].flux[wlbin][obsdirindex].i, &vgrid[n][m].flux[wlbin][obsdirindex].q,
                               &vgrid[n][m].flux[wlbin][obsdirindex].u) == 5);
        }
      }
    }
  }

  fclose(vpkt_grid_file);
}

}  // anonymous namespace

void vpkt_remove_temp_file(const int nts, const int my_rank) {
  std::array<char[MAXFILENAMELENGTH], 3> filenames{};
  snprintf(filenames[0], MAXFILENAMELENGTH, "vspecpol_%.4d_ts%d.tmp", my_rank, nts);
  snprintf(filenames[1], MAXFILENAMELENGTH, "vpkt_grid_%.4d_ts%d.tmp", my_rank, nts);
  snprintf(filenames[2], MAXFILENAMELENGTH, "vpackets_%.4d_ts%d.tmp", my_rank, nts);

  for (const auto *filename : filenames) {
    if (std::filesystem::exists(filename)) {
      std::remove(filename);
      printout("Deleted %s\n", filename);
    }
  }
}

void read_parameterfile_vpkt() {
  FILE *input_file = fopen_required("vpkt.txt", "r");

  // Nobs
  assert_always(fscanf(input_file, "%d", &Nobs) == 1);

  printout("vpkt.txt: Nobs %d directions\n", Nobs);

  // nz_obs_vpkt. Cos(theta) to the observer. A list in the case of many observers
  nz_obs_vpkt.resize(Nobs);
  for (int i = 0; i < Nobs; i++) {
    assert_always(fscanf(input_file, "%lg", &nz_obs_vpkt[i]) == 1);

    if (fabs(nz_obs_vpkt[i]) > 1) {
      printout("Wrong observer direction\n");
      std::abort();
    } else if (nz_obs_vpkt[i] == 1) {
      nz_obs_vpkt[i] = 0.9999;
    } else if (nz_obs_vpkt[i] == -1) {
      nz_obs_vpkt[i] = -0.9999;
    }
  }

  // phi to the observer (degrees). A list in the case of many observers
  phiobs.resize(Nobs);
  for (int i = 0; i < Nobs; i++) {
    double phi_degrees = 0.;
    assert_always(fscanf(input_file, "%lg \n", &phi_degrees) == 1);
    phiobs[i] = phi_degrees * PI / 180.;
    const double theta_degrees = std::acos(nz_obs_vpkt[i]) / PI * 180.;

    printout("vpkt.txt:   direction %d costheta %g (%.1f degrees) phi %g (%.1f degrees)\n", i, nz_obs_vpkt[i],
             theta_degrees, phiobs[i], phi_degrees);
  }

  // Nspectra opacity choices (i.e. Nspectra spectra for each observer)
  int nspectra_customlist_flag = 0;
  assert_always(fscanf(input_file, "%d ", &nspectra_customlist_flag) == 1);

  if (nspectra_customlist_flag != 1) {
    Nspectra = 1;
    exclude.resize(Nspectra, 0);

    exclude[0] = 0;
  } else {
    assert_always(fscanf(input_file, "%d ", &Nspectra) == 1);
    exclude.resize(Nspectra, 0);

    for (int i = 0; i < Nspectra; i++) {
      assert_always(fscanf(input_file, "%d ", &exclude[i]) == 1);

      // The first number should be equal to zero!
      assert_always(exclude[0] == 0);  // The first spectrum should allow for all opacities (exclude[i]=0)
    }
  }

  printout("vpkt.txt: Nspectra %d per observer\n", Nspectra);
  tau_vpkt.resize(Nspectra, 0.);

  // time window. If dum4=1 it restrict vpkt to time windown (dum5,dum6)
  int override_tminmax = 0;
  double vspec_tmin_in_days = 0.;
  double vspec_tmax_in_days = 0.;
  assert_always(fscanf(input_file, "%d %lg %lg \n", &override_tminmax, &vspec_tmin_in_days, &vspec_tmax_in_days) == 3);

  printout("vpkt: compiled with VSPEC_TIMEMIN %.1fd VSPEC_TIMEMAX %1.fd VMTBINS %d\n", VSPEC_TIMEMIN / DAY,
           VSPEC_TIMEMAX / DAY, VMTBINS);
  if (override_tminmax == 1) {
    VSPEC_TIMEMIN_input = vspec_tmin_in_days * DAY;
    VSPEC_TIMEMAX_input = vspec_tmax_in_days * DAY;
    printout("vpkt.txt: VSPEC_TIMEMIN_input %.1fd, VSPEC_TIMEMAX_input %.1fd\n", VSPEC_TIMEMIN_input / DAY,
             VSPEC_TIMEMAX_input / DAY);
  } else {
    VSPEC_TIMEMIN_input = VSPEC_TIMEMIN;
    VSPEC_TIMEMAX_input = VSPEC_TIMEMAX;
    printout(
        "vpkt.txt: VSPEC_TIMEMIN_input %.1fd, VSPEC_TIMEMAX_input %.1fd (inherited from VSPEC_TIMEMIN and "
        "VSPEC_TIMEMAX)\n",
        VSPEC_TIMEMIN_input / DAY, VSPEC_TIMEMAX_input / DAY);
  }

  assert_always(VSPEC_TIMEMIN_input >= VSPEC_TIMEMIN);
  assert_always(VSPEC_TIMEMAX_input <= VSPEC_TIMEMAX);
  assert_always(VSPEC_TIMEMIN_input >= globals::tmin);
  assert_always(VSPEC_TIMEMAX_input <= globals::tmax);

  // frequency window. dum4 restrict vpkt to a frequency range, dum5 indicates the number of ranges,
  // followed by a list of ranges (dum6,dum7)
  int flag_custom_freq_ranges = 0;
  assert_always(fscanf(input_file, "%d ", &flag_custom_freq_ranges) == 1);

  printout("vpkt: compiled with VMNUBINS %d\n", VMNUBINS);
  assert_always(VSPEC_NUMAX > VSPEC_NUMIN);
  printout("vpkt: compiled with VSPEC_NUMAX %g lambda_min %g Å\n", VSPEC_NUMAX, 1e8 * CLIGHT / VSPEC_NUMAX);
  printout("vpkt: compiled with VSPEC_NUMIN %g lambda_max %g Å\n", VSPEC_NUMIN, 1e8 * CLIGHT / VSPEC_NUMIN);

  if (flag_custom_freq_ranges == 1) {
    assert_always(fscanf(input_file, "%d ", &Nrange) == 1);
    VSPEC_NUMIN_input.resize(Nrange, 0.);
    VSPEC_NUMAX_input.resize(Nrange, 0.);

    printout("vpkt.txt: Nrange %d frequency intervals per spectrum per observer\n", Nrange);

    for (int i = 0; i < Nrange; i++) {
      double lmin_vspec_input = 0.;
      double lmax_vspec_input = 0.;
      assert_always(fscanf(input_file, "%lg %lg", &lmin_vspec_input, &lmax_vspec_input) == 2);

      VSPEC_NUMIN_input[i] = CLIGHT / (lmax_vspec_input * 1e-8);
      VSPEC_NUMAX_input[i] = CLIGHT / (lmin_vspec_input * 1e-8);

      assert_always(VSPEC_NUMIN_input[i] >= VSPEC_NUMIN);
      assert_always(VSPEC_NUMAX_input[i] <= VSPEC_NUMAX);
    }
  } else {
    Nrange = 1;

    VSPEC_NUMIN_input.push_back(VSPEC_NUMIN);
    VSPEC_NUMAX_input.push_back(VSPEC_NUMAX);

    printout("vpkt.txt: Nrange 1 frequency interval (inherited from VSPEC_NUMIN and VSPEC_NUMAX)\n");
  }

  for (int i = 0; i < Nrange; i++) {
    printout("vpkt.txt:   range %d lambda [%g, %g] Angstroms\n", i, 1e8 * CLIGHT / VSPEC_NUMAX_input[i],
             1e8 * CLIGHT / VSPEC_NUMIN_input[i]);
  }

  // if dum7=1, vpkt are not created when cell optical depth is larger than cell_is_optically_thick_vpkt
  int overrride_thickcell_tau = 0;
  assert_always(fscanf(input_file, "%d %lg \n", &overrride_thickcell_tau, &cell_is_optically_thick_vpkt) == 2);

  if (overrride_thickcell_tau == 1) {
    printout("vpkt.txt: cell_is_optically_thick_vpkt %lg\n", cell_is_optically_thick_vpkt);
  } else {
    cell_is_optically_thick_vpkt = globals::cell_is_optically_thick;
    printout("vpkt.txt: cell_is_optically_thick_vpkt %lg (inherited from cell_is_optically_thick)\n",
             cell_is_optically_thick_vpkt);
  }

  // Maximum optical depth. If a vpkt reaches dum7 is thrown away
  assert_always(fscanf(input_file, "%lg \n", &tau_max_vpkt) == 1);
  printout("vpkt.txt: tau_max_vpkt %g\n", tau_max_vpkt);

  // Produce velocity grid map if =1
  int in_vgrid_on = 0;
  assert_always(fscanf(input_file, "%d \n", &in_vgrid_on) == 1);
  vgrid_on = in_vgrid_on != 0;
  printout("vpkt.txt: velocity grid map %s\n", (vgrid_on) ? "ENABLED" : "DISABLED");

  if (vgrid_on) {
    double tmin_grid_in_days{NAN};
    double tmax_grid_in_days{NAN};
    // Specify time range for velocity grid map
    assert_always(fscanf(input_file, "%lg %lg \n", &tmin_grid_in_days, &tmax_grid_in_days) == 2);
    tmin_grid = tmin_grid_in_days * DAY;
    tmax_grid = tmax_grid_in_days * DAY;

    printout("vpkt.txt: velocity grid time range tmin_grid %gd tmax_grid %gd\n", tmin_grid / DAY, tmax_grid / DAY);

    // Specify wavelength range: number of intervals (dum9) and limits (dum10,dum11)
    assert_always(fscanf(input_file, "%d ", &Nrange_grid) == 1);

    printout("vpkt.txt: velocity grid frequency intervals %d\n", Nrange_grid);

    nu_grid_max.resize(Nrange_grid, 0.);
    nu_grid_min.resize(Nrange_grid, 0.);
    for (int i = 0; i < Nrange_grid; i++) {
      double range_lambda_min = 0.;
      double range_lambda_max = 0.;
      assert_always(fscanf(input_file, "%lg %lg", &range_lambda_min, &range_lambda_max) == 2);

      nu_grid_max[i] = CLIGHT / (range_lambda_min * 1e-8);
      nu_grid_min[i] = CLIGHT / (range_lambda_max * 1e-8);

      printout("vpkt.txt:   velgrid range %d lambda [%g, %g] Angstroms\n", i, 1e8 * CLIGHT / nu_grid_max[i],
               1e8 * CLIGHT / nu_grid_min[i]);
    }
  }

  fclose(input_file);
}

void vpkt_write_timestep(const int nts, const int my_rank, const bool is_final) {
  if constexpr (!VPKT_ON) {
    return;
  }

  // write specpol of the virtual packets
  char filename_vspecpol[MAXFILENAMELENGTH];

  if (is_final) {
    snprintf(filename_vspecpol, MAXFILENAMELENGTH, "vspecpol_%.4d.out", my_rank);
  } else {
    snprintf(filename_vspecpol, MAXFILENAMELENGTH, "vspecpol_%.4d_ts%d.tmp", my_rank, nts);
  }

  printout("Writing %s\n", filename_vspecpol);
  FILE *vspecpol_file = fopen_required(filename_vspecpol, "w");
  write_vspecpol(vspecpol_file);
  fclose(vspecpol_file);

  if (vgrid_on) {
    char filename_vpktgrid[MAXFILENAMELENGTH];
    if (is_final) {
      snprintf(filename_vpktgrid, MAXFILENAMELENGTH, "vpkt_grid_%.4d.out", my_rank);
    } else {
      snprintf(filename_vpktgrid, MAXFILENAMELENGTH, "vpkt_grid_%.4d_ts%d.tmp", my_rank, nts);
    }

    printout("Writing vpkt grid file %s\n", filename_vpktgrid);
    FILE *vpkt_grid_file = fopen_required(filename_vpktgrid, "w");
    write_vpkt_grid(vpkt_grid_file);
    fclose(vpkt_grid_file);
  }

  if constexpr (VPKT_WRITE_CONTRIBS) {
    vpkt_contrib_file.close();
    char filename_prev[MAXFILENAMELENGTH];
    char filename[MAXFILENAMELENGTH];
    if (is_final) {
      snprintf(filename_prev, MAXFILENAMELENGTH, "vpackets_%.4d_ts%d.tmp", my_rank, nts + 1);
      snprintf(filename, MAXFILENAMELENGTH, "vpackets_%.4d.out", my_rank);
    } else {
      snprintf(filename_prev, MAXFILENAMELENGTH, "vpackets_%.4d_ts%d.tmp", my_rank, nts);
      snprintf(filename, MAXFILENAMELENGTH, "vpackets_%.4d_ts%d.tmp", my_rank, nts + 1);
    }

    std::filesystem::copy_file(filename_prev, filename, std::filesystem::copy_options::overwrite_existing);
    printout("Copying %s to %s\n", filename_prev, filename);

    if (!is_final) {
      vpkt_contrib_file = std::ofstream(filename, std::ios::app);
    }
  }
}

void vpkt_init(const int nts, const int my_rank, const bool continued_from_saved) {
  if constexpr (!VPKT_ON) {
    return;
  }

  init_vspecpol();
  if (vgrid_on) {
    init_vpkt_grid();
  }

  if constexpr (VPKT_WRITE_CONTRIBS) {
    char filename[MAXFILENAMELENGTH];
    snprintf(filename, MAXFILENAMELENGTH, "vpackets_%.4d_ts%d.tmp", my_rank, nts + 1);

    if (continued_from_saved) {
      char filename_prev[MAXFILENAMELENGTH];
      snprintf(filename_prev, MAXFILENAMELENGTH, "vpackets_%.4d_ts%d.tmp", my_rank, nts);
      std::filesystem::copy_file(filename_prev, filename, std::filesystem::copy_options::overwrite_existing);
      printout("Copying %s to %s\n", filename_prev, filename);
    } else {
      // Create new file with header line
      vpkt_contrib_file = std::ofstream(filename, std::ios::trunc);

      vpkt_contrib_file << "#emissiontype trueemissiontype absorption_type absorption_freq";

      for (int obsdirindex = 0; obsdirindex < Nobs; obsdirindex++) {
        vpkt_contrib_file << " dir" << obsdirindex << "_t_arrive_d dir" << obsdirindex << "_nu_rf";
        for (int ind = 0; ind < Nspectra; ind++) {
          vpkt_contrib_file << " dir" << obsdirindex << "_e_rf_" << ind;
        }
      }

      vpkt_contrib_file << "\n";
      vpkt_contrib_file.flush();
      vpkt_contrib_file.close();
    }

    vpkt_contrib_file = std::ofstream(filename, std::ios::app);
  }

  if (continued_from_saved) {
    // Continue simulation: read into temporary files

    read_vspecpol(my_rank, nts);

    if (vgrid_on) {
      read_vpkt_grid(my_rank, nts);
    }
  }
}

auto vpkt_call_estimators(const Packet &pkt, const enum packet_type type_before_rpkt) -> void {
  if constexpr (!VPKT_ON) {
    return;
  }

  // Cut on vpkts
  const int mgi = grid::get_cell_modelgridindex(pkt.where);

  if (grid::modelgrid[mgi].thick != 0) {
    return;
  }

  const double t_current = pkt.prop_time;

  std::stringstream vpkt_contrib_row;

  bool any_dir_escaped = false;
  for (int obsdirindex = 0; obsdirindex < Nobs; obsdirindex++) {
    // loop over different observer directions

    const auto obsdir = std::array<double, 3>{
        sqrt(1 - (nz_obs_vpkt[obsdirindex] * nz_obs_vpkt[obsdirindex])) * cos(phiobs[obsdirindex]),
        sqrt(1 - (nz_obs_vpkt[obsdirindex] * nz_obs_vpkt[obsdirindex])) * sin(phiobs[obsdirindex]),
        nz_obs_vpkt[obsdirindex]};

    const double t_arrive = t_current - (dot(pkt.pos, obsdir) / CLIGHT_PROP);

    bool dir_escaped = false;
    if (t_arrive >= VSPEC_TIMEMIN_input && t_arrive <= VSPEC_TIMEMAX_input) {
      // time selection

      const double doppler = calculate_doppler_nucmf_on_nurf(pkt.pos, obsdir, pkt.prop_time);
      const double nu_rf = pkt.nu_cmf / doppler;
      const double e_rf = pkt.e_cmf / doppler;

      for (int i = 0; i < Nrange; i++) {
        // Loop over frequency ranges

        if ((nu_rf > VSPEC_NUMIN_input[i] && nu_rf < VSPEC_NUMAX_input[i]) ||
            (pkt.absorptionfreq > VSPEC_NUMIN_input[i] && pkt.absorptionfreq < VSPEC_NUMAX_input[i])) {
          // frequency selection
          dir_escaped = rlc_emiss_vpkt(pkt, t_current, t_arrive, nu_rf, e_rf, obsdirindex, obsdir, type_before_rpkt,
                                       vpkt_contrib_row);
          break;  // assume that the frequency ranges do not overlap
        }
      }
    }

    if (dir_escaped) {
      any_dir_escaped = true;
    } else {
      vpkt_contrib_row << " -1. -1.";  // t_arrive_d nu_rf
      for (int ind = 0; ind < Nspectra; ind++) {
        vpkt_contrib_row << " 0.";  // e_rf_diri_j
      }
    }
  }
  if (VPKT_WRITE_CONTRIBS && any_dir_escaped) {
    vpkt_contrib_file << pkt.emissiontype << " " << pkt.trueemissiontype << " " << pkt.absorptiontype << " "
                      << pkt.absorptionfreq;
    vpkt_contrib_file << vpkt_contrib_row.rdbuf() << "\n";
    vpkt_contrib_file.flush();
  }
}

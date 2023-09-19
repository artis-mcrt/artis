#include "vpkt.h"

#include <cmath>
#include <cstring>

#include "atomic.h"
#include "grid.h"
#include "ltepop.h"
#include "rpkt.h"
#include "sn3d.h"
#include "stats.h"
#include "update_grid.h"
#include "vectors.h"

struct stokeparams {
  double i = 0.;
  double q = 0.;
  double u = 0.;
};

struct vspecpol {
  struct stokeparams flux[VMNUBINS];
  float lower_time = NAN;
  float delta_t = NAN;
};

struct vspecpol **vspecpol = nullptr;

float lower_freq_vspec[VMNUBINS];
float delta_freq_vspec[VMNUBINS];

// --------- INPUT PARAMETERS -----------

int Nobs;      // Number of observer directions
int Nspectra;  // Number of virtual packet spectra per observer direction (total + elements switched off)
std::vector<double> nz_obs_vpkt;
std::vector<double> phiobs;
double VSPEC_TIMEMIN_input;
double VSPEC_TIMEMAX_input;
int Nrange;  // Number of wavelength ranges

std::vector<double> VSPEC_NUMIN_input;
std::vector<double> VSPEC_NUMAX_input;
double cell_is_optically_thick_vpkt;
double tau_max_vpkt;

std::vector<int> exclude;  // vector of opacity contribution setups
                           //-1: no line opacity; -2: no bf opacity; -3: no ff opacity; -4: no es opacity,
                           // +ve: exclude element with atomic number's contribution to bound-bound opacity
std::vector<double> tau_vpkt;

// --------- Vstruct packet GRID -----------

struct vgrid {
  std::vector<std::vector<double>> flux;
  double yvel = NAN;
  double zvel = NAN;
};

struct vgrid vgrid_i[VGRID_NY][VGRID_NZ];
struct vgrid vgrid_q[VGRID_NY][VGRID_NZ];
struct vgrid vgrid_u[VGRID_NY][VGRID_NZ];

int Nrange_grid;
double tmin_grid;
double tmax_grid;
std::vector<double> nu_grid_min;
std::vector<double> nu_grid_max;
bool vgrid_on;

double dlogt_vspec = NAN;
double dlognu_vspec = NAN;

// number of virtual packets in a given timestep
int nvpkt;

// number of escaped virtual packet in a given timestep (with tau < tau_max)
int nvpkt_esc1;  // electron scattering event
int nvpkt_esc2;  // kpkt deactivation
int nvpkt_esc3;  // macroatom deactivation

// Virtual packet is killed when tau reaches tau_max_vpkt for ALL the different setups
// E.g. imagine that a packet in the first setup (all elements included) reaches tau = tau_max_vpkt
// because of the element Zi. If we remove Zi, tau now could be lower than tau_max_vpkt and could
// thus contribute to the spectrum.
static auto all_taus_past_taumax(std::vector<double> &tau, const double tau_max) -> bool {
  int count = 0;

  for (int i = 0; i < Nspectra; i++) {
    if (tau[i] > tau_max) {
      count += 1;
    }
  }

  return count == Nspectra;
}

// Routine to add a packet to the outcoming spectrum.
static void add_to_vspecpol(const struct packet &vpkt, const int obsbin, const int ind, const double t_arrive) {
  // Need to decide in which (1) time and (2) frequency bin the vpkt is escaping

  const int ind_comb = Nspectra * obsbin + ind;

  /// Put this into the time grid.
  if (t_arrive > VSPEC_TIMEMIN && t_arrive < VSPEC_TIMEMAX) {
    const int nt = static_cast<int>((log(t_arrive) - log(VSPEC_TIMEMIN)) / dlogt_vspec);
    if (vpkt.nu_rf > VSPEC_NUMIN && vpkt.nu_rf < VSPEC_NUMAX) {
      const int nnu = static_cast<int>((log(vpkt.nu_rf) - log(VSPEC_NUMIN)) / dlognu_vspec);
      const double pktcontrib = vpkt.e_rf / vspecpol[nt][ind_comb].delta_t / delta_freq_vspec[nnu] / 4.e12 / PI /
                                PARSEC / PARSEC / globals::nprocs * 4 * PI;

      safeadd(vspecpol[nt][ind_comb].flux[nnu].i, vpkt.stokes[0] * pktcontrib);
      safeadd(vspecpol[nt][ind_comb].flux[nnu].q, vpkt.stokes[1] * pktcontrib);
      safeadd(vspecpol[nt][ind_comb].flux[nnu].u, vpkt.stokes[2] * pktcontrib);
    }
  }
}

// Routine to add a packet to the outcoming spectrum.
static void add_to_vpkt_grid(const struct packet &vpkt, std::span<const double, 3> vel, const int wlbin,
                             const int obsbin, std::span<const double, 3> obs) {
  double vref1 = NAN;
  double vref2 = NAN;

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
            obs[1] * obs[2] * (1 - obs[0]) / sqrt(1 - obs[0] * obs[0]) * vel[2];
    vref2 = -obs[2] * vel[0] - obs[1] * obs[2] * (1 - obs[0]) / sqrt(1 - obs[0] * obs[0]) * vel[1] +
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
    safeadd(vgrid_i[ny][nz].flux[wlbin][obsbin], vpkt.stokes[0] * vpkt.e_rf);
    safeadd(vgrid_q[ny][nz].flux[wlbin][obsbin], vpkt.stokes[1] * vpkt.e_rf);
    safeadd(vgrid_u[ny][nz].flux[wlbin][obsbin], vpkt.stokes[2] * vpkt.e_rf);
  }
}

static void rlc_emiss_vpkt(const struct packet *const pkt_ptr, const double t_current, const int obsbin,
                           std::span<double, 3> obsdir, const enum packet_type realtype) {
  int snext = 0;
  int mgi = 0;

  struct packet vpkt = *pkt_ptr;

  bool end_packet = false;
  double ldist = 0;
  double t_future = t_current;

  for (int ind = 0; ind < Nspectra; ind++) {
    tau_vpkt[ind] = 0;
  }

  vpkt.dir[0] = obsdir[0];
  vpkt.dir[1] = obsdir[1];
  vpkt.dir[2] = obsdir[2];

  safeincrement(nvpkt);  // increment the number of virtual packet in the given timestep

  double vel_vec[3] = {NAN, NAN, NAN};
  get_velocity(pkt_ptr->pos, vel_vec, t_current);

  // rf frequency and energy
  const double dopplerfactor = doppler_nucmf_on_nurf(vpkt.dir, vel_vec);
  vpkt.nu_rf = vpkt.nu_cmf / dopplerfactor;
  vpkt.e_rf = vpkt.e_cmf / dopplerfactor;

  double Qi = vpkt.stokes[1];
  double Ui = vpkt.stokes[2];

  // ------------ SCATTERING EVENT: dipole function --------------------

  double ref1[3] = {NAN, NAN, NAN};
  double ref2[3] = {NAN, NAN, NAN};
  double pn = NAN;
  double I = NAN;
  double Q = NAN;
  double U = NAN;
  if (realtype == TYPE_RPKT) {
    // Transform Stokes Parameters from the RF to the CMF

    double old_dir_cmf[3] = {NAN, NAN, NAN};
    frame_transform(pkt_ptr->dir, &Qi, &Ui, vel_vec, old_dir_cmf);

    // Need to rotate Stokes Parameters in the scattering plane

    double obs_cmf[3];
    angle_ab(vpkt.dir, vel_vec, obs_cmf);

    meridian(old_dir_cmf, ref1, ref2);

    // This is the i1 angle of Bulla+2015, obtained by computing the angle between the
    // reference axes ref1 and ref2 in the meridian frame and the corresponding axes
    // ref1_sc and ref2_sc in the scattering plane.
    const double i1 = rot_angle(old_dir_cmf, obs_cmf, ref1, ref2);
    const double cos2i1 = cos(2 * i1);
    const double sin2i1 = sin(2 * i1);

    const double Qold = Qi * cos2i1 - Ui * sin2i1;
    const double Uold = Qi * sin2i1 + Ui * cos2i1;

    // Scattering

    const double mu = dot(old_dir_cmf, obs_cmf);

    pn = 3. / (16. * PI) * (1 + pow(mu, 2.) + (pow(mu, 2.) - 1) * Qold);

    const double Inew = 0.75 * ((mu * mu + 1.0) + Qold * (mu * mu - 1.0));
    double Qnew = 0.75 * ((mu * mu - 1.0) + Qold * (mu * mu + 1.0));
    double Unew = 1.5 * mu * Uold;

    Qnew = Qnew / Inew;
    Unew = Unew / Inew;
    I = Inew / Inew;

    // Need to rotate Stokes Parameters out of the scattering plane to the meridian frame

    meridian(obs_cmf, ref1, ref2);

    // This is the i2 angle of Bulla+2015, obtained from the angle THETA between the
    // reference axes ref1_sc and ref2_sc in the scattering plane and ref1 and ref2 in the
    // meridian frame. NB: we need to add PI to transform THETA to i2
    const double i2 = PI + rot_angle(obs_cmf, old_dir_cmf, ref1, ref2);
    const double cos2i2 = cos(2 * i2);
    const double sin2i2 = sin(2 * i2);

    Q = Qnew * cos2i2 + Unew * sin2i2;
    U = -Qnew * sin2i2 + Unew * cos2i2;

    // Transform Stokes Parameters from the CMF to the RF

    const double vel_rev[3] = {-vel_vec[0], -vel_vec[1], -vel_vec[2]};

    frame_transform(obs_cmf, &Q, &U, vel_rev, obsdir);

  } else if (realtype == TYPE_KPKT || realtype == TYPE_MA) {
    // MACROATOM and KPKT: isotropic emission
    I = 1;
    Q = 0;
    U = 0;
    pn = 1 / (4 * PI);
  }

  // compute the optical depth to boundary

  mgi = grid::get_cell_modelgridindex(vpkt.where);
  struct rpkt_cont_opacity kappa_vpkt_cont = {};

  while (!end_packet) {
    // distance to the next cell
    const double sdist = grid::boundary_distance(&vpkt, &snext);
    const double s_cont = sdist * t_current * t_current * t_current / (t_future * t_future * t_future);

    calculate_kappa_rpkt_cont(&vpkt, &kappa_vpkt_cont, false);

    const double kap_cont = kappa_vpkt_cont.total;
    const double kap_cont_nobf = kap_cont - kappa_vpkt_cont.bf;
    const double kap_cont_noff = kap_cont - kappa_vpkt_cont.ff;
    const double kap_cont_noes = kap_cont - kappa_vpkt_cont.es;

    for (int ind = 0; ind < Nspectra; ind++) {
      if (exclude[ind] == -2) {
        tau_vpkt[ind] += kap_cont_nobf * s_cont;
      } else if (exclude[ind] == -3) {
        tau_vpkt[ind] += kap_cont_noff * s_cont;
      } else if (exclude[ind] == -4) {
        tau_vpkt[ind] += kap_cont_noes * s_cont;
      } else {
        tau_vpkt[ind] += kap_cont * s_cont;
      }
    }

    // kill vpkt with high optical depth
    if (all_taus_past_taumax(tau_vpkt, tau_max_vpkt)) {
      return;
    }

    struct packet dummypkt_abort = vpkt;
    move_pkt_withtime(&dummypkt_abort, sdist);
    const double nu_cmf_abort = dummypkt_abort.nu_cmf;
    assert_testmodeonly(nu_cmf_abort <= vpkt.nu_cmf);
    const double d_nu_on_d_l = (nu_cmf_abort - vpkt.nu_cmf) / sdist;

    ldist = 0;
    while (ldist < sdist) {
      const int lineindex = closest_transition(vpkt.nu_cmf, vpkt.next_trans);

      if (lineindex < 0) {
        vpkt.next_trans = globals::nlines + 1;
      } else {
        const double nutrans = globals::linelist[lineindex].nu;

        vpkt.next_trans = lineindex + 1;

        // get distance to the line interaction
        if (vpkt.nu_cmf < nutrans) {
          ldist = 0;
        } else if constexpr (!USE_RELATIVISTIC_DOPPLER_SHIFT) {
          ldist = CLIGHT * t_current * (vpkt.nu_cmf / nutrans - 1);
        } else {
          // use linear interpolation of frequency along the path
          ldist = (nutrans - vpkt.nu_cmf) / d_nu_on_d_l;
        }

        if (ldist < 0.) {
          printout("[warning] get_event: ldist < 0 %g\n", ldist);
        }

        if (ldist > sdist) {
          // exit the while loop if you reach the boundary; go back to the previous transition to start next cell with
          // the excluded line

          vpkt.next_trans -= 1;
          // printout("ldist > sdist : line in the next cell\n");
          break;
        }

        const double t_line = t_current + ldist / CLIGHT;

        const int element = globals::linelist[lineindex].elementindex;
        const int ion = globals::linelist[lineindex].ionindex;
        const int upper = globals::linelist[lineindex].upperlevelindex;
        const int lower = globals::linelist[lineindex].lowerlevelindex;
        const auto A_ul = globals::linelist[lineindex].einstein_A;

        const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nutrans, 3) * A_ul;
        const double B_lu = stat_weight(element, ion, upper) / stat_weight(element, ion, lower) * B_ul;

        const auto n_u = calculate_levelpop(mgi, element, ion, upper);
        const auto n_l = calculate_levelpop(mgi, element, ion, lower);
        const double tau_line = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_line;

        // Check on the element to exclude
        // NB: ldist before need to be computed anyway (I want to move the packets to the
        // line interaction point even if I don't interact)
        const int anumber = get_atomicnumber(element);
        for (int ind = 0; ind < Nspectra; ind++) {
          // If exclude[ind]==-1, I do not include line opacity
          if (exclude[ind] != -1 && (anumber != exclude[ind])) {
            tau_vpkt[ind] += tau_line;
          }
        }

        // kill vpkt with high optical depth
        if (all_taus_past_taumax(tau_vpkt, tau_max_vpkt)) {
          return;
        }
      }
    }

    // virtual packet is still at the starting position
    // move it to cell boundary and go to next cell
    // printf("I'm changing cell. I'm going from nu_cmf = %.e ",dummy_ptr->nu_cmf);

    t_future += (sdist / CLIGHT_PROP);
    vpkt.prop_time = t_future;
    move_pkt(&vpkt, sdist);

    // printout("About to change vpkt cell\n");

    grid::change_cell(&vpkt, snext);
    end_packet = (vpkt.type == TYPE_ESCAPE);
    // printout("Completed change vpkt cell\n");

    // printout("dummy->nu_cmf = %g \n",dummy_ptr->nu_cmf);
    mgi = grid::get_cell_modelgridindex(vpkt.where);
    // break if you reach an empty cell
    if (mgi == grid::get_npts_model()) {
      break;
    }

    // kill vpkt with pass through a thick cell
    if (grid::modelgrid[mgi].thick == 1) {
      return;
    }
  }

  // increment the number of escaped virtual packet in the given timestep
  if (realtype == TYPE_RPKT) {
    safeincrement(nvpkt_esc1);
  } else if (realtype == TYPE_KPKT) {
    safeincrement(nvpkt_esc2);
  } else if (realtype == TYPE_MA) {
    safeincrement(nvpkt_esc3);
  }

  const double t_arrive = t_current - (dot(pkt_ptr->pos, vpkt.dir) / CLIGHT_PROP);
  // -------------- final stokes vector ---------------

  for (int ind = 0; ind < Nspectra; ind++) {
    // printout("obsbin %d spectrum %d tau_vpkt %g\n", obsbin, ind, tau_vpkt[ind]);
    const double prob = pn * exp(-tau_vpkt[ind]);

    assert_always(std::isfinite(prob));

    vpkt.stokes[0] = I * prob;
    vpkt.stokes[1] = Q * prob;
    vpkt.stokes[2] = U * prob;

    for (const auto stokeval : vpkt.stokes) {
      assert_always(std::isfinite(stokeval));
    }

    // bin on fly and produce file with spectrum

    add_to_vspecpol(vpkt, obsbin, ind, t_arrive);
  }

  // vpkt grid

  if (vgrid_on) {
    const double prob = pn * exp(-tau_vpkt[0]);

    vpkt.stokes[0] = I * prob;
    vpkt.stokes[1] = Q * prob;
    vpkt.stokes[2] = U * prob;

    for (int wlbin = 0; wlbin < Nrange_grid; wlbin++) {
      if (vpkt.nu_rf > nu_grid_min[wlbin] && vpkt.nu_rf < nu_grid_max[wlbin]) {  // Frequency selection
        if (t_arrive > tmin_grid && t_arrive < tmax_grid) {                      // Time selection
          add_to_vpkt_grid(vpkt, vel_vec, wlbin, obsbin, obsdir);
        }
      }
    }
  }
}

void init_vspecpol() {
  vspecpol = static_cast<struct vspecpol **>(malloc(VMTBINS * sizeof(struct vspecpol *)));

  const int indexmax = Nspectra * Nobs;
  for (int p = 0; p < VMTBINS; p++) {
    vspecpol[p] = static_cast<struct vspecpol *>(malloc(indexmax * sizeof(struct vspecpol)));
  }

  dlogt_vspec = (log(VSPEC_TIMEMAX) - log(VSPEC_TIMEMIN)) / VMTBINS;
  dlognu_vspec = (log(VSPEC_NUMAX) - log(VSPEC_NUMIN)) / VMNUBINS;

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
        flux.i = 0.0;
        flux.q = 0.0;
        flux.u = 0.0;
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

void read_vspecpol(int my_rank, int nts) {
  char filename[MAXFILENAMELENGTH];

  if (nts % 2 == 0) {
    snprintf(filename, MAXFILENAMELENGTH, "vspecpol_%d_%d_odd.tmp", 0, my_rank);
  } else {
    snprintf(filename, MAXFILENAMELENGTH, "vspecpol_%d_%d_even.tmp", 0, my_rank);
  }

  FILE *vspecpol_file = fopen_required(filename, "rb");

  float a = NAN;
  float b = NAN;
  float c = NAN;

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
      const double yvel = globals::vmax - (n + 0.5) * ybin;
      const double zvel = globals::vmax - (m + 0.5) * zbin;

      vgrid_i[n][m].yvel = yvel;
      vgrid_i[n][m].zvel = zvel;

      vgrid_q[n][m].yvel = yvel;
      vgrid_q[n][m].zvel = zvel;

      vgrid_u[n][m].yvel = yvel;
      vgrid_u[n][m].zvel = zvel;

      vgrid_i[n][m].flux.resize(Nrange_grid, {});
      vgrid_q[n][m].flux.resize(Nrange_grid, {});
      vgrid_u[n][m].flux.resize(Nrange_grid, {});
      for (int wlbin = 0; wlbin < Nrange_grid; wlbin++) {
        vgrid_i[n][m].flux[wlbin] = std::vector<double>(Nobs, 0.);
        vgrid_q[n][m].flux[wlbin] = std::vector<double>(Nobs, 0.);
        vgrid_u[n][m].flux[wlbin] = std::vector<double>(Nobs, 0.);
      }
    }
  }
}

void write_vpkt_grid(FILE *vpkt_grid_file) {
  for (int obsbin = 0; obsbin < Nobs; obsbin++) {
    for (int wlbin = 0; wlbin < Nrange_grid; wlbin++) {
      for (int n = 0; n < VGRID_NY; n++) {
        for (int m = 0; m < VGRID_NZ; m++) {
          fprintf(vpkt_grid_file, "%g ", vgrid_i[n][m].yvel);
          fprintf(vpkt_grid_file, "%g ", vgrid_i[n][m].zvel);

          fprintf(vpkt_grid_file, "%g ", vgrid_i[n][m].flux[wlbin][obsbin]);
          fprintf(vpkt_grid_file, "%g ", vgrid_q[n][m].flux[wlbin][obsbin]);
          fprintf(vpkt_grid_file, "%g ", vgrid_u[n][m].flux[wlbin][obsbin]);

          fprintf(vpkt_grid_file, "\n");
        }
      }
    }
  }
}

void read_vpkt_grid(FILE *vpkt_grid_file) {
  for (int obsbin = 0; obsbin < Nobs; obsbin++) {
    for (int wlbin = 0; wlbin < Nrange_grid; wlbin++) {
      for (int n = 0; n < VGRID_NY; n++) {
        for (int m = 0; m < VGRID_NZ; m++) {
          assert_always(fscanf(vpkt_grid_file, "%lg ", &vgrid_i[n][m].yvel) == 1);
          assert_always(fscanf(vpkt_grid_file, "%lg ", &vgrid_i[n][m].zvel) == 1);

          assert_always(fscanf(vpkt_grid_file, "%lg ", &vgrid_i[n][m].flux[wlbin][obsbin]) == 1);
          assert_always(fscanf(vpkt_grid_file, "%lg ", &vgrid_q[n][m].flux[wlbin][obsbin]) == 1);
          assert_always(fscanf(vpkt_grid_file, "%lg ", &vgrid_u[n][m].flux[wlbin][obsbin]) == 1);

          assert_always(fscanf(vpkt_grid_file, "\n") == 0);
        }
      }
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
      abort();
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
    double tmin_grid_in_days = NAN;
    double tmax_grid_in_days = NAN;
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

auto vpkt_call_estimators(struct packet *pkt_ptr, const enum packet_type realtype) -> void {
  // Cut on vpkts
  const int mgi = grid::get_cell_modelgridindex(pkt_ptr->where);

  if (grid::modelgrid[mgi].thick != 0) {
    return;
  }

  const double t_current = pkt_ptr->prop_time;

  double vel_vec[3];
  get_velocity(pkt_ptr->pos, vel_vec, pkt_ptr->prop_time);

  // this is just to find the next_trans value when is set to 0 (avoid doing that in the vpkt routine for each observer)
  if (pkt_ptr->next_trans == 0) {
    const int lineindex = closest_transition(pkt_ptr->nu_cmf, pkt_ptr->next_trans);  /// returns negative
    if (lineindex < 0) {
      pkt_ptr->next_trans = lineindex + 1;
    }
  }

  for (int obsbin = 0; obsbin < Nobs; obsbin++) {
    // loop over different observer directions

    double obsdir[3] = {sqrt(1 - nz_obs_vpkt[obsbin] * nz_obs_vpkt[obsbin]) * cos(phiobs[obsbin]),
                        sqrt(1 - nz_obs_vpkt[obsbin] * nz_obs_vpkt[obsbin]) * sin(phiobs[obsbin]), nz_obs_vpkt[obsbin]};

    const double t_arrive = t_current - (dot(pkt_ptr->pos, obsdir) / CLIGHT_PROP);

    if (t_arrive >= VSPEC_TIMEMIN_input && t_arrive <= VSPEC_TIMEMAX_input) {
      // time selection

      const double nu_rf = pkt_ptr->nu_cmf / doppler_nucmf_on_nurf(obsdir, vel_vec);

      for (int i = 0; i < Nrange; i++) {
        // Loop over frequency ranges

        if (nu_rf > VSPEC_NUMIN_input[i] && nu_rf < VSPEC_NUMAX_input[i]) {
          // frequency selection

          rlc_emiss_vpkt(pkt_ptr, t_current, obsbin, obsdir, realtype);
        }
      }
    }
  }
}

auto rot_angle(std::span<double, 3> n1, std::span<double, 3> n2, std::span<double, 3> ref1, std::span<double, 3> ref2)
    -> double {
  // Rotation angle from the scattering plane
  // We need to rotate Stokes Parameters to (or from) the scattering plane from (or to)
  // the meridian frame such that Q=1 is in the scattering plane and along ref1

  // ref1_sc is the ref1 axis in the scattering plane ref1 = n1 x ( n1 x n2 )
  const double n1_dot_n2 = dot(n1, n2);
  double ref1_sc[3] = {n1[0] * n1_dot_n2 - n2[0], n1[1] * n1_dot_n2 - n2[1], n1[2] * n1_dot_n2 - n2[2]};
  vec_norm(ref1_sc, ref1_sc);

  double cos_stokes_rot_1 = dot(ref1_sc, ref1);
  const double cos_stokes_rot_2 = dot(ref1_sc, ref2);

  if (cos_stokes_rot_1 < -1) {
    cos_stokes_rot_1 = -1;
  }
  if (cos_stokes_rot_1 > 1) {
    cos_stokes_rot_1 = 1;
  }

  double i = 0;
  if ((cos_stokes_rot_1 > 0) && (cos_stokes_rot_2 > 0)) {
    i = acos(cos_stokes_rot_1);
  } else if ((cos_stokes_rot_1 < 0) && (cos_stokes_rot_2 > 0)) {
    i = M_PI - acos(fabs(cos_stokes_rot_1));
  } else if ((cos_stokes_rot_1 > 0) && (cos_stokes_rot_2 < 0)) {
    i = 2 * M_PI - acos(cos_stokes_rot_1);
  } else if ((cos_stokes_rot_1 < 0) && (cos_stokes_rot_2 < 0)) {
    i = M_PI + acos(fabs(cos_stokes_rot_1));
  }
  if (cos_stokes_rot_1 == 0) {
    i = M_PI / 2.;
  }
  if (cos_stokes_rot_2 == 0) {
    i = 0.0;
  }

  if (!std::isfinite(i)) {
    printout("Warning NaN: %3.6f \t %3.6f \t %3.6f \n", cos_stokes_rot_1, cos_stokes_rot_2, acos(cos_stokes_rot_1));
  }

  return i;
}

// Routine to compute the meridian frame axes ref1 and ref2
void meridian(std::span<const double, 3> n, std::span<double, 3> ref1, std::span<double, 3> ref2) {
  // for ref_1 use (from triple product rule)
  const double n_xylen = sqrt(n[0] * n[0] + n[1] * n[1]);
  ref1[0] = -1. * n[0] * n[2] / n_xylen;
  ref1[1] = -1. * n[1] * n[2] / n_xylen;
  ref1[2] = (1 - (n[2] * n[2])) / n_xylen;

  // for ref_2 use vector product of n_cmf with ref1
  cross_prod(ref1, n, ref2);
}

// Routine to transform the Stokes Parameters from RF to CMF
void frame_transform(std::span<const double, 3> n_rf, double *Q, double *U, std::span<const double, 3> v,
                     std::span<double, 3> n_cmf) {
  double ref1[3] = {NAN, NAN, NAN};
  double ref2[3] = {NAN, NAN, NAN};
  // Meridian frame in the RF
  meridian(n_rf, ref1, ref2);

  const double Q0 = *Q;
  const double U0 = *U;

  // Compute polarisation (which is invariant)
  const double p = sqrt(Q0 * Q0 + U0 * U0);

  // We want to compute the angle between ref1 and the electric field
  double rot_angle = 0;

  if (p > 0) {
    const double cos2rot_angle = Q0 / p;
    const double sin2rot_angle = U0 / p;

    if ((cos2rot_angle > 0) && (sin2rot_angle > 0)) {
      rot_angle = acos(Q0 / p) / 2.;
    } else if ((cos2rot_angle < 0) && (sin2rot_angle > 0)) {
      rot_angle = (M_PI - acos(fabs(cos2rot_angle))) / 2.;
    } else if ((cos2rot_angle < 0) && (sin2rot_angle < 0)) {
      rot_angle = (M_PI + acos(fabs(cos2rot_angle))) / 2.;
    } else if ((cos2rot_angle > 0) && (sin2rot_angle < 0)) {
      rot_angle = (2. * M_PI - acos(fabs(cos2rot_angle))) / 2.;
    } else if (cos2rot_angle == 0) {
      rot_angle = 0.25 * M_PI;
      if (U0 < 0) {
        rot_angle = 0.75 * M_PI;
      }
    }
    if (sin2rot_angle == 0) {
      rot_angle = 0.0;
      if (Q0 < 0) {
        rot_angle = 0.5 * M_PI;
      }
    }
  }

  // Define electric field by linear combination of ref1 and ref2 (using the angle just computed)

  const double elec_rf[3] = {cos(rot_angle) * ref1[0] - sin(rot_angle) * ref2[0],
                             cos(rot_angle) * ref1[1] - sin(rot_angle) * ref2[1],
                             cos(rot_angle) * ref1[2] - sin(rot_angle) * ref2[2]};

  // Aberration
  angle_ab(n_rf, v, n_cmf);

  double elec_cmf[3] = {NAN, NAN, NAN};
  // Lorentz transformation of E
  lorentz(elec_rf, n_rf, v, elec_cmf);

  // Meridian frame in the CMF
  meridian(n_cmf, ref1, ref2);

  // Projection of E onto ref1 and ref2
  const double cosine_elec_ref1 = dot(elec_cmf, ref1);
  const double cosine_elec_ref2 = dot(elec_cmf, ref2);

  // Compute the angle between ref1 and the electric field
  double theta_rot = 0.;
  if ((cosine_elec_ref1 > 0) && (cosine_elec_ref2 < 0)) {
    theta_rot = acos(cosine_elec_ref1);
  } else if ((cosine_elec_ref1 < 0) && (cosine_elec_ref2 > 0)) {
    theta_rot = M_PI + acos(fabs(cosine_elec_ref1));
  } else if ((cosine_elec_ref1 < 0) && (cosine_elec_ref2 < 0)) {
    theta_rot = M_PI - acos(fabs(cosine_elec_ref1));
  } else if ((cosine_elec_ref1 > 0) && (cosine_elec_ref2 > 0)) {
    theta_rot = 2 * M_PI - acos(cosine_elec_ref1);
  }
  if (cosine_elec_ref1 == 0) {
    theta_rot = M_PI / 2.;
  }
  if (cosine_elec_ref2 == 0) {
    theta_rot = 0.0;
  }
  if (cosine_elec_ref1 > 1) {
    theta_rot = 0.0;
  }
  if (cosine_elec_ref1 < -1) {
    theta_rot = M_PI;
  }

  // Compute Stokes Parameters in the CMF
  *Q = cos(2 * theta_rot) * p;
  *U = sin(2 * theta_rot) * p;
}

void lorentz(std::span<const double, 3> e_rf, std::span<const double, 3> n_rf, std::span<const double, 3> v,
             std::span<double, 3> e_cmf) {
  // Lorentz transformations from RF to CMF

  const double beta[3] = {v[0] / CLIGHT, v[1] / CLIGHT, v[2] / CLIGHT};
  double const vsqr = dot(beta, beta);

  const double gamma_rel = 1. / (sqrt(1 - vsqr));

  const double e_par[3] = {dot(e_rf, beta) * beta[0] / (vsqr), dot(e_rf, beta) * beta[1] / (vsqr),
                           dot(e_rf, beta) * beta[2] / (vsqr)};

  const double e_perp[3] = {e_rf[0] - e_par[0], e_rf[1] - e_par[1], e_rf[2] - e_par[2]};

  double b_rf[3] = {NAN, NAN, NAN};
  cross_prod(n_rf, e_rf, b_rf);

  // const double b_par[3] = {dot(b_rf, beta) * beta[0] / (vsqr), dot(b_rf, beta) * beta[1] / (vsqr),
  //                          dot(b_rf, beta) * beta[2] / (vsqr)};

  // const double b_perp[3] = {b_rf[0] - b_par[0], b_rf[1] - b_par[1], b_rf[2] - b_par[2]};

  double v_cr_b[3] = {NAN, NAN, NAN};
  cross_prod(beta, b_rf, v_cr_b);

  // const double v_cr_e[3] = {beta[1] * e_rf[2] - beta[2] * e_rf[1], beta[2] * e_rf[0] - beta[0] * e_rf[2],
  //                           beta[0] * e_rf[1] - beta[1] * e_rf[0]};

  e_cmf[0] = e_par[0] + gamma_rel * (e_perp[0] + v_cr_b[0]);
  e_cmf[1] = e_par[1] + gamma_rel * (e_perp[1] + v_cr_b[1]);
  e_cmf[2] = e_par[2] + gamma_rel * (e_perp[2] + v_cr_b[2]);
  vec_norm(e_cmf, e_cmf);

  // double b_cmf[3];
  // b_cmf[0] = b_par[0] + gamma_rel * (b_perp[0] - v_cr_e[0]);
  // b_cmf[1] = b_par[1] + gamma_rel * (b_perp[1] - v_cr_e[1]);
  // b_cmf[2] = b_par[2] + gamma_rel * (b_perp[2] - v_cr_e[2]);
  // vec_norm(b_cmf, b_cmf);
}
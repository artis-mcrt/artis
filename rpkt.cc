#include "rpkt.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <span>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "ltepop.h"
#include "packet.h"
#include "radfield.h"
#include "sn3d.h"
#include "stats.h"
#include "vectors.h"
#include "vpkt.h"

// kappa in cm^2/g for each bin of each non-empty cell
static float *expansionopacities = nullptr;

// kappa times Planck function for each bin of each non-empty cell
static float *expansionopacity_planck_cumulative = nullptr;
#ifdef MPI_ON
MPI_Win win_expansionopacities = MPI_WIN_NULL;
MPI_Win win_expansionopacity_planck_cumulative = MPI_WIN_NULL;
#endif

void allocate_expansionopacities() {
  if constexpr (!EXPANSIONOPACITIES_ON) {
    return;
  }

  const auto npts_nonempty = grid::get_nonempty_npts_model();
#ifdef MPI_ON
  int my_rank_nonemptycells = npts_nonempty / globals::node_nprocs;
  // rank_in_node 0 gets any remainder
  if (globals::rank_in_node == 0) {
    my_rank_nonemptycells += npts_nonempty - (my_rank_nonemptycells * globals::node_nprocs);
  }
  MPI_Aint size = my_rank_nonemptycells * expopac_nbins * static_cast<MPI_Aint>(sizeof(float));
  int disp_unit = sizeof(float);
  assert_always(MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, globals::mpi_comm_node, &expansionopacities,
                                        &win_expansionopacities) == MPI_SUCCESS);
  assert_always(MPI_Win_shared_query(win_expansionopacities, 0, &size, &disp_unit, &expansionopacities) == MPI_SUCCESS);

  if constexpr (EXPANSION_OPAC_SAMPLE_KAPPAPLANCK) {
    assert_always(MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, globals::mpi_comm_node,
                                          &expansionopacity_planck_cumulative,
                                          &win_expansionopacity_planck_cumulative) == MPI_SUCCESS);
    assert_always(MPI_Win_shared_query(win_expansionopacities, 0, &size, &disp_unit,
                                       &expansionopacity_planck_cumulative) == MPI_SUCCESS);
  }

#else
  expansionopacities = static_cast<float *>(malloc(npts_nonempty * expopac_nbins * sizeof(float)));
  if constexpr (EXPANSION_OPAC_SAMPLE_KAPPAPLANCK) {
    expansionopacity_planck_cumulative = static_cast<float *>(malloc(npts_nonempty * expopac_nbins * sizeof(float)));
  }
#endif
}

auto closest_transition(const double nu_cmf, const int next_trans) -> int
/// for the propagation through non empty cells
// find the next transition lineindex redder than nu_cmf
// return -1 if no transition can be reached
{
  if (next_trans > (globals::nlines - 1)) {
    // packet is tagged as having no more line interactions
    return -1;
  }
  /// if nu_cmf is smaller than the lowest frequency in the linelist,
  /// no line interaction is possible: return negative value as a flag
  if (nu_cmf < globals::linelist[globals::nlines - 1].nu) {
    return -1;
  }

  if (next_trans > 0) {
    /// if left = pkt_ptr->next_trans > 0 we know the next line we should interact with, independent of the packets
    /// current nu_cmf which might be smaller than globals::linelist[left].nu due to propagation errors
    return next_trans;
  }
  // will find the highest frequency (lowest index) line with nu_line <= nu_cmf
  // lower_bound matches the first element where the comparison function is false
  const auto *matchline =
      std::lower_bound(&globals::linelist[next_trans], &globals::linelist[globals::nlines], nu_cmf,
                       [](const auto &line, const double nu_cmf) -> bool { return line.nu > nu_cmf; });
  const int matchindex = std::distance(globals::linelist, matchline);
  if (matchindex >= globals::nlines) {
    return -1;
  }

  return matchindex;
}

static auto get_nu_cmf_abort(struct packet *pkt_ptr, const double abort_dist) -> double {
  // get the frequency change per distance travelled assuming linear change to the abort distance
  struct packet dummypkt_abort = *pkt_ptr;
  // this is done is two parts to get identical results to do_rpkt_step()
  move_pkt_withtime(&dummypkt_abort, abort_dist / 2.);
  move_pkt_withtime(&dummypkt_abort, abort_dist / 2.);
  const double nu_cmf_abort = dummypkt_abort.nu_cmf;
  assert_testmodeonly(nu_cmf_abort <= pkt_ptr->nu_cmf);
  // for USE_RELATIVISTIC_DOPPLER_SHIFT, we will use a linear approximation for
  // the frequency change from start to abort (cell boundary/timestep end)
  return nu_cmf_abort;
}

// wavelength bins are ordered by ascending wavelength (descending frequency)

static constexpr auto get_expopac_bin_nu_upper(const size_t binindex) -> double {
  const auto lambda_lower = expopac_lambdamin + binindex * expopac_deltalambda;
  return 1e8 * CLIGHT / lambda_lower;
}

static constexpr auto get_expopac_bin_nu_lower(const size_t binindex) -> double {
  const auto lambda_upper = expopac_lambdamin + (binindex + 1) * expopac_deltalambda;
  return 1e8 * CLIGHT / lambda_upper;
}

static auto get_event_expansion_opacity(const int modelgridindex, const int nonemptymgi, struct packet *pkt_ptr,
                                        struct rpkt_continuum_absorptioncoeffs &chi_rpkt_cont, const double tau_rnd,
                                        const double abort_dist) -> std::tuple<double, bool> {
  // calculate_chi_rpkt_cont(pkt_ptr->nu_cmf, chi_rpkt_cont, modelgridindex, true);
  const auto doppler = doppler_packet_nucmf_on_nurf(pkt_ptr->pos, pkt_ptr->dir, pkt_ptr->prop_time);

  const auto nu_cmf_abort = get_nu_cmf_abort(pkt_ptr, abort_dist);
  const auto d_nu_on_d_l = (nu_cmf_abort - pkt_ptr->nu_cmf) / abort_dist;

  struct packet dummypkt = *pkt_ptr;
  assert_always(globals::cellcache[tid].cellnumber == modelgridindex);
  double dist = 0.;
  double tau = 0.;
  auto binindex_start =
      static_cast<size_t>(((1e8 * CLIGHT / dummypkt.nu_cmf) - expopac_lambdamin) / expopac_deltalambda);

  for (size_t binindex = binindex_start; binindex < expopac_nbins; binindex++) {
    const auto next_bin_edge_nu = (binindex < 0) ? get_expopac_bin_nu_upper(0) : get_expopac_bin_nu_lower(binindex);
    const auto binedgedist = get_linedistance(dummypkt.prop_time, dummypkt.nu_cmf, next_bin_edge_nu, d_nu_on_d_l);

    // const double chi_cont = chi_rpkt_cont.total * doppler;
    const auto chi_cont = 0.;
    double chi_bb_expansionopac = 0.;
    if (binindex >= 0) {
      const auto kappa = expansionopacities[nonemptymgi * expopac_nbins + binindex];
      // const auto doppler = doppler_packet_nucmf_on_nurf(dummypkt.pos, dummypkt.dir, dummypkt.prop_time);
      // const auto doppler = (pkt_ptr->nu_cmf + d_nu_on_d_l * dist) / pkt_ptr->nu_rf;
      chi_bb_expansionopac = kappa * grid::get_rho(modelgridindex) * doppler;
    }
    double chi_tot = chi_cont + chi_bb_expansionopac;

    if (chi_tot * binedgedist > tau_rnd - tau) {
      // event occurs
      const auto edist = std::max(dist + (tau_rnd - tau) / chi_tot, 0.);
      const bool event_is_boundbound = rng_uniform() <= chi_bb_expansionopac / chi_tot;
      return {edist, event_is_boundbound};
    }

    tau += chi_tot * binedgedist;
    dist += binedgedist;

    if constexpr (!USE_RELATIVISTIC_DOPPLER_SHIFT) {
      move_pkt_withtime(&dummypkt, binedgedist);
    } else {
      // avoid move_pkt_withtime() to skip the standard Doppler shift calculation
      // and use the linear approx instead
      dummypkt.pos[0] += (dummypkt.dir[0] * binedgedist);
      dummypkt.pos[1] += (dummypkt.dir[1] * binedgedist);
      dummypkt.pos[2] += (dummypkt.dir[2] * binedgedist);
      dummypkt.prop_time += binedgedist / CLIGHT_PROP;
      dummypkt.nu_cmf = pkt_ptr->nu_cmf + d_nu_on_d_l * dist;  // should equal nu_trans;
      assert_testmodeonly(dummypkt.nu_cmf <= pkt_ptr->nu_cmf);
    }

    if (dummypkt.nu_cmf <= nu_cmf_abort) {
      // hit edge of cell or timestep limit
      return {std::numeric_limits<double>::max(), false};
    }
  }
  // no more bins, so no opacity below this frequency
  return {std::numeric_limits<double>::max(), false};
}

static auto get_event(const int modelgridindex,
                      struct packet *pkt_ptr,  // pointer to packet object
                      struct rpkt_continuum_absorptioncoeffs &chi_rpkt_cont,
                      const double tau_rnd,    // random optical depth until which the packet travels
                      const double abort_dist  // maximal travel distance before packet leaves cell or time step ends
                      ) -> std::tuple<double, bool>
// returns edist, the distance to the next physical event (continuum or bound-bound) and is_boundbound_event, a
// boolean BE AWARE THAT THIS PROCEDURE SHOULD BE ONLY CALLED FOR NON EMPTY CELLS!!
{
  assert_always(!EXPANSIONOPACITIES_ON);
  assert_testmodeonly(grid::modelgrid[modelgridindex].thick != 1);
  // printout("get_event()\n");
  /// initialize loop variables

  const auto nu_cmf_abort = get_nu_cmf_abort(pkt_ptr, abort_dist);
  const auto d_nu_on_d_l = (nu_cmf_abort - pkt_ptr->nu_cmf) / abort_dist;

  struct packet dummypkt = *pkt_ptr;

  calculate_chi_rpkt_cont(pkt_ptr->nu_cmf, chi_rpkt_cont, modelgridindex, true);
  const double chi_cont =
      chi_rpkt_cont.total * doppler_packet_nucmf_on_nurf(pkt_ptr->pos, pkt_ptr->dir, pkt_ptr->prop_time);
  double tau = 0.;   // optical depth along path
  double dist = 0.;  // position on path
  while (true) {
    /// calculate distance to next line encounter ldist
    /// first select the closest transition in frequency
    /// we need its frequency nu_trans, the element/ion and the corresponding levels
    /// create therefore new variables in packet, which contain next_lowerlevel, ...
    const int lineindex = closest_transition(dummypkt.nu_cmf,
                                             dummypkt.next_trans);  /// returns negative value if nu_cmf > nu_trans
    if (lineindex >= 0) {
      /// line interaction is possible (nu_cmf > nu_trans)

      const double nu_trans = globals::linelist[lineindex].nu;

      // helper variable to overcome numerical problems after line scattering
      // further scattering events should be located at lower frequencies to prevent
      // multiple scattering events of one packet in a single line
      dummypkt.next_trans = lineindex + 1;

      const double ldist = get_linedistance(dummypkt.prop_time, dummypkt.nu_cmf, nu_trans, d_nu_on_d_l);

      const double tau_cont = chi_cont * ldist;

      if (tau_rnd - tau > tau_cont) {
        // got past the continuum optical depth so propagate to the line, and check interaction

        if (nu_trans < nu_cmf_abort) {
          // back up one line, because we didn't reach it before the boundary/timelimit
          pkt_ptr->next_trans = dummypkt.next_trans - 1;

          return {std::numeric_limits<double>::max(), false};
        }

        const int element = globals::linelist[lineindex].elementindex;
        const int ion = globals::linelist[lineindex].ionindex;
        const int upper = globals::linelist[lineindex].upperlevelindex;
        const int lower = globals::linelist[lineindex].lowerlevelindex;
        const double A_ul = einstein_spontaneous_emission(lineindex);
        const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
        const double B_lu = stat_weight(element, ion, upper) / stat_weight(element, ion, lower) * B_ul;

        const double n_u = get_levelpop(modelgridindex, element, ion, upper);
        const double n_l = get_levelpop(modelgridindex, element, ion, lower);

        const double tau_line = std::max(0., (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * dummypkt.prop_time);

        // printout("[debug] get_event:     tau_line %g\n", tau_line);
        // printout("[debug] get_event:       tau_rnd - tau > tau_cont\n");

        if (tau_rnd - tau > tau_cont + tau_line) {
          // total optical depth still below tau_rnd: propagate to the line and continue

          // printout(
          //     "[debug] get_event: tau_rnd - tau > tau_cont + tau_line ... proceed this packets "
          //     "propagation\n");

          dist += ldist;
          tau += tau_cont + tau_line;

          if constexpr (!USE_RELATIVISTIC_DOPPLER_SHIFT) {
            move_pkt_withtime(&dummypkt, ldist);
          } else {
            // avoid move_pkt_withtime() to skip the standard Doppler shift calculation
            // and use the linear approx instead
            dummypkt.pos[0] += (dummypkt.dir[0] * ldist);
            dummypkt.pos[1] += (dummypkt.dir[1] * ldist);
            dummypkt.pos[2] += (dummypkt.dir[2] * ldist);
            dummypkt.prop_time += ldist / CLIGHT_PROP;
            dummypkt.nu_cmf = pkt_ptr->nu_cmf + d_nu_on_d_l * dist;  // should equal nu_trans;
            assert_testmodeonly(dummypkt.nu_cmf <= pkt_ptr->nu_cmf);
          }

          radfield::update_lineestimator(modelgridindex, lineindex,
                                         dummypkt.prop_time * CLIGHT * dummypkt.e_cmf / dummypkt.nu_cmf);

        } else {
          /// bound-bound process occurs
          // printout("[debug] get_event: tau_rnd - tau <= tau_cont + tau_line: bb-process occurs\n");

          pkt_ptr->mastate.element = element;
          pkt_ptr->mastate.ion = ion;
          pkt_ptr->mastate.level = upper;  /// if the MA will be activated it must be in the transitions upper level
          pkt_ptr->mastate.activatingline = lineindex;

          if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
            move_pkt_withtime(&dummypkt, ldist);
            radfield::update_lineestimator(modelgridindex, lineindex,
                                           dummypkt.prop_time * CLIGHT * dummypkt.e_cmf / dummypkt.nu_cmf);
          }

          /// the line and its parameters were already selected by closest_transition!
          // printout("[debug] get_event:         edist %g, abort_dist %g, edist-abort_dist %g, endloop
          // %d\n",edist,abort_dist,edist-abort_dist,endloop);

          pkt_ptr->next_trans = dummypkt.next_trans;

          return {dist + ldist, true};
        }
      } else {
        /// continuum process occurs before reaching the line

        pkt_ptr->next_trans = dummypkt.next_trans - 1;

        return {dist + (tau_rnd - tau) / chi_cont, false};
      }
    } else {
      /// no line interaction possible - check whether continuum process occurs in cell

      const double tau_cont = chi_cont * (abort_dist - dist);

      if (tau_rnd - tau > tau_cont) {
        // no continuum event before abort_dist
        return {std::numeric_limits<double>::max(), false};
      }
      /// continuum process occurs at edist

      pkt_ptr->next_trans = globals::nlines + 1;

      return {dist + (tau_rnd - tau) / chi_cont, false};
    }
  }

  // should have already returned somewhere!
  assert_always(false);
}

static void electron_scatter_rpkt(struct packet *pkt_ptr) {
  /// now make the packet a r-pkt and set further flags
  pkt_ptr->type = TYPE_RPKT;
  pkt_ptr->last_cross = BOUNDARY_NONE;  /// allow all further cell crossings

  const auto vel_vec = get_velocity(pkt_ptr->pos, pkt_ptr->prop_time);

  // Transform Stokes Parameters from the RF to the CMF

  double Qi = pkt_ptr->stokes[1];
  double Ui = pkt_ptr->stokes[2];

  auto old_dir_cmf = frame_transform(pkt_ptr->dir, &Qi, &Ui, vel_vec);

  // Outcoming direction. Compute the new cmf direction from the old direction and the scattering angles (see Kalos &
  // Whitlock 2008)
  double M = 0.;
  double mu = 0.;
  double phisc = 0.;

  if constexpr (DIPOLE) {
    // Assume dipole function (rejecton method, see Code & Whitney 1995)
    double p = 0.;
    double x = 1.;
    while (x > p) {
      const double zrand = rng_uniform();
      const double zrand2 = rng_uniform();
      const double zrand3 = rng_uniform();

      M = 2 * zrand - 1;
      mu = pow(M, 2.);
      phisc = 2 * PI * zrand2;

      // NB: the rotational matrix R here is chosen in the clockwise direction ("+").
      // In Bulla+2015 equation (10) and (12) refer to the specific case shown in Fig.2 where the angle i2
      // is measured in the counter-clockwise direction. Therefore we use the clockwise rotation matrix but
      // with -i1. Here, instead, we calculate the angle in the clockwise direction from 0 to 2PI.
      // For instance, the i1 angle in Fig.2 of Bulla+2015 corresponds to 2PI-i1 here.
      // NB2: the i1 and i2 angles computed in the code (before and after scattering) are instead as in Bulla+2015
      p = (mu + 1) + (mu - 1) * (cos(2 * phisc) * Qi + sin(2 * phisc) * Ui);

      // generate a number between 0 and the maximum of the previous function (2)
      x = 2 * zrand3;
    };
  } else {
    // Assume isotropic scattering
    const double zrand = rng_uniform();
    const double zrand2 = rng_uniform();

    M = 2. * zrand - 1;
    mu = pow(M, 2.);
    phisc = 2 * PI * zrand2;
  }

  const double tsc = acos(M);
  double new_dir_cmf[3];

  if (fabs(old_dir_cmf[2]) < 0.99999) {
    new_dir_cmf[0] = sin(tsc) / sqrt(1. - pow(old_dir_cmf[2], 2.)) *
                         (old_dir_cmf[1] * sin(phisc) - old_dir_cmf[0] * old_dir_cmf[2] * cos(phisc)) +
                     old_dir_cmf[0] * cos(tsc);
    new_dir_cmf[1] = sin(tsc) / sqrt(1 - pow(old_dir_cmf[2], 2.)) *
                         (-old_dir_cmf[0] * sin(phisc) - old_dir_cmf[1] * old_dir_cmf[2] * cos(phisc)) +
                     old_dir_cmf[1] * cos(tsc);
    new_dir_cmf[2] = sin(tsc) * cos(phisc) * sqrt(1 - pow(old_dir_cmf[2], 2.)) + old_dir_cmf[2] * cos(tsc);
  } else {
    new_dir_cmf[0] = sin(tsc) * cos(phisc);
    new_dir_cmf[1] = sin(tsc) * sin(phisc);
    if (old_dir_cmf[2] > 0) {
      new_dir_cmf[2] = cos(tsc);
    } else {
      new_dir_cmf[2] = -cos(tsc);
    }
  }

  // Need to rotate Stokes Parameters in the scattering plane

  std::array<double, 3> ref1{};
  auto ref2 = meridian(old_dir_cmf, ref1);

  // This is the i1 angle of Bulla+2015, obtained by computing the angle between the
  // reference axes ref1 and ref2 in the meridian frame and the corresponding axes
  // ref1_sc and ref2_sc in the scattering plane. It is the supplementary angle of the
  // scatt angle phisc chosen in the rejection technique above (phisc+i1=180 or phisc+i1=540)
  const double i1 = rot_angle(old_dir_cmf, new_dir_cmf, ref1, ref2);
  const double cos2i1 = cos(2 * i1);
  const double sin2i1 = sin(2 * i1);

  const double Qold = Qi * cos2i1 - Ui * sin2i1;
  const double Uold = Qi * sin2i1 + Ui * cos2i1;

  // Scattering

  mu = dot(old_dir_cmf, new_dir_cmf);

  const double Inew = 0.75 * ((mu * mu + 1.0) + Qold * (mu * mu - 1.0));
  double Qnew = 0.75 * ((mu * mu - 1.0) + Qold * (mu * mu + 1.0));
  double Unew = 1.5 * mu * Uold;

  Qnew = Qnew / Inew;
  Unew = Unew / Inew;
  const double I = 1.;  // Inew / Inew

  // Need to rotate Stokes Parameters out of the scattering plane to the meridian frame (Clockwise rotation of PI-i2)

  ref2 = meridian(new_dir_cmf, ref1);

  // This is the i2 angle of Bulla+2015, obtained from the angle THETA between the
  // reference axes ref1_sc and ref2_sc in the scattering plane and ref1 and ref2 in the
  // meridian frame. NB: we need to add PI to transform THETA to i2
  const double i2 = PI + rot_angle(new_dir_cmf, old_dir_cmf, ref1, ref2);
  const double cos2i2 = cos(2 * i2);
  const double sin2i2 = sin(2 * i2);

  double Q = Qnew * cos2i2 + Unew * sin2i2;
  double U = -Qnew * sin2i2 + Unew * cos2i2;

  // Transform Stokes Parameters from the CMF to the RF
  std::array<double, 3> vel_rev{-vel_vec[0], -vel_vec[1], -vel_vec[2]};

  // Update rest frame direction, frequency and energy
  pkt_ptr->dir = frame_transform(new_dir_cmf, &Q, &U, vel_rev);

  pkt_ptr->stokes = {I, Q, U};

  // Check unit vector
  assert_testmodeonly(fabs(vec_len(pkt_ptr->dir) - 1.) < 1.e-6);

  // Finally we want to put in the rest frame energy and frequency.
  // And record that it's now a r-pkt.

  const double dopplerfactor = doppler_packet_nucmf_on_nurf(pkt_ptr->pos, pkt_ptr->dir, pkt_ptr->prop_time);
  pkt_ptr->nu_rf = pkt_ptr->nu_cmf / dopplerfactor;
  pkt_ptr->e_rf = pkt_ptr->e_cmf / dopplerfactor;
}

static void rpkt_event_continuum(struct packet *pkt_ptr, const struct rpkt_continuum_absorptioncoeffs &chi_rpkt_cont,
                                 int modelgridindex) {
  const double nu = pkt_ptr->nu_cmf;

  const double dopplerfactor = doppler_packet_nucmf_on_nurf(pkt_ptr->pos, pkt_ptr->dir, pkt_ptr->prop_time);
  const double chi_cont = chi_rpkt_cont.total * dopplerfactor;
  const double chi_escatter = chi_rpkt_cont.es * dopplerfactor;
  const double chi_ff = chi_rpkt_cont.ff * dopplerfactor;
  const double chi_bf = chi_rpkt_cont.bf * dopplerfactor;

  /// continuum process happens. select due to its probabilities sigma/chi_cont, chi_ff/chi_cont,
  /// chi_bf/chi_cont
  const double zrand = rng_uniform();
  // printout("[debug] rpkt_event:   r-pkt undergoes a continuum transition\n");
  // printout("[debug] rpkt_event:   zrand*chi_cont %g, sigma %g, chi_ff %g, chi_bf %g\n", zrand * chi_cont,
  // sigma, chi_ff, chi_bf);

  if (zrand * chi_cont < chi_escatter) {
    /// electron scattering occurs
    /// in this case the packet stays a R_PKT of same nu_cmf as before (coherent scattering)
    /// but with different direction
    // printout("[debug] rpkt_event:   electron scattering\n");
    pkt_ptr->interactions += 1;
    pkt_ptr->nscatterings += 1;
    pkt_ptr->last_event = LASTEVENT_ELECTRONSCATTERING;
    stats::increment(stats::COUNTER_ESCOUNTER);

    // generate a virtual packet
    vpkt_call_estimators(pkt_ptr, TYPE_RPKT);

    // pkt_ptr->nu_cmf = 3.7474058e+14;
    electron_scatter_rpkt(pkt_ptr);

    /// Electron scattering does not modify the last emission flag
    /// but it updates the last emission position
    pkt_ptr->em_pos = pkt_ptr->pos;
    pkt_ptr->em_time = pkt_ptr->prop_time;

  } else if (zrand * chi_cont < chi_escatter + chi_ff) {
    /// ff: transform to k-pkt
    // printout("[debug] rpkt_event:   free-free transition\n");
    stats::increment(stats::COUNTER_K_STAT_FROM_FF);
    pkt_ptr->interactions += 1;
    pkt_ptr->last_event = 5;
    pkt_ptr->type = TYPE_KPKT;
    pkt_ptr->absorptiontype = -1;
  } else if (zrand * chi_cont < chi_escatter + chi_ff + chi_bf) {
    /// bf: transform to k-pkt or activate macroatom corresponding to probabilities
    // printout("[debug] rpkt_event:   bound-free transition\n");

    pkt_ptr->absorptiontype = -2;

    const double chi_bf_inrest = chi_rpkt_cont.bf;
    assert_always(globals::phixslist[tid].chi_bf_sum[globals::nbfcontinua - 1] == chi_bf_inrest);

    /// Determine in which continuum the bf-absorption occurs
    const double zrand2 = rng_uniform();
    const double chi_bf_rand = zrand2 * chi_bf_inrest;

    // first chi_bf_sum[i] such that chi_bf_sum[i] > chi_bf_rand
    const auto &upperval = std::upper_bound(globals::phixslist[tid].chi_bf_sum,
                                            &globals::phixslist[tid].chi_bf_sum[globals::nbfcontinua - 1], chi_bf_rand);
    const int allcontindex = std::distance(globals::phixslist[tid].chi_bf_sum, upperval);
    assert_always(allcontindex < globals::nbfcontinua);

    const double nu_edge = globals::allcont[allcontindex].nu_edge;
    const int element = globals::allcont[allcontindex].element;
    const int ion = globals::allcont[allcontindex].ion;
    const int level = globals::allcont[allcontindex].level;
    const int phixstargetindex = globals::allcont[allcontindex].phixstargetindex;

    // printout("[debug] rpkt_event:   bound-free: element %d, ion+1 %d, upper %d, ion %d, lower %d\n", element, ion +
    // 1, 0, ion, level); printout("[debug] rpkt_event:   bound-free: nu_edge %g, nu %g\n", nu_edge, nu);

    if constexpr (TRACK_ION_STATS) {
      stats::increment_ion_stats_contabsorption(pkt_ptr, modelgridindex, element, ion);
    }

    /// and decide whether we go to ionisation energy
    const double zrand3 = rng_uniform();
    if (zrand3 < nu_edge / nu) {
      stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_BF);
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 3;

      if constexpr (TRACK_ION_STATS) {
        stats::increment_ion_stats(modelgridindex, element, ion + 1, stats::ION_MACROATOM_ENERGYIN_PHOTOION,
                                   pkt_ptr->e_cmf);
      }

      pkt_ptr->type = TYPE_MA;
      pkt_ptr->mastate.element = element;
      pkt_ptr->mastate.ion = ion + 1;
      const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
      pkt_ptr->mastate.level = upper;
      pkt_ptr->mastate.activatingline = -99;
    }
    /// or to the thermal pool
    else {
      /// transform to k-pkt
      // printout("[debug] rpkt_event:   bound-free: transform to k-pkt\n");
      stats::increment(stats::COUNTER_K_STAT_FROM_BF);
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 4;
      pkt_ptr->type = TYPE_KPKT;
    }
  } else {
    printout("ERROR: could not continuum process\n");
    std::abort();
  }
}

static void rpkt_event_boundbound(struct packet *pkt_ptr, const int mgi) {
  /// bound-bound transition occured
  /// activate macro-atom in corresponding upper-level. Actually all the information
  /// about the macro atoms state has already been set by closest_transition, so
  /// we need here just the activation!
  // printout("[debug] rpkt_event: bound-bound activation of macroatom\n");
  // if (tid == 0) ma_stat_activation_bb++;
  stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_BB);
  pkt_ptr->interactions += 1;
  pkt_ptr->last_event = 1;

  pkt_ptr->absorptiontype = pkt_ptr->mastate.activatingline;
  pkt_ptr->absorptionfreq = pkt_ptr->nu_rf;  // pkt_ptr->nu_cmf;
  pkt_ptr->absorptiondir[0] = pkt_ptr->dir[0];
  pkt_ptr->absorptiondir[1] = pkt_ptr->dir[1];
  pkt_ptr->absorptiondir[2] = pkt_ptr->dir[2];
  pkt_ptr->type = TYPE_MA;

  if constexpr (TRACK_ION_STATS) {
    const int element = pkt_ptr->mastate.element;
    const int ion = pkt_ptr->mastate.ion;
    stats::increment_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_RADEXC, pkt_ptr->e_cmf);

    const int et = pkt_ptr->emissiontype;
    if (et >= 0) {
      const int emissionelement = globals::linelist[et].elementindex;
      const int emissionion = globals::linelist[et].ionindex;
      stats::increment_ion_stats(mgi, emissionelement, emissionion, stats::ION_BOUNDBOUND_ABSORBED,
                                 pkt_ptr->e_cmf / H / pkt_ptr->nu_cmf);
    }
  }

  if constexpr (RECORD_LINESTAT) {
    safeincrement(globals::acounter[pkt_ptr->next_trans - 1]);
  }
}

auto sample_planck_times_expansion_opacity(const int modelgridindex) -> double
// returns a randomly chosen frequency with a distribution of Planck function times the expansion opacity
{
  assert_testmodeonly(EXPANSION_OPAC_SAMPLE_KAPPAPLANCK);
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
  const auto *kappa_planck_bins = &expansionopacity_planck_cumulative[nonemptymgi * expopac_nbins];

  const auto rnd_integral = rng_uniform() * kappa_planck_bins[expopac_nbins - 1];
  const auto *selected_partintegral =
      std::lower_bound(kappa_planck_bins, kappa_planck_bins + expopac_nbins, rnd_integral);
  const auto binindex = std::min(std::distance(kappa_planck_bins, selected_partintegral), expopac_nbins - 1);
  // use a linear interpolation for the frequency within the bin
  const auto bin_nu_lower = get_expopac_bin_nu_lower(binindex);
  const auto delta_nu = get_expopac_bin_nu_upper(binindex) - bin_nu_lower;
  const auto lower_partintegral = (binindex > 0) ? kappa_planck_bins[binindex - 1] : 0.;
  const double nuoffset =
      (rnd_integral - kappa_planck_bins[binindex]) / (kappa_planck_bins[binindex] - lower_partintegral) * delta_nu;

  return bin_nu_lower + nuoffset;
}

static void rpkt_event_thickcell(struct packet *pkt_ptr)
/// Event handling for optically thick cells. Those cells are treated in a grey
/// approximation with electron scattering only.
/// The packet stays an R_PKT of same nu_cmf than before (coherent scattering)
/// but with different direction.
{
  // printout("[debug] rpkt_event_thickcell:   electron scattering\n");
  pkt_ptr->interactions += 1;
  pkt_ptr->nscatterings += 1;
  pkt_ptr->last_event = LASTEVENT_ELECTRONSCATTERING;
  stats::increment(stats::COUNTER_ESCOUNTER);

  emit_rpkt(pkt_ptr);
  /// Electron scattering does not modify the last emission flag
  // pkt_ptr->emissiontype = get_continuumindex(element,ion-1,lower);
  /// but it updates the last emission position
  pkt_ptr->em_pos = pkt_ptr->pos;
  pkt_ptr->em_time = pkt_ptr->prop_time;
}

static void update_estimators(const double e_cmf, const double nu_cmf, const double distance,
                              const double doppler_nucmf_on_nurf, const int modelgridindex, const int nonemptymgi)
/// Update the volume estimators J and nuJ
/// This is done in another routine than move, as we sometimes move dummy
/// packets which do not contribute to the radiation field.
{
  /// Update only non-empty cells
  if (nonemptymgi < 0) {
    return;
  }
  const double distance_e_cmf = distance * e_cmf;

  radfield::update_estimators(nonemptymgi, distance_e_cmf, nu_cmf, doppler_nucmf_on_nurf);

  /// ffheatingestimator does not depend on ion and element, so an array with gridsize is enough.
  safeadd(globals::ffheatingestimator[modelgridindex], distance_e_cmf * globals::chi_rpkt_cont[tid].ffheating);

  if constexpr (USE_LUT_PHOTOION || USE_LUT_BFHEATING) {
    for (int i = 0; i < globals::nbfcontinua_ground; i++) {
      const double nu_edge = globals::groundcont[i].nu_edge;
      if (nu_cmf <= nu_edge) {
        // because groundcont is sorted by nu_edge descending, nu < nu_edge for all remaining items
        return;
      }
      const int element = globals::groundcont[i].element;
      /// Cells with zero abundance for a specific element have zero contribution
      /// (set in calculate_chi_rpkt_cont and therefore do not contribute to
      /// the estimators
      if (grid::get_elem_abundance(modelgridindex, element) > 0) {
        const int ion = globals::groundcont[i].ion;
        const int ionestimindex = get_ionestimindex_nonemptymgi(nonemptymgi, element, ion);

        if constexpr (USE_LUT_PHOTOION) {
          safeadd(globals::gammaestimator[ionestimindex],
                  globals::phixslist[tid].groundcont_gamma_contr[i] * (distance_e_cmf / nu_cmf));

          if (!std::isfinite(globals::gammaestimator[ionestimindex])) {
            printout(
                "[fatal] update_estimators: gamma estimator becomes non finite: mgi %d element %d ion %d gamma_contr "
                "%g, distance_e_cmf_over_nu %g\n",
                modelgridindex, element, ion, globals::phixslist[tid].groundcont_gamma_contr[i],
                distance_e_cmf / nu_cmf);
            std::abort();
          }
        }

        if constexpr (USE_LUT_BFHEATING) {
          safeadd(globals::bfheatingestimator[ionestimindex],
                  globals::phixslist[tid].groundcont_gamma_contr[i] * distance_e_cmf * (1. - nu_edge / nu_cmf));
        }
      }
    }
  }
}

static auto do_rpkt_step(struct packet *pkt_ptr, struct rpkt_continuum_absorptioncoeffs &chi_rpkt_cont, const double t2)
    -> bool
// Routine for moving an r-packet. Similar to do_gamma in objective.
// return value - true if no mgi change, no pkttype change and not reached end of timestep, false otherwise
{
  const int cellindex = pkt_ptr->where;
  const int mgi = grid::get_cell_modelgridindex(cellindex);
  const int nonemptymgi = (mgi != grid::get_npts_model()) ? grid::get_modelcell_nonemptymgi(mgi) : -1;

  // if (pkt_ptr->next_trans > 0) {
  //   printout("[debug] do_rpkt: init: pkt_ptr->nu_cmf %g, nu(pkt_ptr->next_trans=%d)
  //   %g, nu(pkt_ptr->next_trans-1=%d) %g, pkt_ptr->where %d\n", pkt_ptr->nu_cmf, pkt_ptr->next_trans,
  //   globals::linelist[pkt_ptr->next_trans].nu, pkt_ptr->next_trans-1, globals::linelist[pkt_ptr->next_trans-1].nu,
  //   pkt_ptr->where );
  // }

  // Assign optical depth to next physical event. And start counter of
  // optical depth for this path.
  const double zrand = rng_uniform_pos();
  const double tau_next = -1. * log(zrand);

  // Start by finding the distance to the crossing of the grid cell
  // boundaries. sdist is the boundary distance and snext is the
  // grid cell into which we pass.
  auto [sdist, snext] =
      grid::boundary_distance(pkt_ptr->dir, pkt_ptr->pos, pkt_ptr->prop_time, pkt_ptr->where, &pkt_ptr->last_cross);

  if (sdist == 0) {
    grid::change_cell(pkt_ptr, snext);
    const int cellindexnew = pkt_ptr->where;
    const int newmgi = grid::get_cell_modelgridindex(cellindexnew);

    return (pkt_ptr->type == TYPE_RPKT && (newmgi == grid::get_npts_model() || newmgi == mgi));
  }
  const double maxsdist = (GRID_TYPE == GRID_CARTESIAN3D)
                              ? globals::rmax * pkt_ptr->prop_time / globals::tmin
                              : 2 * globals::rmax * (pkt_ptr->prop_time + sdist / CLIGHT_PROP) / globals::tmin;
  if (sdist > maxsdist) {
    printout("[fatal] do_rpkt: Unreasonably large sdist for packet %d. Rpkt. Abort. %g %g %g\n", pkt_ptr->number,
             globals::rmax, pkt_ptr->prop_time / globals::tmin, sdist);
    std::abort();
  }

  if (sdist < 0) {
    const int cellindexnew = pkt_ptr->where;
    printout("[warning] r_pkt: Negative distance (sdist = %g). Abort.\n", sdist);
    printout("[warning] r_pkt: cell %d snext %d\n", cellindexnew, snext);
    printout("[warning] r_pkt: pos %g %g %g\n", pkt_ptr->pos[0], pkt_ptr->pos[1], pkt_ptr->pos[2]);
    printout("[warning] r_pkt: dir %g %g %g\n", pkt_ptr->dir[0], pkt_ptr->dir[1], pkt_ptr->dir[2]);
    printout("[warning] r_pkt: cell corner %g %g %g\n",
             grid::get_cellcoordmin(cellindexnew, 0) * pkt_ptr->prop_time / globals::tmin,
             grid::get_cellcoordmin(cellindexnew, 1) * pkt_ptr->prop_time / globals::tmin,
             grid::get_cellcoordmin(cellindexnew, 2) * pkt_ptr->prop_time / globals::tmin);
    printout("[warning] r_pkt: cell width %g\n", grid::wid_init(cellindexnew, 0) * pkt_ptr->prop_time / globals::tmin);
    assert_always(false);
  }
  if (((snext != -99) && (snext < 0)) || (snext >= grid::ngrid)) {
    printout("[fatal] r_pkt: Heading for inappropriate grid cell. Abort.\n");
    printout("[fatal] r_pkt: Current cell %d, target cell %d.\n", pkt_ptr->where, snext);
    std::abort();
  }

  if (sdist > globals::max_path_step) {
    sdist = globals::max_path_step;
    snext = pkt_ptr->where;
  }

  // At present there is no scattering/destruction process so all that needs to
  // happen is that we determine whether the packet reaches the boundary during the timestep.

  // Find how far it can travel during the time inverval.

  const double tdist = (t2 - pkt_ptr->prop_time) * CLIGHT_PROP;

  assert_always(tdist >= 0);

  const double abort_dist = std::min(tdist, sdist);

  /// Get distance to the next physical event (continuum or bound-bound)
  double edist = -1;
  bool event_is_boundbound = true;
  if (mgi == grid::get_npts_model()) {
    /// for empty cells no physical event occurs. The packets just propagate.
    edist = std::numeric_limits<double>::max();
    pkt_ptr->next_trans = -1;  // skip over lines and search for line list position on the next non-empty cell
  } else if (grid::modelgrid[mgi].thick == 1) {
    /// In the case of optically thick cells, we treat the packets in grey approximation to speed up the calculation

    const double chi_grey = grid::get_kappagrey(mgi) * grid::get_rho(mgi) *
                            doppler_packet_nucmf_on_nurf(pkt_ptr->pos, pkt_ptr->dir, pkt_ptr->prop_time);

    edist = tau_next / chi_grey;
    pkt_ptr->next_trans = -1;
  } else if (EXPANSIONOPACITIES_ON) {
    std::tie(edist, event_is_boundbound) =
        get_event_expansion_opacity(mgi, nonemptymgi, pkt_ptr, chi_rpkt_cont, tau_next, abort_dist);
    pkt_ptr->next_trans = -1;
  } else {
    std::tie(edist, event_is_boundbound) = get_event(mgi, pkt_ptr, chi_rpkt_cont, tau_next, abort_dist);
  }
  assert_always(edist >= 0);

  const auto min_dist = std::min({sdist, abort_dist});

  if (edist == min_dist) {
    // bound-bound or continuum event
    const double doppler_nucmf_on_nurf = move_pkt_withtime(pkt_ptr, edist / 2.);
    update_estimators(pkt_ptr->e_cmf, pkt_ptr->nu_cmf, edist, doppler_nucmf_on_nurf, mgi, nonemptymgi);
    move_pkt_withtime(pkt_ptr, edist / 2.);

    // The previously selected and in pkt_ptr stored event occurs. Handling is done by rpkt_event
    if (grid::modelgrid[mgi].thick == 1) {
      rpkt_event_thickcell(pkt_ptr);
    } else if (event_is_boundbound) {
      if constexpr (EXPANSIONOPACITIES_ON) {
        if constexpr (EXPANSION_OPAC_SAMPLE_KAPPAPLANCK) {
          pkt_ptr->nu_cmf = sample_planck_times_expansion_opacity(mgi);
        }
        rpkt_event_thickcell(pkt_ptr);
      } else {
        rpkt_event_boundbound(pkt_ptr, mgi);
      }
    } else {
      rpkt_event_continuum(pkt_ptr, globals::chi_rpkt_cont[tid], mgi);
    }

    return (pkt_ptr->type == TYPE_RPKT);
  }

  if (sdist == min_dist) {
    // Move it into the new cell.
    const double doppler_nucmf_on_nurf = move_pkt_withtime(pkt_ptr, sdist / 2.);
    update_estimators(pkt_ptr->e_cmf, pkt_ptr->nu_cmf, sdist, doppler_nucmf_on_nurf, mgi, nonemptymgi);
    move_pkt_withtime(pkt_ptr, sdist / 2.);

    int newmgi = mgi;
    if (snext != pkt_ptr->where) {
      grid::change_cell(pkt_ptr, snext);
      const int cellindexnew = pkt_ptr->where;
      newmgi = grid::get_cell_modelgridindex(cellindexnew);
    }

    pkt_ptr->last_event = pkt_ptr->last_event + 100;

    return (pkt_ptr->type == TYPE_RPKT && (newmgi == grid::get_npts_model() || newmgi == mgi));
  }

  if (tdist == min_dist) {
    // reaches end of timestep before cell boundary or interaction
    const double doppler_nucmf_on_nurf = move_pkt_withtime(pkt_ptr, tdist / 2.);
    update_estimators(pkt_ptr->e_cmf, pkt_ptr->nu_cmf, tdist, doppler_nucmf_on_nurf, mgi, nonemptymgi);
    move_pkt_withtime(pkt_ptr, tdist / 2.);
    pkt_ptr->prop_time = t2;
    pkt_ptr->last_event = pkt_ptr->last_event + 1000;

    return false;
  }

  printout("[fatal] do_rpkt: Failed to identify event . Rpkt. edist %g, sdist %g, tdist %g Abort.\n", edist, sdist,
           tdist);
  printout("[fatal] do_rpkt: Trouble was due to packet number %d.\n", pkt_ptr->number);
  std::abort();
}

void do_rpkt(struct packet *pkt_ptr, const double t2, struct rpkt_continuum_absorptioncoeffs &chi_rpkt_cont) {
  while (do_rpkt_step(pkt_ptr, chi_rpkt_cont, t2)) {
    ;
  }
}

void emit_rpkt(struct packet *pkt_ptr) {
  /// now make the packet a r-pkt and set further flags
  pkt_ptr->type = TYPE_RPKT;
  pkt_ptr->last_cross = BOUNDARY_NONE;  /// allow all further cell crossings

  /// Need to assign a new direction. Assume isotropic emission in the cmf

  const auto dir_cmf = get_rand_isotropic_unitvec();

  /// This direction is in the cmf - we want to convert it to the rest
  /// frame - use aberation of angles. We want to convert from cmf to
  /// rest so need -ve velocity.
  const auto vel_vec = get_velocity(pkt_ptr->pos, -1. * pkt_ptr->prop_time);
  /// negative time since we want the backwards transformation here

  pkt_ptr->dir = angle_ab(dir_cmf, vel_vec);
  // printout("[debug] pkt_ptr->dir in RF: %g %g %g\n",pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);

  /// Finally we want to put in the rest frame energy and frequency. And record
  /// that it's now a r-pkt.

  const double dopplerfactor = doppler_packet_nucmf_on_nurf(pkt_ptr->pos, pkt_ptr->dir, pkt_ptr->prop_time);
  pkt_ptr->nu_rf = pkt_ptr->nu_cmf / dopplerfactor;
  pkt_ptr->e_rf = pkt_ptr->e_cmf / dopplerfactor;

  // Reset polarization information
  pkt_ptr->stokes[0] = 1.;
  pkt_ptr->stokes[1] = 0.;
  pkt_ptr->stokes[2] = 0.;

  pkt_ptr->pol_dir = cross_prod(pkt_ptr->dir, std::array<double, 3>{0., 0., 1.});

  if ((dot(pkt_ptr->pol_dir, pkt_ptr->pol_dir)) < 1.e-8) {
    pkt_ptr->pol_dir = cross_prod(pkt_ptr->dir, std::array<double, 3>{0., 1., 0.});
  }

  pkt_ptr->pol_dir = vec_norm(pkt_ptr->pol_dir);
  // printout("initialise pol state of packet %g, %g, %g, %g,
  // %g\n",pkt_ptr->stokes_qu[0],pkt_ptr->stokes_qu[1],pkt_ptr->pol_dir[0],pkt_ptr->pol_dir[1],pkt_ptr->pol_dir[2]);
  // printout("pkt direction %g, %g, %g\n",pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);
}

static auto calculate_chi_freefree(const int modelgridindex, const double nu) -> double
// calculate the free-free absorption coefficient [cm^-1]
// = kappa(free-free) * nne
{
  assert_always(nu > 0.);
  const double g_ff = 1;

  const auto nne = grid::get_nne(modelgridindex);
  const auto T_e = grid::get_Te(modelgridindex);

  double chi_ff = 0.;
  const int nelements = get_nelements();
  for (int element = 0; element < nelements; element++) {
    for (int ion = 0; ion < get_nions(element); ion++) {
      const double nnion = get_nnion(modelgridindex, element, ion);
      const int ioncharge = get_ionstage(element, ion) - 1;
      chi_ff += ioncharge * ioncharge * g_ff * nnion;
    }
  }
  chi_ff *= 3.69255e8 / sqrt(T_e) * pow(nu, -3) * nne * (1 - exp(-HOVERKB * nu / T_e));

  if (!std::isfinite(chi_ff)) {
    printout("ERRORL: chi_ff is non-infinite mgi %d nne %g nu %g T_e %g\n", modelgridindex, nne, nu, T_e);
    std::abort();
  }

  return chi_ff;
}

template <bool usecellhistupdatephixslist>
auto calculate_chi_bf_gammacontr(const int modelgridindex, const double nu) -> double
// bound-free opacity
{
  double chi_bf_sum = 0.;
  if constexpr (usecellhistupdatephixslist && (USE_LUT_PHOTOION || USE_LUT_BFHEATING)) {
    for (int gphixsindex = 0; gphixsindex < globals::nbfcontinua_ground; gphixsindex++) {
      globals::phixslist[tid].groundcont_gamma_contr[gphixsindex] = 0.;
    }
  }

  const auto T_e = grid::get_Te(modelgridindex);
  const auto nne = grid::get_nne(modelgridindex);
  const auto nnetot = grid::get_nnetot(modelgridindex);

  /// The phixslist is sorted by nu_edge in ascending order (longest to shortest wavelength)
  /// If nu < allcont[i].nu_edge no absorption in any of the following continua
  /// is possible, so set their kappas to zero
  // break the list into nu >= nu_edge and the remainder (nu < nu_edge)

  // first element i such that nu < nu_edge[i]
  // const int lastindex = std::upper_bound(globals::allcont_nu_edge, globals::allcont_nu_edge + globals::nbfcontinua,
  // nu,
  //                                        [](const double &nu, const double &nu_edge) { return nu < nu_edge; }) -
  //                       &globals::allcont_nu_edge[0];
  int i = 0;
  const int nbfcontinua = globals::nbfcontinua;
  for (i = 0; i < nbfcontinua; i++) {
    const int element = globals::allcont[i].element;
    const int ion = globals::allcont[i].ion;
    const int level = globals::allcont[i].level;
    /// The bf process happens only if the current cell contains
    /// the involved atomic species

    if ((DETAILED_BF_ESTIMATORS_ON && grid::get_elem_abundance(modelgridindex, element) > 0) ||
        (!DETAILED_BF_ESTIMATORS_ON && ((get_nnion(modelgridindex, element, ion) / nnetot > 1.e-6) || (level == 0)))) {
      const double nu_edge = globals::allcont[i].nu_edge;
      if (nu < nu_edge) {
        break;
      }
      const double nnlevel = usecellhistupdatephixslist ? get_levelpop(modelgridindex, element, ion, level)
                                                        : calculate_levelpop(modelgridindex, element, ion, level);
      const double nu_max_phixs = nu_edge * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table

      if (nu <= nu_max_phixs && nnlevel > 0) {
        const double sigma_bf = photoionization_crosssection_fromtable(globals::allcont[i].photoion_xs, nu_edge, nu);

        const double probability = globals::allcont[i].probability;

        double corrfactor = 1.;  // default to no subtraction of stimulated recombination
        if constexpr (!SEPARATE_STIMRECOMB) {
          double departure_ratio = globals::cellcache[tid].ch_allcont_departureratios[i];
          if (!usecellhistupdatephixslist || departure_ratio < 0) {
            const int upper = globals::allcont[i].upperlevel;
            const double nnupperionlevel = usecellhistupdatephixslist
                                               ? get_levelpop(modelgridindex, element, ion + 1, upper)
                                               : calculate_levelpop(modelgridindex, element, ion + 1, upper);
            const double sf = calculate_sahafact(element, ion, level, upper, T_e, H * nu_edge);
            departure_ratio = nnupperionlevel / nnlevel * nne * sf;  // put that to phixslist
            if (usecellhistupdatephixslist) {
              globals::cellcache[tid].ch_allcont_departureratios[i] = departure_ratio;
            }
          }

          const double stimfactor = departure_ratio * exp(-HOVERKB * nu / T_e);
          corrfactor = std::max(0., 1 - stimfactor);  // photoionisation minus stimulated recombination
        }

        const double sigma_contr = sigma_bf * probability * corrfactor;

        if constexpr (usecellhistupdatephixslist && (USE_LUT_PHOTOION || USE_LUT_BFHEATING)) {
          if (level == 0) {
            const int gphixsindex = globals::allcont[i].index_in_groundphixslist;
            globals::phixslist[tid].groundcont_gamma_contr[gphixsindex] += sigma_contr;
          }
        }

        if constexpr (usecellhistupdatephixslist && DETAILED_BF_ESTIMATORS_ON) {
          globals::phixslist[tid].gamma_contr[i] = sigma_contr;
        }

        const double chi_bf_contr = nnlevel * sigma_contr;
        if (usecellhistupdatephixslist && !std::isfinite(chi_bf_contr)) {
          printout("[fatal] calculate_chi_rpkt_cont: non-finite contribution to chi_bf_contr %g ... abort\n",
                   chi_bf_contr);
          printout("[fatal] phixslist index %d, element %d, ion %d, level %d\n", i, element, ion, level);
          printout("[fatal] Z=%d ionstage %d\n", get_atomicnumber(element), get_ionstage(element, ion));
          printout("[fatal] globals::cell[%d].composition[%d].abundance = %g\n", modelgridindex, element,
                   grid::get_elem_abundance(modelgridindex, element));
          printout("[fatal] nne %g, nnlevel %g, (or %g)\n", grid::get_nne(modelgridindex), nnlevel,
                   get_levelpop(modelgridindex, element, ion, level));
          printout("[fatal] sigma_bf %g, T_e %g, nu %g, nu_edge %g\n", sigma_bf, grid::get_Te(modelgridindex), nu,
                   nu_edge);
          std::abort();
        }

        chi_bf_sum += chi_bf_contr;
        if constexpr (usecellhistupdatephixslist) {
          globals::phixslist[tid].chi_bf_sum[i] = chi_bf_sum;
        }
      } else if constexpr (usecellhistupdatephixslist) {
        // ignore this particular process
        globals::phixslist[tid].chi_bf_sum[i] = chi_bf_sum;
        if constexpr (DETAILED_BF_ESTIMATORS_ON) {
          globals::phixslist[tid].gamma_contr[i] = 0.;
        }
      }
    } else if constexpr (usecellhistupdatephixslist) {
      // no element present or not an important level
      globals::phixslist[tid].chi_bf_sum[i] = chi_bf_sum;
      if constexpr (DETAILED_BF_ESTIMATORS_ON) {
        globals::phixslist[tid].gamma_contr[i] = 0.;
      }
    }
  }

  if constexpr (usecellhistupdatephixslist) {
    for (; i < globals::nbfcontinua; i++) {
      globals::phixslist[tid].chi_bf_sum[i] = chi_bf_sum;
      if constexpr (DETAILED_BF_ESTIMATORS_ON) {
        globals::phixslist[tid].gamma_contr[i] = 0.;
      }
    }
  }

  return chi_bf_sum;
}

void calculate_chi_rpkt_cont(const double nu_cmf, struct rpkt_continuum_absorptioncoeffs &chi_rpkt_cont,
                             const int modelgridindex, const bool usecellhistupdatephixslist) {
  assert_testmodeonly(modelgridindex != grid::get_npts_model());
  assert_testmodeonly(grid::modelgrid[modelgridindex].thick != 1);
  if ((modelgridindex == chi_rpkt_cont.modelgridindex) && (!chi_rpkt_cont.recalculate_required) &&
      (fabs(chi_rpkt_cont.nu / nu_cmf - 1.0) < 1e-4)) {
    // calculated values are a match already
    return;
  }

  const auto nne = grid::get_nne(modelgridindex);

  double sigma = 0.;
  /// free-free absorption
  const double chi_ff = calculate_chi_freefree(modelgridindex, nu_cmf);
  double chi_bf = 0.;
  double chi_ffheating = 0.;

  if (globals::opacity_case == 4) {
    /// First contribution: Thomson scattering on free electrons
    sigma = SIGMA_T * nne;
    // reduced e/s for debugging
    // sigma = 1e-30*sigma;
    // switched off e/s for debugging
    // sigma_cmf = 0. * nne;
    // sigma *= 0.1;

    chi_ffheating = chi_ff;

    /// Third contribution: bound-free absorption
    chi_bf = usecellhistupdatephixslist ? calculate_chi_bf_gammacontr<true>(modelgridindex, nu_cmf)
                                        : calculate_chi_bf_gammacontr<false>(modelgridindex, nu_cmf);

  } else {
    /// in the other cases chi_grey is an mass absorption coefficient
    /// therefore use the mass density
    // sigma = globals::cell[pkt_ptr->where].chi_grey * globals::cell[pkt_ptr->where].rho;
    // sigma = SIGMA_T * nne;

    sigma = 0.;
    // chi_ff = 0.9*sigma;
    // sigma *= 0.1;

    chi_bf = 0.;
  }

  chi_rpkt_cont.nu = nu_cmf;
  chi_rpkt_cont.modelgridindex = modelgridindex;
  chi_rpkt_cont.recalculate_required = false;
  chi_rpkt_cont.total = sigma + chi_bf + chi_ff;
  chi_rpkt_cont.es = sigma;
  chi_rpkt_cont.ff = chi_ff;
  chi_rpkt_cont.bf = chi_bf;
  chi_rpkt_cont.ffheating = chi_ffheating;
  // chi_rpkt_cont_thisthread.bfheating = chi_bfheating;

  if (!std::isfinite(chi_rpkt_cont.total)) {
    printout("[fatal] calculate_chi_rpkt_cont: resulted in non-finite chi_rpkt_cont.total ... abort\n");
    printout("[fatal] es %g, ff %g, bf %g\n", chi_rpkt_cont.es, chi_rpkt_cont.ff, chi_rpkt_cont.bf);
    printout("[fatal] nbfcontinua %d\n", globals::nbfcontinua);
    printout("[fatal] in cell %d with density %g\n", modelgridindex, grid::get_rho(modelgridindex));
    printout("[fatal] pkt_ptr->nu_cmf %g\n", nu_cmf);
    if (std::isfinite(chi_rpkt_cont.es)) {
      chi_rpkt_cont.ff = 0.;
      chi_rpkt_cont.bf = 0.;
      chi_rpkt_cont.total = chi_rpkt_cont.es;
    } else {
      std::abort();
    }
  }
}

void calculate_binned_opacities(const int modelgridindex) {
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
  auto *kappa_bb_bins = &expansionopacities[nonemptymgi * expopac_nbins];
  const auto rho = grid::get_rho(modelgridindex);

  const time_t sys_time_start_calc = time(nullptr);
  const auto temperature = grid::get_TR(modelgridindex);

  printout("calculating binned expansion opacities for cell %d...", modelgridindex);

  const auto t_mid = globals::timesteps[globals::timestep].mid;

  // find the first line with nu below the upper limit of the first bin
  const auto *matchline =
      std::lower_bound(&globals::linelist[0], &globals::linelist[globals::nlines], get_expopac_bin_nu_upper(0),
                       [](const auto &line, const double nu_cmf) -> bool { return line.nu > nu_cmf; });
  int lineindex = std::distance(globals::linelist, matchline);

  for (size_t binindex = 0; binindex < expopac_nbins; binindex++) {
    double bin_linesum = 0.;

    const auto nu_upper = get_expopac_bin_nu_upper(binindex);

    const auto nu_lower = get_expopac_bin_nu_lower(binindex);
    const auto nu_mid = (nu_upper + nu_lower) / 2.;

    const auto delta_nu = nu_upper - nu_lower;

    while (lineindex < globals::nlines && globals::linelist[lineindex].nu >= nu_lower) {
      const float tau_line = get_tau_sobolev(modelgridindex, lineindex, t_mid, false);
      const auto linelambda = 1e8 * CLIGHT / globals::linelist[lineindex].nu;
      bin_linesum += (linelambda / expopac_deltalambda) * -std::expm1(-tau_line);
      lineindex++;
    }

    const float bin_kappa_bb = 1. / (CLIGHT * t_mid * rho) * bin_linesum;
    assert_always(std::isfinite(bin_kappa_bb));
    kappa_bb_bins[binindex] = bin_kappa_bb;

    // calculate_chi_rpkt_cont(nu_mid, globals::chi_rpkt_cont[tid], modelgridindex, false);
    // const auto bin_kappa_cont = globals::chi_rpkt_cont[tid].total / rho;
    const auto bin_kappa_cont = 0.;

    if constexpr (EXPANSION_OPAC_SAMPLE_KAPPAPLANCK) {
      const auto planck_val = radfield::dbb(nu_mid, temperature, 1);
      auto *kappa_planck_bins = &expansionopacity_planck_cumulative[nonemptymgi * expopac_nbins];
      const auto lower_val = ((binindex > 0) ? kappa_planck_bins[binindex - 1] : 0.);
      kappa_planck_bins[binindex] = lower_val + (bin_kappa_bb + bin_kappa_cont) * planck_val * delta_nu;
    }

    // printout("bin %d: lambda %g to %g kappabb %g kappa_cont %g kappa_grey %g kappa_planck_cumulative %g\n", binindex,
    //          1e8 / CLIGHT * nu_upper, 1e8 * CLIGHT / nu_lower, bin_kappa_bb, bin_kappa_cont,
    //          grid::modelgrid[modelgridindex].kappagrey, kappa_planck_cumulative[binindex]);
  }
  printout("took %ld seconds\n", time(nullptr) - sys_time_start_calc);
}
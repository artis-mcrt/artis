#include "rpkt.h"

#ifdef MPI_ON
#include <mpi.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <memory>
#include <span>
#include <tuple>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "ltepop.h"
#include "macroatom.h"
#include "packet.h"
#include "radfield.h"
#include "sn3d.h"
#include "stats.h"
#include "vectors.h"
#include "vpkt.h"

namespace {

constexpr float expopac_lambdamin = 534.5;
constexpr float expopac_lambdamax = 35000.;
constexpr float expopac_deltalambda = 35.5;
constexpr auto expopac_nbins =
    static_cast<std::ptrdiff_t>((expopac_lambdamax - expopac_lambdamin) / expopac_deltalambda);

// kappa in cm^2/g for each bin of each non-empty cell
std::span<float> expansionopacities{};

// kappa times Planck function for each bin of each non-empty cell
std::span<double> expansionopacity_planck_cumulative{};
#ifdef MPI_ON
MPI_Win win_expansionopacities = MPI_WIN_NULL;
MPI_Win win_expansionopacity_planck_cumulative = MPI_WIN_NULL;
#endif

// get the frequency change per distance travelled assuming linear change to the abort distance
// this is done is two parts to get identical results to do_rpkt_step()
auto get_nu_cmf_abort(const std::array<double, 3> &pos, const std::array<double, 3> &dir, const double prop_time,
                      const double nu_rf, const double abort_dist) -> double {
  const auto half_abort_dist = abort_dist / 2.;
  const auto abort_time = prop_time + (half_abort_dist / CLIGHT_PROP) + (half_abort_dist / CLIGHT_PROP);

  const std::array<double, 3> abort_pos{pos[0] + (dir[0] * half_abort_dist) + (dir[0] * half_abort_dist),
                                        pos[1] + (dir[1] * half_abort_dist) + (dir[1] * half_abort_dist),
                                        pos[2] + (dir[2] * half_abort_dist) + (dir[2] * half_abort_dist)};

  const double nu_cmf_abort = nu_rf * calculate_doppler_nucmf_on_nurf(abort_pos, dir, abort_time);

  return nu_cmf_abort;
}

// wavelength bins are ordered by ascending wavelength (descending frequency)

constexpr auto get_expopac_bin_nu_upper(const ptrdiff_t binindex) -> double {
  const auto lambda_lower = expopac_lambdamin + (binindex * expopac_deltalambda);
  return 1e8 * CLIGHT / lambda_lower;
}

constexpr auto get_expopac_bin_nu_lower(const ptrdiff_t binindex) -> double {
  const auto lambda_upper = expopac_lambdamin + ((binindex + 1) * expopac_deltalambda);
  return 1e8 * CLIGHT / lambda_upper;
}

// return edist, the distance to the next physical event (continuum or bound-bound) and is_boundbound_event, a
// boolean BE AWARE THAT THIS PROCEDURE SHOULD BE ONLY CALLED FOR NON EMPTY CELLS!!
auto get_event(const int modelgridindex, const Packet &pkt, const Rpkt_continuum_absorptioncoeffs &chi_rpkt_cont,
               MacroAtomState &mastate,
               const double tau_rnd,     // random optical depth until which the packet travels
               const double abort_dist,  // maximal travel distance before packet leaves cell or time step ends
               const double nu_cmf_abort, const double d_nu_on_d_l, const double doppler, const auto *const linelist,
               const int nlines) -> std::tuple<double, int, bool> {
  assert_testmodeonly(grid::modelgrid[modelgridindex].thick != 1);

  auto pos = pkt.pos;
  auto nu_cmf = pkt.nu_cmf;
  auto e_cmf = pkt.e_cmf;
  auto prop_time = pkt.prop_time;
  int next_trans = pkt.next_trans;

  const double chi_cont = chi_rpkt_cont.total * doppler;
  double tau = 0.;   // optical depth along path
  double dist = 0.;  // position on path
  while (true) {
    // calculate distance to next line encounter ldist
    // first select the closest transition in frequency
    // we need its frequency nu_trans, the element/ion and the corresponding levels
    // create therefore new variables in packet, which contain next_lowerlevel, ...

    // returns negative value if nu_cmf > nu_trans
    if (const int lineindex = closest_transition(nu_cmf, next_trans, nlines, linelist); lineindex >= 0) [[likely]] {
      // line interaction is possible (nu_cmf > nu_trans)

      const double nu_trans = linelist[lineindex].nu;

      // helper variable to overcome numerical problems after line scattering
      // further scattering events should be located at lower frequencies to prevent
      // multiple scattering events of one packet in a single line
      next_trans = lineindex + 1;

      const double ldist = get_linedistance(prop_time, nu_cmf, nu_trans, d_nu_on_d_l);

      const double tau_cont = chi_cont * ldist;

      if (tau_rnd - tau > tau_cont) {
        // got past the continuum optical depth so propagate to the line, and check interaction

        if (nu_trans < nu_cmf_abort) [[unlikely]] {
          // back up one line, because we didn't reach it before the boundary/timelimit

          return {std::numeric_limits<double>::max(), next_trans - 1, false};
        }

        const int element = linelist[lineindex].elementindex;
        const int ion = linelist[lineindex].ionindex;
        const int upper = linelist[lineindex].upperlevelindex;
        const int lower = linelist[lineindex].lowerlevelindex;
        const double A_ul = linelist[lineindex].einstein_A;
        const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
        const double B_lu = stat_weight(element, ion, upper) / stat_weight(element, ion, lower) * B_ul;

        const double n_u = get_levelpop(modelgridindex, element, ion, upper);
        const double n_l = get_levelpop(modelgridindex, element, ion, lower);

        const double tau_line = std::max(0., (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * prop_time);

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
            move_pkt_withtime(pos, pkt.dir, prop_time, pkt.nu_rf, nu_cmf, pkt.e_rf, e_cmf, ldist);
          } else {
            // avoid move_pkt_withtime() to skip the standard Doppler shift calculation
            // and use the linear approx instead
            pos[0] += (pkt.dir[0] * ldist);
            pos[1] += (pkt.dir[1] * ldist);
            pos[2] += (pkt.dir[2] * ldist);
            prop_time += ldist / CLIGHT_PROP;
            nu_cmf = pkt.nu_cmf + d_nu_on_d_l * dist;  // should equal nu_trans;
            assert_testmodeonly(nu_cmf <= pkt.nu_cmf);
          }

          radfield::update_lineestimator(modelgridindex, lineindex, prop_time * CLIGHT * e_cmf / nu_cmf);

        } else {
          // bound-bound process occurs
          // printout("[debug] get_event: tau_rnd - tau <= tau_cont + tau_line: bb-process occurs\n");

          mastate = {.element = element, .ion = ion, .level = upper, .activatingline = lineindex};

          if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
            move_pkt_withtime(pos, pkt.dir, prop_time, pkt.nu_rf, nu_cmf, pkt.e_rf, e_cmf, ldist);
            radfield::update_lineestimator(modelgridindex, lineindex, prop_time * CLIGHT * e_cmf / nu_cmf);
          }

          // the line and its parameters were already selected by closest_transition!
          // printout("[debug] get_event:         edist %g, abort_dist %g, edist-abort_dist %g, endloop
          // %d\n",edist,abort_dist,edist-abort_dist,endloop);

          return {dist + ldist, next_trans, true};
        }
      } else {
        // continuum process occurs before reaching the line

        return {dist + ((tau_rnd - tau) / chi_cont), next_trans - 1, false};
      }
    } else [[unlikely]] {
      // no line interaction possible - check whether continuum process occurs in cell

      const double tau_cont = chi_cont * (abort_dist - dist);

      if (tau_rnd - tau > tau_cont) {
        // no continuum event before abort_dist
        return {std::numeric_limits<double>::max(), next_trans, false};
      }
      // continuum process occurs at edist

      return {dist + ((tau_rnd - tau) / chi_cont), nlines + 1, false};
    }
  }

  // should have already returned somewhere!
  assert_always(false);
}

auto get_event_expansion_opacity(
    const int modelgridindex, const int nonemptymgi, const Packet &pkt,
    const Rpkt_continuum_absorptioncoeffs &chi_rpkt_cont,  // NOLINT(misc-unused-parameters)
    MacroAtomState &mastate, const double tau_rnd, const double nu_cmf_abort, const double d_nu_on_d_l,
    const double doppler) -> std::tuple<double, int, bool> {
  auto pos = pkt.pos;
  const auto nu_rf = pkt.nu_rf;
  auto nu_cmf = pkt.nu_cmf;
  const auto e_rf = pkt.e_rf;
  auto e_cmf = pkt.e_cmf;
  auto prop_time = pkt.prop_time;

  // with thermalisation, we don't keep track of line interactions
  auto next_trans = RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY >= 0. ? -1 : pkt.next_trans;

  assert_always(globals::cellcache[cellcacheslotid].cellnumber == modelgridindex);
  double dist = 0.;
  double tau = 0.;
  auto binindex_start = static_cast<ptrdiff_t>(((1e8 * CLIGHT / nu_cmf) - expopac_lambdamin) / expopac_deltalambda);
  if (binindex_start < 0) {
    binindex_start = -1;
  }

  for (ptrdiff_t binindex = binindex_start; binindex < expopac_nbins; binindex++) {
    const auto next_bin_edge_nu = (binindex < 0) ? get_expopac_bin_nu_upper(0) : get_expopac_bin_nu_lower(binindex);
    const auto binedgedist = get_linedistance(prop_time, nu_cmf, next_bin_edge_nu, d_nu_on_d_l);

    const double chi_cont = chi_rpkt_cont.total * doppler;
    // const auto chi_cont = 0.;
    double chi_bb_expansionopac = 0.;
    if (binindex >= 0) {
      const auto kappa = expansionopacities[(nonemptymgi * expopac_nbins) + binindex];
      chi_bb_expansionopac = kappa * grid::get_rho(modelgridindex) * doppler;
    }

    const double chi_tot = chi_cont + chi_bb_expansionopac;

    if (chi_tot * binedgedist > tau_rnd - tau) {
      // interaction occurs
      if constexpr (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY >= 0.) {
        const auto edist = std::max(dist + ((tau_rnd - tau) / chi_tot), 0.);
        const bool event_is_boundbound = rng_uniform() <= chi_bb_expansionopac / chi_tot;
        return {edist, next_trans, event_is_boundbound};
      } else {
        // re-trace this bin line-by-line
        auto pkt_bin_start = pkt;
        pkt_bin_start.pos = pos;
        pkt_bin_start.nu_rf = nu_rf;
        pkt_bin_start.nu_cmf = nu_cmf;
        pkt_bin_start.e_rf = e_rf;
        pkt_bin_start.e_cmf = e_cmf;
        // expansion opacity was calculated at t_mid, so match it
        pkt_bin_start.prop_time = globals::timesteps[globals::timestep].mid;
        pkt_bin_start.next_trans = -1;
        double edist_after_bin = 0.;
        bool event_is_boundbound = false;
        std::tie(edist_after_bin, next_trans, event_is_boundbound) =
            get_event(modelgridindex, pkt_bin_start, chi_rpkt_cont, mastate, tau_rnd - tau,
                      std::numeric_limits<double>::max(), 0., d_nu_on_d_l, doppler, globals::linelist, globals::nlines);
        // assert_always(edist_after_bin <= 1.1 * binedgedist);
        dist = dist + edist_after_bin;

        return {dist, next_trans, event_is_boundbound};
      }
    }

    tau += chi_tot * binedgedist;
    dist += binedgedist;

    if constexpr (!USE_RELATIVISTIC_DOPPLER_SHIFT) {
      move_pkt_withtime(pos, pkt.dir, prop_time, nu_rf, nu_cmf, e_rf, e_cmf, binedgedist);
    } else {
      // avoid move_pkt_withtime() to skip the standard Doppler shift calculation
      // and use the linear approx instead

      pos[0] += (pkt.dir[0] * binedgedist);
      pos[1] += (pkt.dir[1] * binedgedist);
      pos[2] += (pkt.dir[2] * binedgedist);
      prop_time += binedgedist / CLIGHT_PROP;
      nu_cmf = pkt.nu_cmf + d_nu_on_d_l * dist;  // should equal nu_trans;
      assert_testmodeonly(nu_cmf <= pkt.nu_cmf);
    }

    if (nu_cmf <= nu_cmf_abort) {
      // hit edge of cell or timestep limit
      return {std::numeric_limits<double>::max(), next_trans, false};
    }
  }

  // no more bins, so no opacity and no chance of further interaction below this frequency
  return {std::numeric_limits<double>::max(), next_trans, false};
}

void electron_scatter_rpkt(Packet &pkt) {
  // now make the packet a r-pkt and set further flags
  pkt.type = TYPE_RPKT;
  pkt.last_cross = BOUNDARY_NONE;  // allow all further cell crossings

  const auto vel_vec = get_velocity(pkt.pos, pkt.prop_time);

  // Transform Stokes Parameters from the RF to the CMF

  double Qi = pkt.stokes[1];
  double Ui = pkt.stokes[2];

  const auto old_dir_cmf = frame_transform(pkt.dir, &Qi, &Ui, vel_vec);

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

      M = 2 * zrand - 1;
      mu = pow(M, 2.);
      phisc = 2 * PI * rng_uniform();

      // NB: the rotational matrix R here is chosen in the clockwise direction ("+").
      // In Bulla+2015 equation (10) and (12) refer to the specific case shown in Fig.2 where the angle i2
      // is measured in the counter-clockwise direction. Therefore we use the clockwise rotation matrix but
      // with -i1. Here, instead, we calculate the angle in the clockwise direction from 0 to 2PI.
      // For instance, the i1 angle in Fig.2 of Bulla+2015 corresponds to 2PI-i1 here.
      // NB2: the i1 and i2 angles computed in the code (before and after scattering) are instead as in Bulla+2015
      p = (mu + 1) + (mu - 1) * (cos(2 * phisc) * Qi + sin(2 * phisc) * Ui);

      // generate a number between 0 and the maximum of the previous function (2)
      x = 2. * rng_uniform();
    };
  } else {
    // Assume isotropic scattering
    const double zrand = rng_uniform();

    M = 2. * zrand - 1;
    mu = pow(M, 2.);
    phisc = 2 * PI * rng_uniform();
  }

  const double tsc = acos(M);
  std::array<double, 3> new_dir_cmf{};

  if (fabs(old_dir_cmf[2]) < 0.99999) {
    new_dir_cmf[0] = sin(tsc) / sqrt(1. - pow(old_dir_cmf[2], 2.)) *
                         (old_dir_cmf[1] * sin(phisc) - old_dir_cmf[0] * old_dir_cmf[2] * cos(phisc)) +
                     old_dir_cmf[0] * cos(tsc);
    new_dir_cmf[1] = sin(tsc) / sqrt(1 - pow(old_dir_cmf[2], 2.)) *
                         (-old_dir_cmf[0] * sin(phisc) - old_dir_cmf[1] * old_dir_cmf[2] * cos(phisc)) +
                     old_dir_cmf[1] * cos(tsc);
    new_dir_cmf[2] = sin(tsc) * cos(phisc) * sqrt(1 - pow(old_dir_cmf[2], 2.)) + old_dir_cmf[2] * cos(tsc);
  } else {
    new_dir_cmf = {sin(tsc) * cos(phisc), sin(tsc) * sin(phisc), (old_dir_cmf[2] > 0) ? cos(tsc) : -cos(tsc)};
  }

  // Need to rotate Stokes Parameters in the scattering plane

  const auto [ref1_olddir, ref2_olddir] = meridian(old_dir_cmf);

  // This is the i1 angle of Bulla+2015, obtained by computing the angle between the
  // reference axes ref1 and ref2 in the meridian frame and the corresponding axes
  // ref1_sc and ref2_sc in the scattering plane. It is the supplementary angle of the
  // scatt angle phisc chosen in the rejection technique above (phisc+i1=180 or phisc+i1=540)
  const double i1 = get_rot_angle(old_dir_cmf, new_dir_cmf, ref1_olddir, ref2_olddir);
  const double cos2i1 = cos(2 * i1);
  const double sin2i1 = sin(2 * i1);

  const double Qold = (Qi * cos2i1) - (Ui * sin2i1);
  const double Uold = (Qi * sin2i1) + (Ui * cos2i1);

  // Scattering

  mu = dot(old_dir_cmf, new_dir_cmf);

  const double Inew = 0.75 * ((mu * mu + 1.0) + Qold * (mu * mu - 1.0));
  const double Qnew = (0.75 * ((mu * mu - 1.0) + Qold * (mu * mu + 1.0))) / Inew;
  const double Unew = (1.5 * mu * Uold) / Inew;

  // Need to rotate Stokes Parameters out of the scattering plane to the meridian frame (Clockwise rotation of PI-i2)

  const auto [ref1, ref2] = meridian(new_dir_cmf);

  // This is the i2 angle of Bulla+2015, obtained from the angle THETA between the
  // reference axes ref1_sc and ref2_sc in the scattering plane and ref1 and ref2 in the
  // meridian frame. NB: we need to add PI to transform THETA to i2
  const double i2 = PI + get_rot_angle(new_dir_cmf, old_dir_cmf, ref1, ref2);
  const double cos2i2 = cos(2 * i2);
  const double sin2i2 = sin(2 * i2);

  double Q = (Qnew * cos2i2) + (Unew * sin2i2);
  double U = (-Qnew * sin2i2) + (Unew * cos2i2);

  // Transform Stokes Parameters from the CMF to the RF
  // Update rest frame direction, frequency and energy
  pkt.dir = frame_transform(new_dir_cmf, &Q, &U, std::array<double, 3>{-vel_vec[0], -vel_vec[1], -vel_vec[2]});

  pkt.stokes = {1., Q, U};

  // Check unit vector
  assert_testmodeonly(fabs(vec_len(pkt.dir) - 1.) < 1.e-6);

  // Finally we want to put in the rest frame energy and frequency.
  // And record that it's now a r-pkt.

  const double dopplerfactor = calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);
  pkt.nu_rf = pkt.nu_cmf / dopplerfactor;
  pkt.e_rf = pkt.e_cmf / dopplerfactor;
}

void rpkt_event_continuum(Packet &pkt, const Rpkt_continuum_absorptioncoeffs &chi_rpkt_cont, const int modelgridindex) {
  const double nu = pkt.nu_cmf;

  const double dopplerfactor = calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);
  const double chi_cont = chi_rpkt_cont.total * dopplerfactor;
  const double chi_escatter = chi_rpkt_cont.ffescat * dopplerfactor;
  const double chi_ff = chi_rpkt_cont.ffheat * dopplerfactor;
  const double chi_bf = chi_rpkt_cont.bf * dopplerfactor;

  // continuum process happens. select due to its probabilities sigma/chi_cont, chi_ff/chi_cont,
  // chi_bf/chi_cont
  // printout("[debug] rpkt_event:   r-pkt undergoes a continuum transition\n");
  // printout("[debug] rpkt_event:   zrand*chi_cont %g, sigma %g, chi_ff %g, chi_bf %g\n", zrand * chi_cont,
  // sigma, chi_ff, chi_bf);

  const auto chi_rnd = rng_uniform() * chi_cont;

  if (chi_rnd < chi_escatter) {
    // electron scattering occurs
    // in this case the packet stays a R_PKT of same nu_cmf as before (coherent scattering)
    // but with different direction
    // printout("[debug] rpkt_event:   electron scattering\n");
    stats::increment(stats::COUNTER_INTERACTIONS);
    pkt.nscatterings += 1;
    pkt.last_event = LASTEVENT_ELECTRONSCATTERING;
    stats::increment(stats::COUNTER_ESCOUNTER);

    // generate a virtual packet
    if constexpr (VPKT_ON) {
      vpkt_call_estimators(pkt, TYPE_RPKT);
    }

    // pkt.nu_cmf = 3.7474058e+14;
    electron_scatter_rpkt(pkt);

    // Electron scattering does not modify the last emission flag
    // but it updates the last emission position
    pkt.em_pos = pkt.pos;
    pkt.em_time = pkt.prop_time;

  } else if (chi_rnd < chi_escatter + chi_ff) {
    // ff: transform to k-pkt
    // printout("[debug] rpkt_event:   free-free transition\n");
    stats::increment(stats::COUNTER_K_STAT_FROM_FF);
    stats::increment(stats::COUNTER_INTERACTIONS);
    pkt.last_event = 5;
    pkt.type = TYPE_KPKT;
    pkt.absorptiontype = -1;
  } else if (chi_rnd < chi_escatter + chi_ff + chi_bf) {
    // bf: transform to k-pkt or activate macroatom corresponding to probabilities
    // printout("[debug] rpkt_event:   bound-free transition\n");

    const auto &phixslist = *chi_rpkt_cont.phixslist;

    pkt.absorptiontype = -2;

    const double chi_bf_inrest = chi_rpkt_cont.bf;
    assert_always(phixslist.chi_bf_sum[phixslist.allcontend - 1] == chi_bf_inrest);

    // Determine in which continuum the bf-absorption occurs
    const double chi_bf_rand = rng_uniform() * chi_bf_inrest;

    // first chi_bf_sum[i] such that chi_bf_sum[i] > chi_bf_rand
    const auto allcontindex = std::upper_bound(phixslist.chi_bf_sum.data() + phixslist.allcontbegin,
                                               phixslist.chi_bf_sum.data() + phixslist.allcontend - 1, chi_bf_rand) -
                              phixslist.chi_bf_sum.data();
    assert_always(allcontindex < phixslist.allcontend);

    const double nu_edge = globals::allcont[allcontindex].nu_edge;
    const int element = globals::allcont[allcontindex].element;
    const int ion = globals::allcont[allcontindex].ion;
    const int level = globals::allcont[allcontindex].level;
    const int phixstargetindex = globals::allcont[allcontindex].phixstargetindex;

    // printout("[debug] rpkt_event:   bound-free: element %d, ion+1 %d, upper %d, ion %d, lower %d\n", element, ion +
    // 1, 0, ion, level); printout("[debug] rpkt_event:   bound-free: nu_edge %g, nu %g\n", nu_edge, nu);

    if constexpr (TRACK_ION_STATS) {
      stats::increment_ion_stats_contabsorption(pkt, modelgridindex, element, ion);
    }

    // and decide whether we go to ionisation energy
    if (rng_uniform() < nu_edge / nu) {
      stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_BF);
      stats::increment(stats::COUNTER_INTERACTIONS);
      pkt.last_event = 3;

      if constexpr (TRACK_ION_STATS) {
        stats::increment_ion_stats(modelgridindex, element, ion + 1, stats::ION_MACROATOM_ENERGYIN_PHOTOION, pkt.e_cmf);
      }

      pkt.type = TYPE_MA;
      const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);

      do_macroatom(pkt, {.element = element, .ion = ion + 1, .level = upper, .activatingline = -99});
    }
    // or to the thermal pool
    else {
      // transform to k-pkt
      // printout("[debug] rpkt_event:   bound-free: transform to k-pkt\n");
      stats::increment(stats::COUNTER_K_STAT_FROM_BF);
      stats::increment(stats::COUNTER_INTERACTIONS);
      pkt.last_event = 4;
      pkt.type = TYPE_KPKT;
    }
  } else {
    assert_always(false);
  }
}

// handle bound-bound transition and activate macro-atom in corresponding upper-level
void rpkt_event_boundbound(Packet &pkt, const MacroAtomState &pktmastate, const int mgi) {
  stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_BB);
  stats::increment(stats::COUNTER_INTERACTIONS);
  pkt.last_event = 1;

  pkt.absorptiontype = pktmastate.activatingline;
  pkt.absorptionfreq = pkt.nu_rf;
  pkt.absorptiondir = pkt.dir;
  pkt.type = TYPE_MA;

  if constexpr (TRACK_ION_STATS) {
    const int element = pktmastate.element;
    const int ion = pktmastate.ion;
    stats::increment_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_RADEXC, pkt.e_cmf);

    const int et = pkt.emissiontype;
    if (et >= 0) {
      const int emissionelement = globals::linelist[et].elementindex;
      const int emissionion = globals::linelist[et].ionindex;
      stats::increment_ion_stats(mgi, emissionelement, emissionion, stats::ION_BOUNDBOUND_ABSORBED,
                                 pkt.e_cmf / H / pkt.nu_cmf);
    }
  }

  if constexpr (RECORD_LINESTAT) {
    atomicadd(globals::acounter[pkt.next_trans - 1], 1);
  }

  do_macroatom(pkt, pktmastate);
}

// Handle r-packet interaction in thick cell (grey opacity).
// The packet stays an RPKT of same nu_cmf as before (coherent scattering) but with a different direction.
void rpkt_event_thickcell(Packet &pkt) {
  // printout("[debug] rpkt_event_thickcell:   electron scattering\n");
  stats::increment(stats::COUNTER_INTERACTIONS);
  pkt.nscatterings += 1;
  pkt.last_event = LASTEVENT_ELECTRONSCATTERING;
  stats::increment(stats::COUNTER_ESCOUNTER);

  emit_rpkt(pkt);
  // Electron scattering does not modify the last emission flag but it updates the last emission position
  pkt.em_pos = pkt.pos;
  pkt.em_time = pkt.prop_time;
}

// Update the volume estimators J and nuJ
// This is done in another routine than move, as we sometimes move dummy
// packets which do not contribute to the radiation field.
void update_estimators(const double e_cmf, const double nu_cmf, const double distance,
                       const double doppler_nucmf_on_nurf, const int nonemptymgi,
                       const Rpkt_continuum_absorptioncoeffs &chi_rpkt_cont, const bool thickcell) {
  // Update only non-empty cells
  assert_testmodeonly(nonemptymgi >= 0);
  const double distance_e_cmf = distance * e_cmf;

  radfield::update_estimators(nonemptymgi, distance_e_cmf, nu_cmf, doppler_nucmf_on_nurf, *chi_rpkt_cont.phixslist,
                              thickcell);

  if (thickcell) {
    // chi_rpkt_cont and phixslist are not known for thick cells
    return;
  }

  // ffheatingestimator does not depend on ion and element, so an array with gridsize is enough.
  atomicadd(globals::ffheatingestimator[nonemptymgi], distance_e_cmf * chi_rpkt_cont.ffheat);

  if constexpr (USE_LUT_PHOTOION || USE_LUT_BFHEATING) {
    for (int i = 0; i < globals::nbfcontinua_ground; i++) {
      const double nu_edge = globals::groundcont[i].nu_edge;
      if (nu_cmf <= nu_edge) {
        // because groundcont is sorted by nu_edge descending, nu < nu_edge for all remaining items
        return;
      }
      const int ionestimindex = (nonemptymgi * globals::nbfcontinua_ground) + i;

      if constexpr (USE_LUT_PHOTOION) {
        atomicadd(globals::gammaestimator[ionestimindex],
                  chi_rpkt_cont.phixslist->groundcont_gamma_contr[i] * (distance_e_cmf / nu_cmf));
      }

      if constexpr (USE_LUT_BFHEATING) {
        atomicadd(globals::bfheatingestimator[ionestimindex],
                  chi_rpkt_cont.phixslist->groundcont_gamma_contr[i] * distance_e_cmf * (1. - nu_edge / nu_cmf));
      }
    }
  }
}

// Update an r-packet and return true if no mgi change (or it goes into an empty cell) and no pkttype change and not
// reached end of timestep, otherwise false
auto do_rpkt_step(Packet &pkt, const double t2) -> bool {
  const int cellindex = pkt.where;
  const int mgi = grid::get_cell_modelgridindex(cellindex);
  const int nonemptymgi = (mgi != grid::get_npts_model()) ? grid::get_nonemptymgi_of_mgi(mgi) : -1;

  MacroAtomState pktmastate{};

  THREADLOCALONHOST auto groundcont_gamma_contr = std::make_unique<double[]>(globals::nbfcontinua_ground);
  THREADLOCALONHOST auto chi_bf_sum = std::make_unique<double[]>(globals::nbfcontinua);
  THREADLOCALONHOST auto gamma_contr = std::make_unique<double[]>(globals::bfestimcount);

  THREADLOCALONHOST Phixslist phixslist{
      .groundcont_gamma_contr = std::span(groundcont_gamma_contr.get(), globals::nbfcontinua_ground),
      .chi_bf_sum = std::span(chi_bf_sum.get(), globals::nbfcontinua),
      .gamma_contr = std::span(gamma_contr.get(), globals::bfestimcount),
      .allcontend = 1,
      .allcontbegin = 0,
      .bfestimend = 1,
      .bfestimbegin = 0,
  };

  THREADLOCALONHOST Rpkt_continuum_absorptioncoeffs chi_rpkt_cont{
      .nu = NAN,
      .total = NAN,
      .ffescat = NAN,
      .ffheat = NAN,
      .bf = NAN,
      .modelgridindex = -1,
      .timestep = -1,
      .phixslist = &phixslist,
  };

  // Assign optical depth to next physical event
  const double zrand = rng_uniform_pos();
  const double tau_next = -1. * log(zrand);

  // Start by finding the distance to the crossing of the grid cell
  // boundaries. sdist is the boundary distance and snext is the
  // grid cell into which we pass.
  auto [sdist, snext] = grid::boundary_distance(pkt.dir, pkt.pos, pkt.prop_time, pkt.where, &pkt.last_cross);

  if (sdist == 0) {
    grid::change_cell(pkt, snext);
    const int cellindexnew = pkt.where;
    const int newmgi = grid::get_cell_modelgridindex(cellindexnew);

    return (pkt.type == TYPE_RPKT && (newmgi == grid::get_npts_model() || newmgi == mgi));
  }
  const double maxsdist = (GRID_TYPE == GridType::CARTESIAN3D)
                              ? globals::rmax * pkt.prop_time / globals::tmin
                              : 2 * globals::rmax * (pkt.prop_time + sdist / CLIGHT_PROP) / globals::tmin;
  if (sdist > maxsdist) {
    printout("[fatal] do_rpkt: Unreasonably large sdist for packet %d. Rpkt. Abort. %g %g %g\n", pkt.number,
             globals::rmax, pkt.prop_time / globals::tmin, sdist);
    std::abort();
  }

  if (sdist < 0) {
    const int cellindexnew = pkt.where;
    printout("[warning] r_pkt: Negative distance (sdist = %g). Abort.\n", sdist);
    printout("[warning] r_pkt: cell %d snext %d\n", cellindexnew, snext);
    printout("[warning] r_pkt: pos %g %g %g\n", pkt.pos[0], pkt.pos[1], pkt.pos[2]);
    printout("[warning] r_pkt: dir %g %g %g\n", pkt.dir[0], pkt.dir[1], pkt.dir[2]);
    printout("[warning] r_pkt: cell corner %g %g %g\n",
             grid::get_cellcoordmin(cellindexnew, 0) * pkt.prop_time / globals::tmin,
             grid::get_cellcoordmin(cellindexnew, 1) * pkt.prop_time / globals::tmin,
             grid::get_cellcoordmin(cellindexnew, 2) * pkt.prop_time / globals::tmin);
    printout("[warning] r_pkt: cell width %g\n", grid::wid_init(cellindexnew, 0) * pkt.prop_time / globals::tmin);
    assert_always(false);
  }
  if (((snext != -99) && (snext < 0)) || (snext >= grid::ngrid)) {
    printout("[fatal] r_pkt: Heading for inappropriate grid cell. Abort.\n");
    printout("[fatal] r_pkt: Current cell %d, target cell %d.\n", pkt.where, snext);
    std::abort();
  }

  if (sdist > globals::max_path_step) {
    sdist = globals::max_path_step;
    snext = pkt.where;
  }

  // At present there is no scattering/destruction process so all that needs to
  // happen is that we determine whether the packet reaches the boundary during the timestep.

  // Find how far it can travel during the time interval.

  const double tdist = (t2 - pkt.prop_time) * CLIGHT_PROP;

  assert_always(tdist >= 0);

  const double abort_dist = std::min(tdist, sdist);

  // Get distance to the next physical event (continuum or bound-bound)
  double edist = -1;
  bool event_is_boundbound = true;
  const bool thickcell = grid::modelgrid[mgi].thick == 1;
  if (nonemptymgi < 0) {
    // for empty cells no physical event occurs. The packets just propagate.
    edist = std::numeric_limits<double>::max();
    pkt.next_trans = -1;  // skip over lines and search for line list position on the next non-empty cell
  } else if (thickcell) [[unlikely]] {
    // In the case of optically thick cells, we treat the packets in grey approximation to speed up the calculation

    const double chi_grey = grid::get_kappagrey(mgi) * grid::get_rho(mgi) *
                            calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);

    edist = tau_next / chi_grey;
    pkt.next_trans = -1;
  } else {
    calculate_chi_rpkt_cont(pkt.nu_cmf, chi_rpkt_cont, mgi);

    // for USE_RELATIVISTIC_DOPPLER_SHIFT, we will use a linear approximation for
    // the frequency change from start to abort (cell boundary/timestep end)

    const auto nu_cmf_abort = get_nu_cmf_abort(pkt.pos, pkt.dir, pkt.prop_time, pkt.nu_rf, abort_dist);
    const auto d_nu_on_d_l = (nu_cmf_abort - pkt.nu_cmf) / abort_dist;
    const auto doppler = calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);

    if constexpr (EXPANSIONOPACITIES_ON) {
      std::tie(edist, pkt.next_trans, event_is_boundbound) = get_event_expansion_opacity(
          mgi, nonemptymgi, pkt, chi_rpkt_cont, pktmastate, tau_next, nu_cmf_abort, d_nu_on_d_l, doppler);
    } else {
      std::tie(edist, pkt.next_trans, event_is_boundbound) =
          get_event(mgi, pkt, chi_rpkt_cont, pktmastate, tau_next, abort_dist, nu_cmf_abort, d_nu_on_d_l, doppler,
                    globals::linelist, globals::nlines);
    }
  }
  assert_always(edist >= 0);

  if ((sdist <= tdist) && (sdist <= edist)) {
    // Move it into the new cell.
    const double doppler_nucmf_on_nurf = move_pkt_withtime(pkt, sdist / 2.);
    if (nonemptymgi >= 0) {
      update_estimators(pkt.e_cmf, pkt.nu_cmf, sdist, doppler_nucmf_on_nurf, nonemptymgi, chi_rpkt_cont, thickcell);
    }
    move_pkt_withtime(pkt, sdist / 2.);

    int newmgi = mgi;
    if (snext != pkt.where) {
      grid::change_cell(pkt, snext);
      const int cellindexnew = pkt.where;
      newmgi = grid::get_cell_modelgridindex(cellindexnew);
    }

    pkt.last_event = pkt.last_event + 100;

    return (pkt.type == TYPE_RPKT && (newmgi == grid::get_npts_model() || newmgi == mgi));
  }

  if ((edist <= sdist) && (edist <= tdist)) [[likely]] {
    // bound-bound or continuum event
    const double doppler_nucmf_on_nurf = move_pkt_withtime(pkt, edist / 2.);
    update_estimators(pkt.e_cmf, pkt.nu_cmf, edist, doppler_nucmf_on_nurf, nonemptymgi, chi_rpkt_cont, thickcell);
    move_pkt_withtime(pkt, edist / 2.);

    // The previously selected and in pkt stored event occurs. Handling is done by rpkt_event
    if (thickcell) {
      rpkt_event_thickcell(pkt);
    } else if (event_is_boundbound) {
      if constexpr (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY < 0.) {
        rpkt_event_boundbound(pkt, pktmastate, mgi);
      } else {
        // Probability based thermalisation (i.e. redistibution of the packet frequency) or scattering
        if (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY >= 1. ||
            rng_uniform() < RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY) {
          pkt.nu_cmf = sample_planck_times_expansion_opacity(nonemptymgi);
          // When thermalised, we do not associate the packet with a specific line emission
        }
        rpkt_event_thickcell(pkt);
      }
    } else {
      rpkt_event_continuum(pkt, chi_rpkt_cont, mgi);
    }

    return (pkt.type == TYPE_RPKT);
  }

  if ((tdist <= sdist) && (tdist <= edist)) [[unlikely]] {
    // reaches end of timestep before cell boundary or interaction
    const double doppler_nucmf_on_nurf = move_pkt_withtime(pkt, tdist / 2.);
    if (nonemptymgi >= 0) {
      update_estimators(pkt.e_cmf, pkt.nu_cmf, tdist, doppler_nucmf_on_nurf, nonemptymgi, chi_rpkt_cont, thickcell);
    }
    move_pkt_withtime(pkt, tdist / 2.);
    pkt.prop_time = t2;
    pkt.last_event = pkt.last_event + 1000;

    return false;
  }

  printout("[fatal] do_rpkt: Failed to identify event . Rpkt. edist %g, sdist %g, tdist %g Abort.\n", edist, sdist,
           tdist);
  printout("[fatal] do_rpkt: Trouble was due to packet number %d.\n", pkt.number);
  std::abort();
}

auto calculate_chi_ffheat_nnionpart(const int modelgridindex) -> double {
  const double g_ff = 1;
  double chi_ff_nnionpart = 0.;
  const int nelements = get_nelements();
  for (int element = 0; element < nelements; element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const double nnion = get_nnion(modelgridindex, element, ion);
      const int ioncharge = get_ionstage(element, ion) - 1;
      chi_ff_nnionpart += ioncharge * ioncharge * g_ff * nnion;
    }
  }
  const auto T_e = grid::get_Te(modelgridindex);

  return chi_ff_nnionpart * 3.69255e8 / sqrt(T_e);
}

auto get_chi_ff_nnionpart(const int modelgridindex) -> double {
  if (!use_cellcache || globals::cellcache[cellcacheslotid].cellnumber != modelgridindex) {
    return calculate_chi_ffheat_nnionpart(modelgridindex);
  }

  if (globals::cellcache[cellcacheslotid].chi_ff_nnionpart < 0.) {
    globals::cellcache[cellcacheslotid].chi_ff_nnionpart = calculate_chi_ffheat_nnionpart(modelgridindex);
  }

  return globals::cellcache[cellcacheslotid].chi_ff_nnionpart;
}

// calculate the free-free absorption (to kpkt heating) coefficient [cm^-1]
// = kappa(free-free) * nne
auto calculate_chi_ffheating(const int modelgridindex, const double nu) -> double {
  assert_always(nu > 0.);

  const auto nne = grid::get_nne(modelgridindex);
  const auto T_e = grid::get_Te(modelgridindex);
  const double chi_ff = get_chi_ff_nnionpart(modelgridindex) * pow(nu, -3) * nne * (1 - exp(-HOVERKB * nu / T_e));

  assert_testmodeonly(std::isfinite(chi_ff));

  return chi_ff;
}

// get bound-free opacity
template <bool USECELLHISTANDUPDATEPHIXSLIST>
auto calculate_chi_bf_gammacontr(const int modelgridindex, const double nu, Phixslist *phixslist) -> double {
  assert_always(!USECELLHISTANDUPDATEPHIXSLIST || phixslist != nullptr);

  double chi_bf_sum = 0.;
  if constexpr (USECELLHISTANDUPDATEPHIXSLIST) {
    if constexpr (USE_LUT_PHOTOION || USE_LUT_BFHEATING) {
      std::ranges::fill(phixslist->groundcont_gamma_contr, 0.);
    }
  }

  const auto T_e = grid::get_Te(modelgridindex);
  const auto nne = grid::get_nne(modelgridindex);
  const auto nnetot = grid::get_nnetot(modelgridindex);
  const auto &allcont_nu_edge = globals::allcont_nu_edge;

  // The phixslist is sorted by nu_edge in ascending order (longest to shortest wavelength)
  // If nu < allcont[i].nu_edge no absorption in any of the following continua
  // is possible, so set their kappas to zero
  // break the list into nu >= nu_edge and the remainder (nu < nu_edge)

  int i = 0;
  const int allcontend = static_cast<int>(std::ranges::upper_bound(allcont_nu_edge, nu) - allcont_nu_edge.cbegin());

  const int allcontbegin = std::lower_bound(allcont_nu_edge.data(), allcont_nu_edge.data() + allcontend, nu,
                                            [](const double nu_edge, const double nu_cmf) {
                                              return nu_edge * last_phixs_nuovernuedge < nu_cmf;
                                            }) -
                           allcont_nu_edge.data();

  assert_testmodeonly(allcontbegin >= 0);
  assert_testmodeonly(allcontend <= globals::nbfcontinua);
  assert_testmodeonly(allcontbegin <= allcontend);

  const auto *const allcont = globals::allcont;

  if constexpr (USECELLHISTANDUPDATEPHIXSLIST) {
    phixslist->allcontbegin = allcontbegin;
    phixslist->allcontend = allcontend;

    phixslist->bfestimend =
        static_cast<int>(std::ranges::upper_bound(globals::bfestim_nu_edge, nu) - globals::bfestim_nu_edge.cbegin());

    phixslist->bfestimbegin =
        std::lower_bound(
            globals::bfestim_nu_edge.data(), globals::bfestim_nu_edge.data() + phixslist->bfestimend, nu,
            [](const double nu_edge, const double nu_cmf) { return nu_edge * last_phixs_nuovernuedge < nu_cmf; }) -
        globals::bfestim_nu_edge.data();
  }

  for (i = allcontbegin; i < allcontend; i++) {
    const int element = allcont[i].element;
    const int ion = allcont[i].ion;
    const int level = allcont[i].level;
    const auto bfestimindex =
        (USECELLHISTANDUPDATEPHIXSLIST && DETAILED_BF_ESTIMATORS_ON) ? allcont[i].bfestimindex : -1;
    double sigma_contr = 0.;

    // The bf process happens only if the current cell contains
    // the involved atomic species
    const bool should_keep_this_cont = USECELLHISTANDUPDATEPHIXSLIST
                                           ? globals::cellcache[cellcacheslotid].ch_keep_this_cont[i]
                                           : keep_this_cont(element, ion, level, modelgridindex, nnetot);

    if (should_keep_this_cont) [[likely]] {
      const double nnlevel = USECELLHISTANDUPDATEPHIXSLIST ? globals::cellcache[cellcacheslotid].ch_allcont_nnlevel[i]
                                                           : calculate_levelpop(modelgridindex, element, ion, level);

      if (USECELLHISTANDUPDATEPHIXSLIST || nnlevel > 0) {
        const double nu_edge = allcont[i].nu_edge;
        const double sigma_bf = photoionization_crosssection_fromtable(allcont[i].photoion_xs, nu_edge, nu);

        double corrfactor = 1.;  // default to no subtraction of stimulated recombination
        if constexpr (!SEPARATE_STIMRECOMB) {
          double departure_ratio = globals::cellcache[cellcacheslotid].ch_allcont_departureratios[i];
          if (!USECELLHISTANDUPDATEPHIXSLIST || departure_ratio < 0) {
            const int upper = allcont[i].upperlevel;
            const double nnupperionlevel = USECELLHISTANDUPDATEPHIXSLIST
                                               ? get_levelpop(modelgridindex, element, ion + 1, upper)
                                               : calculate_levelpop(modelgridindex, element, ion + 1, upper);
            const double sf = calculate_sahafact(element, ion, level, upper, T_e, H * nu_edge);
            departure_ratio = nnupperionlevel / nnlevel * nne * sf;  // put that to phixslist
            if (USECELLHISTANDUPDATEPHIXSLIST) {
              globals::cellcache[cellcacheslotid].ch_allcont_departureratios[i] = departure_ratio;
            }
          }

          const double stimfactor = departure_ratio * exp(-HOVERKB * nu / T_e);
          corrfactor = std::max(0., 1 - stimfactor);  // photoionisation minus stimulated recombination
        }

        sigma_contr = sigma_bf * allcont[i].probability * corrfactor;

        if constexpr (USECELLHISTANDUPDATEPHIXSLIST) {
          if ((USE_LUT_PHOTOION || USE_LUT_BFHEATING) && level == 0 && allcont[i].phixstargetindex == 0) {
            phixslist->groundcont_gamma_contr[allcont[i].index_in_groundphixslist] = sigma_contr;
          }
        }

        chi_bf_sum += nnlevel * sigma_contr;
      }
    }
    if constexpr (USECELLHISTANDUPDATEPHIXSLIST) {
      phixslist->chi_bf_sum[i] = chi_bf_sum;
      if constexpr (DETAILED_BF_ESTIMATORS_ON) {
        if (bfestimindex >= 0) {
          phixslist->gamma_contr[bfestimindex] = sigma_contr;
        }
      }
    }
  }

  assert_always(std::isfinite(chi_bf_sum));

  return chi_bf_sum;
}

}  // anonymous namespace

void allocate_expansionopacities() {
  const auto npts_nonempty = grid::get_nonempty_npts_model();
  float *expansionopacities_data{};
  double *expansionopacity_planck_cumulative_data{};

#ifdef MPI_ON
  const auto [_, noderank_nonemptycellcount] =
      get_range_chunk(npts_nonempty, globals::node_nprocs, globals::rank_in_node);

  {
    MPI_Aint size = noderank_nonemptycellcount * expopac_nbins * static_cast<MPI_Aint>(sizeof(float));
    int disp_unit = sizeof(float);
    assert_always(MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, globals::mpi_comm_node,
                                          &expansionopacities_data, &win_expansionopacities) == MPI_SUCCESS);
    assert_always(MPI_Win_shared_query(win_expansionopacities, 0, &size, &disp_unit, &expansionopacities_data) ==
                  MPI_SUCCESS);
  }

  if constexpr (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY >= 0.) {
    MPI_Aint size = noderank_nonemptycellcount * expopac_nbins * static_cast<MPI_Aint>(sizeof(double));
    int disp_unit = sizeof(double);
    assert_always(MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, globals::mpi_comm_node,
                                          &expansionopacity_planck_cumulative_data,
                                          &win_expansionopacity_planck_cumulative) == MPI_SUCCESS);
    assert_always(MPI_Win_shared_query(win_expansionopacity_planck_cumulative, 0, &size, &disp_unit,
                                       &expansionopacity_planck_cumulative_data) == MPI_SUCCESS);
  }

#else
  expansionopacities_data = static_cast<float *>(malloc(npts_nonempty * expopac_nbins * sizeof(float)));
  if constexpr (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY >= 0.) {
    expansionopacity_planck_cumulative_data =
        static_cast<double *>(malloc(npts_nonempty * expopac_nbins * sizeof(double)));
  }
#endif
  expansionopacities = std::span(expansionopacities_data, npts_nonempty * expopac_nbins);
  expansionopacity_planck_cumulative =
      std::span(expansionopacity_planck_cumulative_data,
                expansionopacity_planck_cumulative_data == nullptr ? 0 : npts_nonempty * expopac_nbins);
}

// return a randomly chosen frequency with a distribution of Planck function times the expansion opacity
__host__ __device__ auto sample_planck_times_expansion_opacity(const int nonemptymgi) -> double {
  assert_testmodeonly(RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY >= 0.);

  const auto *kappa_planck_bins = &expansionopacity_planck_cumulative[nonemptymgi * expopac_nbins];

  const auto rnd_integral = rng_uniform() * kappa_planck_bins[expopac_nbins - 1];
  const auto *selected_partintegral =
      std::upper_bound(kappa_planck_bins, kappa_planck_bins + expopac_nbins, rnd_integral);
  const auto binindex = std::min(selected_partintegral - kappa_planck_bins, expopac_nbins - 1);
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < expopac_nbins);

  // use a linear interpolation for the frequency within the bin
  const auto bin_nu_lower = get_expopac_bin_nu_lower(binindex);
  const auto delta_nu = get_expopac_bin_nu_upper(binindex) - bin_nu_lower;
  const double nuoffset = rng_uniform() * delta_nu;
  const double nu = bin_nu_lower + nuoffset;
  return nu;
}

__host__ __device__ void do_rpkt(Packet &pkt, const double t2) {
  while (do_rpkt_step(pkt, t2)) {
    {
    }
  }
}

// make the packet an r-pkt and set further flags
__host__ __device__ void emit_rpkt(Packet &pkt) {
  pkt.type = TYPE_RPKT;
  pkt.last_cross = BOUNDARY_NONE;  // allow all further cell crossings

  // Need to assign a new direction. Assume isotropic emission in the cmf

  const auto dir_cmf = get_rand_isotropic_unitvec();

  // This direction is in the cmf - we want to convert it to the rest
  // frame - use aberration of angles. We want to convert from cmf to
  // rest so need -ve velocity.
  const auto vel_vec = get_velocity(pkt.pos, -1. * pkt.prop_time);
  // negative time since we want the backwards transformation here

  pkt.dir = angle_ab(dir_cmf, vel_vec);
  // printout("[debug] pkt.dir in RF: %g %g %g\n",pkt.dir[0],pkt.dir[1],pkt.dir[2]);

  // Finally we want to put in the rest frame energy and frequency. And record
  // that it's now a r-pkt.

  const double dopplerfactor = calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);
  pkt.nu_rf = pkt.nu_cmf / dopplerfactor;
  pkt.e_rf = pkt.e_cmf / dopplerfactor;

  // Reset polarization information
  pkt.stokes = {1., 0., 0.};

  pkt.pol_dir = cross_prod(pkt.dir, std::array<double, 3>{0., 0., 1.});

  if ((dot(pkt.pol_dir, pkt.pol_dir)) < 1.e-8) {
    pkt.pol_dir = cross_prod(pkt.dir, std::array<double, 3>{0., 1., 0.});
  }

  pkt.pol_dir = vec_norm(pkt.pol_dir);
}

void calculate_chi_rpkt_cont(const double nu_cmf, Rpkt_continuum_absorptioncoeffs &chi_rpkt_cont,
                             const int modelgridindex) {
  assert_testmodeonly(modelgridindex != grid::get_npts_model());
  assert_testmodeonly(grid::modelgrid[modelgridindex].thick != 1);
  if ((modelgridindex == chi_rpkt_cont.modelgridindex) && (globals::timestep == chi_rpkt_cont.timestep) &&
      (fabs((chi_rpkt_cont.nu / nu_cmf) - 1.0) < 1e-4)) {
    // calculated values are a match already
    return;
  }

  const auto nne = grid::get_nne(modelgridindex);

  double chi_escat = 0.;
  // free-free absorption
  const double chi_ff = calculate_chi_ffheating(modelgridindex, nu_cmf);
  double chi_bf = 0.;

  if (globals::opacity_case >= 4) {
    // First contribution: Thomson scattering on free electrons
    chi_escat = SIGMA_T * nne;

    // Third contribution: bound-free absorption
    chi_bf = chi_rpkt_cont.phixslist != nullptr
                 ? calculate_chi_bf_gammacontr<true>(modelgridindex, nu_cmf, chi_rpkt_cont.phixslist)
                 : calculate_chi_bf_gammacontr<false>(modelgridindex, nu_cmf, nullptr);

  } else {
    // in the other cases chi_grey is an mass absorption coefficient
    // therefore use the mass density
    // sigma = SIGMA_T * nne;

    chi_escat = 0.;
    // chi_ff = 0.9*sigma;
    // sigma *= 0.1;

    chi_bf = 0.;
  }

  chi_rpkt_cont.modelgridindex = modelgridindex;
  chi_rpkt_cont.timestep = globals::timestep;
  chi_rpkt_cont.nu = nu_cmf;
  chi_rpkt_cont.ffescat = chi_escat;
  chi_rpkt_cont.bf = chi_bf;
  chi_rpkt_cont.ffheat = chi_ff;
  chi_rpkt_cont.total = chi_rpkt_cont.ffescat + chi_rpkt_cont.bf + chi_rpkt_cont.ffheat;

  if (!std::isfinite(chi_rpkt_cont.total)) {
    printout("[fatal] calculate_chi_rpkt_cont: resulted in non-finite chi_rpkt_cont.total ... abort\n");
    printout("[fatal] es %g, ff %g, bf %g\n", chi_rpkt_cont.ffescat, chi_rpkt_cont.ffheat, chi_rpkt_cont.bf);
    printout("[fatal] nbfcontinua %d\n", globals::nbfcontinua);
    printout("[fatal] in cell %d with density %g\n", modelgridindex, grid::get_rho(modelgridindex));
    printout("[fatal] pkt.nu_cmf %g\n", nu_cmf);
    if (std::isfinite(chi_rpkt_cont.ffescat)) {
      chi_rpkt_cont.ffheat = 0.;
      chi_rpkt_cont.bf = 0.;
      chi_rpkt_cont.total = chi_rpkt_cont.ffescat;
    } else {
      std::abort();
    }
  }
}

#ifdef MPI_ON
void MPI_Bcast_binned_opacities(const ptrdiff_t nonemptymgi, const int root_node_id) {
  if constexpr (EXPANSIONOPACITIES_ON) {
    if (globals::rank_in_node == 0) {
      assert_always(nonemptymgi >= 0);
      MPI_Bcast(&expansionopacities[nonemptymgi * expopac_nbins], expopac_nbins, MPI_FLOAT, root_node_id,
                globals::mpi_comm_internode);

      if constexpr (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY >= 0.) {
        MPI_Bcast(&expansionopacity_planck_cumulative[nonemptymgi * expopac_nbins], expopac_nbins, MPI_DOUBLE,
                  root_node_id, globals::mpi_comm_internode);
      }
    }
  }
}
#endif

void calculate_expansion_opacities(const int modelgridindex) {
  const int nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  const auto rho = grid::get_rho(modelgridindex);

  const auto sys_time_start_calc = std::time(nullptr);
  const auto temperature = grid::get_TR(modelgridindex);

  printout("calculating expansion opacities for cell %d...", modelgridindex);

  const auto t_mid = globals::timesteps[globals::timestep].mid;

  // find the first line with nu below the upper limit of the first bin
  int lineindex = static_cast<int>(
      std::lower_bound(globals::linelist, globals::linelist + globals::nlines, get_expopac_bin_nu_upper(0),
                       [](const auto &line, const double nu_cmf) -> bool { return line.nu > nu_cmf; }) -
      globals::linelist);

  double kappa_planck_cumulative = 0.;

  for (ptrdiff_t binindex = 0; binindex < expopac_nbins; binindex++) {
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
    expansionopacities[(nonemptymgi * expopac_nbins) + binindex] = bin_kappa_bb;

    if constexpr (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY >= 0.) {
      // thread_local Rpkt_continuum_absorptioncoeffs chi_rpkt_cont {};
      // calculate_chi_rpkt_cont(nu_mid, chi_rpkt_cont, nullptr, modelgridindex);
      // const auto bin_kappa_cont = chi_rpkt_cont.total / rho;
      const auto bin_kappa_cont = calculate_chi_ffheating(modelgridindex, nu_mid) / rho;

      const auto planck_val = radfield::dbb(nu_mid, temperature, 1);
      const auto kappa_planck = (bin_kappa_bb + bin_kappa_cont) * planck_val;

      kappa_planck_cumulative += kappa_planck * delta_nu;

      expansionopacity_planck_cumulative[(nonemptymgi * expopac_nbins) + binindex] = kappa_planck_cumulative;
    }
  }
  printout("took %ld seconds\n", std::time(nullptr) - sys_time_start_calc);
}

// r-packet (UVOIR photon packet) propagation: bound-free and free-free continuum opacity,
// distances to the next line or continuum event, and the interaction events themselves
// (electron scattering, continuum absorption, and line absorption into a macro-atom).

#include "rpkt.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <limits>
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
#include "mpi_logging.h"
#include "packet.h"
#include "radfield.h"
#include "random.h"
#include "sn3d.h"
#include "stats.h"
#include "vectors.h"
#include "vpkt.h"

static_assert(!RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value() ||
                  RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.value() >= 0.,
              "RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY must be >= 0. if set");
static_assert(!RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value() ||
                  RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.value() <= 1.,
              "RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY must <= 1.0 if set");

static_assert(!RPKT_USE_EXPANSION_OPACITIES || !VPKT_ON, "VPKT cannot be used with r-packet expansion opacities");

namespace {
// kappa times Planck function for each bin of each non-empty cell
MPI_shared_array<double> expansionopacity_planck_cumulative{};

// get the comoving-frame frequency that the packet will have redshifted to at the abort distance (the cell
// boundary or the end of the timestep, whichever is nearer). The caller turns this into a frequency change
// per unit distance for the linear interpolation along the path.
// The move is done in two halves to give bit-identical results to the two half-steps in do_rpkt_step().
auto get_nu_cmf_abort(const Vec3d& pos, const Vec3d& dir, const double prop_time, const double nu_rf,
                      const double abort_dist) -> double {
  const auto half_abort_dist = abort_dist / 2.;
  const auto abort_time = prop_time + (half_abort_dist / CLIGHT_PROP) + (half_abort_dist / CLIGHT_PROP);

  const Vec3d abort_pos{
      pos[0] + (dir[0] * half_abort_dist) + (dir[0] * half_abort_dist),
      pos[1] + (dir[1] * half_abort_dist) + (dir[1] * half_abort_dist),
      pos[2] + (dir[2] * half_abort_dist) + (dir[2] * half_abort_dist),
  };

  const double nu_cmf_abort = nu_rf * calculate_doppler_nucmf_on_nurf(abort_pos, dir, abort_time);

  return nu_cmf_abort;
}

template <bool USECELLCACHE>
[[nodiscard]] auto get_tau_sobolev(const int nonemptymgi, const int lineindex, const double t_current) -> double {
  const int uniquelevelindex_lower = globals::linelist.uniquelevelindex_lower[lineindex];
  const int uniquelevelindex_upper = globals::linelist.uniquelevelindex_upper[lineindex];

  double n_l{NAN};
  double n_u{NAN};
  if constexpr (USECELLCACHE) {
    n_l = get_cellcache_levelpop(nonemptymgi, uniquelevelindex_lower);
    n_u = get_cellcache_levelpop(nonemptymgi, uniquelevelindex_upper);
  } else {
    const auto element = globals::linelist.elementindex[lineindex];
    const auto ion = globals::linelist.ionindex[lineindex];
    const auto ionuniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);
    const auto lower_uniquelevelindex = globals::linelist.uniquelevelindex_lower[lineindex];
    const auto upper_uniquelevelindex = globals::linelist.uniquelevelindex_upper[lineindex];
    const int lower = lower_uniquelevelindex - ionuniquelevelindexstart;
    const int upper = upper_uniquelevelindex - ionuniquelevelindexstart;
    n_l = calculate_levelpop(nonemptymgi, element, ion, lower);
    n_u = calculate_levelpop(nonemptymgi, element, ion, upper);
  }

  const double B_ul = globals::linelist.B_ul[lineindex];
  const double B_lu = globals::linelist.B_lu[lineindex];

  return std::max(((B_lu * n_l) - (B_ul * n_u)) * HCLIGHTOVERFOURPI * t_current, 0.);
}

// find any line or continuum interaction occurring before frequency decreases to nu_cmf_abort at distance abort_dist
// returns tuple of (distance to event, next transition index for pkt.next_trans, bool for whether line event)
// the next transition index is lineindex + 1 for a line event, may remain the current next_trans if no event occurs,
// and is globals::nlines + 1 for a continuum event
auto get_possible_event(const int nonemptymgi, const Packet& pkt, const ContinuumOpacity& chi_rpkt_cont,
                        MacroAtomState& mastate,
                        const double tau_rnd,  // random optical depth until which the packet travels
                        const double abort_dist,  // maximal travel distance before packet leaves cell or time step ends
                        const double nu_cmf_abort, const double dnu_on_dl, const double doppler,
                        const globals::TransitionLines& linelist) -> std::tuple<double, int, bool> {
  auto pos = pkt.pos;
  auto nu_cmf = pkt.nu_cmf;
  auto e_cmf = pkt.e_cmf;
  auto prop_time = pkt.prop_time;
  auto next_trans = pkt.next_trans;

  const double chi_cont = chi_rpkt_cont.total() * doppler;
  double tau = 0.;  // optical depth along path
  double dist = 0.;  // position on path
  while (true) {
    // step to the next line the packet will redshift onto, accumulating the continuum optical depth over the
    // distance to it. closest_transition() returns a negative index once no line remains at or below nu_cmf, or
    // when the packet has already been tagged as having no further line interactions.
    const int lineindex = closest_transition(nu_cmf, next_trans, linelist.nu);

    if (lineindex < 0) [[unlikely]] {
      // no line interaction possible - check whether continuum process occurs in cell

      const double tau_cont = chi_cont * (abort_dist - dist);

      if (tau_rnd - tau > tau_cont) {
        // no continuum event before abort_dist
        return {std::numeric_limits<double>::max(), next_trans, false};
      }

      // continuum process occurs at edist
      return {dist + ((tau_rnd - tau) / chi_cont), globals::nlines + 1, false};
    }

    // a line interaction is possible, i.e. the packet redshifts onto this line before nu_cmf_abort

    const double nu_trans = linelist.nu[lineindex];

    // advance past this line so that any further event this step is found at a lower frequency. Without this the
    // packet could scatter repeatedly in the same line when rounding leaves nu_cmf marginally above nu_trans.
    next_trans = lineindex + 1;

    const double ldist = get_linedistance(prop_time, nu_cmf, nu_trans, dnu_on_dl);

    const double tau_cont = chi_cont * ldist;

    if (tau_rnd - tau > tau_cont) {
      // got past the continuum optical depth so propagate to the line, and check interaction

      if (nu_trans < nu_cmf_abort) [[unlikely]] {
        // back up one line, because we didn't reach it before the boundary/timelimit

        return {std::numeric_limits<double>::max(), next_trans - 1, false};
      }

      const double tau_line = get_tau_sobolev<true>(nonemptymgi, lineindex, prop_time);

      if ((tau_rnd - tau) <= (tau_cont + tau_line)) {
        // bound-bound process occurs
        const auto element = linelist.elementindex[lineindex];
        const auto ion = linelist.ionindex[lineindex];
        const auto ionuniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);
        const auto upper = globals::linelist.uniquelevelindex_upper[lineindex] - ionuniquelevelindexstart;

        mastate = {.element = element, .ion = ion, .level = upper, .activatingline = lineindex};

        if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
          move_pkt_withtime(pos, pkt.dir, prop_time, pkt.nu_rf, nu_cmf, pkt.e_rf, e_cmf, ldist);
          radfield::update_lineestimator(nonemptymgi, lineindex, prop_time * CLIGHT * e_cmf / nu_cmf);
        }

        // the line and its parameters were already selected by closest_transition!

        return {dist + ldist, next_trans, true};
      }

      // total optical depth still below tau_rnd: propagate to the line and continue

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
        nu_cmf = pkt.nu_cmf + (dnu_on_dl * dist);  // should equal nu_trans;
        assert_testmodeonly(nu_cmf <= pkt.nu_cmf);
        if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
          // keep e_cmf consistent with the linearly-approximated nu_cmf for the line estimator below
          // (e_cmf / nu_cmf = e_rf / nu_rf is invariant along a free path)
          e_cmf = nu_cmf * pkt.e_rf / pkt.nu_rf;
        }
      }

      if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
        radfield::update_lineestimator(nonemptymgi, lineindex, prop_time * CLIGHT * e_cmf / nu_cmf);
      }

    } else {
      // continuum process occurs before reaching the line

      return {dist + ((tau_rnd - tau) / chi_cont), next_trans - 1, false};
    }
  }

  // should have already returned somewhere!
  assert_always(false);
}

auto get_possible_event_expansion_opacity(const int nonemptymgi, Packet& pkt, const ContinuumOpacity& chi_rpkt_cont,
                                          MacroAtomState& mastate, const double tau_rnd, const double nu_cmf_abort,
                                          const double dnu_on_dl, const double doppler) -> std::tuple<double, bool> {
  auto pos = pkt.pos;
  const auto nu_rf = pkt.nu_rf;
  auto nu_cmf = pkt.nu_cmf;
  const auto e_rf = pkt.e_rf;
  auto e_cmf = pkt.e_cmf;
  auto prop_time = pkt.prop_time;

  // with thermalisation or pure scattering, we don't keep track of line interactions

  assert_always(get_cellcache(nonemptymgi).nonemptymgi == nonemptymgi);
  double dist = 0.;
  double tau = 0.;
  // -1 means below the wavelength grid, in which case there is continuum opacity but no expansion opacity
  const auto binindex_start =
      std::max(get_linearbinindex(1e8 * CLIGHT / nu_cmf, expopac_lambdamin, expopac_deltalambda), -1Z);

  for (auto binindex = binindex_start; binindex < expopac_nbins; binindex++) {
    // binindex could be -1, in which case we have only the continuum opacity and no expansion opacity
    const auto next_bin_edge_nu = (binindex < 0) ? get_expopac_bin_nu_upper(0) : get_expopac_bin_nu_lower(binindex);
    const auto binedgedist = get_linedistance(prop_time, nu_cmf, next_bin_edge_nu, dnu_on_dl);

    const double chi_cont = chi_rpkt_cont.total() * doppler;
    double chi_bb_expansionopac = 0.;
    if (binindex >= 0) {
      // opacity in units of [cm^2/g]
      const auto kappa = expansionopacities[(nonemptymgi * expopac_nbins) + binindex];
      // absorption coefficient in units of [1/cm]
      chi_bb_expansionopac = kappa * grid::get_rho(nonemptymgi);
    }

    const double chi_tot = chi_cont + chi_bb_expansionopac;

    if (chi_tot * binedgedist > tau_rnd - tau) {
      // interaction occurs
      if constexpr (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value()) {
        const auto edist = std::max(dist + ((tau_rnd - tau) / chi_tot), 0.);
        // strict comparison, so that a draw of exactly zero cannot select the bound-bound channel
        // when its opacity is zero (rng_uniform() rejects one but can return zero). Otherwise a
        // wavelength bin with no line opacity could take the thermalisation branch, and in a cell
        // with no line opacity at all the Planck sampling that follows has nothing to sample from.
        const bool event_is_boundbound = rng_uniform(get_rngstate(pkt)) < chi_bb_expansionopac / chi_tot;
        return {edist, event_is_boundbound};
      }

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
      auto next_trans = -1;
      std::tie(edist_after_bin, next_trans, event_is_boundbound) =
          get_possible_event(nonemptymgi, pkt_bin_start, chi_rpkt_cont, mastate, tau_rnd - tau,
                             std::numeric_limits<double>::max(), 0., dnu_on_dl, doppler, globals::linelist);
      dist = dist + edist_after_bin;

      return {dist, event_is_boundbound};
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
      nu_cmf = pkt.nu_cmf + (dnu_on_dl * dist);  // should equal nu_trans;
      assert_testmodeonly(nu_cmf <= pkt.nu_cmf);
      if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
        // keep e_cmf consistent with the linearly-approximated nu_cmf, since it seeds the packet copy
        // used for the line-by-line retrace (whose line estimator updates divide e_cmf by nu_cmf)
        e_cmf = nu_cmf * e_rf / nu_rf;
      }
    }

    if (nu_cmf <= nu_cmf_abort) {
      // hit edge of cell or timestep limit
      return {std::numeric_limits<double>::max(), false};
    }
  }

  // no more expansion-opacity bins (the packet is redshifting past the red end of the wavelength
  // grid), but the continuum processes (electron scattering, bound-free, free-free) still provide
  // opacity, so sample the remaining optical depth against them (cf. the vpkt tracer, which adds
  // the continuum opacity independently of the expansion-opacity bin range)
  const double chi_cont = chi_rpkt_cont.total() * doppler;
  if (chi_cont > 0.) {
    return {dist + ((tau_rnd - tau) / chi_cont), false};
  }
  return {std::numeric_limits<double>::max(), false};
}

void electron_scatter_rpkt(Packet& pkt) {
  // now make the packet a r-pkt and set further flags
  pkt.type = TYPE_RPKT;

  const auto vel_vec = get_velocity(pkt.pos, pkt.prop_time);

  // Transform Stokes Parameters from the RF to the CMF

  const auto [old_dir_cmf, q_i_cmf, u_i_cmf] = (POL_ON ? frame_transform(pkt.dir, pkt.stokes_q, pkt.stokes_u, vel_vec)
                                                       : std::make_tuple(angle_ab(pkt.dir, vel_vec), 0., 0.));

  // Outcoming direction. Compute the new cmf direction from the old direction and the scattering angles (see Kalos &
  // Whitlock 2008)
  double M = 0.;
  double phisc = 0.;

  if constexpr (DIPOLE) {
    // Assume dipole function: sample the scattering direction cosine M and azimuth angle phisc by
    // rejection (see Code & Whitney 1995). p is the phase function value for the trial angles and
    // x is a uniform draw from [0, 2] (an upper bound on p); the trial is accepted when x <= p.
    double p = 0.;
    double x = 1.;
    while (x > p) {
      M = (2. * rng_uniform_pos(get_rngstate(pkt))) - 1.;
      const double musquared = pow2(M);
      phisc = 2 * PI * rng_uniform(get_rngstate(pkt));

      // NB: the rotational matrix R here is chosen in the clockwise direction ("+").
      // In Bulla+2015 equation (10) and (12) refer to the specific case shown in Fig.2 where the angle i2
      // is measured in the counter-clockwise direction. Therefore we use the clockwise rotation matrix but
      // with -i1. Here, instead, we calculate the angle in the clockwise direction from 0 to 2PI.
      // For instance, the i1 angle in Fig.2 of Bulla+2015 corresponds to 2PI-i1 here.
      // NB2: the i1 and i2 angles computed in the code (before and after scattering) are instead as in Bulla+2015
      p = (musquared + 1) + ((musquared - 1) * ((cos(2 * phisc) * q_i_cmf) + (sin(2 * phisc) * u_i_cmf)));

      // generate a number between 0 and the maximum of the previous function (2)
      x = 2. * rng_uniform(get_rngstate(pkt));
    }
  } else {
    // Assume isotropic scattering
    M = (2. * rng_uniform(get_rngstate(pkt))) - 1.;
    phisc = 2 * PI * rng_uniform(get_rngstate(pkt));
  }

  Vec3d new_dir_cmf{};

  const double cos_tsc = M;  // M is cos(tsc) by construction
  const double sin_tsc = std::sqrt(1. - pow2(M));

  if (fabs(old_dir_cmf[2]) < 0.99999) {
    const double sin_polar = std::sqrt(1. - pow2(old_dir_cmf[2]));
    const double common_factor = sin_tsc / sin_polar;
    const double cos_phisc = cos(phisc);
    const double sin_phisc = sin(phisc);
    new_dir_cmf = {
        (common_factor * ((old_dir_cmf[1] * sin_phisc) - (old_dir_cmf[0] * old_dir_cmf[2] * cos_phisc))) +
            (old_dir_cmf[0] * cos_tsc),
        (common_factor * ((-old_dir_cmf[0] * sin_phisc) - (old_dir_cmf[1] * old_dir_cmf[2] * cos_phisc))) +
            (old_dir_cmf[1] * cos_tsc),
        (sin_tsc * cos_phisc * sin_polar) + (old_dir_cmf[2] * cos_tsc),
    };
  } else {
    new_dir_cmf = {sin_tsc * cos(phisc), sin_tsc * sin(phisc), (old_dir_cmf[2] > 0) ? cos_tsc : -cos_tsc};
  }

  if constexpr (POL_ON) {
    // Need to rotate Stokes Parameters in the scattering plane
    std::tie(pkt.dir, pkt.stokes_q, pkt.stokes_u, std::ignore) =
        scatter_polarisation_to_rf(old_dir_cmf, new_dir_cmf, q_i_cmf, u_i_cmf, vel_vec);
  } else {
    pkt.dir = angle_ab(new_dir_cmf, Vec3d{-vel_vec[0], -vel_vec[1], -vel_vec[2]});
  }

  // Check unit vector
  assert_testmodeonly(fabs(vec_len(pkt.dir) - 1.) < 1.e-6);

  // Finally we want to put in the rest frame energy and frequency. And record that it's now a r-pkt.

  set_pkt_restframe_from_cmf(pkt);
}

void rpkt_event_continuum(Packet& pkt, const ContinuumOpacity& chi_rpkt_cont) {
  const double nu = pkt.nu_cmf;

  const double dopplerfactor = calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);
  const double chi_cont = chi_rpkt_cont.total() * dopplerfactor;
  const double chi_escatter = chi_rpkt_cont.chi_escatter * dopplerfactor;
  const double chi_ff = chi_rpkt_cont.chi_freefree_heat * dopplerfactor;
  const double chi_bf = chi_rpkt_cont.chi_boundfree * dopplerfactor;

  // continuum process happens. select due to its probabilities sigma/chi_cont, chi_ff/chi_cont, chi_bf/chi_cont

  const auto chi_rnd = rng_uniform(get_rngstate(pkt)) * chi_cont;
  if (chi_rnd < chi_escatter) {
    // electron scattering occurs
    // in this case the packet stays a R_PKT of same nu_cmf as before (coherent scattering) but with different direction
    pkt.nscatterings++;
    stats::increment(stats::Counter::ELECTRON_SCATTERINGS);

    // generate a virtual packet
    if constexpr (VPKT_ON) {
      vpkt::trace_vpkts(pkt, TYPE_RPKT);
    }

    electron_scatter_rpkt(pkt);

    // Electron scattering does not modify the last emission flag but it updates the last emission position
    pkt.em_pos = pkt.pos;
    pkt.em_time = static_cast<float>(pkt.prop_time);

  } else if (chi_rnd < chi_escatter + chi_ff) {
    // ff: transform to k-pkt
    stats::increment(stats::Counter::K_STAT_FROM_FF);
    pkt.type = TYPE_KPKT;
    pkt.absorptiontype = ABSTYPE_FREEFREE;
  } else if (chi_rnd < chi_escatter + chi_ff + chi_bf) {
    // bf: transform to k-pkt or activate macroatom corresponding to probabilities

    const auto& phixslist = chi_rpkt_cont.phixslist;

    pkt.absorptiontype = ABSTYPE_BOUNDFREE;

    const double chi_bf_inrest = chi_rpkt_cont.chi_boundfree;
    assert_testmodeonly(phixslist.chi_bf_sum[phixslist.allcontend - 1] == chi_bf_inrest);

    // Determine in which continuum the bf-absorption occurs
    const double chi_bf_rand = rng_uniform(get_rngstate(pkt)) * chi_bf_inrest;

    // first chi_bf_sum[i] such that chi_bf_sum[i] > chi_bf_rand (or the last one if chi_bf_rand is larger than all
    // chi_bf_sum) gives the index of the continuum
    const auto allcontindex =
        std::ranges::upper_bound(
            phixslist.chi_bf_sum.subspan(phixslist.allcontbegin, phixslist.allcontend - phixslist.allcontbegin - 1),
            chi_bf_rand) -
        phixslist.chi_bf_sum.begin();

    const double nu_edge = globals::allcont.nu_edge[allcontindex];
    const int element = globals::allcont.element[allcontindex];
    const int ion = globals::allcont.ion[allcontindex];
    const int level = globals::allcont.level[allcontindex];
    const int phixstargetindex = globals::allcont.phixstargetindex[allcontindex];

    // decide whether we go to ionisation energy or to the thermal pool
    if (rng_uniform(get_rngstate(pkt)) < nu_edge / nu) {
      stats::increment(stats::Counter::MA_STAT_ACTIVATION_BF);

      do_macroatom(pkt, {
                            .element = element,
                            .ion = ion + 1,
                            .level = get_phixsupperlevel(element, ion, level, phixstargetindex),
                            .activatingline = -99,
                        });
    } else {
      // transform to k-pkt
      stats::increment(stats::Counter::K_STAT_FROM_BF);
      pkt.type = TYPE_KPKT;
    }
  } else {
    assert_always(false);
  }
}

// Update the volume estimators J and nuJ
// This is done in another routine than move, as we sometimes move dummy
// packets which do not contribute to the radiation field.
void update_estimators(const double e_cmf, const double nu_cmf, const double distance, const int nonemptymgi,
                       const ContinuumOpacity& chi_rpkt_cont, const bool thickcell) {
  // Update only non-empty cells
  assert_testmodeonly(nonemptymgi >= 0);
  const double distance_e_cmf = distance * e_cmf;

  radfield::update_estimators(nonemptymgi, distance_e_cmf, nu_cmf, chi_rpkt_cont.phixslist, thickcell);

  if (thickcell) {
    // chi_rpkt_cont and phixslist are not known for thick cells
    return;
  }

  // ffheatingestimator does not depend on ion and element, so an array with gridsize is enough.
  atomicadd(globals::ffheatingestimator[nonemptymgi], distance_e_cmf * chi_rpkt_cont.chi_freefree_heat);

  if constexpr (USE_LUT_PHOTOION || USE_ION_BFHEATING_ESTIMATORS) {
    for (int i = 0; i < globals::nbfcontinua_ground; i++) {
      const double nu_edge = globals::groundcont_nu_edge[i];
      if (nu_cmf <= nu_edge) {
        // because groundcont is sorted by nu_edge ascending, nu_cmf < nu_edge for all remaining items
        return;
      }
      const ptrdiff_t ionestimindex = (static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground) + i;

      if constexpr (USE_LUT_PHOTOION) {
        atomicadd(globals::gammaestimator[ionestimindex],
                  chi_rpkt_cont.phixslist.groundcont_gamma_contr[i] * (distance_e_cmf / nu_cmf));
      }

      if constexpr (USE_ION_BFHEATING_ESTIMATORS) {
        atomicadd(globals::bfheatingestimator[ionestimindex],
                  chi_rpkt_cont.phixslist.groundcont_gamma_contr[i] * distance_e_cmf * (1. - (nu_edge / nu_cmf)));
      }
    }
  }
}

// Update an r-packet and return true if no mgi change (or it goes into an empty cell) and no pkttype change and not
// reached end of timestep, otherwise false
auto do_rpkt_step(Packet& pkt, const double t2, ContinuumOpacity& chi_rpkt_cont) -> bool {
  const auto nonemptymgi = grid::get_propcell_nonemptymgi(pkt.cellindex);

  MacroAtomState pktmastate{};

  // draw random optical depth to next physical event
  const double tau_rnd = -std::log(static_cast<double>(rng_uniform_pos(get_rngstate(pkt))));

  // Finding the distance to the crossing of the grid cell boundaries.
  // boundarydist is the boundary distance to the next grid cell next_cellindex
  const auto [boundarydist, next_cellindex] = grid::boundary_distance(pkt.dir, pkt.pos, pkt.prop_time, pkt.cellindex);

  if (boundarydist == 0) {
    grid::change_cell_or_escape(pkt, next_cellindex);
    const int new_nonemptymgi = grid::get_propcell_nonemptymgi(pkt.cellindex);

    return (pkt.type == TYPE_RPKT && (new_nonemptymgi < 0 || new_nonemptymgi == nonemptymgi));
  }

  // maximum travel distance before t2 (end of timestep) is tdist
  const double tdist = (t2 - pkt.prop_time) * CLIGHT_PROP;

  assert_always(tdist >= 0);

  const double abort_dist = std::min(tdist, boundarydist);

  // Get distance to the next physical event (continuum or bound-bound)
  double edist = -1;
  bool event_is_boundbound = true;
  const bool thickcell = (nonemptymgi >= 0) && (grid::thick_allcells[nonemptymgi] == 1);
  if (nonemptymgi < 0) {
    // for empty cells no physical event occurs. The packets just propagate.
    edist = std::numeric_limits<double>::max();
    pkt.next_trans = -1;  // skip over lines and search for line list position on the next non-empty cell
  } else if (thickcell) {
    // In the case of optically thick cells, we treat the packets in grey approximation to speed up the calculation

    const double chi_grey = grid::kappagrey_allcells[nonemptymgi] * grid::get_rho(nonemptymgi) *
                            calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);

    edist = tau_rnd / chi_grey;
    pkt.next_trans = -1;
  } else {
    calculate_chi_rpkt_cont<true>(pkt.nu_cmf, chi_rpkt_cont, nonemptymgi);

    // for USE_RELATIVISTIC_DOPPLER_SHIFT, we will use a linear approximation for
    // the frequency change from start to abort (cell boundary/timestep end)

    const auto nu_cmf_abort = get_nu_cmf_abort(pkt.pos, pkt.dir, pkt.prop_time, pkt.nu_rf, abort_dist);
    const auto dnu_on_dl = (nu_cmf_abort - pkt.nu_cmf) / abort_dist;
    const auto doppler = calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);

    if constexpr (RPKT_USE_EXPANSION_OPACITIES) {
      std::tie(edist, event_is_boundbound) = get_possible_event_expansion_opacity(
          nonemptymgi, pkt, chi_rpkt_cont, pktmastate, tau_rnd, nu_cmf_abort, dnu_on_dl, doppler);
    } else {
      std::tie(edist, pkt.next_trans, event_is_boundbound) =
          get_possible_event(nonemptymgi, pkt, chi_rpkt_cont, pktmastate, tau_rnd, abort_dist, nu_cmf_abort, dnu_on_dl,
                             doppler, globals::linelist);
    }
  }
  assert_always(edist >= 0);

  if ((edist < boundarydist) && (edist <= tdist)) [[likely]] {
    // bound-bound or continuum event occurs before cell crossing or end of timestep
    move_pkt_withtime(pkt, edist / 2.);
    update_estimators(pkt.e_cmf, pkt.nu_cmf, edist, nonemptymgi, chi_rpkt_cont, thickcell);
    move_pkt_withtime(pkt, edist / 2.);

    // The previously selected event occurs
    stats::increment(stats::Counter::INTERACTIONS);
    if (thickcell) {
      pkt.nscatterings++;
      stats::increment(stats::Counter::ELECTRON_SCATTERINGS);

      emit_rpkt(pkt);
      // Electron scattering does not modify the last emission flag but it updates the last emission position
    } else if (!event_is_boundbound) {
      rpkt_event_continuum(pkt, chi_rpkt_cont);
    } else if constexpr (!RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value()) {
      stats::increment(stats::Counter::MA_STAT_ACTIVATION_BB);

      pkt.absorptiontype = pktmastate.activatingline;
      pkt.absorptionfreq = pkt.nu_rf;

      do_macroatom(pkt, pktmastate);
    } else {
      // Probability based thermalisation (i.e. redistribution of the packet frequency) or scattering
      if (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.value() >= 1. ||
          rng_uniform(get_rngstate(pkt)) < RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.value()) {
        // Thermal redistribution of frequency

        pkt.absorptiontype = pktmastate.activatingline;
        pkt.absorptionfreq = pkt.nu_rf;
        pkt.nu_cmf = sample_planck_times_expansion_opacity(nonemptymgi, get_rngstate(pkt));
        pkt.next_trans = -1;
        // a thermal re-emission at a new frequency, so the packet no longer traces back to the previous emission
        pkt.emissiontype = EMTYPE_NOTSET;
        pkt.trueemissiontype = EMTYPE_NOTSET;
        pkt.trueem_pos = {NAN, NAN, NAN};
        pkt.trueem_time = -1.;

        // re-emit rather than scatter, so that this event is not counted as an electron scattering
        pkt.nscatterings = 0;
      } else {
        // pure scattering, so the packet keeps its co-moving frequency and direction is changed
        pkt.nscatterings++;
        stats::increment(stats::Counter::ELECTRON_SCATTERINGS);
      }
      emit_rpkt(pkt);
    }

    return (pkt.type == TYPE_RPKT);
  }

  if ((boundarydist <= tdist) && (boundarydist <= edist)) {
    // cell crossing event occurs before interaction or end of timestep
    move_pkt_withtime(pkt, boundarydist / 2.);
    if (nonemptymgi >= 0) {
      update_estimators(pkt.e_cmf, pkt.nu_cmf, boundarydist, nonemptymgi, chi_rpkt_cont, thickcell);
    }
    move_pkt_withtime(pkt, boundarydist / 2.);

    if (next_cellindex != pkt.cellindex) {
      grid::change_cell_or_escape(pkt, next_cellindex);
      if (next_cellindex < 0) {
        // we left the grid, so we can stop tracking this packet
        return false;
      }
      const auto new_nonemptymgi = grid::get_propcell_nonemptymgi(pkt.cellindex);
      // if the new cell is empty or the same as the previous one, keep going, otherwise we'll need to change the cell
      // cache
      return ((new_nonemptymgi < 0) || (new_nonemptymgi == nonemptymgi));
    }
    return true;  // if next_cellindex == pkt.cellindex, we reached the maximum path length and are not changing cell
  }

  if ((tdist < boundarydist) && (tdist <= edist)) [[unlikely]] {
    // reaches end of timestep before cell boundary or interaction
    move_pkt_withtime(pkt, tdist / 2.);
    if (nonemptymgi >= 0) {
      update_estimators(pkt.e_cmf, pkt.nu_cmf, tdist, nonemptymgi, chi_rpkt_cont, thickcell);
    }
    move_pkt_withtime(pkt, tdist / 2.);
    pkt.prop_time = t2;

    return false;
  }

  assert_always(false);
  return false;
}

// calculate the free-free absorption (to kpkt heating) coefficient [cm^-1]
// = kappa(free-free) * nne
auto calculate_chi_ffheating(const int nonemptymgi, const double nu, const bool use_cellcache) -> double {
  const auto clumpednne = grid::get_nne(nonemptymgi) * grid::get_clumpfactor(nonemptymgi);
  const auto T_e = grid::Te_allcells[nonemptymgi];
  assert_testmodeonly(!use_cellcache || get_cellcache(nonemptymgi).nonemptymgi == nonemptymgi);
  const auto chi_ff_nnionpart =
      use_cellcache ? get_cellcache(nonemptymgi).chi_ff_nnionpart[0] : calculate_chi_ffheat_nnionpart(nonemptymgi);

  const double chi_ff = chi_ff_nnionpart / pow3(nu) * clumpednne * (1 - exp(-HOVERKB * nu / T_e));

  assert_testmodeonly(std::isfinite(chi_ff));
  assert_testmodeonly(chi_ff >= 0);

  return chi_ff;
}

// sum the bound-free opacity at frequency nu over all photoionisation continua (using the
// binary-searched frequency window of contributing edges). When USECELLHISTANDUPDATEPHIXSLIST is
// true, also record each continuum's contribution in the phixslist and cache the running sum for reuse within the cell.
template <bool USECELLHISTANDUPDATEPHIXSLIST>
auto calculate_chi_bf_gammacontr(const int nonemptymgi, const double nu, Phixslist& phixslist) -> double {
  double chi_bf_sum = 0.;
  if constexpr (USECELLHISTANDUPDATEPHIXSLIST) {
    assert_testmodeonly(std::ssize(phixslist.chi_bf_sum) == globals::nbfcontinua);
    if constexpr (USE_LUT_PHOTOION || USE_ION_BFHEATING_ESTIMATORS) {
      std::ranges::fill(phixslist.groundcont_gamma_contr, 0.);
    }
  }

  const auto T_e = grid::Te_allcells[nonemptymgi];
  const auto clumpednne = grid::get_nne(nonemptymgi) * grid::get_clumpfactor(nonemptymgi);
  const auto nnetot = grid::get_nnetot(nonemptymgi);
  const auto& allcont_nu_edge = globals::allcont.nu_edge;

  // The phixslist is sorted by nu_edge in ascending order, so if nu < allcont[i].nu_edge then no absorption in any of
  // the remaining continua is possible. so set their kappas to zero and break
  const int allcontend = static_cast<int>(std::ranges::upper_bound(allcont_nu_edge, nu) - allcont_nu_edge.begin());

  // require that nu <= nu_edge * last_phixs_nuovernuedge, which can exclude some low-nu edges
  const auto allcontbegin =
      static_cast<int>(std::ranges::lower_bound(allcont_nu_edge.first(allcontend), nu / last_phixs_nuovernuedge) -
                       allcont_nu_edge.begin());

  assert_testmodeonly(allcontbegin >= 0);
  assert_testmodeonly(allcontend <= globals::nbfcontinua);
  assert_testmodeonly(allcontbegin <= allcontend);

  if constexpr (USECELLHISTANDUPDATEPHIXSLIST) {
    phixslist.allcontbegin = allcontbegin;
    phixslist.allcontend = allcontend;

    phixslist.bfestimend =
        static_cast<int>(std::ranges::upper_bound(globals::bfestim_nu_edge, nu) - globals::bfestim_nu_edge.begin());

    phixslist.bfestimbegin = static_cast<int>(
        std::ranges::lower_bound(globals::bfestim_nu_edge.first(phixslist.bfestimend), nu / last_phixs_nuovernuedge) -
        globals::bfestim_nu_edge.begin());
  }

  // const ref these so that the compiler knows they don't change in the loop (and shortens the names)
  const auto& allcont_element = globals::allcont.element;
  const auto& allcont_ion = globals::allcont.ion;
  const auto& allcont_level = globals::allcont.level;
  const auto& allcont_bfestimindex = globals::allcont.bfestimindex;
  const auto& allcont_upperlevel = globals::allcont.upperlevel;
  const auto& allcont_uniquelevelindex = globals::allcont.uniquelevelindex;
  const auto& allcont_index_in_groundphixslist = globals::allcont.index_in_groundphixslist;
  const auto& allcont_probability = globals::allcont.probability;
  const auto& allcont_phixstargetindex = globals::allcont.phixstargetindex;

  for (int i = allcontbegin; i < allcontend; i++) {
    const int element = allcont_element[i];
    const int ion = allcont_ion[i];
    const int level = allcont_level[i];
    double sigma_contr = 0.;

    // The bf process happens only if the current cell contains the involved atomic species
    const bool should_keep_this_cont = USECELLHISTANDUPDATEPHIXSLIST
                                           ? get_cellcache(nonemptymgi).allcont_keep[i]
                                           : keep_this_cont(element, ion, level, nonemptymgi, nnetot);

    if (should_keep_this_cont) [[likely]] {
      const double nnlevel = USECELLHISTANDUPDATEPHIXSLIST ? get_cellcache(nonemptymgi).allcont_nnlevel[i]
                                                           : calculate_levelpop(nonemptymgi, element, ion, level);

      if (USECELLHISTANDUPDATEPHIXSLIST || nnlevel > 0) {
        const double nu_edge = allcont_nu_edge[i];
        const double sigma_bf =
            photoionisation_crosssection_fromtable(get_phixs_table(allcont_uniquelevelindex[i]), nu_edge, nu);

        double corrfactor = 1.;  // default to no subtraction of stimulated recombination
        double modified_departure_ratio = get_cellcache(nonemptymgi).allcont_modified_departureratios[i];

        if (!USECELLHISTANDUPDATEPHIXSLIST || modified_departure_ratio < 0) {
          const int upper = allcont_upperlevel[i];
          const double nnupperionlevel = USECELLHISTANDUPDATEPHIXSLIST
                                             ? get_cellcache_levelpop(nonemptymgi, element, ion + 1, upper)
                                             : calculate_levelpop(nonemptymgi, element, ion + 1, upper);
          const double modified_sahafact =
              SAHACONST * stat_weight(element, ion, level) / stat_weight(element, ion + 1, upper) * std::pow(T_e, -1.5);
          modified_departure_ratio =
              nnupperionlevel / nnlevel * clumpednne * modified_sahafact;  // put that to phixslist
          if (USECELLHISTANDUPDATEPHIXSLIST) {
            get_cellcache(nonemptymgi).allcont_modified_departureratios[i] = modified_departure_ratio;
          }
        }

        const double stimfactor = modified_departure_ratio * exp(-HOVERKB * (nu - nu_edge) / T_e);
        corrfactor = std::max(0., 1 - stimfactor);  // photoionisation minus stimulated recombination

        sigma_contr = sigma_bf * allcont_probability[i] * corrfactor;

        if constexpr (USECELLHISTANDUPDATEPHIXSLIST) {
          if ((USE_LUT_PHOTOION || USE_ION_BFHEATING_ESTIMATORS) && level == 0 && allcont_phixstargetindex[i] == 0) {
            phixslist.groundcont_gamma_contr[allcont_index_in_groundphixslist[i]] = sigma_contr;
          }
        }

        chi_bf_sum += nnlevel * sigma_contr;
      }
    }
    if constexpr (USECELLHISTANDUPDATEPHIXSLIST) {
      phixslist.chi_bf_sum[i] = chi_bf_sum;
      if constexpr (DETAILED_BF_ESTIMATORS_ON) {
        if (allcont_bfestimindex[i] >= 0) {
          phixslist.gamma_contr[allcont_bfestimindex[i]] = sigma_contr;
        }
      }
    }
  }

  assert_always(std::isfinite(chi_bf_sum));

  return chi_bf_sum;
}

}  // anonymous namespace

auto calculate_chi_ffheat_nnionpart(const int nonemptymgi) -> double {
  const double g_ff = 1;
  double chi_ff_nnionpart = 0.;
  const int nelements = get_nelements();
  for (int element = 0; element < nelements; element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const double nnion = get_nnion(nonemptymgi, element, ion);
      const int ioncharge = get_ionstage(element, ion) - 1;
      chi_ff_nnionpart += pow2(ioncharge) * g_ff * nnion;
    }
  }
  const auto T_e = grid::Te_allcells[nonemptymgi];

  return chi_ff_nnionpart * 3.69255e8 / sqrt(T_e);
}

void allocate_expansionopacities() {
  const auto nonempty_npts_model = grid::get_nonempty_npts_model();

  assert_always(expansionopacities.empty());
  if constexpr (RPKT_USE_EXPANSION_OPACITIES || VPKT_USE_EXPANSION_OPACITIES) {
    expansionopacities = MPI_shared_array<float>(nonempty_npts_model * expopac_nbins);
  }

  assert_always(expansionopacity_planck_cumulative.empty());
  if constexpr (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value()) {
    expansionopacity_planck_cumulative = MPI_shared_array<double>(nonempty_npts_model * expopac_nbins);
  }
}

// return a randomly chosen frequency with a distribution of Planck function times the expansion opacity
DEVICE_FUNC auto sample_planck_times_expansion_opacity(const int nonemptymgi, rngstate_type& rngstate) -> double {
  assert_testmodeonly(RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value());

  const std::span<const double> kappa_planck_bins =
      expansionopacity_planck_cumulative.subspan(nonemptymgi * expopac_nbins, expopac_nbins);
  assert_always(kappa_planck_bins.back() > 0);
  const auto rnd_integral = rng_uniform(rngstate) * kappa_planck_bins[expopac_nbins - 1];
  const auto binindex = std::min(index_upperbound(kappa_planck_bins, rnd_integral), expopac_nbins - 1);
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < expopac_nbins);

  // use a linear interpolation for the frequency within the bin
  const auto bin_nu_lower = get_expopac_bin_nu_lower(binindex);
  const auto delta_nu = get_expopac_bin_nu_upper(binindex) - bin_nu_lower;
  const double nuoffset = rng_uniform(rngstate) * delta_nu;
  const double nu = bin_nu_lower + nuoffset;
  return nu;
}

DEVICE_FUNC void do_rpkt(Packet& pkt, const double t2, ContinuumOpacity& chi_rpkt_cont) {
  while (do_rpkt_step(pkt, t2, chi_rpkt_cont)) {
    {
    }
  }
}

// make the packet an r-pkt and set further flags
DEVICE_FUNC void emit_rpkt(Packet& pkt) {
  pkt.type = TYPE_RPKT;

  // Need to assign a new direction. Assume isotropic emission in the cmf

  const auto dir_cmf = get_rand_isotropic_unitvec(get_rngstate(pkt));

  // This direction is in the cmf - we want to convert it to the rest
  // frame - use aberration of angles. We want to convert from cmf to
  // rest so need -ve velocity.
  const auto vel_vec = get_velocity(pkt.pos, -pkt.prop_time);
  // negative time since we want the backwards transformation here

  pkt.dir = angle_ab(dir_cmf, vel_vec);

  // Finally we want to put in the rest frame energy and frequency. And record that it's now a r-pkt.

  set_pkt_restframe_from_cmf(pkt);

  if constexpr (POL_ON) {
    // Reset to unpolarised
    pkt.stokes_u = 0.;
    pkt.stokes_q = 0.;
  }

  pkt.em_pos = pkt.pos;
  pkt.em_time = static_cast<float>(pkt.prop_time);
}

template <bool USECELLHISTANDUPDATEPHIXSLIST>
void calculate_chi_rpkt_cont(const double nu_cmf, ContinuumOpacity& chi_rpkt_cont, const int nonemptymgi) {
  assert_testmodeonly(grid::thick_allcells[nonemptymgi] != 1);
  if ((nonemptymgi == chi_rpkt_cont.nonemptymgi) && (globals::timestep == chi_rpkt_cont.timestep) &&
      (fabs((chi_rpkt_cont.nu / nu_cmf) - 1.0) < 1e-4)) {
    // calculated values are a match already
    return;
  }

  const auto nne = grid::get_nne(nonemptymgi);

  // free-free absorption
  chi_rpkt_cont.chi_freefree_heat = calculate_chi_ffheating(nonemptymgi, nu_cmf, USECELLHISTANDUPDATEPHIXSLIST);

  // Thomson scattering on free electrons
  chi_rpkt_cont.chi_escatter = SIGMA_T * nne;

  // bound-free absorption
  chi_rpkt_cont.chi_boundfree =
      calculate_chi_bf_gammacontr<USECELLHISTANDUPDATEPHIXSLIST>(nonemptymgi, nu_cmf, chi_rpkt_cont.phixslist);

  chi_rpkt_cont.nonemptymgi = nonemptymgi;
  chi_rpkt_cont.timestep = globals::timestep;
  chi_rpkt_cont.nu = nu_cmf;
}

// specialize calculate_chi_rpkt_cont templates with true and false:
template void calculate_chi_rpkt_cont<true>(const double nu_cmf, ContinuumOpacity& chi_rpkt_cont,
                                            const int nonemptymgi);
template void calculate_chi_rpkt_cont<false>(const double nu_cmf, ContinuumOpacity& chi_rpkt_cont,
                                             const int nonemptymgi);

void MPI_Bcast_binned_opacities(const ptrdiff_t nonemptymgi, const int root_node_id) {
  if (globals::rank_in_node == 0) {
    assert_always(nonemptymgi >= 0);
    if constexpr (RPKT_USE_EXPANSION_OPACITIES || VPKT_USE_EXPANSION_OPACITIES) {
      MPI_Bcast_safe(expansionopacities.subspan(nonemptymgi * expopac_nbins, expopac_nbins), root_node_id,
                     globals::mpi_comm_internode);
    }

    if constexpr (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value()) {
      MPI_Bcast_safe(expansionopacity_planck_cumulative.subspan(nonemptymgi * expopac_nbins, expopac_nbins),
                     root_node_id, globals::mpi_comm_internode);
    }
  }
}

void calculate_expansion_opacities(const int nonemptymgi) {
  const auto rho = grid::get_rho(nonemptymgi);

  const auto sys_time_start_calc = std::chrono::steady_clock::now();
  const auto temperature = grid::Te_allcells[nonemptymgi];

  printlog("calculating expansion opacities for cell {}...", grid::get_mgi_of_nonemptymgi(nonemptymgi));

  const auto t_mid = globals::timesteps[globals::timestep].mid;

  // find the first line with nu below the upper limit of the first bin
  int lineindex = static_cast<int>(
      std::ranges::lower_bound(globals::linelist.nu, get_expopac_bin_nu_upper(0), std::ranges::greater{}) -
      globals::linelist.nu.begin());

  double kappa_planck_cumulative = 0.;

  for (auto binindex = 0Z; binindex < expopac_nbins; binindex++) {
    double bin_linesum = 0.;
    const auto nu_lower = get_expopac_bin_nu_lower(binindex);

    while (lineindex < globals::nlines && globals::linelist.nu[lineindex] >= nu_lower) {
      const auto tau_line = static_cast<float>(get_tau_sobolev<false>(nonemptymgi, lineindex, t_mid));
      const auto linelambda = 1e8 * CLIGHT / globals::linelist.nu[lineindex];
      bin_linesum += (linelambda / expopac_deltalambda) * -std::expm1(-tau_line);
      lineindex++;
    }
    // opacity in units of [cm^2/g]
    const auto bin_kappa_bb = static_cast<float>(1. / (CLIGHT * t_mid * rho) * bin_linesum);
    assert_always(std::isfinite(bin_kappa_bb));

    if constexpr (RPKT_USE_EXPANSION_OPACITIES || VPKT_USE_EXPANSION_OPACITIES) {
      expansionopacities[(nonemptymgi * expopac_nbins) + binindex] = bin_kappa_bb;
    }

    if constexpr (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value()) {
      const auto nu_upper = get_expopac_bin_nu_upper(binindex);
      const auto nu_mid = (nu_upper + nu_lower) / 2.;
      const auto bin_kappa_cont = calculate_chi_ffheating(nonemptymgi, nu_mid, false) / rho;

      const auto planck_val = radfield::planck(nu_mid, temperature);
      const auto kappa_planck = (bin_kappa_bb + bin_kappa_cont) * planck_val;

      const auto delta_nu = nu_upper - nu_lower;
      kappa_planck_cumulative += kappa_planck * delta_nu;

      expansionopacity_planck_cumulative[(nonemptymgi * expopac_nbins) + binindex] = kappa_planck_cumulative;
    }
  }
  const auto expansion_opacity_duration =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_calc).count();
  printlnlog("took {:.1f} seconds", expansion_opacity_duration);
}

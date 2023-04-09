#include "macroatom.h"

#include <gsl/gsl_integration.h>

#include <algorithm>
#include <cmath>

#include "artisoptions.h"
#include "globals.h"
#include "grid.h"
#include "ltepop.h"
#include "nonthermal.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "rpkt.h"
#include "sn3d.h"
#include "stats.h"
#include "vpkt.h"

// save to the macroatom_*.out file
constexpr bool LOG_MACROATOM = false;

static FILE *macroatom_file = nullptr;

static void calculate_macroatom_transitionrates(const int modelgridindex, const int element, const int ion,
                                                const int level, const double t_mid, struct chlevels *const chlevel) {
  // printout("Calculating transition rates for element %d ion %d level %d\n", element, ion, level);
  double *processrates = chlevel->processrates;
  const auto T_e = grid::get_Te(modelgridindex);
  const auto nne = grid::get_nne(modelgridindex);
  const double epsilon_current = epsilon(element, ion, level);
  const double statweight = stat_weight(element, ion, level);

  /// Downward transitions within the current ionisation stage:
  /// radiative/collisional deexcitation and internal downward jumps
  processrates[MA_ACTION_RADDEEXC] = 0.;
  processrates[MA_ACTION_COLDEEXC] = 0.;
  processrates[MA_ACTION_INTERNALDOWNSAME] = 0.;
  const int ndowntrans = get_ndowntrans(element, ion, level);
  for (int i = 0; i < ndowntrans; i++) {
    const int lower = globals::elements[element].ions[ion].levels[level].downtrans[i].targetlevelindex;
    const auto A_ul = globals::elements[element].ions[ion].levels[level].downtrans[i].einstein_A;
    const double epsilon_target = epsilon(element, ion, lower);
    const double epsilon_trans = epsilon_current - epsilon_target;

    const double R =
        rad_deexcitation_ratecoeff(modelgridindex, element, ion, level, lower, epsilon_trans, A_ul, statweight, t_mid);
    const double C = col_deexcitation_ratecoeff(T_e, nne, epsilon_trans, element, ion, level, i);

    const double individ_internal_down_same = (R + C) * epsilon_target;

    const double sum_epstrans_rad_deexc = R * epsilon_trans;
    const double individ_col_deexc = C * epsilon_trans;

    processrates[MA_ACTION_RADDEEXC] += sum_epstrans_rad_deexc;
    processrates[MA_ACTION_COLDEEXC] += individ_col_deexc;
    processrates[MA_ACTION_INTERNALDOWNSAME] += individ_internal_down_same;

    chlevel->sum_epstrans_rad_deexc[i] = processrates[MA_ACTION_RADDEEXC];
    chlevel->sum_internal_down_same[i] = processrates[MA_ACTION_INTERNALDOWNSAME];

    // printout("checking downtrans %d to level %d: R %g, C %g, epsilon_trans %g\n",i,lower,R,C,epsilon_trans);
  }

  /// Downward transitions to lower ionisation stages:
  /// radiative/collisional recombination and internal downward jumps
  processrates[MA_ACTION_RADRECOMB] = 0.;
  processrates[MA_ACTION_COLRECOMB] = 0.;
  processrates[MA_ACTION_INTERNALDOWNLOWER] = 0.;
  // checks only if there is a lower ion, doesn't make sure that Z(ion)=Z(ion-1)+1
  if (ion > 0 && level <= get_maxrecombininglevel(element, ion)) {
    // nlevels = get_nlevels(element,ion-1);
    const int nlevels = get_ionisinglevels(element, ion - 1);
    // nlevels = get_ionisinglevels(element,ion-1);
    for (int lower = 0; lower < nlevels; lower++) {
      const double epsilon_target = epsilon(element, ion - 1, lower);
      const double epsilon_trans = epsilon_current - epsilon_target;

      const double R = rad_recombination_ratecoeff(T_e, nne, element, ion, level, lower, modelgridindex);
      const double C = col_recombination_ratecoeff(modelgridindex, element, ion, level, lower, epsilon_trans);

      processrates[MA_ACTION_INTERNALDOWNLOWER] += (R + C) * epsilon_target;

      processrates[MA_ACTION_RADRECOMB] += R * epsilon_trans;
      processrates[MA_ACTION_COLRECOMB] += C * epsilon_trans;
    }
  }

  /// Calculate sum for upward internal transitions
  /// transitions within the current ionisation stage
  processrates[MA_ACTION_INTERNALUPSAME] = 0.;
  const int nuptrans = get_nuptrans(element, ion, level);
  for (int i = 0; i < nuptrans; i++) {
    const int upper = globals::elements[element].ions[ion].levels[level].uptrans[i].targetlevelindex;
    const int lineindex = globals::elements[element].ions[ion].levels[level].uptrans[i].lineindex;
    const double epsilon_trans = epsilon(element, ion, upper) - epsilon_current;

    const double R = rad_excitation_ratecoeff(modelgridindex, element, ion, level, i, epsilon_trans, lineindex, t_mid);
    const double C = col_excitation_ratecoeff(T_e, nne, element, ion, level, i, epsilon_trans, statweight);
    const double NT =
        nonthermal::nt_excitation_ratecoeff(modelgridindex, element, ion, level, upper, epsilon_trans, lineindex);

    const double individ_internal_up_same = (R + C + NT) * epsilon_current;

    processrates[MA_ACTION_INTERNALUPSAME] += individ_internal_up_same;
    chlevel->sum_internal_up_same[i] = processrates[MA_ACTION_INTERNALUPSAME];
  }

  assert_always(std::isfinite(processrates[MA_ACTION_INTERNALUPSAME]));

  /// Transitions to higher ionisation stages
  processrates[MA_ACTION_INTERNALUPHIGHERNT] = 0.;
  processrates[MA_ACTION_INTERNALUPHIGHER] = 0.;
  const int ionisinglevels = get_ionisinglevels(element, ion);
  if (ion < get_nions(element) - 1 && level < ionisinglevels) {
    if (NT_ON) {
      processrates[MA_ACTION_INTERNALUPHIGHERNT] =
          nonthermal::nt_ionization_ratecoeff(modelgridindex, element, ion) * epsilon_current;
    }

    for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++) {
      // const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
      // const double epsilon_trans = epsilon(element, ion + 1, upper) - epsilon(element, ion, level);
      const double epsilon_trans = get_phixs_threshold(element, ion, level, phixstargetindex);

      const double R = get_corrphotoioncoeff(element, ion, level, phixstargetindex, modelgridindex);
      const double C = col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);

      processrates[MA_ACTION_INTERNALUPHIGHER] += (R + C) * epsilon_current;
    }
  }
}

static int do_macroatom_internal_down_same(int modelgridindex, int element, int ion, int level, double t_mid,
                                           double total_internal_down_same) {
  const int ndowntrans = get_ndowntrans(element, ion, level);

  // printout("[debug] do_ma:   internal downward jump within current ionstage\n");

  /// Randomly select the occuring transition
  const double zrand = rng_uniform();
  const double targetval = zrand * total_internal_down_same;

  const double *sum_internal_down_same =
      globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].sum_internal_down_same;

  // first sum_internal_down_same[i] such that sum_internal_down_same[i] > targetval
  const double *const upperval =
      std::upper_bound(&sum_internal_down_same[0], &sum_internal_down_same[ndowntrans], targetval);
  const int downtransindex = upperval - &sum_internal_down_same[0];

  assert_always(downtransindex < ndowntrans);
  const int lower = globals::elements[element].ions[ion].levels[level].downtrans[downtransindex].targetlevelindex;

  return lower;
}

static void do_macroatom_raddeexcitation(struct packet *pkt_ptr, const int modelgridindex, const int element,
                                         const int ion, const int level, const double rad_deexc,
                                         const int activatingline, const double t_mid) {
  /// radiative deexcitation of MA: emitt rpkt
  /// randomly select which line transitions occurs
  int linelistindex = -99;
  const int ndowntrans = get_ndowntrans(element, ion, level);

  double *sum_epstrans_rad_deexc =
      globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].sum_epstrans_rad_deexc;

  const double zrand = rng_uniform();
  const double targetval = zrand * rad_deexc;

  // first sum_epstrans_rad_deexc[i] such that sum_epstrans_rad_deexc[i] > targetval
  const double *const upperval =
      std::upper_bound(&sum_epstrans_rad_deexc[0], &sum_epstrans_rad_deexc[ndowntrans], targetval);
  const int downtransindex = upperval - &sum_epstrans_rad_deexc[0];

  assert_always(downtransindex < ndowntrans);
  linelistindex = globals::elements[element].ions[ion].levels[level].downtrans[downtransindex].lineindex;

  if (linelistindex == activatingline) {
    stats::increment(stats::COUNTER_RESONANCESCATTERINGS);
  }

  if constexpr (RECORD_LINESTAT) {
    safeincrement(globals::ecounter[linelistindex]);
  }

  const int lower = globals::elements[element].ions[ion].levels[level].downtrans[downtransindex].targetlevelindex;

  // printout("[debug] do_ma:   jump to level %d\n", lower);

  const double epsilon_trans = epsilon(element, ion, level) - epsilon(element, ion, lower);

  double oldnucmf;
  if (pkt_ptr->last_event == 1) oldnucmf = pkt_ptr->nu_cmf;
  pkt_ptr->nu_cmf = epsilon_trans / H;

  if (pkt_ptr->last_event == 1) {
    if (oldnucmf < pkt_ptr->nu_cmf) {
      stats::increment(stats::COUNTER_UPSCATTER);
    } else {
      stats::increment(stats::COUNTER_DOWNSCATTER);
    }
  }

  stats::increment(stats::COUNTER_MA_STAT_DEACTIVATION_BB);
  pkt_ptr->interactions += 1;
  pkt_ptr->last_event = 0;

  // emit the rpkt in a random direction
  emitt_rpkt(pkt_ptr);

  // the r-pkt can only interact with lines redder than the current one
  pkt_ptr->next_trans = linelistindex + 1;
  pkt_ptr->emissiontype = linelistindex;
  vec_copy(pkt_ptr->em_pos, pkt_ptr->pos);
  pkt_ptr->em_time = pkt_ptr->prop_time;
  pkt_ptr->nscatterings = 0;

#ifdef VPKT_ON
  const int realtype = 3;
  vpkt_call_estimators(pkt_ptr, pkt_ptr->prop_time, realtype);
#endif
}

static void do_macroatom_radrecomb(struct packet *pkt_ptr, const int modelgridindex, const int element, int *ion,
                                   int *level, const double rad_recomb) {
  const auto T_e = grid::get_Te(modelgridindex);
  const float nne = grid::get_nne(modelgridindex);
  const double epsilon_current = epsilon(element, *ion, *level);
  const int upperion = *ion;
  const int upperionlevel = *level;
  /// Randomly select a continuum
  double zrand = rng_uniform();
  double rate = 0;
  const int nlevels = get_ionisinglevels(element, upperion - 1);
  int lower = 0;
  for (lower = 0; lower < nlevels; lower++) {
    const double epsilon_trans = epsilon_current - epsilon(element, upperion - 1, lower);
    const double R = rad_recombination_ratecoeff(T_e, nne, element, upperion, upperionlevel, lower, modelgridindex);

    rate += R * epsilon_trans;

    // printout("[debug] do_ma:   R %g, deltae %g\n",R,(epsilon(element, upperion, upperionlevel) - epsilon(element,
    // upperion - 1, lower))); printout("[debug] do_ma:   rate to level %d of ion %d = %g\n", lower, upperion - 1,
    // rate); printout("[debug] do_ma:   zrand*rad_recomb = %g\n", zrand * rad_recomb);

    if (zrand * rad_recomb < rate) {
      break;
    }
  }
  if (zrand * rad_recomb >= rate) {
    printout(
        "%s: From Z=%d ionstage %d level %d, could not select lower level to recombine to. zrand %g * rad_recomb %g >= "
        "rate %g",
        __func__, get_element(element), get_ionstage(element, *ion), *level, zrand, rad_recomb, rate);
    abort();
  }

  /// set the new state
  *ion = upperion - 1;
  *level = lower;

  pkt_ptr->nu_cmf = select_continuum_nu(element, upperion - 1, lower, upperionlevel, T_e);

#ifndef FORCE_LTE
  // mabfcount[pkt_ptr->where] += pkt_ptr->e_cmf;
  // mabfcount_thermal[pkt_ptr->where] += pkt_ptr->e_cmf*(1-nu_threshold/pkt_ptr->nu_cmf);
  // matotem[pkt_ptr->where] += pkt_ptr->e_cmf;
#endif

  // printout("%s: From Z=%d ionstage %d, recombining to ionstage %d level %d\n",
  //          __func__, get_element(element), get_ionstage(element, *ion + 1), get_ionstage(element, *ion), lower);
  // printout("[debug] do_ma:   pkt_ptr->nu_cmf %g\n",pkt_ptr->nu_cmf);

  if (!std::isfinite(pkt_ptr->nu_cmf)) {
    printout("[fatal] rad recombination of MA: selected frequency not finite ... abort\n");
    abort();
  }
  stats::increment(stats::COUNTER_MA_STAT_DEACTIVATION_FB);
  pkt_ptr->interactions += 1;
  pkt_ptr->last_event = 2;

  /// Finally emit the packet into a randomly chosen direction, update the continuum opacity and set some flags
  emitt_rpkt(pkt_ptr);

  if constexpr (TRACK_ION_STATS) {
    stats::increment_ion_stats(modelgridindex, element, upperion, stats::ION_RADRECOMB_MACROATOM,
                               pkt_ptr->e_cmf / H / pkt_ptr->nu_cmf);

    const double escape_prob = get_rpkt_escape_prob(pkt_ptr, pkt_ptr->prop_time);

    stats::increment_ion_stats(modelgridindex, element, upperion, stats::ION_RADRECOMB_ESCAPED,
                               pkt_ptr->e_cmf / H / pkt_ptr->nu_cmf * escape_prob);
  }

  pkt_ptr->next_trans = 0;  /// continuum transition, no restrictions for further line interactions
  pkt_ptr->emissiontype = get_continuumindex(element, *ion, lower, upperionlevel);
  vec_copy(pkt_ptr->em_pos, pkt_ptr->pos);
  pkt_ptr->em_time = pkt_ptr->prop_time;
  pkt_ptr->nscatterings = 0;

#ifdef VPKT_ON
  const int realtype = 3;
  vpkt_call_estimators(pkt_ptr, pkt_ptr->prop_time, realtype);
#endif
}

static void do_macroatom_ionisation(const int modelgridindex, const int element, int *ion, int *level,
                                    const double epsilon_current, const double internal_up_higher) {
  const auto T_e = grid::get_Te(modelgridindex);
  const float nne = grid::get_nne(modelgridindex);

  int upper = -1;
  /// Randomly select the occuring transition
  const double zrand = rng_uniform();
  double rate = 0.;
  for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, *ion, *level); phixstargetindex++) {
    upper = get_phixsupperlevel(element, *ion, *level, phixstargetindex);
    // const double epsilon_trans = epsilon(element, *ion + 1, upper) - epsilon_current;
    const double epsilon_trans = get_phixs_threshold(element, *ion, *level, phixstargetindex);
    const double R = get_corrphotoioncoeff(element, *ion, *level, phixstargetindex, modelgridindex);
    const double C = col_ionization_ratecoeff(T_e, nne, element, *ion, *level, phixstargetindex, epsilon_trans);
    rate += (R + C) * epsilon_current;
    if (zrand * internal_up_higher < rate) break;
  }
  if (zrand * internal_up_higher >= rate) {
    printout(
        "%s: From Z=%d ionstage %d level %d, could not select upper level to ionise to. zrand %g * internal_up_higher "
        "%g >= rate %g\n",
        __func__, get_element(element), get_ionstage(element, *ion), *level, zrand, internal_up_higher, rate);
    abort();
  }

  assert_always(upper >= 0);

  /// and set the macroatom's new state
  *ion += 1;
  *level = upper;
}

void do_macroatom(struct packet *pkt_ptr, const int timestep)
/// Material for handling activated macro atoms.
{
  const int tid = get_thread_num();
  const double t_mid = globals::time_step[timestep].mid;

  // printout("[debug] do MA\n");

  const int cellindex = pkt_ptr->where;
  const int modelgridindex = grid::get_cell_modelgridindex(cellindex);
  const auto T_e = grid::get_Te(modelgridindex);
  const float nne = grid::get_nne(modelgridindex);

  assert_always(grid::modelgrid[modelgridindex].thick != 1);  // macroatom should not be used in thick cells

  /// calculate occupation number for active MA level ////////////////////////////////////
  /// general QUESTION: is it better to calculate the n_1 (later the n_ionstage and
  /// U_ionstage) here where we need them or once in update_grid for each grid cell
  /// not sure whether this reduces the number of calculations, as number of grid cells
  /// is much larger than number of pellets (next question: connection to number of
  /// photons)
  const int element = pkt_ptr->mastate.element;
  int ion = pkt_ptr->mastate.ion;
  int level = pkt_ptr->mastate.level;

  const int activatingline = pkt_ptr->mastate.activatingline;
  if (pkt_ptr->absorptiontype > 0 && activatingline > 0 && activatingline != pkt_ptr->absorptiontype) {
    printout("error: mismatched absorptiontype %d != activatingline = %d pkt last_event %d emissiontype %d\n",
             pkt_ptr->absorptiontype, activatingline, pkt_ptr->last_event, pkt_ptr->emissiontype);
  }

  const int ion_in = ion;
  const int level_in = level;
  const double nu_cmf_in = pkt_ptr->nu_cmf;
  const double nu_rf_in = pkt_ptr->nu_rf;

  if constexpr (TRACK_ION_STATS) {
    stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYIN_TOTAL, pkt_ptr->e_cmf);
  }

  int jumps = 0;

  bool end_packet = false;
  while (!end_packet) {
    // ionisinglevels = get_ionisinglevels(element,ion);

    /// Set this here to 1 to overcome problems in cells which have zero population
    /// in some ionisation stage. This is possible because the dependence on the
    /// originating levels population cancels out in the macroatom transition probabilities
    /// which are based on detailed balance.

    // printout("[debug] %s Z=%d ionstage %d level %d, jumps %d\n", __func__, get_element(element),
    // get_ionstage(element,ion), level, jumps);

    assert_always(ion >= 0);
    assert_always(ion < get_nions(element));

    const double epsilon_current = epsilon(element, ion, level);
    // const int ndowntrans = get_ndowntrans(element, ion, level);
    const int nuptrans = get_nuptrans(element, ion, level);

    assert_always(globals::cellhistory[tid].cellnumber == modelgridindex);

    struct chlevels *chlevel = &globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level];
    const double *processrates = chlevel->processrates;
    /// If there are no precalculated rates available then calculate them
    if (processrates[MA_ACTION_COLDEEXC] < 0) {
      calculate_macroatom_transitionrates(modelgridindex, element, ion, level, t_mid, chlevel);
    }

    // for debugging the transition rates:
    // {
    //   printout("macroatom element %d ion %d level %d\n", element, ion, level);
    //
    // const char *actionlabel[MA_ACTION_COUNT] = {
    //   "MA_ACTION_RADDEEXC", "MA_ACTION_COLDEEXC", "MA_ACTION_RADRECOMB",
    //   "MA_ACTION_COLRECOMB", "MA_ACTION_INTERNALDOWNSAME", "MA_ACTION_INTERNALDOWNLOWER",
    //   "MA_ACTION_INTERNALUPSAME", "MA_ACTION_INTERNALUPHIGHER", "MA_ACTION_INTERNALUPHIGHERNT"};

    //   for (enum ma_action action = 0; action < MA_ACTION_COUNT; action++)
    //     printout("actions: %30s %g\n", actionlabel[action], processrates[action]);
    // }

    // select transition according to probabilities
    double total_transitions = 0.;
    for (int action = 0; action < MA_ACTION_COUNT; action++) {
      total_transitions += processrates[action];
    }
    assert_always(total_transitions > 0);

    enum ma_action selected_action = MA_ACTION_COUNT;
    double zrand = rng_uniform();
    const double randomrate = zrand * total_transitions;
    double rate = 0.;
    for (int action = 0; action < MA_ACTION_COUNT; action++) {
      rate += processrates[action];
      if (rate > randomrate) {
        selected_action = static_cast<enum ma_action>(action);
        break;
      }
    }

    assert_always(rate > randomrate);

    switch (selected_action) {
      case MA_ACTION_RADDEEXC: {
        // printout("[debug] do_ma:   radiative deexcitation\n");
        // printout("[debug] do_ma:   jumps = %d\n", jumps);

        do_macroatom_raddeexcitation(pkt_ptr, modelgridindex, element, ion, level, processrates[MA_ACTION_RADDEEXC],
                                     activatingline, t_mid);

        if constexpr (TRACK_ION_STATS) {
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYOUT_RADDEEXC,
                                     pkt_ptr->e_cmf);

          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_BOUNDBOUND_MACROATOM,
                                     pkt_ptr->e_cmf / H / pkt_ptr->nu_cmf);

          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYOUT_TOTAL,
                                     pkt_ptr->e_cmf);
        }

        if constexpr (LOG_MACROATOM) {
          fprintf(macroatom_file, "%8d %14d %2d %12d %12d %9d %9d %9d %11.5e %11.5e %11.5e %11.5e %9d\n", timestep,
                  modelgridindex, get_element(element), get_ionstage(element, ion_in), get_ionstage(element, ion),
                  level_in, level, activatingline, nu_cmf_in, pkt_ptr->nu_cmf, nu_rf_in, pkt_ptr->nu_rf, jumps);
        }

        end_packet = true;
        break;
      }

      case MA_ACTION_COLDEEXC: {
        /// collisional deexcitation of macro atom => convert the packet into a k-packet
        // printout("[debug] do_ma:   collisional deexcitation\n");
        // printout("[debug] do_ma: jumps = %d\n", jumps);

        stats::increment(stats::COUNTER_MA_STAT_DEACTIVATION_COLLDEEXC);
        pkt_ptr->interactions += 1;
        pkt_ptr->last_event = 10;

        if constexpr (TRACK_ION_STATS) {
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYOUT_COLLDEEXC,
                                     pkt_ptr->e_cmf);
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYOUT_TOTAL,
                                     pkt_ptr->e_cmf);
        }

        pkt_ptr->type = TYPE_KPKT;
        end_packet = true;
#ifndef FORCE_LTE
        safeadd(globals::colheatingestimator[modelgridindex], pkt_ptr->e_cmf);
#endif
        break;
      }

      case MA_ACTION_INTERNALDOWNSAME: {
        pkt_ptr->interactions += 1;
        jumps++;
        level = do_macroatom_internal_down_same(modelgridindex, element, ion, level, t_mid,
                                                processrates[MA_ACTION_INTERNALDOWNSAME]);

        break;
      }

      case MA_ACTION_RADRECOMB: {
        /// Radiative recombination of MA: emitt a continuum-rpkt
        // printout("[debug] do_ma:   radiative recombination\n");
        // printout("[debug] do_ma:   jumps = %d\n", jumps);
        // printout("[debug] do_ma:   element %d, ion %d, level %d\n", element, ion, level);

        if constexpr (TRACK_ION_STATS) {
          // stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYOUT_RADRECOMB,
          // pkt_ptr->e_cmf); stats::increment_ion_stats(modelgridindex, element, ion,
          // stats::ION_MACROATOM_ENERGYOUT_TOTAL, pkt_ptr->e_cmf);
        }

        do_macroatom_radrecomb(pkt_ptr, modelgridindex, element, &ion, &level, processrates[MA_ACTION_RADRECOMB]);
        end_packet = true;
        break;
      }

      case MA_ACTION_COLRECOMB: {
        /// collisional recombination of macro atom => convert the packet into a k-packet
        // printout("[debug] do_ma:   collisonal recombination\n");
        // printout("[debug] do_ma: jumps = %d\n",jumps);
        stats::increment(stats::COUNTER_MA_STAT_DEACTIVATION_COLLRECOMB);
        pkt_ptr->interactions += 1;
        pkt_ptr->last_event = 11;

        if constexpr (TRACK_ION_STATS) {
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYOUT_COLLRECOMB,
                                     pkt_ptr->e_cmf);
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYOUT_TOTAL,
                                     pkt_ptr->e_cmf);
        }

        pkt_ptr->type = TYPE_KPKT;
        end_packet = true;
#ifndef FORCE_LTE
        safeadd(globals::colheatingestimator[modelgridindex], pkt_ptr->e_cmf);
#endif
        break;
      }

      case MA_ACTION_INTERNALDOWNLOWER: {
        // printout("[debug] do_ma:   internal downward jump to lower ionstage\n");
        pkt_ptr->interactions += 1;
        jumps++;

        stats::increment(stats::COUNTER_MA_STAT_INTERNALDOWNLOWER);

        /// Randomly select the occuring transition
        zrand = rng_uniform();
        // zrand = 1. - 1e-14;
        rate = 0.;
        // nlevels = get_nlevels(element,ion-1);

        const int nlevels = get_ionisinglevels(element, ion - 1);
        // nlevels = get_ionisinglevels(element,ion-1);
        int lower;
        for (lower = 0; lower < nlevels; lower++) {
          const double epsilon_target = epsilon(element, ion - 1, lower);
          const double epsilon_trans = epsilon_current - epsilon_target;
          const double R = rad_recombination_ratecoeff(T_e, nne, element, ion, level, lower, modelgridindex);
          const double C = col_recombination_ratecoeff(modelgridindex, element, ion, level, lower, epsilon_trans);
          rate += (R + C) * epsilon_target;
          if (zrand * processrates[MA_ACTION_INTERNALDOWNLOWER] < rate) break;
        }
        /// and set the macroatom's new state

        if constexpr (TRACK_ION_STATS) {
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYOUT_INTERNAL,
                                     pkt_ptr->e_cmf);
        }

        ion -= 1;
        level = lower;

        if constexpr (TRACK_ION_STATS) {
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYIN_INTERNAL,
                                     pkt_ptr->e_cmf);
        }

        if (lower >= nlevels) {
          printout("internal_down_lower  %g\n", processrates[MA_ACTION_INTERNALDOWNLOWER]);
          printout("abort at rate %g, zrand %g\n", rate, zrand);
          abort();
        }
        if (get_ionstage(element, ion) == 0 && lower == 0) {
          printout("internal downward transition to ground level occured ... abort\n");
          printout("element %d, ion %d, level %d, lower %d\n", element, ion, level, lower);
          printout("Z %d, ionstage %d, energy %g\n", get_element(element), get_ionstage(element, ion - 1),
                   globals::elements[element].ions[ion - 1].levels[lower].epsilon);
          printout("[debug] do_ma:   internal downward jump to lower ionstage\n");
          abort();
        }
        break;
      }

      case MA_ACTION_INTERNALUPSAME: {
        // printout("[debug] do_ma:   internal upward jump within current ionstage\n");
        pkt_ptr->interactions += 1;
        jumps++;

        /// randomly select the occuring transition
        zrand = rng_uniform();
        const double *sum_internal_up_same =
            globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].sum_internal_up_same;

        const double targetval = zrand * processrates[MA_ACTION_INTERNALUPSAME];

        // first sum_internal_up_same[i] such that sum_internal_up_same[i] > targetval
        const double *const upperval =
            std::upper_bound(&sum_internal_up_same[0], &sum_internal_up_same[nuptrans], targetval);
        const int uptransindex = upperval - &sum_internal_up_same[0];

        assert_always(uptransindex < nuptrans);
        const int upper = globals::elements[element].ions[ion].levels[level].uptrans[uptransindex].targetlevelindex;

        level = upper;
        break;
      }

      case MA_ACTION_INTERNALUPHIGHER: {
        // printout("[debug] do_ma:   internal upward jump to next ionstage\n");
        pkt_ptr->interactions += 1;
        jumps++;

        stats::increment(stats::COUNTER_MA_STAT_INTERNALUPHIGHER);

        if constexpr (TRACK_ION_STATS) {
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYOUT_INTERNAL,
                                     pkt_ptr->e_cmf);
        }

        do_macroatom_ionisation(modelgridindex, element, &ion, &level, epsilon_current,
                                processrates[MA_ACTION_INTERNALUPHIGHER]);

        if constexpr (TRACK_ION_STATS) {
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYIN_INTERNAL,
                                     pkt_ptr->e_cmf);
        }

        break;
      }

      case MA_ACTION_INTERNALUPHIGHERNT: {
        pkt_ptr->interactions += 1;
        // ion += 1;
        if constexpr (TRACK_ION_STATS) {
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYOUT_INTERNAL,
                                     pkt_ptr->e_cmf);
        }

        ion = nonthermal::nt_random_upperion(modelgridindex, element, ion, false);
        level = 0;
        stats::increment(stats::COUNTER_MA_STAT_INTERNALUPHIGHERNT);

        if constexpr (TRACK_ION_STATS) {
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYIN_INTERNAL,
                                     pkt_ptr->e_cmf);
        }
        // printout("Macroatom non-thermal ionisation to Z=%d ionstage %d level %d\n", get_element(element), ion,
        // level);
        break;
      }

      case MA_ACTION_COUNT: {
        printout("ERROR: Problem selecting MA_ACTION\n");
        abort();
      }
    }
  }  /// endwhile

  if (pkt_ptr->trueemissiontype == EMTYPE_NOTSET) {
    pkt_ptr->trueemissiontype = pkt_ptr->emissiontype;
    pkt_ptr->trueemissionvelocity = vec_len(pkt_ptr->em_pos) / pkt_ptr->em_time;
    pkt_ptr->trueem_time = pkt_ptr->em_time;
  }

  /// procedure ends only after a change to r or k packets has taken place
}

/// Calculation of radiative rates ///////////////////////////////////////////////////////

void macroatom_open_file(const int my_rank) {
  if constexpr (!LOG_MACROATOM) return;
  char filename[MAXFILENAMELENGTH];
  snprintf(filename, MAXFILENAMELENGTH, "macroatom_%.4d.out", my_rank);
  assert_always(macroatom_file == nullptr);
  macroatom_file = fopen_required(filename, "w");
  fprintf(macroatom_file, "%8s %14s %2s %12s %12s %9s %9s %9s %11s %11s %11s %11s %9s\n", "timestep", "modelgridindex",
          "Z", "ionstage_in", "ionstage_out", "level_in", "level_out", "activline", "nu_cmf_in", "nu_cmf_out",
          "nu_rf_in", "nu_rf_out", "jumps");
}

void macroatom_close_file() {
  if (macroatom_file != nullptr) {
    fclose(macroatom_file);
  }
}

double rad_deexcitation_ratecoeff(const int modelgridindex, const int element, const int ion, const int upper,
                                  const int lower, const double epsilon_trans, const float A_ul,
                                  const double upperstatweight, const double t_current)
/// radiative deexcitation rate: paperII 3.5.2
// multiply by upper level population to get a rate per second
{
  assert_always(upper > lower);

  const double n_u = get_levelpop(modelgridindex, element, ion, upper);
  const double n_l = get_levelpop(modelgridindex, element, ion, lower);

  double R = 0.0;

  // if ((n_u > 1.1 * MINPOP) && (n_l > 1.1 * MINPOP))
  {
    const double nu_trans = epsilon_trans / H;

    // const double A_ul = einstein_spontaneous_emission(lineindex);
    const double B_ul = CLIGHTSQUAREDOVERTWOH / std::pow(nu_trans, 3) * A_ul;
    const double B_lu = upperstatweight / stat_weight(element, ion, lower) * B_ul;

    const double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current;

    if (tau_sobolev > 1e-100) {
      const double beta = 1.0 / tau_sobolev * (-std::expm1(-tau_sobolev));
      // const double beta = 1.0;
      R = A_ul * beta;
    } else {
      // printout("[warning] rad_deexcitation: tau_sobolev %g <= 0, set beta=1\n",tau_sobolev);
      // printout("[warning] rad_deexcitation: element %d, ion %d, upper %d, lower %d\n",element,ion,upper,lower);
      // printout("[warning] rad_deexcitation: n_l %g, n_u %g, B_lu %g, B_ul %g\n",n_l,n_u,B_lu,B_ul);
      // printout("[warning] rad_deexcitation: T_e %g, T_R %g, W %g in model cell
      // %d\n",grid::get_Te(modelgridindex),get_TR(modelgridindex),get_W(modelgridindex),modelgridindex);
      R = 0.0;
      // printout("[fatal] rad_excitation: tau_sobolev <= 0 ... %g abort",tau_sobolev);
      // abort();
    }

    // printout("[debug] rad_rates_down: Z=%d, ionstage %d, upper %d, lower %d\n", get_element(element),
    // get_ionstage(element, ion), upper, lower); printout("[debug] rad_deexc: A_ul %g, tau_sobolev %g, n_u %g\n", A_ul,
    // tau_sobolev, n_u);
    assert_testmodeonly(std::isfinite(R));
  }

  return R;
}

double rad_excitation_ratecoeff(const int modelgridindex, const int element, const int ion, const int lower,
                                const int i, const double epsilon_trans, int lineindex, const double t_current)
/// radiative excitation rate: paperII 3.5.2
// multiply by lower level population to get a rate per second
{
  const int upper = globals::elements[element].ions[ion].levels[lower].uptrans[i].targetlevelindex;

  const double n_u = get_levelpop(modelgridindex, element, ion, upper);
  const double n_l = get_levelpop(modelgridindex, element, ion, lower);
  double R = 0.0;
  // if ((n_u >= 1.1 * MINPOP) && (n_l >= 1.1 * MINPOP))
  {
    const double nu_trans = epsilon_trans / H;
    const double A_ul = globals::elements[element].ions[ion].levels[lower].uptrans[i].einstein_A;
    const double B_ul = CLIGHTSQUAREDOVERTWOH / std::pow(nu_trans, 3) * A_ul;
    const double B_lu = stat_weight(element, ion, upper) / stat_weight(element, ion, lower) * B_ul;

    const double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current;

    if (tau_sobolev > 1e-100) {
      double beta = 1.0 / tau_sobolev * (-std::expm1(-tau_sobolev));

      const double R_over_J_nu = n_l > 0. ? (B_lu - B_ul * n_u / n_l) * beta : B_lu * beta;

      R = R_over_J_nu * radfield::radfield(nu_trans, modelgridindex);

      if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
        if (!globals::initial_iteration) {
          // check for a detailed line flux estimator to replace the binned/blackbody radiation field estimate
          const int jblueindex = radfield::get_Jblueindex(lineindex);
          if (jblueindex >= 0) {
            const double Jb_lu = radfield::get_Jb_lu(modelgridindex, jblueindex);
            const double R_Jb = R_over_J_nu * Jb_lu;
            // const int contribcount = radfield_get_Jb_lu_contribcount(modelgridindex, jblueindex);
            // const double R_radfield = R_over_J_nu * radfield::radfield(nu_trans, modelgridindex);
            // const double linelambda = 1e8 * CLIGHT / nu_trans;
            // printout("Using detailed rad excitation lambda %5.1f contribcont %d R(Jblue) %g R(radfield) %g R_Jb/R
            // %g\n",
            //          linelambda, contribcount, R_Jb, R_radfield, R_Jb / R_radfield);
            // printout("  (for transition Z=%02d ionstage %d lower %d upper %d)\n",
            //          get_element(element), get_ionstage(element, ion), lower, upper);
            R = R_Jb;
          }
        }
      }
    } else {
      // printout("[warning] rad_excitation: tau_sobolev %g <= 0, set beta=1\n",tau_sobolev);
      // printout("[warning] rad_excitation: element %d, ion %d, upper %d, lower %d\n",element,ion,upper,lower);
      // printout("[warning] rad_excitation: n_l %g, n_u %g, B_lu %g, B_ul %g\n",n_l,n_u,B_lu,B_ul);
      // printout("[warning] rad_excitation: T_e %g, T_R %g, W
      // %g\n",grid::get_Te(modelgridindex),get_TR(modelgridindex),get_W(modelgridindex));
      R = 0.;
    }

    assert_testmodeonly(R >= 0.);
    assert_testmodeonly(std::isfinite(R));
  }

  return R;
}

double rad_recombination_ratecoeff(const float T_e, const float nne, const int element, const int upperion,
                                   const int upper, const int lower, const int modelgridindex)
/// radiative recombination rate: paperII 3.5.2
// multiply by upper level population to get a rate per second
{
  // it's probably faster to only check this condition outside this function
  // in a case where this wasn't checked, the function will return zero anyway
  // if (upper > get_maxrecombininglevel(element, upperion))
  //   return 0.;

  double R = 0.0;
  const int nphixstargets = get_nphixstargets(element, upperion - 1, lower);
  for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
    if (get_phixsupperlevel(element, upperion - 1, lower, phixstargetindex) == upper) {
      R = nne *
          get_spontrecombcoeff(
              element, upperion - 1, lower, phixstargetindex,
              T_e);  // +
                     // stimrecombestimator_save[pkt_ptr->where*get_nelements()*get_max_nions()+element*get_max_nions()+(ion-1)]);
      // printout("calculate rad_recombination: element %d, ion %d, upper %d, -> lower %d, nne %g, spontrecombcoeff
      // %g\n",element,ion,upper,lower,nne,get_spontrecombcoeff(element, ion-1, lower, T_e));

      if (modelgridindex >= 0 && SEPARATE_STIMRECOMB) {
        R += nne * get_stimrecombcoeff(element, upperion - 1, lower, phixstargetindex, modelgridindex);
      }
      break;
    }
  }

  assert_testmodeonly(std::isfinite(R));

  return R;
}

double stim_recombination_ratecoeff(const float nne, const int element, const int upperion, const int upper,
                                    const int lower, const int modelgridindex)
/// radiative recombination rate: paperII 3.5.2
// multiply by upper level population to get a rate per second
{
  // it's probably faster to only check this condition outside this function
  // in a case where this wasn't checked, the function will return zero anyway
  // if (upper > get_maxrecombininglevel(element, upperion))
  //   return 0.;

  double R = 0.0;
  const int nphixstargets = get_nphixstargets(element, upperion - 1, lower);
  for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
    if (get_phixsupperlevel(element, upperion - 1, lower, phixstargetindex) == upper) {
      R = nne * get_stimrecombcoeff(element, upperion - 1, lower, phixstargetindex, modelgridindex);
      break;
    }
  }

  assert_always(std::isfinite(R));

  return R;
}

double col_recombination_ratecoeff(const int modelgridindex, const int element, const int upperion, const int upper,
                                   const int lower, const double epsilon_trans)
// multiply by upper level population to get a rate per second
{
  // it's probably faster to only check this condition outside this function
  // in a case where this wasn't checked, the function will return zero anyway
  // if (upper > get_maxrecombininglevel(element, upperion))
  //   return 0.;

  const int nphixstargets = get_nphixstargets(element, upperion - 1, lower);
  for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
    if (get_phixsupperlevel(element, upperion - 1, lower, phixstargetindex) == upper) {
      const float nne = grid::get_nne(modelgridindex);
      const auto T_e = grid::get_Te(modelgridindex);
      const double fac1 = epsilon_trans / KB / T_e;
      const int ionstage = get_ionstage(element, upperion);

      /// Seaton approximation: Mihalas (1978), eq.5-79, p.134
      /// select gaunt factor according to ionic charge
      double g;
      if (ionstage - 1 == 1) {
        g = 0.1;
      } else if (ionstage - 1 == 2) {
        g = 0.2;
      } else {
        g = 0.3;
      }

      const double sigma_bf = (globals::elements[element].ions[upperion - 1].levels[lower].photoion_xs[0] *
                               get_phixsprobability(element, upperion - 1, lower, phixstargetindex));

      const double sf = calculate_sahafact(element, upperion - 1, lower, upper, T_e, epsilon_trans);

      double C = nne * nne * sf * 1.55e13 * pow(T_e, -0.5) * g * sigma_bf * exp(-fac1) / fac1;

      return C;
    }
  }

  return 0.;
}

double col_ionization_ratecoeff(const float T_e, const float nne, const int element, const int ion, const int lower,
                                const int phixstargetindex, const double epsilon_trans)
/// collisional ionization rate: paperII 3.5.1
// multiply by lower level population to get a rate per second
{
  assert_testmodeonly(phixstargetindex >= 0);
  assert_testmodeonly(phixstargetindex < get_nphixstargets(element, ion, lower));

  /// Seaton approximation: Mihalas (1978), eq.5-79, p.134
  /// select gaunt factor according to ionic charge
  double g;
  const int ionstage = get_ionstage(element, ion);
  if (ionstage == 1) {
    g = 0.1;
  } else if (ionstage == 2) {
    g = 0.2;
  } else {
    g = 0.3;
  }

  const double fac1 = epsilon_trans / KB / T_e;

  const double sigma_bf = globals::elements[element].ions[ion].levels[lower].photoion_xs[0] *
                          get_phixsprobability(element, ion, lower, phixstargetindex);
  const double C = nne * 1.55e13 * pow(T_e, -0.5) * g * sigma_bf * exp(-fac1) / fac1;  /// photoionization at the edge

  // printout("[debug] col_ion: nne %g, T_e %g, g %g, epsilon_trans %g, sigma_bf %g\n",
  // nne,T_e,g,epsilon_trans,sigma_bf);
  assert_testmodeonly(std::isfinite(C));

  return C;
}

double col_deexcitation_ratecoeff(const float T_e, const float nne, const double epsilon_trans, int element, int ion,
                                  int upper, int i)
// multiply by upper level population to get a rate per second
{
  const int lower = globals::elements[element].ions[ion].levels[upper].downtrans[i].targetlevelindex;
  const double upperstatweight = stat_weight(element, ion, upper);
  const double lowerstatweight = stat_weight(element, ion, lower);
  const double coll_str_thisline = globals::elements[element].ions[ion].levels[upper].downtrans[i].coll_str;
  double C = 0.;
  if (coll_str_thisline < 0) {
    const double forbidden = globals::elements[element].ions[ion].levels[upper].downtrans[i].forbidden;
    if (!forbidden)  // alternative: (coll_strength > -1.5) i.e. to catch -1
    {
      /// permitted E1 electric dipole transitions
      /// collisional deexcitation: formula valid only for atoms!!!!!!!!!!!
      /// Rutten script eq. 3.33. p.50
      // f = osc_strength(element,ion,upper,lower);
      // C = n_u * 2.16 * pow(fac1,-1.68) * pow(T_e,-1.5) *
      // stat_weight(element,ion,lower)/stat_weight(element,ion,upper)  * nne * f;
      const double trans_osc_strength = globals::elements[element].ions[ion].levels[upper].downtrans[i].osc_strength;

      const double eoverkt = epsilon_trans / (KB * T_e);
      /// Van-Regemorter formula, Mihalas (1978), eq.5-75, p.133
      const double g_bar = 0.2;  /// this should be read in from transitions data: it is 0.2 for transitions nl -> n'l'
                                 /// and 0.7 for transitions nl -> nl'
      // test = 0.276 * exp(fac1) * gsl_sf_expint_E1(fac1);
      /// crude approximation to the already crude Van-Regemorter formula

      // double test = 0.276 * exp(fac1) * (-0.5772156649 - log(fac1));
      // double Gamma = (g_bar > test) ? g_bar : test;

      // optimisation
      const double gauntfac =
          (eoverkt > 0.33421) ? g_bar : 0.276 * std::exp(eoverkt) * (-0.5772156649 - std::log(eoverkt));

      const double g_ratio = lowerstatweight / upperstatweight;

      C = C_0 * 14.51039491 * nne * std::sqrt(T_e) * trans_osc_strength * std::pow(H_ionpot / epsilon_trans, 2) *
          eoverkt * g_ratio * gauntfac;
    } else  // alterative: (coll_strength > -3.5) to catch -2 or -3
    {
      // forbidden transitions: magnetic dipole, electric quadropole...
      // could be Axelrod? or Maurer
      C = nne * 8.629e-6 * 0.01 * lowerstatweight / std::sqrt(T_e);
    }
  } else  // positive values are treated as effective collision strengths
  {
    // from Osterbrock and Ferland, p51
    // statweight_target is LOWER LEVEL stat weight
    C = nne * 8.629e-6 * coll_str_thisline / upperstatweight / std::sqrt(T_e);
    // test test
    // C = n_u * nne * 8.629e-6 * pow(T_e,-0.5) * 0.01 * statweight_target;
  }

  return C;
}

double col_excitation_ratecoeff(const float T_e, const float nne, int element, int ion, int lower, int i,
                                const double epsilon_trans, const double lowerstatweight)
// multiply by lower level population to get a rate per second
{
  // assert_testmodeonly(i < get_nuptrans(element, ion, lower));
  double C = 0.;
  const double coll_strength = globals::elements[element].ions[ion].levels[lower].uptrans[i].coll_str;
  const double eoverkt = epsilon_trans / (KB * T_e);

  if (coll_strength < 0) {
    const bool forbidden = globals::elements[element].ions[ion].levels[lower].uptrans[i].forbidden;
    if (!forbidden)  // alternative: (coll_strength > -1.5) i.e. to catch -1
    {
      const double trans_osc_strength = globals::elements[element].ions[ion].levels[lower].uptrans[i].osc_strength;
      /// permitted E1 electric dipole transitions
      /// collisional excitation: formula valid only for atoms!!!!!!!!!!!
      /// Rutten script eq. 3.32. p.50
      // C = n_l * 2.16 * pow(eoverkt,-1.68) * pow(T_e,-1.5) * exp(-eoverkt) * nne *
      // osc_strength(element,ion,upper,lower);

      // Van-Regemorter formula, Mihalas (1978), eq.5-75, p.133
      const double g_bar = 0.2;  // this should be read in from transitions data: it is 0.2 for transitions nl -> n'l'
                                 // and 0.7 for transitions nl -> nl'
      // test = 0.276 * exp(eoverkt) * gsl_sf_expint_E1(eoverkt);
      /// crude approximation to the already crude Van-Regemorter formula
      const double exp_eoverkt = exp(eoverkt);

      const double test = 0.276 * exp_eoverkt * (-0.5772156649 - std::log(eoverkt));
      const double Gamma = g_bar > test ? g_bar : test;
      C = C_0 * nne * std::sqrt(T_e) * 14.51039491 * trans_osc_strength * pow(H_ionpot / epsilon_trans, 2) * eoverkt /
          exp_eoverkt * Gamma;
    } else  // alterative: (coll_strength > -3.5) to catch -2 or -3
    {
      // forbidden transitions: magnetic dipole, electric quadropole...
      // Axelrod's approximation (thesis 1980)
      const int upper = globals::elements[element].ions[ion].levels[lower].uptrans[i].targetlevelindex;
      const double upperstatweight = stat_weight(element, ion, upper);
      C = nne * 8.629e-6 * 0.01 * std::exp(-eoverkt) * upperstatweight / std::sqrt(T_e);
    }
  } else {
    // from Osterbrock and Ferland, p51
    C = nne * 8.629e-6 * coll_strength * std::exp(-eoverkt) / lowerstatweight / std::sqrt(T_e);
  }

  return C;
}

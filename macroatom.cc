#include "macroatom.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <format>
#include <fstream>
#include <ios>
#include <span>

#if defined(STDPAR_ON) || defined(_OPENMP_ON)
#include <mutex>
#endif
#include <numeric>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "ltepop.h"
#include "nonthermal.h"
#include "packet.h"
#include "radfield.h"
#include "random.h"
#include "ratecoeff.h"
#include "rpkt.h"
#include "sn3d.h"
#include "stats.h"
#include "vpkt.h"

namespace {

// save to the macroatom_*.out file
constexpr bool LOG_MACROATOM = false;

std::fstream macroatom_file;

[[nodiscard]] auto get_sum_internal_down_same_exceptlast(const int uniquelevelindex) -> std::span<const double> {
  const auto ndowntrans = get_ndowntrans(uniquelevelindex);
  return std::span{globals::cellcache[cellcacheslotid].allmacroatomictransitions}.subspan(
      globals::alllevels.matransblock_start[uniquelevelindex] + get_ndowntrans(uniquelevelindex), ndowntrans - 1);
}

[[nodiscard]] auto get_sum_internal_up_same_exceptlast(const int uniquelevelindex) -> std::span<const double> {
  return std::span{globals::cellcache[cellcacheslotid].allmacroatomictransitions}.subspan(
      globals::alllevels.matransblock_start[uniquelevelindex] + (2 * get_ndowntrans(uniquelevelindex)),
      get_nuptrans(uniquelevelindex) - 1);
}

[[nodiscard]] auto get_sum_epstrans_rad_deexc_exceptlast(const int uniquelevelindex) -> std::span<const double> {
  return std::span{globals::cellcache[cellcacheslotid].allmacroatomictransitions}.subspan(
      globals::alllevels.matransblock_start[uniquelevelindex], get_ndowntrans(uniquelevelindex) - 1);
}

auto calculate_macroatom_transitionrates(const int nonemptymgi, const int element, const int ion, const int level,
                                         const int ionuniquelevelindexstart, const double t_mid,
                                         const globals::AllTransitions& alltrans) {
  // printout("Calculating transition rates for element %d ion %d level %d\n", element, ion, level);
  auto processrates = std::array<double, MA_ACTION_COUNT>{};

  const auto uniquelevelindex = ionuniquelevelindexstart + level;

  const auto transblock = std::span{globals::cellcache[cellcacheslotid].allmacroatomictransitions};
  const auto alllevels_matransblock_start = globals::alllevels.matransblock_start[uniquelevelindex];

  const auto T_e = grid::get_Te(nonemptymgi);
  const auto nne = grid::get_nne(nonemptymgi);
  const double epsilon_current = epsilon(uniquelevelindex);
  const double statweight = stat_weight(uniquelevelindex);
  const auto nnlevel = get_levelpop(nonemptymgi, uniquelevelindex);

  // Downward transitions within the current ionisation stage:
  // radiative/collisional deexcitation and internal downward jumps
  double sum_internal_down_same = 0.;
  double sum_raddeexc = 0.;
  double sum_coldeexc = 0.;
  const auto alltrans_startdown = get_alltrans_startdown(uniquelevelindex);
  const auto ndowntrans = get_ndowntrans(uniquelevelindex);

  const auto arr_sum_epstrans_rad_deexc = transblock.subspan(alllevels_matransblock_start, ndowntrans);
  const auto arr_sum_internal_down_same = transblock.subspan(alllevels_matransblock_start + ndowntrans, ndowntrans);
  for (int i = 0; i < ndowntrans; i++) {
    const auto alltransindex = alltrans_startdown + i;
    const int lower = alltrans.targetlevelindex[alltransindex];
    const auto A_ul = alltrans.einstein_A[alltransindex];
    const auto lower_uniquelevelindex = ionuniquelevelindexstart + lower;
    const double epsilon_target = epsilon(lower_uniquelevelindex);
    const double epsilon_trans = epsilon_current - epsilon_target;
    const auto lower_statweight = stat_weight(lower_uniquelevelindex);

    const double R = rad_deexcitation_ratecoeff(nonemptymgi, lower_uniquelevelindex, epsilon_trans, A_ul, statweight,
                                                lower_statweight, nnlevel, t_mid);
    const double C = col_deexcitation_ratecoeff(T_e, nne, epsilon_trans, statweight, lower_statweight, alltransindex);

    sum_raddeexc += R * epsilon_trans;
    sum_coldeexc += C * epsilon_trans;
    sum_internal_down_same += (R + C) * epsilon_target;

    arr_sum_epstrans_rad_deexc[i] = sum_raddeexc;
    arr_sum_internal_down_same[i] = sum_internal_down_same;

    // printout("checking downtrans %d to level %d: R %g, C %g, epsilon_trans %g\n",i,lower,R,C,epsilon_trans);
  }
  processrates[MA_ACTION_RADDEEXC] = sum_raddeexc;
  processrates[MA_ACTION_COLDEEXC] = sum_coldeexc;
  processrates[MA_ACTION_INTERNALDOWNSAME] = sum_internal_down_same;

  // Calculate sum for upward internal transitions
  // transitions within the current ionisation stage
  double sum_internal_up_same = 0.;
  const int nuptrans = get_nuptrans(uniquelevelindex);
  const auto arr_sum_internal_up_same = transblock.subspan(alllevels_matransblock_start + (2 * ndowntrans), nuptrans);
  const auto alltrans_startup = get_alltrans_startup(uniquelevelindex);
  for (int ii = 0; ii < nuptrans; ii++) {
    const auto alltransindex = alltrans_startup + ii;
    const auto upper = alltrans.targetlevelindex[alltransindex];
    const auto upper_uniquelevelindex = ionuniquelevelindexstart + upper;
    const double epsilon_trans = epsilon(upper_uniquelevelindex) - epsilon_current;
    const auto upper_statweight = stat_weight(upper_uniquelevelindex);

    const double R = rad_excitation_ratecoeff(nonemptymgi, upper_uniquelevelindex, upper_statweight,
                                              alltrans.einstein_A[alltransindex], epsilon_trans, nnlevel, statweight,
                                              alltransindex, t_mid);
    const double C = col_excitation_ratecoeff(T_e, nne, upper_statweight, alltransindex, epsilon_trans, statweight);
    const double NT = nonthermal::nt_excitation_ratecoeff(nonemptymgi, level, upper, alltransindex);

    sum_internal_up_same += (R + C + NT) * epsilon_current;
    arr_sum_internal_up_same[ii] = sum_internal_up_same;
  }
  processrates[MA_ACTION_INTERNALUPSAME] = sum_internal_up_same;

  assert_always(std::isfinite(processrates[MA_ACTION_INTERNALUPSAME]));

  // Downward transitions to lower ionisation stages:
  // radiative/collisional recombination and internal downward jumps
  // checks only if there is a lower ion, doesn't make sure that Z(ion)=Z(ion-1)+1
  double sum_internal_down_lower = 0.;
  double sum_radrecomb = 0.;
  double sum_colrecomb = 0.;
  if (ion > 0 && level <= get_maxrecombininglevel(element, ion)) {
    const int nlevels = get_nlevels_ionising(element, ion - 1);
    const auto lowerionuniquelevelindexstart = globals::elements[element].ions[ion - 1].uniquelevelindexstart;
    for (int lower = 0; lower < nlevels; lower++) {
      const double epsilon_target = epsilon(lowerionuniquelevelindexstart + lower);
      const double epsilon_trans = epsilon_current - epsilon_target;

      const double R = rad_recombination_ratecoeff(T_e, nne, element, ion, level, lower, nonemptymgi);
      const double C = col_recombination_ratecoeff(T_e, nne, element, ion, level, lower, epsilon_trans);

      sum_internal_down_lower += (R + C) * epsilon_target;

      sum_radrecomb += R * epsilon_trans;
      sum_colrecomb += C * epsilon_trans;
    }
  }
  processrates[MA_ACTION_INTERNALDOWNLOWER] = sum_internal_down_lower;
  processrates[MA_ACTION_RADRECOMB] = sum_radrecomb;
  processrates[MA_ACTION_COLRECOMB] = sum_colrecomb;

  // Transitions to higher ionisation stages
  double sum_up_highernt = 0.;
  double sum_up_higher = 0.;
  const int ionisinglevels = get_nlevels_ionising(element, ion);
  if (ion < get_nions(element) - 1 && level < ionisinglevels) {
    if (NT_ON) {
      sum_up_highernt = nonthermal::nt_ionisation_ratecoeff(nonemptymgi, element, ion) * epsilon_current;
    }

    const auto nphixstargets = get_nphixstargets(uniquelevelindex);
    for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
      const double epsilon_trans = get_phixs_threshold(uniquelevelindex, phixstargetindex);

      const double R = get_corrphotoioncoeff(element, ion, level, phixstargetindex, nonemptymgi);
      const double C = col_ionisation_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);

      sum_up_higher += (R + C) * epsilon_current;
    }
  }
  processrates[MA_ACTION_INTERNALUPHIGHERNT] = sum_up_highernt;
  processrates[MA_ACTION_INTERNALUPHIGHER] = sum_up_higher;

  return processrates;
}

// radiative deexcitation
void do_macroatom_raddeexcitation(Packet& pkt, const int ionuniquelevelindexstart, const int uniquelevelindex,
                                  const int activatingline, const double epsilon_current, const double totalrate) {
  // randomly select which line transitions occurs

  const double targetval = rng_uniform() * totalrate;
  const auto sum_epstrans_rad_deexc_exceptlast = get_sum_epstrans_rad_deexc_exceptlast(uniquelevelindex);

  // first sum_epstrans_rad_deexc[i] such that sum_epstrans_rad_deexc[i] > targetval
  const auto downtransindex = std::ranges::upper_bound(sum_epstrans_rad_deexc_exceptlast, targetval) -
                              sum_epstrans_rad_deexc_exceptlast.begin();
  const auto alltrans_startdown = get_alltrans_startdown(uniquelevelindex);

  const auto lineindex = globals::alltrans.lineindex[alltrans_startdown + downtransindex];

  if (lineindex == activatingline) {
    stats::increment(stats::COUNTER_RESONANCESCATTERINGS);
  }

  if constexpr (RECORD_LINESTAT) {
    atomicadd(globals::ecounter[lineindex], 1);
  }

  const auto uniquelevelindexlower =
      ionuniquelevelindexstart + globals::alltrans.targetlevelindex[alltrans_startdown + downtransindex];

  const double epsilon_trans = epsilon_current - epsilon(uniquelevelindexlower);

  const double oldnucmf = pkt.nu_cmf;
  pkt.nu_cmf = epsilon_trans / H;

  if (activatingline >= 0) {
    stats::increment((oldnucmf < pkt.nu_cmf) ? stats::COUNTER_UPSCATTER : stats::COUNTER_DOWNSCATTER);
  }

  stats::increment(stats::COUNTER_MA_STAT_DEACTIVATION_BB);
  stats::increment(stats::COUNTER_INTERACTIONS);

  // emit the rpkt in a random direction
  emit_rpkt(pkt);

  // the r-pkt can only interact with lines redder than the current one
  pkt.next_trans = lineindex + 1;
  pkt.emissiontype = lineindex;
  pkt.em_pos = pkt.pos;
  pkt.em_time = static_cast<float>(pkt.prop_time);
  pkt.nscatterings = 0;
}

// get the level index of the lower ionisation stage after a randomly selected radiative recombination and update
// counters
[[nodiscard]] auto do_macroatom_radrecomb(Packet& pkt, const int nonemptymgi, const int element, const int upperion,
                                          const int upperionlevel, const double rad_recomb) -> int {
  const auto T_e = grid::get_Te(nonemptymgi);
  const auto nne = grid::get_nne(nonemptymgi);
  const double epsilon_current = epsilon(element, upperion, upperionlevel);
  // Randomly select a continuum
  const double targetval = rng_uniform() * rad_recomb;
  double rate = 0;
  const int nlevels = get_nlevels_ionising(element, upperion - 1);
  int lowerionlevel = -1;
  for (int tmp_lowerionlevel = 0; tmp_lowerionlevel < nlevels; tmp_lowerionlevel++) {
    const double epsilon_trans = epsilon_current - epsilon(element, upperion - 1, tmp_lowerionlevel);
    const double R =
        rad_recombination_ratecoeff(T_e, nne, element, upperion, upperionlevel, tmp_lowerionlevel, nonemptymgi);

    rate += R * epsilon_trans;

    if (targetval < rate) {
      lowerionlevel = tmp_lowerionlevel;
      break;
    }
  }
  assert_always(lowerionlevel >= 0);

  // set the new state
  const int lowerion = upperion - 1;

  pkt.nu_cmf = select_continuum_nu(element, upperion - 1, lowerionlevel, upperionlevel, T_e);

  stats::increment(stats::COUNTER_MA_STAT_DEACTIVATION_FB);
  stats::increment(stats::COUNTER_INTERACTIONS);

  // Finally emit the packet into a randomly chosen direction, update the continuum opacity and set some flags
  emit_rpkt(pkt);

  pkt.next_trans = -1;  // continuum transition, no restrictions for further line interactions
  pkt.emissiontype = get_emtype_continuum(element, lowerion, lowerionlevel, upperionlevel);
  pkt.em_pos = pkt.pos;
  pkt.em_time = static_cast<float>(pkt.prop_time);
  pkt.nscatterings = 0;
  return lowerionlevel;
}

// get the level index of the upper ionisation stage after randomly-selected photoionisation or thermal collisional
// ionisation and update counters
[[nodiscard]] auto do_macroatom_ionisation(const int nonemptymgi, const int element, const int ion, const int level,
                                           const double epsilon_current, const double internal_up_higher) -> int {
  const auto T_e = grid::get_Te(nonemptymgi);
  const auto nne = grid::get_nne(nonemptymgi);

  // Randomly select the occurring transition
  const double targetrate = rng_uniform() * internal_up_higher;
  double rate = 0.;
  const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
  const int nphixstargets = get_nphixstargets(uniquelevelindex);
  for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
    const double epsilon_trans = get_phixs_threshold(element, ion, level, phixstargetindex);
    const double R = get_corrphotoioncoeff(element, ion, level, phixstargetindex, nonemptymgi);
    const double C = col_ionisation_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);
    rate += (R + C) * epsilon_current;
    if (rate > targetrate) {
      // set the macroatom's new state
      return get_phixsupperlevel(uniquelevelindex, phixstargetindex);
    }
  }

  assert_always(false);
  return -1;
}

constexpr auto gaunt_factor(const int ionstage) -> double {
  if (ionstage == 1) {
    return 0.1;
  }
  if (ionstage == 2) {
    return 0.2;
  }
  return 0.3;
}

}  // anonymous namespace

// handle activated macro atoms
__host__ __device__ void do_macroatom(Packet& pkt, const MacroAtomState& pktmastate) {
  const auto nonemptymgi = grid::get_propcell_nonemptymgi(pkt.where);
  assert_testmodeonly(nonemptymgi >= 0);
  const auto T_e = grid::get_Te(nonemptymgi);

  const double t_mid = globals::timesteps[globals::timestep].mid;

  // printout("[debug] do MA\n");

  const auto nne = grid::get_nne(nonemptymgi);

  assert_testmodeonly(grid::modelgrid[nonemptymgi].thick != 1);  // macroatom should not be used in thick cells

  // calculate occupation number for active MA level ////////////////////////////////////
  // general QUESTION: is it better to calculate the n_1 (later the n_ionstage and
  // U_ionstage) here where we need them or once in update_grid for each grid cell
  // not sure whether this reduces the number of calculations, as number of grid cells
  // is much larger than number of pellets (next question: connection to number of
  // photons)
  const int element = pktmastate.element;
  int ion = pktmastate.ion;
  int level = pktmastate.level;

  const int activatingline = pktmastate.activatingline;
  assert_testmodeonly(pkt.absorptiontype < 0 || activatingline < 0 || activatingline == pkt.absorptiontype);

  const int ion_in = ion;
  const int level_in = level;
  const double nu_cmf_in = pkt.nu_cmf;
  const double nu_rf_in = pkt.nu_rf;

  bool end_packet = false;
  while (!end_packet) {
    // Set this here to 1 to overcome problems in cells which have zero population
    // in some ionisation stage. This is possible because the dependence on the
    // originating levels population cancels out in the macroatom transition probabilities
    // which are based on detailed balance.

    assert_testmodeonly(ion >= 0);
    assert_testmodeonly(ion < get_nions(element));
    const auto ionuniquelevelindexstart = globals::elements[element].ions[ion].uniquelevelindexstart;
    const int uniquelevelindex = ionuniquelevelindexstart + level;

    const double epsilon_current = epsilon(uniquelevelindex);
    auto& processrates = globals::cellcache[cellcacheslotid].alllevels_maprocessrates[uniquelevelindex];
    {
#if (defined(STDPAR_ON) || defined(_OPENMP_ON)) && !defined(GPU_ON)
      const auto lock = std::lock_guard<std::mutex>(globals::mutex_cellcachemacroatom[uniquelevelindex]);
#endif

      assert_testmodeonly(globals::cellcache[cellcacheslotid].nonemptymgi == nonemptymgi);

      // If there are no precalculated rates available then calculate them
      if (processrates[0] < 0.) {
        processrates = calculate_macroatom_transitionrates(nonemptymgi, element, ion, level, ionuniquelevelindexstart,
                                                           t_mid, globals::alltrans);
      }
    }

    // select transition according to probabilities
    std::array<double, MA_ACTION_COUNT> cumulative_transitions{};
    std::partial_sum(processrates.cbegin(), processrates.cend(), cumulative_transitions.begin());

    const double randomrate = rng_uniform() * cumulative_transitions[MA_ACTION_COUNT - 1];

    // first cumulative_transitions[i] such that cumulative_transitions[i] > randomrate
    const auto selected_action = static_cast<int>(std::ranges::upper_bound(cumulative_transitions, randomrate) -
                                                  cumulative_transitions.cbegin());

    switch (selected_action) {
      case MA_ACTION_RADDEEXC: {
        // printout("[debug] do_ma:   radiative deexcitation\n");
        do_macroatom_raddeexcitation(pkt, ionuniquelevelindexstart, uniquelevelindex, activatingline, epsilon_current,
                                     processrates[MA_ACTION_RADDEEXC]);

        if constexpr (LOG_MACROATOM) {
          macroatom_file << std::format(
              "{:8d} {:14d} {:2d} {:12d} {:12d} {:9d} {:9d} {:9d} {:11.5e} {:11.5e} {:11.5e} {:11.5e}\n",
              globals::timestep, grid::get_mgi_of_nonemptymgi(nonemptymgi), get_atomicnumber(element),
              get_ionstage(element, ion_in), get_ionstage(element, ion), level_in, level, activatingline, nu_cmf_in,
              pkt.nu_cmf, nu_rf_in, pkt.nu_rf);
        }

        end_packet = true;
        break;
      }

      case MA_ACTION_COLDEEXC: {
        // collisional deexcitation of macro atom => convert the packet into a k-packet
        // printout("[debug] do_ma:   collisional deexcitation\n");

        stats::increment(stats::COUNTER_MA_STAT_DEACTIVATION_COLLDEEXC);
        stats::increment(stats::COUNTER_INTERACTIONS);

        pkt.type = TYPE_KPKT;
        end_packet = true;
        if constexpr (!DIRECT_COL_HEAT) {
          atomicadd(globals::colheatingestimator[nonemptymgi], pkt.e_cmf);
        }
        break;
      }

      case MA_ACTION_INTERNALDOWNSAME: {
        stats::increment(stats::COUNTER_INTERACTIONS);
        // Randomly select the occurring transition
        const double targetval = rng_uniform() * processrates[MA_ACTION_INTERNALDOWNSAME];

        // first sum_internal_down_same[i] such that sum_internal_down_same[i] > targetval
        const auto sum_internal_down_same_exceptlast = get_sum_internal_down_same_exceptlast(uniquelevelindex);
        const auto downtransindex = std::ranges::upper_bound(sum_internal_down_same_exceptlast, targetval) -
                                    sum_internal_down_same_exceptlast.begin();

        level = globals::alltrans.targetlevelindex[get_alltrans_startdown(uniquelevelindex) + downtransindex];

        break;
      }

      case MA_ACTION_RADRECOMB: {
        // Radiative recombination of MA: emitt a continuum-rpkt
        // printout("[debug] do_ma:   radiative recombination\n");
        // printout("[debug] do_ma:   element %d, ion %d, level %d\n", element, ion, level);

        level = do_macroatom_radrecomb(pkt, nonemptymgi, element, ion, level, processrates[MA_ACTION_RADRECOMB]);
        ion -= 1;
        end_packet = true;
        break;
      }

      case MA_ACTION_COLRECOMB: {
        // collisional recombination of macro atom => convert the packet into a k-packet
        // printout("[debug] do_ma:   collisonal recombination\n");
        stats::increment(stats::COUNTER_MA_STAT_DEACTIVATION_COLLRECOMB);
        stats::increment(stats::COUNTER_INTERACTIONS);

        pkt.type = TYPE_KPKT;
        end_packet = true;
        if constexpr (!DIRECT_COL_HEAT) {
          atomicadd(globals::colheatingestimator[nonemptymgi], pkt.e_cmf);
        }
        break;
      }

      case MA_ACTION_INTERNALDOWNLOWER: {
        // printout("[debug] do_ma:   internal downward jump to lower ionstage\n");
        stats::increment(stats::COUNTER_INTERACTIONS);
        stats::increment(stats::COUNTER_MA_STAT_INTERNALDOWNLOWER);

        // Randomly select the occurring transition
        const double targetrate = rng_uniform() * processrates[MA_ACTION_INTERNALDOWNLOWER];
        double rate = 0.;

        const int nlevels = get_nlevels_ionising(element, ion - 1);
        int lower = -1;
        const auto lowerionuniquelevelindexstart = globals::elements[element].ions[ion - 1].uniquelevelindexstart;
        for (int tmp_lower = 0; tmp_lower < nlevels; tmp_lower++) {
          const double epsilon_target = epsilon(lowerionuniquelevelindexstart + tmp_lower);
          const double epsilon_trans = epsilon_current - epsilon_target;
          const double R = rad_recombination_ratecoeff(T_e, nne, element, ion, level, tmp_lower, nonemptymgi);
          const double C = col_recombination_ratecoeff(T_e, nne, element, ion, level, tmp_lower, epsilon_trans);
          rate += (R + C) * epsilon_target;
          if (rate > targetrate) {
            lower = tmp_lower;
            break;
          }
        }
        assert_always(lower >= 0);

        ion--;
        level = lower;

        break;
      }

      case MA_ACTION_INTERNALUPSAME: {
        // printout("[debug] do_ma:   internal upward jump within current ionstage\n");
        stats::increment(stats::COUNTER_INTERACTIONS);

        // randomly select the occurring transition
        const auto sum_internal_up_same_exceptlast = get_sum_internal_up_same_exceptlast(uniquelevelindex);

        const double targetval = rng_uniform() * processrates[MA_ACTION_INTERNALUPSAME];

        // first sum_internal_up_same[i] such that sum_internal_up_same[i] > targetval
        const auto uptransindex = std::ranges::upper_bound(sum_internal_up_same_exceptlast, targetval) -
                                  sum_internal_up_same_exceptlast.begin();
        const auto alltrans_startup = get_alltrans_startup(uniquelevelindex);

        const int upper = globals::alltrans.targetlevelindex[alltrans_startup + uptransindex];

        level = upper;
        break;
      }

      case MA_ACTION_INTERNALUPHIGHER: {
        // printout("[debug] do_ma:   internal upward jump to next ionstage\n");
        stats::increment(stats::COUNTER_INTERACTIONS);

        stats::increment(stats::COUNTER_MA_STAT_INTERNALUPHIGHER);

        level = do_macroatom_ionisation(nonemptymgi, element, ion, level, epsilon_current,
                                        processrates[MA_ACTION_INTERNALUPHIGHER]);
        ion += 1;

        break;
      }

      case MA_ACTION_INTERNALUPHIGHERNT: {
        stats::increment(stats::COUNTER_INTERACTIONS);

        ion = nonthermal::nt_random_upperion(nonemptymgi, element, ion, false);
        level = 0;
        stats::increment(stats::COUNTER_MA_STAT_INTERNALUPHIGHERNT);

        break;
      }

      case MA_ACTION_COUNT: {
        printout("ERROR: Problem selecting MA_ACTION\n");
        std::abort();
      }

      default:
        if constexpr (TESTMODE) {
          printout("ERROR: Unknown macroatom selected_action type %d\n", selected_action);
          assert_testmodeonly(false);
        } else {
          __builtin_unreachable();
        }
    }
  }

  // TODO Luke: we should probably only do this if the packet has become a r-packet, otherwise we should set
  // trueemissiontype to EM_TYPE_NOTSET, but this method has already been published. If the difference is small for
  // nebular Type Ias then just fix it.
  if (pkt.trueemissiontype == EMTYPE_NOTSET) {
    pkt.trueemissiontype = pkt.emissiontype;
    pkt.trueem_pos = pkt.em_pos;
    pkt.trueem_time = pkt.em_time;
  }

  if (pkt.type == TYPE_RPKT) {
    if constexpr (VPKT_ON) {
      vpkt::call_estimators(pkt, TYPE_MA);
    }
  }
}

void macroatom_open_file(const int my_rank) {
  if constexpr (!LOG_MACROATOM) {
    return;
  }

  macroatom_file = fstream_required(std::format("macroatom_{:04d}.out", my_rank), std::ios::out | std::ios::trunc);

  macroatom_file << std::format("{:8s} {:14s} {:2s} {:12s} {:12s} {:9s} {:9s} {:9s} {:11s} {:11s} {:11s} {:11s}\n",
                                "timestep", "modelgridindex", "Z", "ionstage_in", "ionstage_out", "level_in",
                                "level_out", "activline", "nu_cmf_in", "nu_cmf_out", "nu_rf_in", "nu_rf_out");
}

void macroatom_close_file() {
  if (macroatom_file.is_open()) {
    macroatom_file.close();
  }
}

// radiative deexcitation rate: paperII 3.5.2
// multiply by upper level population to get a rate per second
auto rad_deexcitation_ratecoeff(const int nonemptymgi, const int lower_uniquelevelindex, const double epsilon_trans,
                                const float A_ul, const double upperstatweight, const double lowerstatweight,
                                const double nnlevelupper, const double t_current) -> double {
  const double n_l = get_levelpop(nonemptymgi, lower_uniquelevelindex);

  const double nu_trans = epsilon_trans / H;

  const double B_ul = CLIGHTSQUAREDOVERTWOH / std::pow(nu_trans, 3) * A_ul;
  const double B_lu = upperstatweight / lowerstatweight * B_ul;

  const double tau_sobolev = (B_lu * n_l - B_ul * nnlevelupper) * HCLIGHTOVERFOURPI * t_current;

  if (tau_sobolev > 1e-100) {
    const double beta = 1.0 / tau_sobolev * (-std::expm1(-tau_sobolev));
    const auto R = A_ul * beta;
    assert_testmodeonly(std::isfinite(R));
    return R;
  }
  return 0.;
}

// radiative excitation rate: paperII 3.5.2
// multiply by lower level population to get a rate per second
auto rad_excitation_ratecoeff(const int nonemptymgi, const int upper_uniquelevelindex, const double upper_statweight,
                              const double einstein_A, const double epsilon_trans, const double nnlevel_lower,
                              const double statweight_lower, const int alltransindex, const double t_current)
    -> double {
  const double n_u = get_levelpop(nonemptymgi, upper_uniquelevelindex);
  const double nu_trans = epsilon_trans / H;
  const double B_ul = CLIGHTSQUAREDOVERTWOH / std::pow(nu_trans, 3) * einstein_A;
  const double B_lu = upper_statweight / statweight_lower * B_ul;

  const double tau_sobolev = (B_lu * nnlevel_lower - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current;

  if (tau_sobolev > 1e-100) {
    const double beta = 1.0 / tau_sobolev * (-std::expm1(-tau_sobolev));

    const double R_over_J_nu = nnlevel_lower > 0. ? (B_lu - B_ul * n_u / nnlevel_lower) * beta : B_lu * beta;

    if (DETAILED_LINE_ESTIMATORS_ON && !globals::lte_iteration) {
      // check for a detailed line flux estimator to replace the binned/blackbody radiation field estimate
      if (const int jblueindex = radfield::get_Jblueindex(globals::alltrans.lineindex[alltransindex]);
          jblueindex >= 0) {
        return R_over_J_nu * radfield::get_Jb_lu(nonemptymgi, jblueindex);
      }
    }

    const double R = R_over_J_nu * radfield::radfield(nu_trans, nonemptymgi);

    assert_testmodeonly(R >= 0.);
    assert_testmodeonly(std::isfinite(R));
    return R;
  }

  return 0.;
}

// radiative recombination rate: paperII 3.5.2
// multiply by upper level population to get a rate per second
auto rad_recombination_ratecoeff(const float T_e, const float nne, const int element, const int upperion,
                                 const int upperionlevel, const int lowerionlevel, const int nonemptymgi) -> double {
  // it's faster to only check this condition outside this function than to check it for every level
  assert_testmodeonly(upperionlevel <= get_maxrecombininglevel(element, upperion));

  const int lowerion = upperion - 1;
  const auto lowerionlower_uniquelevelindex = get_uniquelevelindex(element, lowerion, lowerionlevel);
  const int nphixstargets = get_nphixstargets(lowerionlower_uniquelevelindex);
  for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
    if (get_phixsupperlevel(lowerionlower_uniquelevelindex, phixstargetindex) == upperionlevel) {
      const auto R_spont = nne * get_spontrecombcoeff(lowerionlower_uniquelevelindex, phixstargetindex, T_e);
      assert_testmodeonly(std::isfinite(R_spont));

      if constexpr (SEPARATE_STIMRECOMB) {
        const auto R_stimb = nne * get_stimrecombcoeff(element, lowerion, lowerionlevel, phixstargetindex, nonemptymgi);
        assert_testmodeonly(std::isfinite(R_stimb));
        return R_spont + R_stimb;
      }

      return R_spont;
    }
  }

  return 0.;
}

auto stim_recombination_ratecoeff(const float nne, const int element, const int upperion, const int upper,
                                  const int lower, const int nonemptymgi) -> double {
  if constexpr (SEPARATE_STIMRECOMB) {
    const int nphixstargets = get_nphixstargets(element, upperion - 1, lower);
    for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
      if (get_phixsupperlevel(element, upperion - 1, lower, phixstargetindex) == upper) {
        const double R = nne * get_stimrecombcoeff(element, upperion - 1, lower, phixstargetindex, nonemptymgi);
        return R;
      }
    }
  }

  return 0.;
}

// multiply by upper level population to get a rate per second
auto col_recombination_ratecoeff(const float T_e, const float nne, const int element, const int upperion,
                                 const int upper, const int lower, const double epsilon_trans) -> double {
  const auto lowerionlower_uniquelevelindex = get_uniquelevelindex(element, upperion - 1, lower);
  const int nphixstargets = get_nphixstargets(lowerionlower_uniquelevelindex);
  for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
    if (get_phixsupperlevel(lowerionlower_uniquelevelindex, phixstargetindex) == upper) {
      // Seaton approximation: Mihalas (1978), eq.5-79, p.134

      const double fac1 = epsilon_trans / KB / T_e;

      // select gaunt factor according to ionic charge
      const double g = gaunt_factor(get_ionstage(element, upperion - 1));

      const double sigma_bf = (get_phixs_table(lowerionlower_uniquelevelindex)[0] *
                               get_phixsprobability(lowerionlower_uniquelevelindex, phixstargetindex));

      const double sf = calculate_sahafact(element, upperion - 1, lower, upper, T_e, epsilon_trans);

      const double C = nne * nne * sf * 1.55e13 * pow(T_e, -0.5) * g * sigma_bf * exp(-fac1) / fac1;

      return C;
    }
  }

  return 0.;
}

// collisional ionisation rate: paperII 3.5.1
// multiply by lower level population to get a rate per second
auto col_ionisation_ratecoeff(const float T_e, const float nne, const int element, const int ion, const int lower,
                              const int phixstargetindex, const double epsilon_trans) -> double {
  assert_testmodeonly(phixstargetindex >= 0);
  assert_testmodeonly(phixstargetindex < get_nphixstargets(element, ion, lower));

  // Seaton approximation: Mihalas (1978), eq.5-79, p.134

  const double g = gaunt_factor(get_ionstage(element, ion));

  const double fac1 = epsilon_trans / KB / T_e;

  const double sigma_bf =
      get_phixs_table(element, ion, lower)[0] * get_phixsprobability(element, ion, lower, phixstargetindex);
  const double C = nne * 1.55e13 * pow(T_e, -0.5) * g * sigma_bf * exp(-fac1) / fac1;  // photoionisation at the edge

  // printout("[debug] col_ion: nne %g, T_e %g, g %g, epsilon_trans %g, sigma_bf %g\n",
  // nne,T_e,g,epsilon_trans,sigma_bf);
  assert_testmodeonly(std::isfinite(C));

  return C;
}

// multiply by upper level population to get a rate per second
auto col_deexcitation_ratecoeff(const float T_e, const float nne, const double epsilon_trans,
                                const double upperstatweight, const double lowerstatweight, const int alltransindex)
    -> double {
  const auto coll_strength = globals::alltrans.coll_str[alltransindex];
  if (coll_strength < 0) {
    if (!globals::alltrans.forbidden[alltransindex])  // alternative: (coll_strength > -1.5) i.e. to catch -1
    {
      // permitted E1 electric dipole transitions
      // collisional deexcitation: formula valid only for atoms!!!!!!!!!!!
      // Rutten script eq. 3.33. p.50
      // f = osc_strength(element,ion,upper,lower);
      // C = n_u * 2.16 * pow(fac1,-1.68) * pow(T_e,-1.5) *
      // stat_weight(element,ion,lower)/stat_weight(element,ion,upper)  * nne * f;
      const double trans_osc_strength = globals::alltrans.osc_strength[alltransindex];

      const double eoverkt = epsilon_trans / (KB * T_e);
      // Van-Regemorter formula, Mihalas (1978), eq.5-75, p.133
      constexpr double g_bar = 0.2;  // this should be read in from transitions data: it is 0.2 for transitions nl ->
                                     // n'l' and 0.7 for transitions nl -> nl'
      // test = 0.276 * exp(fac1) * gsl_sf_expint_E1(fac1);
      // crude approximation to the already crude Van-Regemorter formula

      // double test = 0.276 * exp(fac1) * (-EULERGAMMA - log(fac1));
      // double Gamma = (g_bar > test) ? g_bar : test;

      // optimisation
      const double gauntfac =
          (eoverkt > 0.33421) ? g_bar : 0.276 * std::exp(eoverkt) * (-EULERGAMMA - std::log(eoverkt));

      const double g_ratio = lowerstatweight / upperstatweight;

      return C_0 * 14.51039491 * nne * std::sqrt(T_e) * trans_osc_strength * std::pow(H_ionpot / epsilon_trans, 2) *
             eoverkt * g_ratio * gauntfac;
    }

    // forbidden transitions: magnetic dipole, electric quadropole...
    // could be Axelrod? or Maurer
    return nne * 8.629e-6 * 0.01 * lowerstatweight / std::sqrt(T_e);
  }

  // positive coll_str is treated as effective collision strength

  // from Osterbrock and Ferland, p51
  return nne * 8.629e-6 * static_cast<double>(coll_strength) / upperstatweight / std::sqrt(T_e);
}

// multiply by lower level population to get a rate per second
auto col_excitation_ratecoeff(const float T_e, const float nne, const double upperstatweight, const int alltransindex,
                              const double epsilon_trans, const double lowerstatweight) -> double {
  const auto coll_strength = globals::alltrans.coll_str[alltransindex];
  const double eoverkt = epsilon_trans / (KB * T_e);

  if (coll_strength < 0) {
    const bool forbidden = globals::alltrans.forbidden[alltransindex];
    if (!forbidden) {
      // alternative condition: (coll_strength > -1.5) i.e. to catch -1
      const double trans_osc_strength = globals::alltrans.osc_strength[alltransindex];
      // permitted E1 electric dipole transitions
      // collisional excitation: formula valid only for atoms!!!!!!!!!!!
      // Rutten script eq. 3.32. p.50
      // C = n_l * 2.16 * pow(eoverkt,-1.68) * pow(T_e,-1.5) * exp(-eoverkt) * nne *
      // osc_strength(element,ion,upper,lower);

      // Van-Regemorter formula, Mihalas (1978), eq.5-75, p.133
      constexpr double g_bar = 0.2;  // this should be read in from transitions data: it is 0.2 for transitions nl ->
                                     // n'l' and 0.7 for transitions nl -> nl'
      // test = 0.276 * std::exp(eoverkt) * gsl_sf_expint_E1(eoverkt);
      // crude approximation to the already crude Van-Regemorter formula
      const double exp_eoverkt = std::exp(eoverkt);

      const double test = 0.276 * exp_eoverkt * (-EULERGAMMA - std::log(eoverkt));
      const double Gamma = g_bar > test ? g_bar : test;
      return C_0 * nne * std::sqrt(T_e) * 14.51039491 * trans_osc_strength * pow(H_ionpot / epsilon_trans, 2) *
             eoverkt / exp_eoverkt * Gamma;
    }

    // alternative condition: (coll_strength > -3.5) to catch -2 or -3

    // forbidden transitions: magnetic dipole, electric quadropole...
    // Axelrod's approximation (thesis 1980)

    return nne * 8.629e-6 * 0.01 * std::exp(-eoverkt) * upperstatweight / std::sqrt(T_e);
  }

  // from Osterbrock and Ferland, p51
  return nne * 8.629e-6 * static_cast<double>(coll_strength) * std::exp(-eoverkt) / lowerstatweight / std::sqrt(T_e);
}

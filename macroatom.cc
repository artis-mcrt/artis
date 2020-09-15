#include <gsl/gsl_integration.h>
#include "sn3d.h"
#include "atomic.h"
#include "grid_init.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "rpkt.h"
#include "vectors.h"

// constant for van-Regemorter approximation.
#define C_0 (5.465e-11)

// save to the macroatom_*.out file
#define LOG_MACROATOM false

static FILE *macroatom_file;


static inline double get_individ_rad_deexc(int element, int ion, int level, int i)
{
  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_rad_deexc[i];
}


static inline double get_individ_internal_down_same(int element, int ion, int level, int i)
{
  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_down_same[i];
}


static inline double get_individ_internal_up_same(int element, int ion, int level, int i)
{
  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_up_same[i];
}


static void get_macroatom_transitionrates(
  const int modelgridindex, const int element, const int ion, const int level, const double t_mid, double *processrates,
  chlevels_struct *const chlevel)
{
  const float T_e = get_Te(modelgridindex);
  const float nne = get_nne(modelgridindex);
  const double epsilon_current = epsilon(element, ion, level);

  /// Downward transitions within the current ionisation stage:
  /// radiative/collisional deexcitation and internal downward jumps
  processrates[MA_ACTION_RADDEEXC] = 0.;
  processrates[MA_ACTION_COLDEEXC] = 0.;
  processrates[MA_ACTION_INTERNALDOWNSAME] = 0.;
  const int ndowntrans = get_ndowntrans(element, ion, level);
  for (int i = 0; i < ndowntrans; i++)
  {
    const int lineindex = elements[element].ions[ion].levels[level].downtrans_lineindicies[i];
    const int lower = linelist[lineindex].lowerlevelindex;
    const double epsilon_target = epsilon(element, ion, lower);
    const double epsilon_trans = epsilon_current - epsilon_target;

    const double R = rad_deexcitation_ratecoeff(modelgridindex, element, ion, level, lower, epsilon_trans, lineindex, t_mid);
    const double C = col_deexcitation_ratecoeff(T_e, nne, epsilon_trans, lineindex);

    const double individ_internal_down_same = (R + C) * epsilon_target;

    const double individ_rad_deexc = R * epsilon_trans;
    const double individ_col_deexc = C * epsilon_trans;

    chlevel->individ_rad_deexc[i] = individ_rad_deexc;
    chlevel->individ_internal_down_same[i] = individ_internal_down_same;

    processrates[MA_ACTION_RADDEEXC] += individ_rad_deexc;
    processrates[MA_ACTION_COLDEEXC] += individ_col_deexc;
    processrates[MA_ACTION_INTERNALDOWNSAME] += individ_internal_down_same;

    #ifdef DEBUG_ON
      if (debuglevel == 2)
        printout("checking downtrans %d to level %d: R %g, C %g, epsilon_trans %g\n",i,lower,R,C,epsilon_trans);
    #endif
  }

  /// Downward transitions to lower ionisation stages:
  /// radiative/collisional recombination and internal downward jumps
  processrates[MA_ACTION_RADRECOMB] = 0.;
  processrates[MA_ACTION_COLRECOMB] = 0.;
  processrates[MA_ACTION_INTERNALDOWNLOWER] = 0.;
  if (ion > 0 && level <= get_maxrecombininglevel(element, ion)) ///checks only if there is a lower ion, doesn't make sure that Z(ion)=Z(ion-1)+1
  {
    //nlevels = get_nlevels(element,ion-1);
    const int nlevels = get_ionisinglevels(element, ion - 1);
    //nlevels = get_ionisinglevels(element,ion-1);
    for (int lower = 0; lower < nlevels; lower++)
    {
      const double epsilon_target = epsilon(element, ion - 1, lower);
      const double epsilon_trans = epsilon_current - epsilon_target;

      const double R = rad_recombination_ratecoeff(T_e, nne, element, ion, level, lower, modelgridindex);
      //printout("rad recombination of element %d, ion %d, level %d, to lower level %d has rate %g\n",element,ion,level,lower,R);
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
  for (int i = 0; i < nuptrans; i++)
  {
    const int lineindex = elements[element].ions[ion].levels[level].uptrans_lineindicies[i];
    const int upper = linelist[lineindex].upperlevelindex;
    const double epsilon_trans = epsilon(element, ion, upper) - epsilon_current;

    double individ_internal_up_same = 0.;
    const double R = rad_excitation_ratecoeff(modelgridindex, element, ion, level, upper, epsilon_trans, lineindex, t_mid);
    const double C = col_excitation_ratecoeff(T_e, nne, lineindex, epsilon_trans);
    // const double NT = nt_excitation_ratecoeff(modelgridindex, element, ion, level, upper, epsilon_trans, lineindex);
    const double NT = 0.;

    individ_internal_up_same = (R + C + NT) * epsilon_current;

    chlevel->individ_internal_up_same[i] = individ_internal_up_same;

    processrates[MA_ACTION_INTERNALUPSAME] += individ_internal_up_same;

    #ifdef DEBUG_ON
      if (debuglevel == 2)
        printout("checking uptrans %d to level %d: R %g, C %g, epsilon_trans %g\n",i,upper,R,C,epsilon_trans);
    #endif
  }
  if (!isfinite(processrates[MA_ACTION_INTERNALUPSAME]))
    printout("fatal: internal_up_same has nan contribution\n");

  /// Transitions to higher ionisation stages
  processrates[MA_ACTION_INTERNALUPHIGHERNT] = 0.;
  processrates[MA_ACTION_INTERNALUPHIGHER] = 0.;
  const int ionisinglevels = get_ionisinglevels(element,ion);
  if (ion < get_nions(element) - 1 && level < ionisinglevels)
  {
    if (NT_ON)
    {
      processrates[MA_ACTION_INTERNALUPHIGHERNT] = nt_ionization_ratecoeff(modelgridindex, element, ion) * epsilon_current;
    }

    for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++)
    {
      const double epsilon_trans = get_phixs_threshold(element, ion, level, phixstargetindex);

      const double R = get_corrphotoioncoeff(element, ion, level, phixstargetindex, modelgridindex);
      const double C = col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);

      processrates[MA_ACTION_INTERNALUPHIGHER] += (R + C) * epsilon_current;
    }
  }
}


static void do_macroatom_raddeexcitation(
  PKT *pkt_ptr, const int modelgridindex, const int element, const int ion, const int level, const double rad_deexc,
  const double total_transitions, const int activatingline)
{
  ///radiative deexcitation of MA: emitt rpkt
  ///randomly select which line transitions occurs
  const double zrand = gsl_rng_uniform(rng);
  double rate = 0.;
  int linelistindex = -99;
  const int ndowntrans = get_ndowntrans(element, ion, level);
  for (int i = 0; i < ndowntrans; i++)
  {
    rate += get_individ_rad_deexc(element, ion, level, i);
    if (zrand * rad_deexc < rate)
    {
      linelistindex = elements[element].ions[ion].levels[level].downtrans_lineindicies[i];
      break;
    }
  }
  assert(linelistindex >= 0);
  #ifdef RECORD_LINESTAT
    if (tid == 0) ecounter[linelistindex]++;    /// This way we will only record line statistics from OMP-thread 0
                                                /// With an atomic pragma or a thread-private structure with subsequent
                                                /// reduction this could be extended to all threads. However, I'm not
                                                /// sure if this is worth the additional computational expenses.
  #endif
  const int lower = linelist[linelistindex].lowerlevelindex;
  #ifdef DEBUG_ON
    if (debuglevel == 2) printout("[debug] do_ma:   jump to level %d\n", lower);
  #endif
  const double epsilon_trans = epsilon(element, ion, level) - epsilon(element, ion, lower);

  double oldnucmf;
  if (pkt_ptr->last_event == 1)
    oldnucmf = pkt_ptr->nu_cmf;
  pkt_ptr->nu_cmf = epsilon_trans / H;
  //if (tid == 0)
  //{
  if (pkt_ptr->last_event == 1)
  {
    if (oldnucmf < pkt_ptr->nu_cmf)
      upscatter++;
    else
      downscatter++;
  }
  //}

  #ifdef DEBUG_ON
    if (!isfinite(pkt_ptr->nu_cmf))
    {
      printout("[fatal] rad deexcitation of MA: selected frequency not finite ... abort\n");
      abort();
    }
    if (linelistindex < 0)
    {
      printout("[fatal] problem in selecting radiative downward transition of MA zrand %g, rate %g, rad_deexc %g, ndowntrans %d\n", zrand, rate, rad_deexc, ndowntrans);
      printout("[fatal] total_transitions %g, element %d, ion %d, level %d\n", total_transitions, element, ion, level);
      abort();
    }
    if (debuglevel == 2)
      printout("[debug] do_ma: calculate_kappa_rpkt_cont after MA deactivation\n");
    //if (tid == 0) ma_stat_deactivation_bb++;
    ma_stat_deactivation_bb++;
    pkt_ptr->interactions += 1;
    pkt_ptr->last_event = 0;
  #endif

  /// Emit the rpkt in a random direction
  emitt_rpkt(pkt_ptr);

  if (linelistindex == activatingline)
    resonancescatterings++;
  else
  {
    calculate_kappa_rpkt_cont(pkt_ptr, modelgridindex);
  }

  /// NB: the r-pkt can only interact with lines redder than the current one
  pkt_ptr->next_trans = linelistindex + 1;
  pkt_ptr->emissiontype = linelistindex;
  vec_copy(pkt_ptr->em_pos, pkt_ptr->pos);
  pkt_ptr->em_time = pkt_ptr->prop_time;
  pkt_ptr->nscatterings = 0;
  //printout("next possible line encounter %d\n",pkt_ptr->next_trans);
  #ifndef FORCE_LTE
    //matotem[pkt_ptr->where] += pkt_ptr->e_cmf;
  #endif
}


static void do_macroatom_radrecomb(
  PKT *pkt_ptr, const int modelgridindex, const int element, int *ion, int *level,
  const double rad_recomb)
{
  const float T_e = get_Te(modelgridindex);
  const float nne = get_nne(modelgridindex);
  const double epsilon_current = epsilon(element, *ion, *level);
  const int upperion = *ion;
  const int upperionlevel = *level;
  /// Randomly select a continuum
  double zrand = gsl_rng_uniform(rng);
  double rate = 0;
  const int nlevels = get_ionisinglevels(element, upperion - 1);
  int lower = 0;
  for (lower = 0; lower < nlevels; lower++)
  {
    const double epsilon_trans = epsilon_current - epsilon(element, upperion - 1, lower);
    const double R = rad_recombination_ratecoeff(T_e, nne, element, upperion, upperionlevel, lower, modelgridindex);

    rate += R * epsilon_trans;
    #ifdef DEBUG_ON
      if (debuglevel == 2)
      {
        printout("[debug] do_ma:   R %g, deltae %g\n",R,(epsilon(element, upperion, upperionlevel) - epsilon(element, upperion - 1, lower)));
        printout("[debug] do_ma:   rate to level %d of ion %d = %g\n", lower, upperion - 1, rate);
        printout("[debug] do_ma:   zrand*rad_recomb = %g\n", zrand * rad_recomb);
      }
    #endif
    if (zrand * rad_recomb < rate)
    {
      break;
    }
  }
  if (zrand * rad_recomb >= rate)
  {
    printout("%s: From Z=%d ionstage %d level %d, could not select lower level to recombine to. zrand %g * rad_recomb %g >= rate %g",
             __func__, get_element(element), get_ionstage(element, *ion), *level, zrand, rad_recomb, rate);
    abort();
  }

  /// set the new state
  *ion = upperion - 1;
  *level = lower;

  pkt_ptr->nu_cmf = select_continuum_nu(element, upperion - 1, lower, upperionlevel, T_e);

  #ifndef FORCE_LTE
    //mabfcount[pkt_ptr->where] += pkt_ptr->e_cmf;
    //mabfcount_thermal[pkt_ptr->where] += pkt_ptr->e_cmf*(1-nu_threshold/pkt_ptr->nu_cmf);
    //matotem[pkt_ptr->where] += pkt_ptr->e_cmf;
  #endif

/*      if (element == 6)
  {
    //printout("%g, %g, %g\n",pkt_ptr->e_cmf,nu_threshold,pkt_ptr->e_cmf/nu_threshold/H);
    cell[pkt_ptr->where].radrecomb[ion-1] += pkt_ptr->e_cmf/pkt_ptr->nu_cmf/H;
  }*/

  #ifdef DEBUG_ON
    if (debuglevel == 2)
    {
      printout("%s: From Z=%d ionstage %d, recombining to ionstage %d level %d\n",
               __func__, get_element(element), get_ionstage(element, *ion + 1), get_ionstage(element, *ion), lower);
      printout("[debug] do_ma:   pkt_ptr->nu_cmf %g\n",pkt_ptr->nu_cmf);
    }
    if (!isfinite(pkt_ptr->nu_cmf))
    {
      printout("[fatal] rad recombination of MA: selected frequency not finite ... abort\n");
      abort();
    }
    //if (tid == 0) ma_stat_deactivation_fb++;
    ma_stat_deactivation_fb++;
    pkt_ptr->interactions += 1;
    pkt_ptr->last_event = 2;
    if (debuglevel == 2) printout("[debug] do_ma: calculate_kappa_rpkt_cont after MA recombination\n");
  #endif

  /// Finally emit the packet into a randomly chosen direction, update the continuum opacity and set some flags
  emitt_rpkt(pkt_ptr);

  #if (TRACK_ION_STATS)
  increment_ion_stats(modelgridindex, element, upperion, ION_COUNTER_RADRECOMB_MACROATOM, pkt_ptr->e_cmf / H / pkt_ptr->nu_cmf);

  const double escape_prob = get_rpkt_escape_prob(pkt_ptr, t_current);

  increment_ion_stats(modelgridindex, element, upperion, ION_COUNTER_RADRECOMB_ESCAPED, pkt_ptr->e_cmf / H / pkt_ptr->nu_cmf * escape_prob);
  #endif

  calculate_kappa_rpkt_cont(pkt_ptr, modelgridindex);

  pkt_ptr->next_trans = 0;       /// continuum transition, no restrictions for further line interactions
  pkt_ptr->emissiontype = get_continuumindex(element, *ion, lower, upperionlevel);
  vec_copy(pkt_ptr->em_pos, pkt_ptr->pos);
  pkt_ptr->em_time = pkt_ptr->prop_time;
  pkt_ptr->nscatterings = 0;
}


static void do_macroatom_ionisation(
  const int modelgridindex, const int element, int *ion, int *level,
  const double epsilon_current, const double internal_up_higher)
{
  const float T_e = get_Te(modelgridindex);
  const float nne = get_nne(modelgridindex);

  int upper;
  /// Randomly select the occuring transition
  const double zrand = gsl_rng_uniform(rng);
  double rate = 0.;
  for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, *ion, *level); phixstargetindex++)
  {
    upper = get_phixsupperlevel(element, *ion, *level, phixstargetindex);
    const double epsilon_trans = epsilon(element, *ion + 1, upper) - epsilon_current;
    const double R = get_corrphotoioncoeff(element, *ion, *level, phixstargetindex, modelgridindex);
    const double C = col_ionization_ratecoeff(T_e, nne, element, *ion, *level, phixstargetindex, epsilon_trans);
    rate += (R + C) * epsilon_current;
    if (zrand * internal_up_higher < rate)
      break;
  }
  if (zrand * internal_up_higher >= rate)
  {
    printout("%s: From Z=%d ionstage %d level %d, could not select upper level to ionise to. zrand %g * internal_up_higher %g >= rate %g\n",
             __func__, get_element(element), get_ionstage(element, *ion), *level, zrand, internal_up_higher, rate);
    abort();
  }

  /// and set the macroatom's new state
  *ion += 1;
  *level = upper;
}


void do_macroatom(PKT *pkt_ptr, const int timestep)
/// Material for handling activated macro atoms.
{
  const double t_mid = time_step[timestep].mid;

  //printout("[debug] do MA\n");

  const int cellindex = pkt_ptr->where;
  const int modelgridindex = cell[cellindex].modelgridindex;
  const float T_e = get_Te(modelgridindex);
  const float nne = get_nne(modelgridindex);

  assert(modelgrid[modelgridindex].thick != 1); // macroatom should not be used in thick cells

  /// calculate occupation number for active MA level ////////////////////////////////////
  /// general QUESTION: is it better to calculate the n_1 (later the n_ionstage and
  /// U_ionstage) here where we need them or once in update_grid for each grid cell
  /// not sure whether this reduces the number of calculations, as number of grid cells
  /// is much larger than number of pellets (next question: connection to number of
  /// photons)
  const int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;

  const int activatingline = mastate[tid].activatingline;
  if (pkt_ptr->absorptiontype > 0 && activatingline > 0 && activatingline != pkt_ptr->absorptiontype)
  {
    printout("error: mismatched absorptiontype %d != activatingline = %d pkt last_event %d emissiontype %d\n",
             pkt_ptr->absorptiontype, activatingline, pkt_ptr->last_event, pkt_ptr->emissiontype);
  }

  const int ion_in = ion;
  const int level_in = level;
  const double nu_cmf_in = pkt_ptr->nu_cmf;
  const double nu_rf_in = pkt_ptr->nu_rf;

  #if (TRACK_ION_STATS)
  increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYIN_TOTAL, pkt_ptr->e_cmf);
  #endif

  /// dummy-initialize these to nonsense values, if something goes wrong with the real
  /// initialization we should see errors

  // debuglevel = 2;
  #ifdef DEBUG_ON
    if (debuglevel == 2)
      printout("[debug] =============entering do_ma\n");
    int jumps = 0;
    int jump = -99;
  #endif

  bool end_packet = false;
  while (!end_packet)
  {
    //ionisinglevels = get_ionisinglevels(element,ion);

    /// Set this here to 1 to overcome problems in cells which have zero population
    /// in some ionisation stage. This is possible because the dependence on the
    /// originating levels population cancels out in the macroatom transition probabilities
    /// which are based on detailed balance.

    #ifdef DEBUG_ON
      const int ionstage = get_ionstage(element,ion);
      // if (Z == 26 && ionstage == 1) debuglevel = 2000;

      if (debuglevel == 2000)
        printout("[debug] %s Z=%d ionstage %d level %d, jumps %d\n", __func__, get_element(element), ionstage, level, jumps);

      assert(ion >= 0);
      assert(ion < get_nions(element));
    #endif

    const double epsilon_current = epsilon(element, ion, level);
    const int ndowntrans = get_ndowntrans(element, ion, level);
    const int nuptrans = get_nuptrans(element, ion, level);
    //nlevels_nextion  ///not needed as long we only ionise to the ground state
    //nlevels_lowerion ///only needed if level = 0, this won't happen too often

    #ifdef DEBUG_ON
      if (debuglevel == 2) printout("[debug] do_ma: ndowntrans %d, nuptrans %d\n",ndowntrans,nuptrans);
    #endif

    assert(cellhistory[tid].cellnumber == modelgridindex);

    double *processrates = cellhistory[tid].chelements[element].chions[ion].chlevels[level].processrates;

    /// If there are no precalculated rates available then calculate them
    if (processrates[MA_ACTION_COLDEEXC] < 0)
      get_macroatom_transitionrates(
        modelgridindex, element, ion, level, t_mid, processrates,
        &cellhistory[tid].chelements[element].chions[ion].chlevels[level]);

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
    for (int action = 0; action < MA_ACTION_COUNT; action++)
    {
      total_transitions += processrates[action];
    }

    enum ma_action selected_action;
    double zrand = gsl_rng_uniform(rng);
    //printout("zrand %g\n",zrand);
    const double randomrate = zrand * total_transitions;
    double rate = 0.;
    for (int action = 0; action < MA_ACTION_COUNT; action++)
    {
      rate += processrates[action];
      if (rate > randomrate)
      {
        selected_action = (enum ma_action)(action);
        break;
      }
    }

    if (rate <= randomrate)
    {
      printout("[fatal] do_ma: problem with random numbers .. abort\n");
      //printout("[debug]    rad_down %g, col_down %g, internal_down %g, internal_up %g\n",rad_down,col_down,internal_down,internal_up);
      printout("[debug] do_ma: element %d, ion %d, level %d, total_transitions %g\n",element,ion,level,total_transitions);
      printout("[debug]    col_deexc %g, col_recomb %g, rad_deexc %g, rad_recomb %g\n",processrates[MA_ACTION_COLDEEXC],processrates[MA_ACTION_COLRECOMB],processrates[MA_ACTION_RADDEEXC],processrates[MA_ACTION_RADRECOMB]);
      printout("[debug]    internal_down_same %g, internal_down_lower %g\n",processrates[MA_ACTION_INTERNALDOWNSAME],processrates[MA_ACTION_INTERNALDOWNLOWER]);
      printout("[debug]    internal_up_same %g, internal_up_higher %g\n",processrates[MA_ACTION_INTERNALUPSAME],processrates[MA_ACTION_INTERNALUPHIGHER]);
      printout("[debug]    zrand %g\n",zrand);
      printout("[debug]    jumps %d\n",jumps);
      printout("[debug]    pkt_ptr->number %d\n",pkt_ptr->number);

      debuglevel = 777;

      if (ion > 0) ///checks only if there is a lower ion, doesn't make sure that Z(ion)=Z(ion-1)+1
      {
        printout("[debug]    check recombination\n");
        //nlevels = get_nlevels(element,ion-1);
        const int nlevels = get_ionisinglevels(element, ion - 1);
        //nlevels = get_ionisinglevels(element,ion-1);
        for (int lower = 0; lower < nlevels; lower++)
        {
          const double epsilon_target = epsilon(element, ion - 1,lower);
          const double epsilon_trans = epsilon_current - epsilon_target;
          const double R = rad_recombination_ratecoeff(T_e, nne, element, ion, level, lower, modelgridindex);
          const double C = col_recombination_ratecoeff(modelgridindex, element, ion, level, lower, epsilon_trans);
          printout("[debug]    recombination to ion %d, level %d, epsilon_target %g, epsilon_trans %g, R %g, C %g\n",ion-1,lower,epsilon_target,epsilon_trans,R,C);
        }
      }

      printout("[debug]    check deexcitation\n");
      printout("[debug]    ndowntrans %d %d\n", ndowntrans, get_ndowntrans(element, ion, level));
      for (int i = 0; i < ndowntrans; i++)
      {
        const int lineindex = elements[element].ions[ion].levels[level].downtrans_lineindicies[i];
        const int lower = linelist[lineindex].lowerlevelindex;
        const double epsilon_trans = epsilon_current - epsilon(element, ion, lower);
        const double R = rad_deexcitation_ratecoeff(modelgridindex, element, ion, level, lower, epsilon_trans, lineindex, t_mid);
        const double C = col_deexcitation_ratecoeff(T_e, nne, epsilon_trans, lineindex);
        printout("[debug]    deexcitation to level %d, epsilon_trans %g, epsilon_trans %g, R %g, C %g\n",lower,epsilon_trans,epsilon_trans,R,C);
      }

      printout("[debug]    check excitation\n");
      printout("[debug]    nuptrans %d %d\n",nuptrans,get_nuptrans(element, ion, level));
      for (int i = 0; i < nuptrans; i++)
      {
        const int lineindex = elements[element].ions[ion].levels[level].uptrans_lineindicies[i];
        const int upper = linelist[lineindex].upperlevelindex;
        const double epsilon_trans = epsilon(element, ion, upper) - epsilon_current;
        const double R = rad_excitation_ratecoeff(modelgridindex, element, ion, level, upper, epsilon_trans, lineindex, t_mid);
        const double C = col_excitation_ratecoeff(T_e, nne, lineindex, epsilon_trans);
        printout("[debug]    excitation to level %d, epsilon_trans %g, R %g, C %g\n",upper,epsilon_trans,R,C);
      }

      if (ion < get_nions(element) - 1)  //&& get_ionstage(element,ion) < get_element(element)+1)
      {
        printout("[debug]    check ionisation\n");
        for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++)
        {
          const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
          const double epsilon_trans = epsilon(element, ion + 1, upper) - epsilon_current;
          const double R = get_corrphotoioncoeff(element, ion, level, phixstargetindex, modelgridindex);
          const double C = col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);
          printout("[debug]    ionisation to ion %d, level %d, epsilon_trans %g, R %g, C %g\n", ion + 1, upper, epsilon_trans, R, C);
          break;
        }
      }

      abort();
    }

    #ifdef DEBUG_ON
      if (debuglevel == 2)
      {
        printout("[debug] do_ma: element %d, ion %d, level %d\n",element,ion,level);
        printout(
          "[debug] do_ma:   rad_deexc %g, col_deexc %g, internal_down_same %g, rad_recomb %g, col_recomb %g, internal_down_lower %g, internal_up_same %g, internal_up_higher %g\n",
          processrates[MA_ACTION_RADDEEXC],processrates[MA_ACTION_COLDEEXC],processrates[MA_ACTION_INTERNALDOWNSAME],processrates[MA_ACTION_RADRECOMB],processrates[MA_ACTION_COLRECOMB],
          processrates[MA_ACTION_INTERNALDOWNLOWER], processrates[MA_ACTION_INTERNALUPSAME],processrates[MA_ACTION_INTERNALUPHIGHER]);
      }
      if (total_transitions <= 0.)
      {
        printout("[debug] do_ma: element %d, ion %d, level %d, total_transitions = %g\n",element,ion,level,total_transitions);
        printout("[debug]    col_deexc %g, col_recomb %g, rad_deexc %g, rad_recomb %g\n",processrates[MA_ACTION_COLDEEXC],processrates[MA_ACTION_COLRECOMB],processrates[MA_ACTION_RADDEEXC],processrates[MA_ACTION_RADRECOMB]);
        printout("[debug]    internal_down_same %g, internal_down_lower %g\n",processrates[MA_ACTION_INTERNALDOWNSAME],processrates[MA_ACTION_INTERNALDOWNLOWER]);
        printout("[debug]    internal_up_same %g, internal_up_higher %g\n",processrates[MA_ACTION_INTERNALUPSAME],processrates[MA_ACTION_INTERNALUPHIGHER]);
        printout("[debug]    zrand %g\n",zrand);
        printout("[debug]    jumps %d, jump %d\n",jumps,jump);
        printout("[debug]    pkt_ptr->number %d, pkt_ptr->where %d\n",pkt_ptr->number,cellindex);
        printout("[debug]    groundlevelpop of current ion in current cell %g\n",modelgrid[cell[cellindex].modelgridindex].composition[element].groundlevelpop[ion]);

        double R = 0.0;
        double C = 0.0;
        for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++)
        {
          const int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
          const double epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
          R += get_corrphotoioncoeff(element,ion,level,phixstargetindex,modelgridindex);
          C += col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);
          printout("epsilon_current %g, epsilon_trans %g, photion %g, colion %g, internal_up_higher %g\n",epsilon_current,epsilon_trans,R,C,(R + C) * epsilon_current);
        }

        const float T_R = get_TR(modelgridindex);
        const float W = get_W(modelgridindex);
        printout("modelgridindex %d, T_R %g, T_e %g, W %g, T_J %g\n",modelgridindex,T_R,get_Te(modelgridindex),W,get_TJ(modelgridindex));
        #if (!NO_LUT_PHOTOION)
          const double gammacorr = W * interpolate_corrphotoioncoeff(element,ion,level,0,T_R);
          const int index_in_groundlevelcontestimor = elements[element].ions[ion].levels[level].closestgroundlevelcont;
          const double renorm  = corrphotoionrenorm[modelgridindex * nelements * maxion + index_in_groundlevelcontestimor];
          printout("gammacorr %g, index %d, renorm %g, total %g\n",gammacorr,index_in_groundlevelcontestimor,renorm,gammacorr*renorm);
        #endif

        //abort();
      }
    #endif

    switch (selected_action)
    {
      case MA_ACTION_RADDEEXC:
      {
        #ifdef DEBUG_ON
        // if (debuglevel == 2)
        // {
        //   printout("[debug] do_ma:   radiative deexcitation\n");
        //   printout("[debug] do_ma:   jumps = %d\n",jumps);
        // }
        #endif

        do_macroatom_raddeexcitation(pkt_ptr, modelgridindex, element, ion, level, processrates[MA_ACTION_RADDEEXC], total_transitions, activatingline);

        #if (TRACK_ION_STATS)
        increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_RADDEEXC, pkt_ptr->e_cmf);

        increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_BOUNDBOUND_MACROATOM, pkt_ptr->e_cmf / H / pkt_ptr->nu_cmf);

        increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_TOTAL, pkt_ptr->e_cmf);
        #endif

        if (LOG_MACROATOM)
        {
          fprintf(macroatom_file, "%8d %14d %2d %12d %12d %9d %9d %9d %11.5e %11.5e %11.5e %11.5e %9d\n",
                  timestep, modelgridindex, get_element(element), get_ionstage(element, ion_in), get_ionstage(element, ion),
                  level_in, level, activatingline, nu_cmf_in, pkt_ptr->nu_cmf, nu_rf_in, pkt_ptr->nu_rf, jumps);
        }

        end_packet = true;
        break;
      }

      case MA_ACTION_COLDEEXC:
      {
        ///collisional deexcitation of macro atom => convert the packet into a k-packet
        #ifdef DEBUG_ON
          if (debuglevel == 2)
          {
            printout("[debug] do_ma:   collisonal deexcitation\n");
            printout("[debug] do_ma: jumps = %d\n", jumps);
          }
          //if (tid == 0) ma_stat_deactivation_colldeexc++;
          ma_stat_deactivation_colldeexc++;
          pkt_ptr->interactions += 1;
          pkt_ptr->last_event = 10;
        #endif

        #if (TRACK_ION_STATS)
        increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_COLLDEEXC, pkt_ptr->e_cmf);
        increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_TOTAL, pkt_ptr->e_cmf);
        #endif

        pkt_ptr->type = TYPE_KPKT;
        end_packet = true;
        #ifndef FORCE_LTE
          //matotem[pkt_ptr->where] += pkt_ptr->e_cmf;
          #ifdef _OPENMP
            #pragma omp atomic
          #endif
          colheatingestimator[modelgridindex] += pkt_ptr->e_cmf;
        #endif
        break;
      }

      case MA_ACTION_INTERNALDOWNSAME:
      {
        #ifdef DEBUG_ON
          if (debuglevel == 2)
            printout("[debug] do_ma:   internal downward jump within current ionstage\n");
          pkt_ptr->interactions += 1;
          jumps++;
          jump = 0;
        #endif

        /// Randomly select the occuring transition
        zrand = gsl_rng_uniform(rng);
        int lower = -99;
        double rate = 0.;
        for (int i = 0; i < ndowntrans; i++)
        {
          rate += get_individ_internal_down_same(element, ion, level, i);
          if (zrand * processrates[MA_ACTION_INTERNALDOWNSAME] < rate)
          {
            const int lineindex = elements[element].ions[ion].levels[level].downtrans_lineindicies[i];
            lower = linelist[lineindex].lowerlevelindex;
            break;
          }
        }

        level = lower;

        #ifdef DEBUG_ON
          if (debuglevel == 2)
            printout("[debug] do_ma:   to level %d\n", lower);
          if (get_ionstage(element,ion) == 0 && lower == 0)
          {
            printout("internal downward transition to ground level occured ... abort\n");
            printout("element %d, ion %d, level %d, lower %d\n", element, ion, level, lower);
            printout("Z %d, ionstage %d, energy %g\n",
                     get_element(element), get_ionstage(element,ion), elements[element].ions[ion].levels[lower].epsilon);
            printout("[debug] do_ma:   internal downward jump within current ionstage\n");
            abort();
          }
        #endif
        break;
      }

      case MA_ACTION_RADRECOMB:
      {
        /// Radiative recombination of MA: emitt a continuum-rpkt
        #ifdef DEBUG_ON
          // if (debuglevel == 2)
          // {
          //   printout("[debug] do_ma:   radiative recombination\n");
          //   printout("[debug] do_ma:   jumps = %d\n",jumps);
          //   printout("[debug] do_ma:   element %d, ion %d, level %d\n",element,ion,level);
          // }
        #endif

        #if (TRACK_ION_STATS)
        // increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_RADRECOMB, pkt_ptr->e_cmf);
        // increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_TOTAL, pkt_ptr->e_cmf);
        #endif

        do_macroatom_radrecomb(pkt_ptr, modelgridindex, element, &ion, &level, processrates[MA_ACTION_RADRECOMB]);
        end_packet = true;
        break;
      }

      case MA_ACTION_COLRECOMB:
      {
        ///collisional recombination of macro atom => convert the packet into a k-packet
        #ifdef DEBUG_ON
          if (debuglevel == 2)
          {
            printout("[debug] do_ma:   collisonal recombination\n");
            printout("[debug] do_ma: jumps = %d\n",jumps);
          }
          //if (tid == 0) ma_stat_deactivation_collrecomb++;
          ma_stat_deactivation_collrecomb++;
          pkt_ptr->interactions += 1;
          pkt_ptr->last_event = 11;
        #endif

        #if (TRACK_ION_STATS)
        increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_COLLRECOMB, pkt_ptr->e_cmf);
        increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_TOTAL, pkt_ptr->e_cmf);
        #endif

        pkt_ptr->type = TYPE_KPKT;
        end_packet = true;
        #ifndef FORCE_LTE
          //matotem[pkt_ptr->where] += pkt_ptr->e_cmf;
          #ifdef _OPENMP
            #pragma omp atomic
          #endif
          colheatingestimator[modelgridindex] += pkt_ptr->e_cmf;
        #endif
        break;
      }

      case MA_ACTION_INTERNALDOWNLOWER:
      {
        #ifdef DEBUG_ON
          if (debuglevel == 2)
            printout("[debug] do_ma:   internal downward jump to lower ionstage\n");
          pkt_ptr->interactions += 1;
          jumps++;
          jump = 1;
        #endif

        ma_stat_internaldownlower++;

        /// Randomly select the occuring transition
        zrand = gsl_rng_uniform(rng);
        //zrand = 1. - 1e-14;
        rate = 0.;
        //nlevels = get_nlevels(element,ion-1);

        const int nlevels = get_ionisinglevels(element, ion - 1);
        //nlevels = get_ionisinglevels(element,ion-1);
        int lower;
        for (lower = 0; lower < nlevels; lower++)
        {
          const double epsilon_target = epsilon(element, ion - 1, lower);
          const double epsilon_trans = epsilon_current - epsilon_target;
          const double R = rad_recombination_ratecoeff(T_e, nne, element, ion, level, lower, modelgridindex);
          const double C = col_recombination_ratecoeff(modelgridindex, element, ion, level, lower, epsilon_trans);
          rate += (R + C) * epsilon_target;
          if (zrand * processrates[MA_ACTION_INTERNALDOWNLOWER] < rate)
            break;
        }
        /// and set the macroatom's new state

        #if (TRACK_ION_STATS)
        increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_INTERNAL, pkt_ptr->e_cmf);
        #endif

        ion -= 1;
        level = lower;

        #if (TRACK_ION_STATS)
        increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYIN_INTERNAL, pkt_ptr->e_cmf);
        #endif

        #ifdef DEBUG_ON
          if (lower >= nlevels)
          {
            printout("internal_down_lower  %g\n", processrates[MA_ACTION_INTERNALDOWNLOWER]);
            printout("abort at rate %g, zrand %g\n",rate,zrand);
            abort();
          }
          if (get_ionstage(element,ion-1) == 0 && lower == 0)
  	  //        if (ion-1 == 0 && lower == 0)
          {
            printout("internal downward transition to ground level occured ... abort\n");
            printout("element %d, ion %d, level %d, lower %d\n",element,ion,level,lower);
            printout("Z %d, ionstage %d, energy %g\n",get_element(element),get_ionstage(element,ion-1),elements[element].ions[ion-1].levels[lower].epsilon);
            printout("[debug] do_ma:   internal downward jump to lower ionstage\n");
            abort();
          }
        #endif
        break;
      }

      case MA_ACTION_INTERNALUPSAME:
      {
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] do_ma:   internal upward jump within current ionstage\n");
          pkt_ptr->interactions += 1;
          jumps++;
          jump = 2;
        #endif

        ///randomly select the occuring transition
        zrand = gsl_rng_uniform(rng);
        int upper = -99;
        rate = 0.;
        for (int i = 0; i < nuptrans; i++)
        {
          rate += get_individ_internal_up_same(element, ion, level, i);
          if (zrand * processrates[MA_ACTION_INTERNALUPSAME] < rate)
          {
            const int lineindex = elements[element].ions[ion].levels[level].uptrans_lineindicies[i];
            upper = linelist[lineindex].upperlevelindex;
            break;
          }
        }
        ///and set the macroatom's new state
        level = upper;
        break;
      }

      case MA_ACTION_INTERNALUPHIGHER:
      {
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] do_ma:   internal upward jump to next ionstage\n");
          pkt_ptr->interactions += 1;
          jumps++;
          jump = 3;
        #endif

        ma_stat_internaluphigher++;

        #if (TRACK_ION_STATS)
        increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_INTERNAL, pkt_ptr->e_cmf);
        #endif

        do_macroatom_ionisation(modelgridindex, element, &ion, &level, epsilon_current, processrates[MA_ACTION_INTERNALUPHIGHER]);

        #if (TRACK_ION_STATS)
        increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYIN_INTERNAL, pkt_ptr->e_cmf);
        #endif

        break;
      }

      case MA_ACTION_INTERNALUPHIGHERNT:
      {
        pkt_ptr->interactions += 1;
        // ion += 1;
        #if (TRACK_ION_STATS)
        increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_INTERNAL, pkt_ptr->e_cmf);
        #endif

        ion = nt_random_upperion(modelgridindex, element, ion, false);
        level = 0;
        ma_stat_internaluphighernt++;

        #if (TRACK_ION_STATS)
        increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYIN_INTERNAL, pkt_ptr->e_cmf);
        #endif
        // printout("Macroatom non-thermal ionisation to Z=%d ionstage %d level %d\n", get_element(element), ion, level);
        break;
      }

      case MA_ACTION_COUNT:
      {
        printout("ERROR: Somehow selected MA_ACTION_COUNT\n");
        abort();
      }
    }
  }///endwhile

  #ifdef DEBUG_ON
    //debuglevel = 4;
    //debuglevel = 2;
  #endif

  if (pkt_ptr->trueemissiontype < 0)
  {
    pkt_ptr->trueemissiontype = pkt_ptr->emissiontype;
    pkt_ptr->trueemissionvelocity = vec_len(pkt_ptr->em_pos) / pkt_ptr->em_time;
    pkt_ptr->trueem_time = pkt_ptr->em_time;
  }

  /// procedure ends only after a change to r or k packets has taken place
}


/// Calculation of radiative rates ///////////////////////////////////////////////////////


void macroatom_open_file(const int my_rank)
{
  if (!LOG_MACROATOM)
    return;
  char filename[100];
  sprintf(filename, "macroatom_%.4d.out",my_rank);
  macroatom_file = fopen_required(filename, "w");
  fprintf(macroatom_file, "%8s %14s %2s %12s %12s %9s %9s %9s %11s %11s %11s %11s %9s\n",
          "timestep", "modelgridindex", "Z", "ionstage_in", "ionstage_out", "level_in", "level_out", "activline",
          "nu_cmf_in", "nu_cmf_out", "nu_rf_in", "nu_rf_out", "jumps");
}


void macroatom_close_file(void)
{
  if (LOG_MACROATOM)
  {
    fclose(macroatom_file);
  }
}


double rad_deexcitation_ratecoeff(
  const int modelgridindex, const int element, const int ion, const int upper, const int lower,
  const double epsilon_trans, const int lineindex, const double t_current)
///radiative deexcitation rate: paperII 3.5.2
// multiply by upper level population to get a rate per second
{
  #ifdef DEBUG_ON
  if (upper <= lower)
  {
    printout("[fatal] rad_deexcitation: tried to calculate upward transition ... abort\n");
    abort();
  }
  #endif

  const double n_u = calculate_exclevelpop(modelgridindex, element, ion, upper);
  const double n_l = calculate_exclevelpop(modelgridindex, element, ion, lower);

  double R = 0.0;

  // if ((n_u > 1.1 * MINPOP) && (n_l > 1.1 * MINPOP))
  {
    const double nu_trans = epsilon_trans / H;

    const double A_ul = einstein_spontaneous_emission(lineindex);
    const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
    const double B_lu = stat_weight(element, ion, upper) / stat_weight(element, ion, lower) * B_ul;

    const double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current;

    if (tau_sobolev > 1e-100)
    {
      const double beta = 1.0 / tau_sobolev * (-expm1(-tau_sobolev));
      // const double beta = 1.0;
      R = A_ul * beta;
    }
    else
    {
      //printout("[warning] rad_deexcitation: tau_sobolev %g <= 0, set beta=1\n",tau_sobolev);
      //printout("[warning] rad_deexcitation: element %d, ion %d, upper %d, lower %d\n",element,ion,upper,lower);
      //printout("[warning] rad_deexcitation: n_l %g, n_u %g, B_lu %g, B_ul %g\n",n_l,n_u,B_lu,B_ul);
      //printout("[warning] rad_deexcitation: T_e %g, T_R %g, W %g in model cell %d\n",get_Te(modelgridindex),get_TR(modelgridindex),get_W(modelgridindex),modelgridindex);
      R = 0.0;
      //printout("[fatal] rad_excitation: tau_sobolev <= 0 ... %g abort",tau_sobolev);
      //abort();
    }

    #ifdef DEBUG_ON
      if (debuglevel == 2)
      {
        // printout("[debug] rad_rates_down: Z=%d, ionstage %d, upper %d, lower %d\n", get_element(element), get_ionstage(element, ion), upper, lower);
        // printout("[debug] rad_deexc: A_ul %g, tau_sobolev %g, n_u %g\n", A_ul, tau_sobolev, n_u);
      }
      else if (debuglevel == 777)
        printout("[debug] rad_deexc: A_ul %g, tau_sobolev %g, n_u %g\n", A_ul, tau_sobolev, n_u);
      if (!isfinite(R))
      {
        printout("fatal a1: abort\n");
        abort();
      }
    #endif
  }

  return R;
}


double rad_excitation_ratecoeff(
  const int modelgridindex, const int element, const int ion, const int lower,
  const int upper, const double epsilon_trans, const int lineindex, const double t_current)
///radiative excitation rate: paperII 3.5.2
// multiply by lower level population to get a rate per second
{
  #ifdef DEBUG_ON
  if (upper <= lower)
  {
    printout("[fatal] rad_excitation: tried to calculate downward transition ... abort\n");
    abort();
  }
  #endif

  const double n_u = calculate_exclevelpop(modelgridindex, element, ion, upper);
  const double n_l = calculate_exclevelpop(modelgridindex, element, ion, lower);
  double R = 0.0;
  // if ((n_u >= 1.1 * MINPOP) && (n_l >= 1.1 * MINPOP))
  {
    const double nu_trans = epsilon_trans / H;
    const double A_ul = einstein_spontaneous_emission(lineindex);
    const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
    const double B_lu = stat_weight(element, ion, upper) / stat_weight(element, ion, lower) * B_ul;

    const double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current;

    if (tau_sobolev > 1e-100)
    {
      double beta = 1.0 / tau_sobolev * (-expm1(-tau_sobolev));
      //printout("[check] rad_excitation: %g, %g, %g\n",1.0/tau_sobolev,exp(-tau_sobolev),1.0/tau_sobolev * (1. - exp(-tau_sobolev)));
      //n_u2 = calculate_levelpop_fromreflevel(pkt_ptr->where,element,ion,upper,lower,mastate[tid].nnlevel);
      //R = (B_lu*mastate[tid].nnlevel - B_ul * n_u2) * beta * radfield(nu_trans,pkt_ptr->where);

      const double R_over_J_nu = n_l > 0. ? (B_lu - B_ul * n_u / n_l) * beta : B_lu * beta;

      // const double linelambda = 1e8 * CLIGHT / nu_trans;
      // if (use_cellhist && false) //  || (linelambda > 7000)
      {
        R = R_over_J_nu * radfield(nu_trans, modelgridindex);
      }
      // else
      {
        // R = R_over_J_nu * radfield_dbb_mgi(nu_trans, modelgridindex);
        // R = 0.;
      }
      if (!initial_iteration)
      {
        // check for a detailed line flux estimator to replace the binned/blackbody radiation field estimate
        const int jblueindex = radfield_get_Jblueindex(lineindex);
        if (jblueindex >= 0)
        {
          const double Jb_lu = radfield_get_Jb_lu(modelgridindex, jblueindex);
          const double R_Jb = R_over_J_nu * Jb_lu;
          // const int contribcount = radfield_get_Jb_lu_contribcount(modelgridindex, jblueindex);
          // const double R_radfield = R_over_J_nu * radfield(nu_trans, modelgridindex);
          // const double linelambda = 1e8 * CLIGHT / nu_trans;
          // printout("Using detailed rad excitation lambda %5.1f contribcont %d R(Jblue) %g R(radfield) %g R_Jb/R %g\n",
          //          linelambda, contribcount, R_Jb, R_radfield, R_Jb / R_radfield);
          // printout("  (for transition Z=%02d ionstage %d lower %d upper %d)\n",
          //          get_element(element), get_ionstage(element, ion), lower, upper);
          R = R_Jb;
        }
      }
    }
    else
    {
      //printout("[warning] rad_excitation: tau_sobolev %g <= 0, set beta=1\n",tau_sobolev);
      //printout("[warning] rad_excitation: element %d, ion %d, upper %d, lower %d\n",element,ion,upper,lower);
      //printout("[warning] rad_excitation: n_l %g, n_u %g, B_lu %g, B_ul %g\n",n_l,n_u,B_lu,B_ul);
      //printout("[warning] rad_excitation: T_e %g, T_R %g, W %g\n",get_Te(modelgridindex),get_TR(modelgridindex),get_W(modelgridindex));
      R = 0.;

      //printout("[fatal] rad_excitation: tau_sobolev <= 0 ... %g abort",tau_sobolev);
      //abort();
    }

    #ifdef DEBUG_ON
    if (R < 0)
    {
      const double g_u = statw_upper(lineindex);
      const double g_u2 = stat_weight(element, ion, upper);
      const double g_l = statw_lower(lineindex);
      const double g_l2 = stat_weight(element, ion, lower);
      printout("Negative excitation rate from level %d to %d\n", lower, upper);
      printout("n_l %g, n_u %g, g_l %g (?=%g), g_u %g (?=%g)\n", n_l, n_u, g_l, g_l2, g_u, g_u2);
      printout("n_u/n_l %g, g_u/g_l %g\n", n_u / n_l, g_u / g_l);
      printout("radfield(nutrans=%g) = %g\n", nu_trans, radfield(nu_trans, modelgridindex));
      abort();
    }
    if (debuglevel == 2)
    {
      // printout("[debug] rad_rates_up: Z=%d, ionstage %d, upper, lower, A_ul, n_u: %d, %d, %d, %d, %g, %g\n",
      //          get_element(element), get_ionstage(element, ion), upper, lower, A_ul, n_l);
      // printout("[debug] rad_exc: A_ul %g, tau_sobolev %g, n_u %g, n_l %g, radfield %g\n",
      //          A_ul, tau_sobolev, n_u, n_l, radfield(nu_trans,modelgridindex));
    }
    else if (debuglevel == 777)
      printout("[debug] rad_exc: A_ul %g, tau_sobolev %g, n_u %g, n_l %g, radfield %g\n", A_ul, tau_sobolev, n_u, n_l, radfield(nu_trans, modelgridindex));
    if (!isfinite(R))
    {
      printout("[fatal] rad_excitation: abort\n");
      printout("[fatal] rad_excitation: R %g, B_lu %g, B_ul %g, n_u %g, n_l %g, radfield %g,tau_sobolev %g, t_current %g\n",
               R,B_lu,B_ul,n_u,n_l,radfield(nu_trans,modelgridindex),tau_sobolev,t_current);
      printout("[fatal] rad_excitation: %g, %g, %g\n", 1.0 / tau_sobolev, exp(-tau_sobolev), 1.0 / tau_sobolev * (1. - exp(-tau_sobolev)));
      abort();
    }
    #endif
  }

  return R;
}


double rad_recombination_ratecoeff(
  const float T_e, const float nne, const int element, const int upperion, const int upper, const int lower, const int modelgridindex)
///radiative recombination rate: paperII 3.5.2
// multiply by upper level population to get a rate per second
{
  // it's probably faster to only check this condition outside this function
  // in a case where this wasn't checked, the function will return zero anyway
  // if (upper > get_maxrecombininglevel(element, upperion))
  //   return 0.;

  double R = 0.0;
  const int nphixstargets = get_nphixstargets(element, upperion - 1,lower);
  for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++)
  {
    if (get_phixsupperlevel(element,upperion-1,lower,phixstargetindex) == upper)
    {
      R = nne * get_spontrecombcoeff(element, upperion - 1, lower, phixstargetindex, T_e);// + stimrecombestimator_save[pkt_ptr->where*nelements*maxion+element*maxion+(ion-1)]);
      //printout("calculate rad_recombination: element %d, ion %d, upper %d, -> lower %d, nne %g, spontrecombcoeff %g\n",element,ion,upper,lower,nne,get_spontrecombcoeff(element, ion-1, lower, T_e));

      if (modelgridindex >= 0 && SEPARATE_STIMRECOMB)
      {
        R += nne * get_stimrecombcoeff(element, upperion - 1, lower, phixstargetindex, modelgridindex);
      }
      break;
    }
  }

  #ifdef DEBUG_ON
  assert(isfinite(R));
  #endif

  return R;
}


double stim_recombination_ratecoeff(
  const float nne, const int element, const int upperion, const int upper, const int lower, const int modelgridindex)
///radiative recombination rate: paperII 3.5.2
// multiply by upper level population to get a rate per second
{
  // it's probably faster to only check this condition outside this function
  // in a case where this wasn't checked, the function will return zero anyway
  // if (upper > get_maxrecombininglevel(element, upperion))
  //   return 0.;

  double R = 0.0;
  const int nphixstargets = get_nphixstargets(element, upperion - 1,lower);
  for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++)
  {
    if (get_phixsupperlevel(element, upperion - 1, lower, phixstargetindex) == upper)
    {
      R = nne * get_stimrecombcoeff(element, upperion - 1, lower, phixstargetindex, modelgridindex);
    }
    break;
  }

  assert(isfinite(R));

  return R;
}


/// Calculation of collisional rates /////////////////////////////////////////////////////
///***************************************************************************/


double col_deexcitation_ratecoeff(const float T_e, const float nne, const double epsilon_trans, const int lineindex)
// multiply by upper level population to get a rate per second
{
  double C;
  const double coll_str_thisline = get_coll_str(lineindex);
  const double upperstatweight = statw_upper(lineindex);
  if (coll_str_thisline < 0)
  {
    const double statweight_target = statw_lower(lineindex);
    const bool forbidden = linelist[lineindex].forbidden;
    if (!forbidden) // alternative: (coll_strength > -1.5) i.e. to catch -1
    {
      ///permitted E1 electric dipole transitions
      ///collisional deexcitation: formula valid only for atoms!!!!!!!!!!!
      ///Rutten script eq. 3.33. p.50
      //f = osc_strength(element,ion,upper,lower);
      //C = n_u * 2.16 * pow(fac1,-1.68) * pow(T_e,-1.5) * stat_weight(element,ion,lower)/stat_weight(element,ion,upper)  * nne * f;

      const double eoverkt = epsilon_trans / (KB * T_e);
      ///Van-Regemorter formula, Mihalas (1978), eq.5-75, p.133
      const double g_bar = 0.2; ///this should be read in from transitions data: it is 0.2 for transitions nl -> n'l' and 0.7 for transitions nl -> nl'
      //test = 0.276 * exp(fac1) * gsl_sf_expint_E1(fac1);
      /// crude approximation to the already crude Van-Regemorter formula

      //double test = 0.276 * exp(fac1) * (-0.5772156649 - log(fac1));
      //double Gamma = (g_bar > test) ? g_bar : test;

      //optimisation
      const double gauntfac = (eoverkt > 0.33421) ? g_bar : 0.276 * exp(eoverkt) * (-0.5772156649 - log(eoverkt));

      const double g_ratio = statweight_target / upperstatweight;

      C = C_0 * 14.51039491 * nne * sqrt(T_e) * osc_strength(lineindex) * pow(H_ionpot / epsilon_trans, 2) * eoverkt * g_ratio * gauntfac;
    }
    else // alterative: (coll_strength > -3.5) to catch -2 or -3
    {
      //forbidden transitions: magnetic dipole, electric quadropole...
      //could be Axelrod? or Maurer
      C = nne * 8.629e-6 * 0.01 * statweight_target / sqrt(T_e);
    }
  }
  else //positive values are treated as effective collision strengths
  {
    //from Osterbrock and Ferland, p51
    //statweight_target is LOWER LEVEL stat weight
    C = nne * 8.629e-6 * coll_str_thisline / upperstatweight / sqrt(T_e);
    // test test
    //C = n_u * nne * 8.629e-6 * pow(T_e,-0.5) * 0.01 * statweight_target;
  }

  #ifdef DEBUG_ON
    /*if (debuglevel == 2)
    {
      //printout("[debug] col_deexc: element %d, ion %d, upper %d, lower %d\n",element,ion,upper,lower);
      printout("[debug] col_deexc: n_u %g, nne %g, T_e %g, f_ul %g, epsilon_trans %g, Gamma %g, g_ratio %g\n",n_u,nne,T_e,osc_strength(lineindex),epsilon_trans,Gamma,g_ratio);
    }*/
    //printout("col_deexc(%d,%d,%d,%d) %g\n",element,ion,upper,lower,C);
    assert(isfinite(C));
  #endif

  return C;
}


double col_excitation_ratecoeff(const float T_e, const float nne, const int lineindex, const double epsilon_trans)
// multiply by lower level population to get a rate per second
{
  double C;
  const double coll_strength = get_coll_str(lineindex);
  const double eoverkt = epsilon_trans / (KB * T_e);

  #ifdef DEBUG_ON
  // if (upper <= lower)
  // {
  //   printout("[fatal] col_excitation: tried to calculate downward transition ... abort");
  //   abort();
  // }
  #endif

  if (coll_strength < 0)
  {
    const bool forbidden = linelist[lineindex].forbidden;
    if (!forbidden) // alternative: (coll_strength > -1.5) i.e. to catch -1
    {
      /// permitted E1 electric dipole transitions
      /// collisional excitation: formula valid only for atoms!!!!!!!!!!!
      /// Rutten script eq. 3.32. p.50
      //C = n_l * 2.16 * pow(eoverkt,-1.68) * pow(T_e,-1.5) * exp(-eoverkt) * nne * osc_strength(element,ion,upper,lower);

      // Van-Regemorter formula, Mihalas (1978), eq.5-75, p.133
      const double g_bar = 0.2; // this should be read in from transitions data: it is 0.2 for transitions nl -> n'l' and 0.7 for transitions nl -> nl'
      // test = 0.276 * exp(eoverkt) * gsl_sf_expint_E1(eoverkt);
      /// crude approximation to the already crude Van-Regemorter formula
      const double exp_eoverkt = exp(eoverkt);

      const double test = 0.276 * exp_eoverkt * (-0.5772156649 - log(eoverkt));
      const double Gamma = g_bar > test ? g_bar : test;
      C = C_0 * nne * sqrt(T_e) * 14.51039491 * osc_strength(lineindex) * pow(H_ionpot/epsilon_trans, 2) * eoverkt / exp_eoverkt * Gamma;
    }
    else // alterative: (coll_strength > -3.5) to catch -2 or -3
    {
      // forbidden transitions: magnetic dipole, electric quadropole...
      // Axelrod's approximation (thesis 1980)
      C = nne * 8.629e-6 * 0.01 * exp(-eoverkt) * statw_upper(lineindex) / sqrt(T_e);
    }
  }
  else
  {
    //from Osterbrock and Ferland, p51
    C = nne * 8.629e-6 * coll_strength * exp(-eoverkt) / statw_lower(lineindex) / sqrt(T_e);
    //test test
    //C = n_l * nne * 8.629e-6 * pow(T_e,-0.5) * 0.01 * exp(-fac1) * statw_upper(lineindex);
  }

  #ifdef DEBUG_ON
    // if (debuglevel == 2)
    // {
      //printout("[debug] col_exc: element %d, ion %d, lower %d, upper %d\n",element,ion,lower,upper);
      //printout("[debug] col_exc: n_l %g, nne %g, T_e %g, f_ul %g, epsilon_trans %g, Gamma %g\n",n_l, nne,T_e,osc_strength(lineindex),epsilon_trans,Gamma);
    // }

    //if (!isfinite(C))
    //{
    //  printout("fatal a5: abort\n");
      //printout("[debug] col_exc: element %d, ion %d, lower %d, upper %d\n",element,ion,lower,upper);
      //printout("[debug] col_exc: n_l %g, nne %g, T_e %g, f_ul %g, epsilon_trans %g, Gamma %g\n",n_l, nne,T_e,osc_strength(lineindex),epsilon_trans,Gamma);
      //printout("[debug] col_exc: g_bar %g, fac1 %g, test %g, %g, %g, %g\n",g_bar,fac1,test,0.276 * exp(fac1),-0.5772156649 - log(fac1),0.276 * exp(fac1) * (-0.5772156649 - log(fac1)));
    //  printout("[debug] col_exc: get_coll_str(lineindex) %g statw_upper(lineindex) %g mastate[tid].statweight %g\n", get_coll_str(lineindex),statw_upper(lineindex),mastate[tid].statweight);
    //  abort();
    //}
  #endif

  return C;
}


double col_recombination_ratecoeff(
  const int modelgridindex, const int element, const int upperion, const int upper, const int lower, const double epsilon_trans)
// multiply by upper level population to get a rate per second
{
  // it's probably faster to only check this condition outside this function
  // in a case where this wasn't checked, the function will return zero anyway
  // if (upper > get_maxrecombininglevel(element, upperion))
  //   return 0.;

  const int nphixstargets = get_nphixstargets(element, upperion - 1, lower);
  for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++)
  {
    if (get_phixsupperlevel(element, upperion - 1, lower, phixstargetindex) == upper)
    {
      const float nne = get_nne(modelgridindex);
      const float T_e = get_Te(modelgridindex);
      const double fac1 = epsilon_trans / KB / T_e;
      const int ionstage = get_ionstage(element, upperion);

      ///Seaton approximation: Mihalas (1978), eq.5-79, p.134
      ///select gaunt factor according to ionic charge
      double g;
      if (ionstage - 1 == 1)
        g = 0.1;
      else if (ionstage - 1 == 2)
        g = 0.2;
      else
        g = 0.3;

      const double sigma_bf = (
        elements[element].ions[upperion - 1].levels[lower].photoion_xs[0] *
        get_phixsprobability(element, upperion - 1, lower, phixstargetindex));

      double C = nne * nne * calculate_sahafact(element, upperion - 1, lower, upper, T_e, epsilon_trans) *
                       1.55e13 * pow(T_e, -0.5) * g * sigma_bf * exp(-fac1) / fac1;

      return C;
    }
  }

  return 0.;
}


double col_ionization_ratecoeff(
  const float T_e, const float nne, const int element, const int ion, const int lower,
  const int phixstargetindex, const double epsilon_trans)
/// collisional ionization rate: paperII 3.5.1
// multiply by lower level population to get a rate per second
{
  #ifdef DEBUG_ON
  if (phixstargetindex > get_nphixstargets(element,ion,lower))
  {
    printout("[fatal] col_ionization called with phixstargetindex %d > nphixstargets %d",
             phixstargetindex, get_nphixstargets(element, ion, lower));
    abort();
  }
  #endif

  ///Seaton approximation: Mihalas (1978), eq.5-79, p.134
  ///select gaunt factor according to ionic charge
  double g;
  const int ionstage = get_ionstage(element,ion);
  if (ionstage == 1)
    g = 0.1;
  else if (ionstage == 2)
    g = 0.2;
  else
    g = 0.3;

  const double fac1 = epsilon_trans / KB / T_e;

  const double sigma_bf = elements[element].ions[ion].levels[lower].photoion_xs[0] * get_phixsprobability(element,ion,lower,phixstargetindex);
  const double C = nne * 1.55e13 * pow(T_e,-0.5) * g * sigma_bf * exp(-fac1) / fac1; ///photoionization at the edge

  #ifdef DEBUG_ON
    if (debuglevel == 777)
      printout("[debug] col_ion: nne %g, T_e %g, g %g, epsilon_trans %g, sigma_bf %g\n", nne,T_e,g,epsilon_trans,sigma_bf);
    assert(isfinite(C));
  #endif

  return C;
}


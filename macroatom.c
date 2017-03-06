#include "sn3d.h"
#include <gsl/gsl_integration.h>
#include "atomic.h"
#include "grid_init.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "rpkt.h"

// constant for van-Regemorter approximation.
#define C_0 (5.465e-11)


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


static void set_cellhistory_transitionrates(
  const int modelgridindex, const int element, const int ion, const int level, const double t_mid)
{
  const float T_e = get_Te(modelgridindex);
  const float nne = get_nne(modelgridindex);
  const double epsilon_current = epsilon(element, ion, level);

  /// Downward transitions within the current ionisation stage:
  /// radiative/collisional deexcitation and internal downward jumps
  double rad_deexc = 0.;
  double col_deexc = 0.;
  double internal_down_same = 0.;
  const int ndowntrans = elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
  for (int i = 1; i <= ndowntrans; i++)
  {
    const int lower = elements[element].ions[ion].levels[level].downtrans[i].targetlevel;
    const double epsilon_target = epsilon(element, ion, lower);
    const double epsilon_trans = elements[element].ions[ion].levels[level].downtrans[i].epsilon_trans;
    const int lineindex = elements[element].ions[ion].levels[level].downtrans[i].lineindex;

    const double R = rad_deexcitation_ratecoeff(modelgridindex, element, ion, level, lower, epsilon_trans, lineindex, t_mid);
    const double C = col_deexcitation_ratecoeff(T_e, nne, epsilon_trans, lineindex);

    const double individ_rad_deexc = R * epsilon_trans;
    const double individ_col_deexc = C * epsilon_trans;
    const double individ_internal_down_same = (R + C) * epsilon_target;

    cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_rad_deexc[i] = individ_rad_deexc;
    cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_down_same[i] = individ_internal_down_same;

    rad_deexc += individ_rad_deexc;
    col_deexc += individ_col_deexc;
    internal_down_same += individ_internal_down_same;

    #ifdef DEBUG_ON
      if (debuglevel == 2) printout("checking downtrans %d to level %d: R %g, C %g, epsilon_trans %g\n",i,lower,R,C,epsilon_trans);
    #endif
  }
  cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_deexc = rad_deexc;
  cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_deexc = col_deexc;
  cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_same = internal_down_same;

  /// Downward transitions to lower ionisation stages:
  /// radiative/collisional recombination and internal downward jumps
  double rad_recomb = 0.;
  double col_recomb = 0.;
  double internal_down_lower = 0.;
  if (ion > 0) ///checks only if there is a lower ion, doesn't make sure that Z(ion)=Z(ion-1)+1
  {
    //nlevels = get_nlevels(element,ion-1);
    const int nlevels = get_bfcontinua(element,ion-1);
    //nlevels = get_ionisinglevels(element,ion-1);
    for (int lower = 0; lower < nlevels; lower++)
    {
      const double epsilon_target = epsilon(element, ion - 1, lower);
      const double epsilon_trans = epsilon_current - epsilon_target;

      const double R = rad_recombination_ratecoeff(modelgridindex, element, ion, level, lower);
      //printout("rad recombination of element %d, ion %d, level %d, to lower level %d has rate %g\n",element,ion,level,lower,R);
      const double C = col_recombination_ratecoeff(T_e, nne, element, ion, level, lower, epsilon_trans);

      rad_recomb += R * epsilon_trans;
      col_recomb += C * epsilon_trans;

      internal_down_lower += (R + C) * epsilon_target;
    }
  }
  cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_recomb = rad_recomb;
  cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_recomb = col_recomb;
  cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_lower = internal_down_lower;

  /// Calculate sum for upward internal transitions
  /// transitions within the current ionisation stage
  double internal_up_same = 0.;
  const int nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
  for (int i = 1; i <= nuptrans; i++)
  {
    const int upper = elements[element].ions[ion].levels[level].uptrans[i].targetlevel;
    const double epsilon_trans = elements[element].ions[ion].levels[level].uptrans[i].epsilon_trans;
    const int lineindex = elements[element].ions[ion].levels[level].uptrans[i].lineindex;

    const double R = rad_excitation_ratecoeff(modelgridindex, element, ion, level, upper, epsilon_trans, lineindex, t_mid);
    const double C = col_excitation_ratecoeff(T_e, nne, lineindex, epsilon_trans);

    const double individ_internal_up_same = (R + C) * epsilon_current;

    cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_up_same[i] = individ_internal_up_same;

    internal_up_same += individ_internal_up_same;

    #ifdef DEBUG_ON
      if (debuglevel == 2) printout("checking uptrans %d to level %d: R %g, C %g, epsilon_trans %g\n",i,upper,R,C,epsilon_trans);
      if (!isfinite(internal_up_same)) {printout("fatal: internal_up_same has nan contribution\n");}
    #endif
  }
  cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_same = internal_up_same;

  /// Transitions to higher ionisation stages
  double internal_up_higher = 0.;

  const int ionisinglevels = get_bfcontinua(element,ion);
  if (ion < get_nions(element) - 1 && level < ionisinglevels)  //&& get_ionstage(element,ion) < get_element(element)+1)
  {
    // if (NT_ON)
    //   internal_up_higher += nt_ionization_ratecoeff(modelgridindex, element, ion);
    for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++)
    {
      const int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
      const double epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;

      const double R = get_corrphotoioncoeff(element,ion,level,phixstargetindex,modelgridindex);
      const double C = col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);

      internal_up_higher += (R + C) * epsilon_current;
    }
  }
  cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_higher = internal_up_higher;
}


static void do_macroatom_raddeexcitation(
  PKT *restrict pkt_ptr, const int element, const int ion, const int level, const double rad_deexc,
  const double total_transitions, const double t_current)
{

  ///radiative deexcitation of MA: emitt rpkt
  ///randomly select which line transitions occurs
  const double zrand = gsl_rng_uniform(rng);
  //zrand = 1. - 1e-14; /// ONLY FOR DEBUG!!!
  double rate = 0.;
  int linelistindex = -99;
  int i;
  const int ndowntrans = elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
  double epsilon_trans = -100;
  for (i = 1; i <= ndowntrans; i++)
  {
    rate += get_individ_rad_deexc(element, ion, level, i);
    if (zrand * rad_deexc < rate)
    {
      const int lower = elements[element].ions[ion].levels[level].downtrans[i].targetlevel;
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] do_ma:   jump to level %d\n", lower);
      #endif
      epsilon_trans = elements[element].ions[ion].levels[level].downtrans[i].epsilon_trans;
      linelistindex = elements[element].ions[ion].levels[level].downtrans[i].lineindex;
      //linelistindex = elements[element].ions[ion].levels[level].transitions[level-lower-1].linelistindex;
      #ifdef RECORD_LINESTAT
        if (tid == 0) ecounter[linelistindex]++;    /// This way we will only record line statistics from OMP-thread 0
                                                    /// With an atomic pragma or a thread-private structure with subsequent
                                                    /// reduction this could be extended to all threads. However, I'm not
                                                    /// sure if this is worth the additional computational expenses.
      #endif
      break;
    }
  }
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
    if (i > ndowntrans)
    {
      printout("[fatal] problem in selecting radiative downward transition of MA zrand %g, rate %g, rad_deexc %g, i %d, ndowntrans %d\n", zrand, rate, rad_deexc, i, ndowntrans);
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
  emitt_rpkt(pkt_ptr, t_current);

  if (linelistindex == mastate[tid].activatingline)
    resonancescatterings++;
  else
    calculate_kappa_rpkt_cont(pkt_ptr, t_current);

  /// NB: the r-pkt can only interact with lines redder than the current one
  pkt_ptr->next_trans = linelistindex + 1;
  pkt_ptr->emissiontype = linelistindex;
  pkt_ptr->em_pos[0] = pkt_ptr->pos[0];
  pkt_ptr->em_pos[1] = pkt_ptr->pos[1];
  pkt_ptr->em_pos[2] = pkt_ptr->pos[2];
  pkt_ptr->em_time = t_current;
  pkt_ptr->nscatterings = 0;
  //printout("next possible line encounter %d\n",pkt_ptr->next_trans);
  #ifndef FORCE_LTE
    //matotem[pkt_ptr->where] += pkt_ptr->e_cmf;
  #endif
}


static void do_macroatom_radrecomb(
  PKT *restrict pkt_ptr, const int modelgridindex, const int element, const int ion, const int level,
  const double rad_recomb, const double t_current)
{
  const float T_e = get_Te(modelgridindex);

  const double epsilon_current = epsilon(element, ion, level);
  /// Randomly select a continuum
  double zrand = gsl_rng_uniform(rng);
  //zrand = 1. - 1e-14;
  double rate = 0;
  //nlevels = get_nlevels(element,ion-1);
  const int nlevels = get_bfcontinua(element, ion - 1);
  //nlevels = get_ionisinglevels(element,ion-1);
  int lower;
  double epsilon_trans;
  for (lower = 0; lower < nlevels; lower++)
  {
    epsilon_trans = epsilon_current - epsilon(element,ion - 1,lower);
    const double R = rad_recombination_ratecoeff(modelgridindex, element, ion, level, lower);
    rate += R * epsilon_trans;
    #ifdef DEBUG_ON
      if (debuglevel == 2)
      {
        printout("[debug] do_ma:   R %g, deltae %g\n",R,(epsilon(element, ion, level) - epsilon(element, ion - 1, lower)));
        printout("[debug] do_ma:   rate to level %d of ion %d = %g\n", lower, ion - 1, rate);
        printout("[debug] do_ma:   zrand*rad_recomb = %g\n", zrand * rad_recomb);
      }
    #endif
    if (zrand * rad_recomb < rate) break;
  }
  /// and set its threshold frequency
  const double nu_threshold = epsilon_trans / H;
  #ifdef DEBUG_ON
    if (debuglevel == 2) printout("[debug] do_ma:   going to level %d of ion %d of element %d\n", lower, ion - 1, element);
  #endif

  /// Then randomly sample the packets frequency according to the continuums
  /// energy distribution and set some flags
  //zrand = gsl_rng_uniform(rng);
  //zrand = 1. - zrand;  /// Make sure that 0 < zrand <= 1
  //pkt_ptr->nu_cmf = nu_threshold * (1 - KB*T_e/H/nu_threshold*log(zrand));
  //pkt_ptr->nu_cmf = nu_threshold;


  zrand = gsl_rng_uniform(rng);
  zrand = 1. - zrand;  /// Make sure that 0 < zrand <= 1
  mastate[tid].element = element;
  mastate[tid].ion = ion - 1;
  mastate[tid].level = lower;
  const double intaccuracy = 1e-2;        /// Fractional accuracy of the integrator
  gsl_integration_workspace *wsp = gsl_integration_workspace_alloc(1024);
  gslintegration_paras intparas;
  intparas.T = T_e;
  intparas.nu_edge = nu_threshold;   /// Global variable which passes the threshold to the integrator
  gsl_function F_alpha_sp;
  //F_alpha_sp.function = &alpha_sp_integrand_gsl;
  F_alpha_sp.function = &alpha_sp_E_integrand_gsl;
  F_alpha_sp.params = &intparas;
  const double deltanu = nu_threshold * NPHIXSNUINCREMENT;
  const double nu_max_phixs = nu_threshold * last_phixs_nuovernuedge; //nu of the uppermost point in the phixs table
  double error;
  double total_alpha_sp;
  gsl_integration_qag(&F_alpha_sp, nu_threshold, nu_max_phixs, 0, intaccuracy, 1024, GSL_INTEG_GAUSS61, wsp, &total_alpha_sp, &error);
  double alpha_sp_old = total_alpha_sp;
  double nu_lower = nu_threshold;
  for (int i = 1; i < NPHIXSPOINTS; i++)
  // LJS: this loop could probably be made a bit faster
  // use the overlap with the previous integral and add on a piece each time instead of recalculating the
  // integral over the entire region
  {
    // the reason the lower limit of integration is incremented is that most of the probability distribution is at the low
    // frequency end, so this minimizes the number of iterations needed
    double alpha_sp;
    nu_lower += deltanu;
    /// Spontaneous recombination and bf-cooling coefficient don't depend on the cutted radiation field
    gsl_integration_qag(&F_alpha_sp, nu_lower, nu_max_phixs, 0, intaccuracy, 1024, GSL_INTEG_GAUSS61, wsp, &alpha_sp, &error);
    //alpha_sp *= FOURPI * sf;
    //if (zrand > alpha_sp/get_spontrecombcoeff(element,ion-1,lower,pkt_ptr->where)) break;
    if (zrand >= alpha_sp / total_alpha_sp)
    {
      const double nuoffset = (total_alpha_sp * zrand - alpha_sp_old) / (alpha_sp - alpha_sp_old) * deltanu;
      nu_lower = nu_threshold + (i - 1) * deltanu + nuoffset;
      break;
    }
    //printout("[debug] macroatom: zrand %g, step %d, alpha_sp %g, total_alpha_sp %g, alpha_sp/total_alpha_sp %g, nu_lower %g\n",zrand,i,alpha_sp,total_alpha_sp,alpha_sp/total_alpha_sp,nu_lower);
    alpha_sp_old = alpha_sp;
  }
  gsl_integration_workspace_free(wsp);
  if (nu_lower == nu_threshold)
  {
    nu_lower = nu_max_phixs;
  }

  pkt_ptr->nu_cmf = nu_lower;
  //printout("nu_lower %g, nu_threshold %g, nu_left %g, nu_right %g\n",nu_lower,nu_threshold,nu_threshold+(i-1)*deltanu,nu_threshold+(i)*deltanu);


  /*
  ///k-pkt emission rule
  double bfcooling_coeff,total_bfcooling_coeff,bfcooling_coeff_old;
  int ii;
  double bfcooling_integrand_gsl_2(double nu, void *paras);
  gsl_function F_bfcooling;
  F_bfcooling.function = &bfcooling_integrand_gsl_2;
  zrand = gsl_rng_uniform(rng);
  zrand = 1. - zrand;  /// Make sure that 0 < zrand <= 1
  mastate[tid].element = element;
  mastate[tid].ion = ion-1;
  mastate[tid].level = lower;
  intparas.T = T_e;
  intparas.nu_edge = nu_threshold;   /// Global variable which passes the threshold to the integrator
  F_bfcooling.params = &intparas;
  deltanu = nu_threshold * NPHIXSNUINCREMENT;
  gsl_integration_qag(&F_bfcooling, nu_threshold, nu_max_phixs, 0, intaccuracy, 1024, 6, wsp, &total_bfcooling_coeff, &error);
  bfcooling_coeff = total_bfcooling_coeff;
  for (ii= 0; ii < NPHIXSPOINTS; ii++)
  {
    bfcooling_coeff_old = bfcooling_coeff;
    if (ii > 0)
    {
      nu_lower = nu_threshold + ii*deltanu;
      /// Spontaneous recombination and bf-cooling coefficient don't depend on the cutted radiation field
      gsl_integration_qag(&F_bfcooling, nu_lower, nu_max_phixs, 0, intaccuracy, 1024, 6, wsp, &bfcooling_coeff, &error);
      //bfcooling_coeff *= FOURPI * sf;
      //if (zrand > bfcooling_coeff/get_bfcooling(element,ion,level,pkt_ptr->where)) break;
    }
    //printout("zrand %g, bfcooling_coeff %g, total_bfcooling_coeff %g, nu_lower %g\n",zrand,bfcooling_coeff,total_bfcooling_coeff,nu_lower);
    if (zrand >= bfcooling_coeff/total_bfcooling_coeff) break;
  }
  if (ii==NPHIXSPOINTS)
  {
    printout("kpkt emitts bf-photon at upper limit\n");
    nu_lower = nu_threshold * last_phixs_nuovernuedge; // + ii*deltanu;
  }
  else if (ii > 0)
  {
    nuoffset = (total_bfcooling_coeff*zrand - bfcooling_coeff_old) / (bfcooling_coeff-bfcooling_coeff_old) * deltanu;
    nu_lower = nu_threshold + (ii-1)*deltanu + nuoffset;
  }
  else
    nu_lower = nu_threshold; // + ii*deltanu;
  //printout("emitt at nu %g\n",nu_lower);
  pkt_ptr->nu_cmf = nu_lower;
    */

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
    if (debuglevel == 2) printout("[debug] do_ma:   pkt_ptr->nu_cmf %g\n",pkt_ptr->nu_cmf);
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
  emitt_rpkt(pkt_ptr, t_current);
  calculate_kappa_rpkt_cont(pkt_ptr, t_current);
  pkt_ptr->next_trans = 0;       /// continuum transition, no restrictions for further line interactions
  pkt_ptr->emissiontype = get_continuumindex(element, ion - 1, lower);
  pkt_ptr->em_pos[0] = pkt_ptr->pos[0];
  pkt_ptr->em_pos[1] = pkt_ptr->pos[1];
  pkt_ptr->em_pos[2] = pkt_ptr->pos[2];
  pkt_ptr->em_time = t_current;
  pkt_ptr->nscatterings = 0;
}


double do_macroatom(PKT *restrict pkt_ptr, const double t1, const double t2, const int timestep)
/// Material for handling activated macro atoms.
{
  double t_current = t1; // this will keep track of time in the calculation
  const double t_mid = time_step[timestep].mid;

  //printout("[debug] do MA\n");

  const int cellindex = pkt_ptr->where;
  const int modelgridindex = cell[cellindex].modelgridindex;
  const float T_e = get_Te(modelgridindex);
  const float nne = get_nne(modelgridindex);

  /// calculate occupation number for active MA level ////////////////////////////////////
  /// general QUESTION: is it better to calculate the n_1 (later the n_ionstage and
  /// U_ionstage) here where we need them or once in update_grid for each grid cell
  /// not sure whether this reduces the number of calculations, as number of grid cells
  /// is much larger than number of pellets (next question: connection to number of
  /// photons)
  const int element = mastate[tid].element;

  /// dummy-initialize these to nonsense values, if something goes wrong with the real
  /// initialization we should see errors
  double epsilon_trans = -100.;
  int upper = -100;
  int lower = -100;

  //debuglevel = 2;
  #ifdef DEBUG_ON
    if (debuglevel == 2)
      printout("[debug] =============entering do_ma\n");
    int jumps = 0;
    int jump = -99;
  #endif

  //debuglevel = 2;
  bool end_packet = false;
  while (!end_packet)
  {
    const int ion = mastate[tid].ion;
    const int level = mastate[tid].level;
    // mastate[tid].statweight = stat_weight(element, ion, level);
    //ionisinglevels = get_ionisinglevels(element,ion);

    /// Set this here to 1 to overcome problems in cells which have zero population
    /// in some ionisation stage. This is possible because the dependence on the
    /// originating levels population cancels out in the macroatom transition probabilities
    /// which are based on detailed balance.
    mastate[tid].nnlevel = 1.;
    //mastate[tid].nnlevel = get_levelpop(modelgridindex,element,ion,level);

    #ifdef DEBUG_ON
      //if (element == 1 && ion == 1 && level == 149) debuglevel = 2;
      //if (element == 1 && ion == 1 && level == 150) debuglevel = 2;
      //if (element == 1 && ion == 1 && level == 151) debuglevel = 2;
      //if (element == 1 && ion == 1 && level == 61) debuglevel = 2;
      if (get_ionstage(element,ion) == 0 && level == 0)
      {
        printout("element %d, ion %d, level %d, ionstage %d\n",element,ion,level,get_ionstage(element,ion));
        debuglevel = 2;
      }
      if (ion >= get_nions(element))
      {
        printout("[fatal] do_ma: ion does not exist ... abort\n");
        abort();
      }
      if (debuglevel == 2) printout("[debug] do_ma: element %d, ion %d, level %d:\n",element,ion,level);
      if (debuglevel == 2) printout("[debug] do_ma: jumps = %d\n",jumps);
    #endif

    const double epsilon_current = epsilon(element, ion, level);
    const int ndowntrans = elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
    const int nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
    //nlevels_nextion  ///not needed as long we only ionise to the ground state
    //nlevels_lowerion ///only needed if level = 0, this won't happen too often

    #ifdef DEBUG_ON
      if (debuglevel == 2) printout("[debug] do_ma: ndowntrans %d, nuptrans %d\n",ndowntrans,nuptrans);
      if (debuglevel == 2) printout("[debug] do_ma: col_deexc_stored %g\n",cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_deexc);
    #endif
    if (cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_deexc < 0 || level == 0)
    {
      /// If there are no precalculated rates available we must calculate them
      set_cellhistory_transitionrates(modelgridindex, element, ion, level, t_mid);
    }

    const double rad_deexc = cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_deexc;
    const double col_deexc = cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_deexc;
    const double rad_recomb = cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_recomb;
    const double col_recomb = cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_recomb;
    const double internal_down_same = cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_same;
    const double internal_up_same = cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_same;
    const double internal_down_lower = cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_lower;
    const double internal_up_higher = cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_higher;

    /// select transition according to probabilities /////////////////////////////////////

    const double total_transitions = rad_deexc + col_deexc + internal_down_same + rad_recomb + col_recomb + internal_down_lower + internal_up_same + internal_up_higher;
    double processrates[8] = {rad_deexc, col_deexc, internal_down_same, rad_recomb, col_recomb,
                              internal_down_lower, internal_up_same, internal_up_higher};
    enum ma_action actions[8] = {MA_ACTION_RADDEEXC, MA_ACTION_COLDEEXC, MA_ACTION_INTERNALDOWNSAME, MA_ACTION_RADRECOMB, MA_ACTION_COLRECOMB,
                                 MA_ACTION_INTERNALDOWNLOWER, MA_ACTION_INTERNALUPSAME, MA_ACTION_INTERNALUPHIGHER};
    double zrand = gsl_rng_uniform(rng);
    //printout("zrand %g\n",zrand);
    const double randomrate = zrand * total_transitions;
    double rate = 0.;
    for (int p = 0; p < 8; p++)
    {
      rate += processrates[p];
      if (randomrate < rate)
      {
        mastate[tid].lastaction = actions[p];
        break;
      }
      else if (p == 8)
      {
        printout("[fatal] do_ma: problem with random numbers .. abort\n");
        //printout("[debug]    rad_down %g, col_down %g, internal_down %g, internal_up %g\n",rad_down,col_down,internal_down,internal_up);
        printout("[debug] do_ma: element %d, ion %d, level %d, total_transitions %g\n",element,ion,level,total_transitions);
        printout("[debug]    col_deexc %g, col_recomb %g, rad_deexc %g, rad_recomb %g\n",col_deexc,col_recomb,rad_deexc,rad_recomb);
        printout("[debug]    internal_down_same %g, internal_down_lower %g\n",internal_down_same,internal_down_lower);
        printout("[debug]    internal_up_same %g, internal_up_higher %g\n",internal_up_same,internal_up_higher);
        printout("[debug]    zrand %g\n",zrand);
        printout("[debug]    jumps %d\n",jumps);
        printout("[debug]    pkt_ptr->number %d\n",pkt_ptr->number);

        debuglevel = 777;

        if (ion > 0) ///checks only if there is a lower ion, doesn't make sure that Z(ion)=Z(ion-1)+1
        {
          printout("[debug]    check recombination\n");
          //nlevels = get_nlevels(element,ion-1);
          const int nlevels = get_bfcontinua(element, ion - 1);
          //nlevels = get_ionisinglevels(element,ion-1);
          for (lower = 0; lower < nlevels; lower++)
          {
            const double epsilon_target = epsilon(element, ion - 1,lower);
            epsilon_trans = epsilon_current - epsilon_target;
            const double R = rad_recombination_ratecoeff(modelgridindex, element, ion, level, lower);
            const double C = col_recombination_ratecoeff(T_e, nne, element, ion, level, lower, epsilon_trans);
            printout("[debug]    recombination to ion %d, level %d, epsilon_target %g, epsilon_trans %g, R %g, C %g\n",ion-1,lower,epsilon_target,epsilon_trans,R,C);
          }
        }

        printout("[debug]    check deexcitation\n");
        printout("[debug]    ndowntrans %d %d\n",ndowntrans,elements[element].ions[ion].levels[level].downtrans[0].targetlevel);
        for (int i = 1; i <= ndowntrans; i++)
        {
          lower = elements[element].ions[ion].levels[level].downtrans[i].targetlevel;
          epsilon_trans = elements[element].ions[ion].levels[level].downtrans[i].epsilon_trans;
          const int lineindex = elements[element].ions[ion].levels[level].downtrans[i].lineindex;
          const double R = rad_deexcitation_ratecoeff(modelgridindex, element, ion, level, lower, epsilon_trans, lineindex, t_mid);
          const double C = col_deexcitation_ratecoeff(T_e, nne, epsilon_trans, lineindex);
          printout("[debug]    deexcitation to level %d, epsilon_trans %g, epsilon_trans %g, R %g, C %g\n",lower,epsilon_trans,epsilon_trans,R,C);
        }

        printout("[debug]    check excitation\n");
        printout("[debug]    nuptrans %d %d\n",nuptrans,elements[element].ions[ion].levels[level].uptrans[0].targetlevel);
        for (int i = 1; i <= nuptrans; i++)
        {
          upper = elements[element].ions[ion].levels[level].uptrans[i].targetlevel;
          epsilon_trans = elements[element].ions[ion].levels[level].uptrans[i].epsilon_trans;
          const int lineindex = elements[element].ions[ion].levels[level].uptrans[i].lineindex;
          const double R = rad_excitation_ratecoeff(modelgridindex, element, ion, lower, upper, epsilon_trans, lineindex, t_mid);
          const double C = col_excitation_ratecoeff(T_e, nne, lineindex, epsilon_trans);
          printout("[debug]    excitation to level %d, epsilon_trans %g, R %g, C %g\n",upper,epsilon_trans,R,C);
        }

        if (ion < get_nions(element) - 1)  //&& get_ionstage(element,ion) < get_element(element)+1)
        {
          printout("[debug]    check ionisation\n");
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++)
          {
            upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
            epsilon_trans = epsilon(element, ion + 1, upper) - epsilon_current;
            const double R = get_corrphotoioncoeff(element, ion, level, phixstargetindex, modelgridindex);
            const double C = col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);
            printout("[debug]    ionisation to ion %d, level %d, epsilon_trans %g, R %g, C %g\n", ion + 1, upper, epsilon_trans, R, C);
            break;
          }
        }

        abort();
      }
    }

    #ifdef DEBUG_ON
      if (debuglevel == 2)
      {
        printout("[debug] do_ma: element %d, ion %d, level %d\n",element,ion,level);
        printout("[debug] do_ma:   rad_deexc %g, col_deexc %g, internal_down_same %g, rad_recomb %g, col_recomb %g, internal_down_lower %g, internal_up_same %g, internal_up_higher %g\n",rad_deexc,col_deexc,internal_down_same,rad_recomb,col_recomb, internal_down_lower, internal_up_same,internal_up_higher);
      }
      if (total_transitions <= 0.)
      {
        printout("[debug] do_ma: element %d, ion %d, level %d, total_transitions = %g\n",element,ion,level,total_transitions);
        printout("[debug]    col_deexc %g, col_recomb %g, rad_deexc %g, rad_recomb %g\n",col_deexc,col_recomb,rad_deexc,rad_recomb);
        printout("[debug]    internal_down_same %g, internal_down_lower %g\n",internal_down_same,internal_down_lower);
        printout("[debug]    internal_up_same %g, internal_up_higher %g\n",internal_up_same,internal_up_higher);
        printout("[debug]    zrand %g\n",zrand);
        printout("[debug]    jumps %d, jump %d\n",jumps,jump);
        printout("[debug]    pkt_ptr->number %d, pkt_ptr->where %d\n",pkt_ptr->number,cellindex);
        printout("[debug]    groundlevelpop of current ion in current cell %g\n",modelgrid[cell[cellindex].modelgridindex].composition[element].groundlevelpop[ion]);
        printout("[debug]    levelpop %g\n",mastate[tid].nnlevel);

        double R = 0.0;
        double C = 0.0;
        for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++)
        {
          upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
          epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
          R += get_corrphotoioncoeff(element,ion,level,phixstargetindex,modelgridindex);
          C += col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);
          printout("epsilon_current %g, epsilon_trans %g, photion %g, colion %g, internal_up_higher %g, saved_internal_up_higher %g\n",epsilon_current,epsilon_trans,R,C,(R + C) * epsilon_current,cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_higher);
        }

        float T_R = get_TR(modelgridindex);
        float W = get_W(modelgridindex);
        printout("modelgridindex %d, T_R %g, T_e %g, W %g, T_J %g\n",modelgridindex,T_R,get_Te(modelgridindex),W,get_TJ(modelgridindex));
        #if (!NO_LUT_PHOTOION)
          double gammacorr = W * interpolate_corrphotoioncoeff(element,ion,level,0,T_R);
          int index_in_groundlevelcontestimor = elements[element].ions[ion].levels[level].closestgroundlevelcont;
          double renorm  = corrphotoionrenorm[modelgridindex * nelements * maxion + index_in_groundlevelcontestimor];
          printout("gammacorr %g, index %d, renorm %g, total %g\n",gammacorr,index_in_groundlevelcontestimor,renorm,gammacorr*renorm);
        #endif

        //abort();
      }
    #endif

    switch (mastate[tid].lastaction)
    {
      case MA_ACTION_RADDEEXC:
        #ifdef DEBUG_ON
          if (debuglevel == 2)
          {
            printout("[debug] do_ma:   radiative deexcitation\n");
            printout("[debug] do_ma:   jumps = %d\n",jumps);
          }
        #endif
        do_macroatom_raddeexcitation(pkt_ptr, element, ion, level, rad_deexc, total_transitions, t_current);
        end_packet = true;
        break;

      case MA_ACTION_COLDEEXC:
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

      case MA_ACTION_INTERNALDOWNSAME:
        #ifdef DEBUG_ON
          if (debuglevel == 2)
            printout("[debug] do_ma:   internal downward jump within current ionstage\n");
          pkt_ptr->interactions += 1;
          jumps++;
          jump = 0;
        #endif

        /// Randomly select the occuring transition
        zrand = gsl_rng_uniform(rng);
        lower = -99;
        rate = 0.;
        for (int i = 1; i <= ndowntrans; i++)
        {
          rate += get_individ_internal_down_same(element, ion, level, i);
          if (zrand * internal_down_same < rate)
          {
            lower = elements[element].ions[ion].levels[level].downtrans[i].targetlevel;
            break;
          }
        }
        /// and set the macroatom's new state
        mastate[tid].ion = ion;
        mastate[tid].level = lower;

        #ifdef DEBUG_ON
          if (debuglevel == 2)
            printout("[debug] do_ma:   to level %d\n", lower);
          if (get_ionstage(element,ion) == 0 && lower == 0)
          {
            printout("internal downward transition to ground level occured ... abort\n");
            printout("element %d, ion %d, level %d, lower %d\n", element, ion, level, lower);
            printout("Z %d, ionstage %d, energy %g\n",
                     elements[element].anumber, get_ionstage(element,ion), elements[element].ions[ion].levels[lower].epsilon);
            printout("[debug] do_ma:   internal downward jump within current ionstage\n");
            abort();
          }
        #endif
        break;

      case MA_ACTION_RADRECOMB:
        /// Radiative recombination of MA: emitt a continuum-rpkt
        #ifdef DEBUG_ON
          if (debuglevel == 2)
          {
            printout("[debug] do_ma:   radiative recombination\n");
            printout("[debug] do_ma:   jumps = %d\n",jumps);
            printout("[debug] do_ma:   element %d, ion %d, level %d\n",element,ion,level);
          }
        #endif
        do_macroatom_radrecomb(pkt_ptr, modelgridindex, element, ion, level, rad_recomb, t_current);
        end_packet = true;
        break;

      case MA_ACTION_COLRECOMB:
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

      case MA_ACTION_INTERNALDOWNLOWER:
        #ifdef DEBUG_ON
          if (debuglevel == 2)
            printout("[debug] do_ma:   internal downward jump to lower ionstage\n");
          pkt_ptr->interactions += 1;
          jumps++;
          jump = 1;
        #endif

        /// Randomly select the occuring transition
        zrand = gsl_rng_uniform(rng);
        //zrand = 1. - 1e-14;
        rate = 0.;
        //nlevels = get_nlevels(element,ion-1);
        const int nlevels = get_bfcontinua(element, ion - 1);
        //nlevels = get_ionisinglevels(element,ion-1);
        for (lower = 0; lower < nlevels; lower++)
        {
          const double epsilon_target = epsilon(element, ion - 1, lower);
          epsilon_trans = epsilon_current - epsilon_target;
          const double R = rad_recombination_ratecoeff(modelgridindex, element, ion, level, lower);
          const double C = col_recombination_ratecoeff(T_e, nne, element, ion, level, lower, epsilon_trans);
          rate += (R + C) * epsilon_target;
          if (zrand * internal_down_lower < rate)
            break;
        }
        /// and set the macroatom's new state
        mastate[tid].ion = ion - 1;
        mastate[tid].level = lower;

        #ifdef DEBUG_ON
          if (lower >= nlevels)
          {
            printout("internal_down_lower  %g\n",internal_down_lower);
            printout("abort at rate %g, zrand %g\n",rate,zrand);
            abort();
          }
          if (get_ionstage(element,ion-1) == 0 && lower == 0)
  	  //        if (ion-1 == 0 && lower == 0)
          {
            printout("internal downward transition to ground level occured ... abort\n");
            printout("element %d, ion %d, level %d, lower %d\n",element,ion,level,lower);
            printout("Z %d, ionstage %d, energy %g\n",elements[element].anumber,get_ionstage(element,ion-1),elements[element].ions[ion-1].levels[lower].epsilon);
            printout("[debug] do_ma:   internal downward jump to lower ionstage\n");
            abort();
          }
        #endif
        break;

      case MA_ACTION_INTERNALUPSAME:
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] do_ma:   internal upward jump within current ionstage\n");
          pkt_ptr->interactions += 1;
          jumps++;
          jump = 2;
        #endif

        ///randomly select the occuring transition
        zrand = gsl_rng_uniform(rng);
        upper = -99;
        rate = 0.;
        for (int i = 1; i <= nuptrans; i++)
        {
          rate += get_individ_internal_up_same(element, ion, level, i);
          if (zrand*internal_up_same < rate)
          {
            upper = elements[element].ions[ion].levels[level].uptrans[i].targetlevel;
            break;
          }
        }
        ///and set the macroatom's new state
        mastate[tid].ion = ion;
        mastate[tid].level = upper;
        break;

      case MA_ACTION_INTERNALUPHIGHER:
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] do_ma:   internal upward jump to next ionstage\n");
          pkt_ptr->interactions += 1;
          jumps++;
          jump = 3;
        #endif

        /// Randomly select the occuring transition
        zrand = gsl_rng_uniform(rng);
        rate = 0.;
        // if (NT_ON)
        //   rate += nt_ionization_ratecoeff(modelgridindex, element, ion);
        if (zrand * internal_up_higher < rate)
        {
          upper = 0;
          epsilon_trans = epsilon(element, ion + 1, upper) - epsilon_current;
        }
        else
        {
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++)
          {
            upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
            epsilon_trans = epsilon(element, ion + 1, upper) - epsilon_current;
            const double R = get_corrphotoioncoeff(element, ion, level, phixstargetindex, modelgridindex);
            const double C = col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);
            rate += (R + C) * epsilon_current;
            if (zrand * internal_up_higher < rate)
              break;
          }
        }

        /// and set the macroatom's new state
        mastate[tid].ion = ion + 1;
        mastate[tid].level = upper;
        break;

      default:
        break;
    }
  }///endwhile

  #ifdef DEBUG_ON
    //debuglevel = 4;
    //debuglevel = 2;
  #endif

  /// procedure ends only after a change to r or k packets has taken place and
  /// returns then the actual time, which is the same as the input t1
  /// internal transitions are carried out until a type change occurs
  return t_current;
}


/// Calculation of radiative rates ///////////////////////////////////////////////////////


double rad_deexcitation_ratecoeff(int modelgridindex, int element, int ion, int upper, int lower, double epsilon_trans, int lineindex, double t_current)
///radiative deexcitation rate: paperII 3.5.2
{
  #ifdef DEBUG_ON
  if (upper <= lower)
  {
    printout("[fatal] rad_deexcitation: tried to calculate upward transition ... abort\n");
    abort();
  }
  #endif

  const double statweight_target = statw_lower(lineindex);
  const double statweight = statw_upper(lineindex);

  const double n_u = get_levelpop(modelgridindex,element,ion,upper);
  const double n_l = get_levelpop(modelgridindex,element,ion,lower);

  double R = 0.0;

  if ((n_u > 1.1 * MINPOP) && (n_l > 1.1 * MINPOP))
  {
    const double nu_trans = epsilon_trans / H;

    const double A_ul = einstein_spontaneous_emission(lineindex);
    const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
    const double B_lu = statweight / statweight_target * B_ul;

    const double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current;

    if (tau_sobolev > 0)
    {
      const double beta = 1.0 / tau_sobolev * (1 - exp(-tau_sobolev));
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
        printout("[debug] rad_rates_down: element, ion, upper, lower %d, %d, %d, %d\n",element,ion,upper,lower);
        printout("[debug] rad_deexc: A_ul %g, tau_sobolev %g, n_u %g\n",A_ul,tau_sobolev,n_u);
      }
      if (debuglevel == 777)
        printout("[debug] rad_deexc: A_ul %g, tau_sobolev %g, n_u %g\n",A_ul,tau_sobolev,n_u);
      if (!isfinite(R))
      {
        printout("fatal a1: abort\n");
        abort();
      }
    #endif
  }
  else
  {
    R = 0.0;
  }

  return R;
}


double rad_excitation_ratecoeff(int modelgridindex, int element, int ion, int lower, int upper, double epsilon_trans, int lineindex, double t_current)
///radiative excitation rate: paperII 3.5.2
{
  #ifdef DEBUG_ON
  if (upper <= lower)
  {
    printout("[fatal] rad_excitation: tried to calculate downward transition ... abort\n");
    abort();
  }
  #endif

  const double statweight = statw_lower(lineindex);
  const double statweight_target = statw_upper(lineindex);

  const double n_u = get_levelpop(modelgridindex,element,ion,upper);
  const double n_l = get_levelpop(modelgridindex,element,ion,lower);
  double R = 0.0;
  if ((n_u >= 1.1 * MINPOP) && (n_l >= 1.1 * MINPOP))
  {
    const double nu_trans = epsilon_trans / H;
    const double A_ul = einstein_spontaneous_emission(lineindex);
    const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans,3) * A_ul;
    const double B_lu = statweight_target / statweight * B_ul;

    const double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current;

    if (tau_sobolev > 0)
    {
      double beta = 1.0 / tau_sobolev * (1. - exp(-tau_sobolev));
      //printout("[check] rad_excitation: %g, %g, %g\n",1.0/tau_sobolev,exp(-tau_sobolev),1.0/tau_sobolev * (1. - exp(-tau_sobolev)));
      //n_u2 = calculate_levelpop_fromreflevel(pkt_ptr->where,element,ion,upper,lower,mastate[tid].nnlevel);
      //R = (B_lu*mastate[tid].nnlevel - B_ul * n_u2) * beta * radfield(nu_trans,pkt_ptr->where);
      R = (B_lu - B_ul * n_u / n_l) * beta * radfield(nu_trans, modelgridindex);
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
      const double g_u2 = stat_weight(element,ion,upper);
      const double g_l = statw_lower(lineindex);
      const double g_l2 = stat_weight(element,ion,lower);
      printout("Negative excitation rate from level %d to %d\n",lower,upper);
      printout("n_l %g, n_u %g, g_l %g (?=%g), g_u %g (?=%g)\n",n_l,n_u,g_l,g_l2,g_u,g_u2);
      printout("n_u/n_l %g, g_u/g_l %g\n",n_u/n_l,g_u/g_l);
      printout("radfield(nutrans=%g) = %g\n",nu_trans,radfield(nu_trans,modelgridindex));
      abort();
    }
    if (debuglevel == 2)
    {
      printout("[debug] rad_rates_up: element, ion, upper, lower, A_ul, n_u: %d, %d, %d, %d, %g, %g\n",element,ion,upper,lower,A_ul,n_l);
      printout("[debug] rad_exc: A_ul %g, tau_sobolev %g, n_u %g, n_l %g, radfield %g\n",A_ul,tau_sobolev,n_u,n_l,radfield(nu_trans,modelgridindex));
    }
    if (debuglevel == 777)
      printout("[debug] rad_exc: A_ul %g, tau_sobolev %g, n_u %g, n_l %g, radfield %g\n",A_ul,tau_sobolev,n_u,n_l,radfield(nu_trans,modelgridindex));
    if (!isfinite(R))
    {
      printout("[fatal] rad_excitation: abort\n");
      printout("[fatal] rad_excitation: R %g, mastate[tid].nnlevel %g, B_lu %g, B_ul %g, n_u %g, n_l %g, radfield %g,tau_sobolev %g, t_current %g\n",R,mastate[tid].nnlevel,B_lu,B_ul,n_u,n_l,radfield(nu_trans,modelgridindex),tau_sobolev,t_current);
      printout("[fatal] rad_excitation: %g, %g, %g\n",1.0/tau_sobolev,exp(-tau_sobolev),1.0/tau_sobolev * (1. - exp(-tau_sobolev)));
      abort();
    }
    #endif
  }
  else
  {
    R = 0.0;
  }

  return R;
}


double rad_recombination_ratecoeff(int modelgridindex, int element, int upperion, int upper, int lower)
///radiative recombination rate: paperII 3.5.2
{
  const float nne = get_nne(modelgridindex);
  double R = 0.0;
  const int nphixstargets = get_nphixstargets(element,upperion-1,lower);
  for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++)
  {
    if (get_phixsupperlevel(element,upperion-1,lower,phixstargetindex) == upper)
    {
      R = nne * get_spontrecombcoeff(element,upperion-1,lower,phixstargetindex,modelgridindex);// + stimrecombestimator_save[pkt_ptr->where*nelements*maxion+element*maxion+(ion-1)]);
      //printout("calculate rad_recombination: element %d, ion %d, upper %d, -> lower %d, n_u %g, nne %g, spontrecombcoeff %g\n",element,ion,upper,lower,mastate[tid].nnlevel,nne,get_spontrecombcoeff(element, ion-1, lower, pkt_ptr->where));
      break;
    }
  }

  #ifdef DEBUG_ON
    //printout("[debug]    rad_recombiantion: R %g\n",R);
    if (!isfinite(R))
    {
      printout("fatal a2: abort\n");
      abort();
    }
  #endif

  return R;
}


/// Calculation of collisional rates /////////////////////////////////////////////////////
///***************************************************************************/


double col_deexcitation_ratecoeff(float T_e, float nne, double epsilon_trans, int lineindex)
{
  double C;
  const double coll_str_thisline = get_coll_str(lineindex);
  const double statweight = statw_upper(lineindex);
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

      const double fac1 = epsilon_trans / KB / T_e;
      ///Van-Regemorter formula, Mihalas (1978), eq.5-75, p.133
      const double g_bar = 0.2; ///this should be read in from transitions data: it is 0.2 for transitions nl -> n'l' and 0.7 for transitions nl -> nl'
      //test = 0.276 * exp(fac1) * gsl_sf_expint_E1(fac1);
      /// crude approximation to the already crude Van-Regemorter formula

      //double test = 0.276 * exp(fac1) * (-0.5772156649 - log(fac1));
      //double Gamma = (g_bar > test) ? g_bar : test;

      //optimisation
      const double gauntfac = (fac1 > 0.33421) ? g_bar : 0.276 * exp(fac1) * (-0.5772156649 - log(fac1));

      const double g_ratio = statweight_target / statweight;

      C = C_0 * 14.51039491 * nne * sqrt(T_e) * osc_strength(lineindex) * pow(H_ionpot/epsilon_trans,2) * fac1 * g_ratio * gauntfac;
    }
    else // alterative: (coll_strength > -3.5) to catch -2 or -3
    {
      //forbidden transitions: magnetic dipole, electric quadropole...
      //could be Axelrod? or Maurer
      C = nne * 8.629e-6 * pow(T_e,-0.5) * 0.01 * statweight_target;
    }
  }
  else //positive values are treated as effective collision strengths
  {
    //from Osterbrock and Ferland, p51
    //mastate[tid].statweight is UPPER LEVEL stat weight
    //statweight_target is LOWER LEVEL stat weight
    C = nne * 8.629e-6 * pow(T_e, -0.5) * coll_str_thisline / statweight;
    // test test
    //C = n_u * nne * 8.629e-6 * pow(T_e,-0.5) * 0.01 * statweight_target;
  }

  #ifdef DEBUG_ON
    //int element = mastate[tid].element;
    //int ion = mastate[tid].ion;
    /*if (debuglevel == 2)
    {
      //printout("[debug] col_deexc: element %d, ion %d, upper %d, lower %d\n",element,ion,upper,lower);
      printout("[debug] col_deexc: n_u %g, nne %g, T_e %g, f_ul %g, epsilon_trans %g, Gamma %g, g_ratio %g\n",n_u,nne,T_e,osc_strength(lineindex),epsilon_trans,Gamma,g_ratio);
    }*/
    //printout("col_deexc(%d,%d,%d,%d) %g\n",element,ion,upper,lower,C);
    if (!isfinite(C))
    {
      printout("fatal a7: abort\n");
      abort();
    }
  #endif

  return C;
}


double col_excitation_ratecoeff(float T_e, float nne, int lineindex, double epsilon_trans)
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

      const double test = 0.276 * exp(eoverkt) * (-0.5772156649 - log(eoverkt));
      const double Gamma = g_bar > test ? g_bar : test;
      C = C_0 * nne * sqrt(T_e) * 14.51039491 * osc_strength(lineindex) * pow(H_ionpot/epsilon_trans, 2) * eoverkt * exp(-eoverkt) * Gamma;
    }
    else // alterative: (coll_strength > -3.5) to catch -2 or -3
    {
      // forbidden transitions: magnetic dipole, electric quadropole...
      C = nne * 8.629e-6 * 0.01 * pow(T_e,-0.5) * exp(-eoverkt) * statw_upper(lineindex);
    }
  }
  else
  {
    //from Osterbrock and Ferland, p51
    C = nne * 8.629e-6 * pow(T_e, -0.5) * coll_strength * exp(-eoverkt) / statw_lower(lineindex);
    //test test
    //C = n_l * nne * 8.629e-6 * pow(T_e,-0.5) * 0.01 * exp(-fac1) * statw_upper(lineindex);
  }

  #ifdef DEBUG_ON
    //int element = mastate[tid].element;
    //int ion = mastate[tid].ion
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


double col_recombination_ratecoeff(float T_e, float nne, int element, int upperion, int upper, int lower, double epsilon_trans)
{
  const double fac1 = epsilon_trans / KB / T_e;
  const int ionstage = get_ionstage(element,upperion);
  const float sigma_bf_all_targets = elements[element].ions[upperion-1].levels[lower].photoion_xs[0];

  const int nphixstargets = get_nphixstargets(element,upperion-1,lower);
  for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++)
  {
    if (get_phixsupperlevel(element,upperion-1,lower,phixstargetindex) == upper)
    {
      ///Seaton approximation: Mihalas (1978), eq.5-79, p.134
      ///select gaunt factor according to ionic charge
      double g;
      if (ionstage-1 == 1)
        g = 0.1;
      else if (ionstage-1 == 2)
        g = 0.2;
      else
        g = 0.3;

      const double sigma_bf = sigma_bf_all_targets * get_phixsprobability(element,upperion-1,lower,phixstargetindex);

      const double C = nne * nne * get_sahafact(element,upperion-1,lower,phixstargetindex,T_e,epsilon_trans) *
                       1.55e13 * pow(T_e,-0.5) * g * sigma_bf * exp(-fac1) / fac1;

      #ifdef DEBUG_ON
        /*if (debuglevel == 777)
        {
          printout("get_sahafact %g, fac1 %g, C %g\n",get_sahafact(element,ion-1,lower,phixstargetindex,T_e,epsilon_trans),fac1,C);
          ///means n_u*nne * detailed_balancing of c_ikappa
          printout("[debug] col_recomb: n_u %g, nne %g, T_e %g, g %g, epsilon_trans %g, sigma_bf %g\n",n_u, nne,T_e,g,epsilon_trans,sigma_bf);
        }
        if (!isfinite(C))
        {
          printout("fatal a8: abort\n");
          abort();
        }*/
      #endif

      return C;
    }
  }

  return 0.0;
}


double col_ionization_ratecoeff(float T_e, float nne, int element, int ion, int lower, int phixstargetindex, double epsilon_trans)
/// collisional ionization rate: paperII 3.5.1
{
  #ifdef DEBUG_ON
  if (phixstargetindex > get_nphixstargets(element,ion,lower))
  {
    printout("[fatal] col_ionization called with phixstargetindex %g > nphixstargets %g",phixstargetindex,get_nphixstargets(element,ion,lower));
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

  const float sigma_bf_all_targets = elements[element].ions[ion].levels[lower].photoion_xs[0];
  const double sigma_bf = sigma_bf_all_targets * get_phixsprobability(element,ion,lower,phixstargetindex);
  const double C = nne * 1.55e13 * pow(T_e,-0.5) * g * sigma_bf * exp(-fac1) / fac1; ///photoionization at the edge

  #ifdef DEBUG_ON
    if (debuglevel == 777)
      printout("[debug] col_ion: nne %g, T_e %g, g %g, epsilon_trans %g, sigma_bf %g\n", nne,T_e,g,epsilon_trans,sigma_bf);
    if (!isfinite(C))
    {
      printout("fatal a6: abort\n");
      abort();
    }
  #endif

  return C;
}


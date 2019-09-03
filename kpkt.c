#include "sn3d.h"
#include <gsl/gsl_integration.h>
#include "atomic.h"
#include "grid_init.h"
#include "kpkt.h"
#include "ltepop.h"
#include "macroatom.h"
#include "ratecoeff.h"
#include "rpkt.h"
#include "vectors.h"

extern inline int get_coolinglistoffset(int element, int ion);


static int get_ncoolingterms(int element, int ion)
{
  return elements[element].ions[ion].ncoolingterms;
}


static double get_cooling_ion_coll_exc(const int modelgridindex, const int element, const int ion, const double T_e, const double nne)
{
  double C_exc = 0.;
  const int nlevels = get_nlevels(element, ion);

  /// excitation to same ionization stage
  /// -----------------------------------
  for (int level = 0; level < nlevels; level++)
  {
    const double nnlevel = calculate_exclevelpop(modelgridindex, element, ion, level);
    const double epsilon_current = epsilon(element,ion,level);
    const int nuptrans = get_nuptrans(element, ion, level);
    for (int ii = 0; ii < nuptrans; ii++)
    {
      const int lineindex = elements[element].ions[ion].levels[level].uptrans_lineindicies[ii];
      const int upper = linelist[lineindex].upperlevelindex;
      //printout("    excitation to level %d possible\n",upper);
      const double epsilon_trans = epsilon(element,ion,upper) - epsilon_current;
      const double C = nnlevel * col_excitation_ratecoeff(T_e, nne, lineindex, epsilon_trans) * epsilon_trans;
      C_exc += C;
    }
  }

  return C_exc;
}


static double get_bfcoolingcoeff(int element, int ion, int level, int phixstargetindex, float T_e)
{
  const int lowerindex = floor(log(T_e / MINTEMP) / T_step_log);
  if (lowerindex < TABLESIZE-1)
  {
    const int upperindex = lowerindex + 1;
    const double T_lower =  MINTEMP * exp(lowerindex * T_step_log);
    const double T_upper =  MINTEMP * exp(upperindex * T_step_log);

    const double f_upper = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfcooling_coeff[upperindex];
    const double f_lower = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfcooling_coeff[lowerindex];

    return (f_lower + (f_upper - f_lower) / (T_upper - T_lower) * (T_e - T_lower));
  }
  else
    return elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfcooling_coeff[TABLESIZE-1];
}


void calculate_cooling_rates(const int modelgridindex, heatingcoolingrates_t *heatingcoolingrates)
// Calculate the cooling rates for a given cell and store them for each ion
// optionally store components (ff, bf, collisional) in heatingcoolingrates struct
{
  const float nne = get_nne(modelgridindex);
  const float T_e = get_Te(modelgridindex);

  double C_total = 0.;
  double C_ff_all = 0.;   /// free-free creation of rpkts
  double C_fb_all = 0.;   /// free-bound creation of rpkt
  double C_exc_all = 0.;  /// collisional excitation of macroatoms
  double C_ionization_all = 0.;  /// collisional ionisation of macroatoms
  for (int element = 0; element < nelements; element++)
  {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++)
    {
      double C_ion = 0.;   /// all cooling for an ion
      const int nionisinglevels = get_ionisinglevels(element, ion);
      const double nncurrention = ionstagepop(modelgridindex, element, ion);

      /// ff creation of rpkt
      const int ioncharge = get_ionstage(element,ion) - 1;
      if (ioncharge > 0)
      {
        const double C_ff_ion = 1.426e-27 * sqrt(T_e) * pow(ioncharge, 2) * nncurrention * nne;
        C_ff_all += C_ff_ion;
        C_ion += C_ff_ion;
      }

      const double C_exc_ion = get_cooling_ion_coll_exc(modelgridindex, element, ion, T_e, nne);
      C_exc_all += C_exc_ion;
      C_ion += C_exc_ion;

      if (ion < nions - 1)
      {
        for (int level = 0; level < nionisinglevels; level++)
        {
          //printout("[debug] do_kpkt: element %d, ion %d, level %d\n",element,ion,level);
          const double epsilon_current = epsilon(element,ion,level);
          const double nnlevel = calculate_exclevelpop(modelgridindex, element, ion, level);
          //printout("    ionisation possible\n");
          /// ionization to higher ionization stage
          /// -------------------------------------
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            const int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            const double epsilon_trans = epsilon(element, ion + 1, upper) - epsilon_current;
            //printout("cooling list: col_ionization\n");
            const double C_ionization_ion_thistarget = nnlevel * col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans) * epsilon_trans;
            C_ionization_all += C_ionization_ion_thistarget;
            C_ion += C_ionization_ion_thistarget;
          }

          /// fb creation of r-pkt
          /// free bound rates are calculated from the lower ion, but associated to the higher ion
          /// --------------------
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            const int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            const double nnupperlevel = calculate_exclevelpop(modelgridindex, element, ion + 1, upper);

            const double C_fb_ion_thistarget = get_bfcoolingcoeff(element, ion, level, phixstargetindex, T_e) * nnupperlevel * nne;
            C_fb_all += C_fb_ion_thistarget;
            C_ion += C_fb_ion_thistarget;
          }
        }
      }

      C_total += C_ion;
      modelgrid[modelgridindex].cooling[element].contrib[ion] = C_ion;
    }
  }
  modelgrid[modelgridindex].totalcooling = C_total;

  // only used in the T_e solver and write_to_estimators file
  if (heatingcoolingrates != NULL)
  {
    heatingcoolingrates->cooling_collisional = C_exc_all + C_ionization_all;
    heatingcoolingrates->cooling_fb = C_fb_all;
    heatingcoolingrates->cooling_ff = C_ff_all;
  }
}


static void calculate_kpkt_rates_ion(int modelgridindex, int element, int ion, int indexionstart, double oldcoolingsum)
/// Set up the global cooling list and determine the important entries
/// by applying a cut to the total cooling rate. Then sort the global
/// cooling list by the strength of the individual process contribution.
{
  const float nne = get_nne(modelgridindex);
  const float T_e = get_Te(modelgridindex);
  //double T_R = get_TR(modelgridindex);
  //double W = get_W(modelgridindex);

  /// calculate rates for
  //C_ff = 0.;   /// free-free creation of rpkts
  //C_fb = 0.;   /// free-bound creation of rpkt
  //C_exc = 0.;  /// collisional excitation of macroatoms
  //C_ion = 0.;  /// collisional ionisation of macroatoms
  double contrib = oldcoolingsum;
  int i = indexionstart;

  //printout("[debug] do_kpkt: element %d\n",element);
  const int nions = get_nions(element);
  const int nlevels_currention = get_nlevels(element,ion);

  const int ionisinglevels = get_ionisinglevels(element,ion);
  const double nncurrention = ionstagepop(modelgridindex,element,ion);

  /// ff creation of rpkt
  /// -------------------
  const int ioncharge = get_ionstage(element,ion) - 1;
  //printout("[debug] ioncharge %d, nncurrention %g, nne %g\n",ion,nncurrention,nne);
  if (ioncharge > 0)
  {
    const double C = 1.426e-27 * sqrt(T_e) * pow(ioncharge,2) * nncurrention * nne;
    //C_ff += C;
    //C_ion += C;
    //cellhistory[tid].coolinglist[i].contribution = C;
    contrib += C;
    cellhistory[tid].coolinglist[i].contribution = contrib;
    cellhistory[tid].coolinglist[i].type = COOLINGTYPE_FF;
    cellhistory[tid].coolinglist[i].element = element;
    cellhistory[tid].coolinglist[i].ion = ion;
    cellhistory[tid].coolinglist[i].level = -99;
    cellhistory[tid].coolinglist[i].upperlevel = -99;
    //if (contrib < oldcoolingsum) printout("contrib %g < oldcoolingsum %g: C%g, element %d, ion %d, level %d, coolingtype %d, i %d, low %d\n",contrib,oldcoolingsum,C,element,ion,-99,cellhistory[tid].coolinglist[i].type,i,low);
    i++;
  }

  for (int level = 0; level < nlevels_currention; level++)
  {
    //printout("[debug] do_kpkt: element %d, ion %d, level %d\n",element,ion,level);
    const double epsilon_current = epsilon(element,ion,level);
    double nnlevel = calculate_exclevelpop(modelgridindex, element, ion, level);

    /// excitation to same ionization stage
    /// -----------------------------------
    int nuptrans = get_nuptrans(element, ion, level);
    for (int ii = 0; ii < nuptrans; ii++)
    {
      const int lineindex = elements[element].ions[ion].levels[level].uptrans_lineindicies[ii];
      const int upper = linelist[lineindex].upperlevelindex;
      //printout("    excitation to level %d possible\n",upper);
      const double epsilon_trans = epsilon(element, ion, upper) - epsilon_current;
      const double C = nnlevel * col_excitation_ratecoeff(T_e, nne, lineindex, epsilon_trans) * epsilon_trans;
      contrib += C;
      cellhistory[tid].coolinglist[i].contribution = contrib;
      cellhistory[tid].coolinglist[i].type = COOLINGTYPE_COLLEXC;
      cellhistory[tid].coolinglist[i].element = element;
      cellhistory[tid].coolinglist[i].ion = ion;
      cellhistory[tid].coolinglist[i].level = level;
      cellhistory[tid].coolinglist[i].upperlevel = upper;
      cellhistory[tid].coolinglist[i].lineindex = lineindex;
      //if (contrib < oldcoolingsum) printout("contrib %g < oldcoolingsum %g: C%g, element %d, ion %d, level %d, coolingtype %d, i %d, low %d\n",contrib,oldcoolingsum,C,element,ion,level,cellhistory[tid].coolinglist[i].type,i,low);
      i++;
      //linecounter++;
    }

    if (ion < (nions - 1) && level < ionisinglevels) ///check whether further ionisation stage available
    {
      //printout("    ionisation possible\n");
      /// ionization to higher ionization stage
      /// -------------------------------------
      /// for the moment we deal only with ionisations to the next ions groundlevel
      /// to allow an easy generalisation this was only made sure by col_ionization
      /// for speed issues (reduced number of calls to epsilon) it is now done also
      /// here explicitly
      const int nphixstargets = get_nphixstargets(element, ion, level);
      for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++)
      {
        const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
        const double epsilon_upper = epsilon(element, ion + 1, upper);
        const double epsilon_trans = epsilon_upper - epsilon_current;
        //printout("cooling list: col_ionization\n");
        const double C = nnlevel * col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans) * epsilon_trans;

        contrib += C;
        cellhistory[tid].coolinglist[i].contribution = contrib;
        cellhistory[tid].coolinglist[i].type = COOLINGTYPE_COLLION;
        cellhistory[tid].coolinglist[i].element = element;
        cellhistory[tid].coolinglist[i].ion = ion;
        cellhistory[tid].coolinglist[i].level = level;
        cellhistory[tid].coolinglist[i].upperlevel = upper;
        //if (contrib < oldcoolingsum) printout("contrib %g < oldcoolingsum %g: C%g, element %d, ion %d, level %d, coolingtype %d, i %d, low %d\n",contrib,oldcoolingsum,C,element,ion,level,cellhistory[tid].coolinglist[i].type,i,low);
        i++;
      }
      //C_ion += C;
      //C = 0.;
      //cellhistory[tid].coolinglist[i].contribution = C;
      //}


      /// fb creation of r-pkt
      /// free bound rates are calculated from the lower ion, but associated to the higher ion
      /// --------------------
      for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
      {
        const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
        const double nnupperlevel = calculate_exclevelpop(modelgridindex,element,ion + 1, upper);

        const double C = get_bfcoolingcoeff(element, ion, level, phixstargetindex, T_e) * nnupperlevel * nne;
        contrib += C;
        cellhistory[tid].coolinglist[i].contribution = contrib;
        cellhistory[tid].coolinglist[i].type = COOLINGTYPE_FB;
        cellhistory[tid].coolinglist[i].element = element;
        cellhistory[tid].coolinglist[i].ion = ion;
        cellhistory[tid].coolinglist[i].level = level;
        cellhistory[tid].coolinglist[i].upperlevel = upper;
        i++;
      }
      //printout("element %d, ion %d, level %d, T_e %g, alpha_E - alpha %g, bfcoolingcoeff %g\n",element,ion,level,T_e,C/E_threshold,C);
      //double interpolate_stimulated_recomb(int element, int ion, int level, double T);
      //C = get_bfcooling(element,ion,level,modelgridindex) + (nnnextionlevel*nne * W * interpolate_stimulated_bfcoolingcoeff(element,ion,level,T_R));
      //printout("nonE %g, E %g \n",interpolate_bfcoolingcoeff(element,ion,level,T_e), interpolate_stimulated_bfcoolingcoeff(element,ion,level,T_R));
      //C = get_bfcooling(element,ion,level,modelgridindex) + (nnnextionlevel*nne * (stimrecombestimator_E_save[pkt_ptr->where*nelements*maxion+element*maxion+ion]-stimrecombestimator_save[pkt_ptr->where*nelements*maxion+element*maxion+ion]) * E_threshold);
      //printout("element %d, ion %d, modified %g, usual %g, diff %g, nonE %g\n",element,ion,stimrecombestimator_E_save[pkt_ptr->where*nelements*maxion+element*maxion+ion],stimrecombestimator_save[pkt_ptr->where*nelements*maxion+element*maxion+ion],(stimrecombestimator_E_save[pkt_ptr->where*nelements*maxion+element*maxion+ion]-stimrecombestimator_save[pkt_ptr->where*nelements*maxion+element*maxion+ion])*E_threshold,interpolate_bfcoolingcoeff(element,ion,level,T_e));
      //printout("alpha_sp %g , alpha_st %g\n",get_spontrecombcoeff(element,ion,level,T_e),interpolate_stimulated_recomb(element,ion,level,T_R));
      //epsilon_upper = epsilon(element,ion+1,0);
      //E_threshold = epsilon_upper - epsilon_current;
      //C = interpolate_spontrecombcoeff(element,ion,level,T_e) * E_threshold * nnnextionlevel*nne;
      //C_ion += C;
      //C_fb += C;
      //cellhistory[tid].coolinglist[i].contribution = C;
      //if (contrib < oldcoolingsum) printout("contrib %g < oldcoolingsum %g: C%g, element %d, ion %d, level %d, coolingtype %d, i %d, low %d\n",contrib,oldcoolingsum,C,element,ion,level,cellhistory[tid].coolinglist[i].type,i,low);
    }
  }

  assert(indexionstart == get_coolinglistoffset(element, ion));
  assert(i == indexionstart + get_ncoolingterms(element, ion));
  assert(fabs((modelgrid[modelgridindex].cooling[element].contrib[ion] + oldcoolingsum - contrib) / contrib) < 1e-3);

  //importantcoolingterms = ncoolingterms;
  //cellhistory[tid].totalcooling = contrib;
}


static double planck(const double nu, const double T)
/// returns intensity for frequency nu and temperature T according
/// to the Planck distribution
{
  return TWOHOVERCLIGHTSQUARED * pow(nu, 3) / expm1(HOVERKB * nu / T);
}


static double sample_planck(const double T)
/// returns a randomly chosen frequency according to the Planck
/// distribution of temperature T
{
  const double nu_peak = 5.879e10 * T;
  if (nu_peak > nu_max_r || nu_peak < nu_min_r)
    printout("[warning] sample_planck: intensity peaks outside frequency range\n");

  const double B_peak = planck(nu_peak, T);

  double nu;
  bool endloop = false;
  // int i = 0;
  while (!endloop)
  {
    // i++;
    double zrand = gsl_rng_uniform(rng);
    double zrand2 = gsl_rng_uniform(rng);
    nu = nu_min_r + zrand * (nu_max_r - nu_min_r);
    if (zrand2 * B_peak <= planck(nu, T))
      endloop = true;
    //printout("[debug] sample_planck: planck_sampling %d\n", i);
  }

  return nu;
}


double do_kpkt_bb(PKT *restrict pkt_ptr, const double t1)
/// Now routine to deal with a k-packet. Similar idea to do_gamma.
{
  //double nne = cell[pkt_ptr->where].nne ;
  int cellindex = pkt_ptr->where;
  const int modelgridindex = cell[cellindex].modelgridindex;
  const float T_e = get_Te(modelgridindex);
  double t_current = t1;

  pkt_ptr->nu_cmf = sample_planck(T_e);
  if (!isfinite(pkt_ptr->nu_cmf))
  {
    printout("[fatal] do_kpkt_bb: selected frequency not finite ... abort\n");
    abort();
  }
  /// and then emitt the packet randomly in the comoving frame
  emitt_rpkt(pkt_ptr, t_current);
  if (debuglevel == 2)
    printout("[debug] calculate_kappa_rpkt after kpkt to rpkt by ff\n");
  cellindex = pkt_ptr->where;
  pkt_ptr->next_trans = 0;      ///FLAG: transition history here not important, cont. process
  //if (tid == 0) k_stat_to_r_bb++;
  k_stat_to_r_bb++;
  pkt_ptr->interactions++;
  pkt_ptr->last_event = 6;
  pkt_ptr->emissiontype = -9999999;
  vec_copy(pkt_ptr->em_pos, pkt_ptr->pos);
  pkt_ptr->em_time = t_current;
  pkt_ptr->nscatterings = 0;

  return t_current;
}


double do_kpkt(PKT *restrict pkt_ptr, double t1, double t2, int nts)
/// Now routine to deal with a k-packet. Similar idea to do_gamma.
//{
//  double do_kpkt_bb(PKT *pkt_ptr, double t1, double t2);
//  return do_kpkt_bb(pkt_ptr, t1, t2);
//}
{
  const int cellindex = pkt_ptr->where;
  const int modelgridindex = cell[cellindex].modelgridindex;

  /// don't calculate cooling rates after each cell crossings anylonger
  /// but only if we really get a kpkt and they hadn't been calculated already
  //if (cellhistory[tid].totalcooling == COOLING_UNDEFINED)
/*  int ondemand = 1;
  if (modelgrid[modelgridindex].totalcooling == COOLING_UNDEFINED)
  {
    //printout("initial calculation of all cooling rates\n");
    coolingratecalccounter++;
    ondemand = 0;
    calculate_kpkt_rates(modelgridindex);
  }*/
  //printout("totalcooling %g\n",modelgrid[modelgridindex].totalcooling );

  //printout("[debug] do_kpkt: propagate k-pkt\n");
  //double cut = 0.99; //1.;


  const float T_e = get_Te(modelgridindex);
  double deltat = 0.;
  if (nts < n_kpktdiffusion_timesteps)
    deltat = kpktdiffusion_timescale * time_step[nts].width;
  //double deltat = 1./(nne*1.02e-12*pow(T_e/1e4,0.843));
  //printout("kpkt diffusion time simple %g, advanced %g\n",deltat,1/(nne*1.02e-12*pow(T_e/1e4,0.843)));
  double t_current = t1 + deltat;

  if (t_current <= t2)
  {
    vec_scale(pkt_ptr->pos, t_current / t1);

    /// Randomly select the occuring cooling process out of the important ones
    double coolingsum = 0.;
    double zrand = gsl_rng_uniform(rng);
    //if (debuglevel == 2) printout("do_kpkt: totalcooling %g, zrand %g, cut %g\n",cellhistory[tid].totalcooling,zrand,COOLINGCUT);
    //printout("do_kpkt: totalcooling %g, zrand %g, cut %g\n",cellhistory[tid].totalcooling,zrand,COOLINGCUT);
    /*for (i = 0; i < importantcoolingterms; i++)
    {
      coolingsum += cellhistory[tid].coolinglist[i].contribution;
      if (debuglevel == 2) printout("do_kpkt: loop i %d, coolingsum %g\n",i,coolingsum);
      //printout("do_kpkt: loop i %d, coolingsum %g, contr %g\n",i,coolingsum,cellhistory[tid].coolinglist[i].contribution);
      if (zrand*cellhistory[tid].totalcooling < coolingsum) break;
    }*/


    const double rndcool = zrand * modelgrid[modelgridindex].totalcooling;
    //printout("rndcool %g\n",rndcool);
    double oldcoolingsum;
    int element;
    int ion;
    for (element = 0; element < nelements; element++)
    {
      const int nions = get_nions(element);
      for (ion = 0; ion < nions; ion++)
      {
        oldcoolingsum = coolingsum;
        coolingsum += modelgrid[modelgridindex].cooling[element].contrib[ion];
        // printout("Z=%d, ionstage %d, coolingsum %g\n", get_element(element), get_ionstage(element, ion), coolingsum);
        if (coolingsum > rndcool) break;
      }
      if (coolingsum > rndcool) break;
    }
    // printout("kpkt selected Z=%d ionstage %d\n", get_element(element), get_ionstage(element, ion));

  //  #ifdef DEBUG_ON
      if (element >= nelements || ion >= get_nions(element))
      {
        printout("do_kpkt: problem selecting a cooling process ... abort\n");
        printout("do_kpkt: tried to select element %d, ion %d (Z=%d ionstage %d)\n",element, ion, get_element(element), get_ionstage(element, ion));
        printout("do_kpkt: totalcooling %g, coolingsum %g, rndcool %g\n",modelgrid[modelgridindex].totalcooling,coolingsum,rndcool);
        printout("do_kpkt: modelgridindex %d, cellno %d, nne %g\n",modelgridindex,pkt_ptr->where,get_nne(modelgridindex));
        for (element = 0; element < nelements; element++)
        {
          const int nions = get_nions(element);
          for (ion = 0; ion < nions; ion++)
          {
            printout("do_kpkt: element %d, ion %d, coolingcontr %g\n",element,ion,modelgrid[modelgridindex].cooling[element].contrib[ion]);
          }
        }
        abort();
      }
  //  #endif

    // debuglevel = 2;
    //printout("element %d, ion %d, coolingsum %g\n",element,ion,coolingsum);
    int ilow = get_coolinglistoffset(element,ion);
    int low = ilow;
    int high = low + get_ncoolingterms(element,ion)-1;
    //printout("element %d, ion %d, low %d, high %d\n",element,ion,low,high);
    if (cellhistory[tid].coolinglist[low].contribution == COOLING_UNDEFINED)
    {
      //printout("calculate kpkt rates on demand\n");
      calculate_kpkt_rates_ion(modelgridindex,element,ion,low,oldcoolingsum);
    }
    //if (cellhistory[tid].coolinglist[high].contribution != coolingsum);
    //{
    //  printout("ondemand %d\n",ondemand);
    //  printout("element %d, ion %d, oldcoolingsum %g, coolingsum %g\n",element,ion,oldcoolingsum,coolingsum);
    //  printout("contrib[ilow] %g, contrib[high] %g\n",cellhistory[tid].coolinglist[ilow].contribution,cellhistory[tid].coolinglist[high].contribution);
    //}
    //low = 0;
    //high = importantcoolingterms - 1;
    //rndcool = zrand*cellhistory[tid].totalcooling;
    int i = -1;
    while (low <= high)
    {
      i = (low + high) / 2;
      if (cellhistory[tid].coolinglist[i].contribution > rndcool)
      {
        if (i == ilow || cellhistory[tid].coolinglist[i-1].contribution < rndcool)
          break; /// found (1)
        else
          high = i - 1;
      }
      else
        low = i + 1;
      //else if (cellhistory[tid].coolinglist[i].contribution < rndcool)
      //  low = i + 1;
      //else
      //  break; /// found (2)
    }
    if (low > high)
    {
      printout("do_kpkt: error occured while selecting a cooling channel: low %d, high %d, i %d, rndcool %g\n",low,high,i,rndcool);
      printout("element %d, ion %d, offset %d, terms %d, coolingsum %g\n",element,ion,get_coolinglistoffset(element,ion),get_ncoolingterms(element,ion),coolingsum);
      printout("oldcoolingsum %g, coolingsum %g\n",oldcoolingsum,coolingsum);
      //for (i=ilow-1; i<=ilow+get_ncoolingterms(element,ion); i++)
      //  printout("  i %d, contrib %g\n",i,cellhistory[tid].coolinglist[i].contribution);

      //printout("ondemand %d\n",ondemand);

      printout("lower %g, %g, %g\n",cellhistory[tid].coolinglist[get_coolinglistoffset(element,ion)-1].contribution,cellhistory[tid].coolinglist[get_coolinglistoffset(element,ion)].contribution,cellhistory[tid].coolinglist[get_coolinglistoffset(element,ion)+1].contribution);
      int finalpos = get_coolinglistoffset(element,ion)+get_ncoolingterms(element,ion)-1;
      printout("upper %g, %g, %g\n",cellhistory[tid].coolinglist[finalpos-1].contribution,cellhistory[tid].coolinglist[finalpos].contribution,cellhistory[tid].coolinglist[finalpos+1]);
    }


    if (debuglevel == 2) printout("do_kpkt: selected process %d, coolingsum %g, importantcoolingterms %d\n",i,coolingsum,importantcoolingterms);
    //printout("do_kpkt: selected process %d, coolingsum %g, importantcoolingterms %d, type %d\n",i,coolingsum,importantcoolingterms,cellhistory[tid].coolinglist[i].type);

    // printout("element Z=%d, ion_stage %d, leve %d upper %d offset %d, terms %d, coolingsum %g\n",
    //   get_element(cellhistory[tid].coolinglist[i].element),
    //   get_ionstage(cellhistory[tid].coolinglist[i].element, cellhistory[tid].coolinglist[i].ion),
    //   cellhistory[tid].coolinglist[i].level,
    //   cellhistory[tid].coolinglist[i].upperlevel,
    //   get_coolinglistoffset(element,ion), get_ncoolingterms(element,ion), coolingsum);

    if (cellhistory[tid].coolinglist[i].type == COOLINGTYPE_FF)
    {
      /// The k-packet converts directly into a r-packet by free-free-emission.
      /// Need to select the r-packets frequency and a random direction in the
      /// co-moving frame.
      // printout("[debug] do_kpkt: k-pkt -> free-free\n");
      //kdecay.to_r++;

      /// Sample the packets comoving frame frequency according to paperII 5.4.3 eq.41
      //zrand = gsl_rng_uniform(rng);   /// delivers zrand in [0,1[
      //zrand = 1. - zrand;             /// make sure that log gets a zrand in ]0,1]
      zrand = gsl_rng_uniform_pos(rng);   /// delivers zrand in ]0,1[
      pkt_ptr->nu_cmf = -KB * T_e / H * log(zrand);
      //pkt_ptr->nu_cmf = 3.7474058e+14;

      if (!isfinite(pkt_ptr->nu_cmf))
      {
        printout("[fatal] ff cooling: selected frequency not finite ... abort\n");
        abort();
      }
      /// and then emitt the packet randomly in the comoving frame
      emitt_rpkt(pkt_ptr,t_current);
      pkt_ptr->next_trans = 0;      ///FLAG: transition history here not important, cont. process
      //if (tid == 0) k_stat_to_r_ff++;
      k_stat_to_r_ff++;
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 6;
      pkt_ptr->emissiontype = -9999999;
      vec_copy(pkt_ptr->em_pos, pkt_ptr->pos);
      pkt_ptr->em_time = t_current;
      pkt_ptr->nscatterings = 0;
      #ifndef FORCE_LTE
        //kffcount[pkt_ptr->where] += pkt_ptr->e_cmf;
      #endif
    }
    else if (cellhistory[tid].coolinglist[i].type == COOLINGTYPE_FB)
    {
      /// The k-packet converts directly into a r-packet by free-bound-emission.
      /// Need to select the r-packets frequency and a random direction in the
      /// co-moving frame.
      const int element = cellhistory[tid].coolinglist[i].element;
      const int lowerion = cellhistory[tid].coolinglist[i].ion;
      const int level = cellhistory[tid].coolinglist[i].level;
      const int upper = cellhistory[tid].coolinglist[i].upperlevel;
      // const double nu_threshold = get_phixs_threshold(element, ion, level, phixstargetindex)

      #ifdef DEBUG_ON
        // printout("[debug] do_kpkt: k-pkt -> free-bound\n");
        // printout("[debug] do_kpkt: element  %d, ion %d, level %d, upper %d, nu_threshold %g\n",element,ion,level,upper,nu_threshold);
      #endif

      /// then randomly sample the packets frequency according to the continuums
      /// energy distribution and set some flags
      //zrand = gsl_rng_uniform(rng);   /// delivers zrand in [0,1[
      //zrand = 1. - zrand;   /// convert it to ]0,1]
      //pkt_ptr->nu_cmf = nu_threshold * (1 - KB*T_e/H/nu_threshold*log(zrand));
      //pkt_ptr->nu_cmf = nu_threshold * (1+sqrt(1+(4*KB*T_e/H/nu_threshold)))/2 * (1 - KB*T_e/H/nu_threshold*log(zrand));
      //pkt_ptr->nu_cmf = nu_threshold;

      // Sample the packets comoving frame frequency according to paperII 4.2.2
      //zrand = gsl_rng_uniform(rng);
      //if (zrand < 0.5)
      {
        pkt_ptr->nu_cmf = select_continuum_nu(element, lowerion, level, upper, T_e);
      }
      // else
      // {
      //   ///Emitt like a BB
      //   pkt_ptr->nu_cmf = sample_planck(T_e);
      // }

      // printout("[debug] do_kpkt: pkt_ptr->nu_cmf %g\n",pkt_ptr->nu_cmf);

      // and then emitt the packet randomly in the comoving frame
      emitt_rpkt(pkt_ptr, t_current);

      #if (TRACK_ION_STATS)
      increment_ion_stats(modelgridindex, element, lowerion + 1, ION_COUNTER_RADRECOMB_KPKT, pkt_ptr->e_cmf / H / pkt_ptr->nu_cmf);
      const double escape_prob = get_rpkt_escape_prob(pkt_ptr, t_current);
      increment_ion_stats(modelgridindex, element, lowerion + 1, ION_COUNTER_RADRECOMB_ESCAPED, pkt_ptr->e_cmf / H / pkt_ptr->nu_cmf * escape_prob);
      #endif

      pkt_ptr->next_trans = 0;      ///FLAG: transition history here not important, cont. process
      //if (tid == 0) k_stat_to_r_fb++;
      k_stat_to_r_fb++;
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 7;
      pkt_ptr->emissiontype = get_continuumindex(element, lowerion, level, upper);
      pkt_ptr->trueemissiontype = pkt_ptr->emissiontype;
      vec_copy(pkt_ptr->em_pos, pkt_ptr->pos);
      pkt_ptr->em_time = t_current;
      pkt_ptr->nscatterings = 0;
    }
    else if (cellhistory[tid].coolinglist[i].type == COOLINGTYPE_COLLEXC)
    {
      /// the k-packet activates a macro-atom due to collisional excitation
      if (debuglevel == 2) printout("[debug] do_kpkt: k-pkt -> collisional excitation of MA\n");
      const int element = cellhistory[tid].coolinglist[i].element;
      const int ion = cellhistory[tid].coolinglist[i].ion;
      const int upper = cellhistory[tid].coolinglist[i].upperlevel;
      mastate[tid].element = element;
      mastate[tid].ion = ion;
      mastate[tid].level = upper;
      mastate[tid].activatingline = -99;

      #if (TRACK_ION_STATS)
      increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYIN_COLLEXC, pkt_ptr->e_cmf);
      #endif

      pkt_ptr->type = TYPE_MA;
      //if (tid == 0) ma_stat_activation_collexc++;
      ma_stat_activation_collexc++;
      //if (tid == 0) k_stat_to_ma_collexc++;
      k_stat_to_ma_collexc++;
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 8;
      pkt_ptr->trueemissiontype = -1; // since this is below zero, macroatom will set it
      pkt_ptr->trueemissionvelocity = -1;
      #ifndef FORCE_LTE
        //maabs[pkt_ptr->where] += pkt_ptr->e_cmf;
        //kffcount[pkt_ptr->where] += pkt_ptr->e_cmf;
      #endif

      return do_macroatom(pkt_ptr, t_current, t2, nts);
    }
    else if (cellhistory[tid].coolinglist[i].type == COOLINGTYPE_COLLION)
    {
      /// the k-packet activates a macro-atom due to collisional ionisation
      if (debuglevel == 2) printout("[debug] do_kpkt: k-pkt -> collisional ionisation of MA\n");
      const int element = cellhistory[tid].coolinglist[i].element;
      const int ion = cellhistory[tid].coolinglist[i].ion + 1;
      const int upper = cellhistory[tid].coolinglist[i].upperlevel;
      mastate[tid].element = element;
      mastate[tid].ion = ion;
      mastate[tid].level = upper;
      mastate[tid].activatingline = -99;

      #if (TRACK_ION_STATS)
      increment_ion_stats(modelgridindex, element, ion, ION_COUNTER_MACROATOM_ENERGYIN_COLLION, pkt_ptr->e_cmf);
      #endif

      pkt_ptr->type = TYPE_MA;
      //if (tid == 0) ma_stat_activation_collion++;
      ma_stat_activation_collion++;
      //if (tid == 0) k_stat_to_ma_collion++;
      k_stat_to_ma_collion++;
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 9;
      pkt_ptr->trueemissiontype = -1; // since this is below zero, macroatom will set it
      pkt_ptr->trueemissionvelocity = -1;
      #ifndef FORCE_LTE
        //maabs[pkt_ptr->where] += pkt_ptr->e_cmf;
        //kffcount[pkt_ptr->where] += pkt_ptr->e_cmf;
      #endif

      return do_macroatom(pkt_ptr, t_current, t2, nts);
    }
    else
    {
      printout("[fatal] do_kpkt: coolinglist.type mismatch\n");
      printout("[fatal] do_kpkt: zrand %g, modelgrid[modelgridindex].totalcooling %g, coolingsum %g, i %d\n",zrand,modelgrid[modelgridindex].totalcooling,  coolingsum,i);
      printout("[fatal] do_kpkt: COOLINGCUT %d, importantcoolingterms %d, tid %d\n",COOLINGCUT,importantcoolingterms,tid);
      printout("[fatal] do_kpkt: cellhistory[tid].coolinglist[i].type %d\n",cellhistory[tid].coolinglist[i].type);
      printout("[fatal] do_kpkt: pkt_ptr->where %d, mgi %d\n",pkt_ptr->where,modelgridindex);
      abort();
    }

    return t_current;
  }
  else
  {
    vec_scale(pkt_ptr->pos, t2 / t1);
    return PACKET_SAME;
  }
}


/*static int compare_coolinglistentry(const void *p1, const void *p2)
/// Helper function to sort the coolinglist by the strength of the
/// individual cooling contributions.
{
  ionscoolinglist_t *a1, *a2;
  a1 = (ionscoolinglist_t *)(p1);
  a2 = (ionscoolinglist_t *)(p2);
  //printf("%d %d %d %d %g\n",a1->elementindex,a1->ionindex,a1->lowerlevelindex,a1->upperlevelindex,a1->nu);
  //printf("%d %d %d %d %g\n",a2->elementindex,a2->ionindex,a2->lowerlevelindex,a2->upperlevelindex,a2->nu);
  //printf("%g\n",a2->nu - a1->nu);
  if (a1->contribution - a2->contribution < 0)
    return 1;
  else if (a1->contribution - a2->contribution > 0)
    return -1;
  else
    return 0;
}*/


/*double get_bfcooling_direct(int element, int ion, int level, int cellnumber)
/// Returns the rate for bfheating. This can be called during packet propagation
/// or update_grid. Therefore we need to decide whether a cell history is
/// known or not.
{
  double bfcooling;

  gsl_integration_workspace *wsp;
  gslintegration_paras intparas;
  double bfcooling_integrand_gsl(double nu, void *paras);
  gsl_function F_bfcooling;
  F_bfcooling.function = &bfcooling_integrand_gsl;
  double intaccuracy = 1e-2;        /// Fractional accuracy of the integrator
  double error;
  double nu_max_phixs;

  float T_e = cell[cellnumber].T_e;
  float nne = cell[cellnumber].nne;
  double nnionlevel = get_groundlevelpop(cellnumber,element,ion+1);
  //upper = cellhistory[tid].coolinglist[i].upperlevel;
  double nu_threshold = (epsilon(element,ion+1,0) - epsilon(element,ion,level)) / H;
  nu_max_phixs = nu_threshold * last_phixs_nuovernuedge; //nu of the uppermost point in the phixs table

  mastate[tid].element = element;
  mastate[tid].ion = ion;
  mastate[tid].level = level;
  intparas.T = T_e;
  intparas.nu_edge = nu_threshold;   /// Global variable which passes the threshold to the integrator
  F_bfcooling.params = &intparas;
  gsl_integration_qag(&F_bfcooling, nu_threshold, nu_max_phixs, 0, intaccuracy, 1024, 6, wsp, &bfcooling, &error);
  bfcooling *= nnionlevel*nne*4*PI*calculate_sahafact(element,ion,level,upperionlevel,T_e,nu_threshold*H);

  return bfcooling;
}*/


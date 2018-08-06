#include <gsl/gsl_roots.h>
#include "sn3d.h"
#include "atomic.h"
#include "assert.h"
#include "grid_init.h"
#include "kpkt.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "nltepop.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "thermalbalance.h"
#include "update_grid.h"


extern inline double get_abundance(int modelgridindex, int element);


void precalculate_partfuncts(int modelgridindex)
/// The partition functions depend only on T_R and W. This means they don't
/// change during any iteration on T_e. Therefore their precalculation was
/// taken out of calculate_populations to save runtime.
{
  /// Precalculate partition functions for each ion in every cell
  /// this saves a factor 10 in calculation time of Saha-Boltzman populations
  for (int element = 0; element < nelements; element++)
  {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++)
    {
      //printout("precalc element %d, ion %d, mgi %d\n",element,ion,modelgridindex);
      //cell[cellnumber].composition[element].ltepartfunct[ion] = calculate_ltepartfunct(element,ion,T_R);
      modelgrid[modelgridindex].composition[element].partfunct[ion] = calculate_partfunct(element,ion,modelgridindex);
    }

  }
}

static void calculate_double_decay_chain(
  double initabund1, double meanlife1,
  double initabund2, double meanlife2,
  double t_current,
  double *abund1, double *abund2, double *abund3)
{
  // calculate abundances from double decay, e.g., Ni56 -> Co56 -> Fe56
  // initabund3 is assumed to be zero, so the abundance of species 3 is only from decays of species 2
  const double lambda1 = 1 / meanlife1;
  const double lambda2 = 1 / meanlife2;

  *abund1 = initabund1 * exp(-lambda1 * t_current);

  *abund2 = (
    initabund2 * exp(-lambda2 * t_current) +
    initabund1 * lambda1 / (lambda1 - lambda2) * (exp(-lambda2 * t_current) - exp(-lambda1 * t_current)));

  *abund3 = (
    (initabund2 + initabund1) * (lambda1 - lambda2) -
    initabund2 * lambda1 * exp(-lambda2 * t_current) +
    initabund2 * lambda2 * exp(-lambda2 * t_current) -
    initabund1 * lambda1 * exp(-lambda2 * t_current) +
    initabund1 * lambda2 * exp(-lambda1 * t_current)) / (lambda1 - lambda2);

  // printout("calculate_double_decay_chain: initabund1 %g, initabund2 %g\n", initabund1, initabund2);
  // printout("calculate_double_decay_chain: abund1 %g, abund2 %g abund3 %g\n", abund1, abund2, abund3);

  // ensure that the decays haven't altered the total abundance of all three species
  assert(fabs((initabund1 + initabund2) - (*abund1 + *abund2 + *abund3)) < 0.001);
}


static void update_abundances(const int modelgridindex, const int timestep, const double t_current)
/// Updates the mass fractions of elements associated with the decay sequence
/// (56)Ni -> (56)Co -> (56)Fe at the onset of each timestep
/// Parameters: - modelgridindex: the grid cell for which to update the abundances
///             - t_current: current time (here mid of current timestep)
{
  const double timediff = t_current - t_model;
  if (homogeneous_abundances)
  {
    const double ni56_in = elements[get_elementindex(28)].abundance; // assume all ni56
    const double co56_in = elements[get_elementindex(27)].abundance; // assume all co56
    const double fe_in = elements[get_elementindex(26)].abundance;

    double ni56frac = 0.;
    double co56frac = 0.;
    double fe56frac_fromdecay = 0.;
    calculate_double_decay_chain(ni56_in, T56NI, co56_in, T56CO, timediff, &ni56frac, &co56frac, &fe56frac_fromdecay);

    //fe_in = cell[modelgridindex].f_fe_init;
    for (int element = nelements - 1; element >= 0; element--)
    {
      const int atomic_number = get_element(element);
      if (atomic_number == 28)
      {
        const double nifrac = get_fnistable(modelgridindex) + ni56frac;
        modelgrid[modelgridindex].composition[element].abundance = nifrac;
      }
      else if (atomic_number == 27)
      {
        const double cofrac = get_fcostable(modelgridindex) + co56frac;
        modelgrid[modelgridindex].composition[element].abundance = cofrac;
      }
      else if (atomic_number == 26)
      {
        const double fefrac = fe_in + fe56frac_fromdecay;
        modelgrid[modelgridindex].composition[element].abundance = fefrac;
      }
    }
  }
  else
  {
    // Ni56 -> Co56 -> Fe56
    // abundances from the input model
    const double ni56_in = get_modelradioabund(modelgridindex, NUCLIDE_NI56);
    const double co56_in = get_modelradioabund(modelgridindex, NUCLIDE_CO56);
    double ni56frac = 0.;
    double co56frac = 0.;
    double fe56frac_fromdecay = 0.;
    calculate_double_decay_chain(ni56_in, T56NI, co56_in, T56CO, timediff, &ni56frac, &co56frac, &fe56frac_fromdecay);

    // Ni57 -> Co57 -> Fe57
    const double ni57_in = get_modelradioabund(modelgridindex, NUCLIDE_NI57);
    const double co57_in = get_modelradioabund(modelgridindex, NUCLIDE_CO57);
    double ni57frac = 0.;
    double co57frac = 0.;
    double fe57frac_fromdecay = 0.;
    calculate_double_decay_chain(ni57_in, T57NI, co57_in, T57CO, timediff, &ni57frac, &co57frac, &fe57frac_fromdecay);

    // Fe52 -> Mn52 -> Cr52
    const double fe52_in = get_modelradioabund(modelgridindex, NUCLIDE_FE52);
    double fe52frac = 0.;
    double mn52frac = 0.;
    double cr52frac_fromdecay = 0.;
    calculate_double_decay_chain(fe52_in, T52FE, 0., T52MN, timediff, &fe52frac, &mn52frac, &cr52frac_fromdecay);

    // Cr48 -> V48 -> Ti48
    const double cr48_in = get_modelradioabund(modelgridindex, NUCLIDE_CR48);
    double cr48frac = 0.;
    double v48frac = 0.;
    double ti48frac_fromdecay = 0.;
    calculate_double_decay_chain(cr48_in, T48CR, 0., T48V, timediff, &cr48frac, &v48frac, &ti48frac_fromdecay);

    // printout("model cell %d, has input radioactive ni56_in %g, co56_in %g, fe52_in %g\n",modelgridindex,ni56_in,co56_in,fe52_in);

    for (int element = nelements-1; element >= 0; element--)
    {
      const int atomic_number = get_element(element);
      if (atomic_number == 28)
      {
        const double nifrac = get_fnistable(modelgridindex) + ni56frac + ni57frac;
        modelgrid[modelgridindex].composition[element].abundance = nifrac;
      }
      else if (atomic_number == 27)
      {
        const double cofrac = get_fcostable(modelgridindex) + co56frac + co57frac;
        modelgrid[modelgridindex].composition[element].abundance = cofrac;
      }
      else if (atomic_number == 26)
      {
        const double fefrac = get_ffestable(modelgridindex) + fe52frac + fe56frac_fromdecay + fe57frac_fromdecay;
        modelgrid[modelgridindex].composition[element].abundance = fefrac;
      }
      else if (atomic_number == 25)
      {
        const double mnfrac = get_fmnstable(modelgridindex) + mn52frac;
        modelgrid[modelgridindex].composition[element].abundance = mnfrac;
      }
      else if (atomic_number == 24)
      {
        const double crfrac = get_fcrstable(modelgridindex) + cr48frac + cr52frac_fromdecay;
        modelgrid[modelgridindex].composition[element].abundance = crfrac;
      }
      else if (atomic_number == 23)
      {
        const double vfrac = get_fvstable(modelgridindex) + v48frac;
        modelgrid[modelgridindex].composition[element].abundance = vfrac;
      }
      else if (atomic_number == 22)
      {
        const double tifrac = get_ftistable(modelgridindex) + ti48frac_fromdecay;
        modelgrid[modelgridindex].composition[element].abundance = tifrac;
      }
    }
    // printout("model cell %d at t_current %g has frac: Ni %g Co %g Fe %g, stable: Ni %g Co %g Fe %g\n",
    //          modelgridindex, t_current,
    //          modelgrid[modelgridindex].composition[get_elementindex(28)].abundance,
    //          modelgrid[modelgridindex].composition[get_elementindex(27)].abundance,
    //          modelgrid[modelgridindex].composition[get_elementindex(26)].abundance,
    //          get_fnistabel(modelgridindex), get_fcostable(modelgridindex), get_ffestable(modelgridindex));
  }

  calculate_deposition_rate_density(modelgridindex, timestep);
}


static void write_to_estimators_file(const int n, const int timestep)
{
  if (mg_associated_cells[n] > 0)
  {
    //fprintf(estimators_file,"%d %g %g %g %g %d ",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),modelgrid[n].thick);
    //fprintf(estimators_file,"%d %g %g %g %g %g ",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),grey_optical_depth);
    fprintf(estimators_file, "timestep %d modelgridindex %d TR %g Te %g W %g TJ %g grey_depth %g nne %g\n",
            timestep, n, get_TR(n), get_Te(n), get_W(n), get_TJ(n), modelgrid[n].grey_depth, get_nne(n));
    //fprintf(estimators_file,"%d %g %g %g %g %g %g %g
    //",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),grey_optical_depth,grey_optical_deptha,compton_optical_depth);

    if (NLTE_POPS_ON)  //  && timestep % 2 == 0
      nltepop_write_to_file(n, timestep);

    const float T_e = get_Te(n);
    const float nne = get_nne(n);

    for (int element = 0; element < nelements; element++)
    {
      if (get_abundance(n, element) <= 0.) // skip elements with no abundance
        continue;

      fprintf(estimators_file, "populations    Z=%2d", get_element(element));
      const int nions = get_nions(element);
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions; ion++)
      {
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ionstagepop(n, element, ion));
      }
      fprintf(estimators_file, "\n");

      const bool printdebug = false;

      bool assume_lte = true;
      bool per_gmpop = true;
      const bool lower_superlevel_only = false;

      fprintf(estimators_file, "RRC_LTE_Nahar  Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions; ion++)
      {
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                calculate_ionrecombcoeff(-1, T_e, element, ion, assume_lte, false, printdebug, lower_superlevel_only, per_gmpop));
      }
      fprintf(estimators_file, "\n");

      per_gmpop = false;

      // fprintf(estimators_file, "AlphaLTE_R*nne Z=%2d", get_element(element));
      // for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //   fprintf(estimators_file, "              ");
      // for (int ion = 0; ion < nions; ion++)
      // {
      //   fprintf(estimators_file, "  %d: %9.3e",
      //           get_ionstage(element, ion),
      //           calculate_ionrecombcoeff(n, T_e, element, ion, assume_lte, false, printdebug, lower_superlevel_only, per_gmpop) * nne);
      // }
      // fprintf(estimators_file, "\n");

      assume_lte = false;

      fprintf(estimators_file, "Alpha_R*nne    Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions; ion++)
      {
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                calculate_ionrecombcoeff(n, T_e, element, ion, assume_lte, false, printdebug, lower_superlevel_only, per_gmpop) * nne);
      }
      fprintf(estimators_file, "\n");

      // fprintf(estimators_file, "Alpha_C*nne    Z=%2d", get_element(element));
      // for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //   fprintf(estimators_file, "              ");
      // for (int ion = 0; ion < nions; ion++)
      // {
      //   fprintf(estimators_file, "  %d: %9.3e",
      //           get_ionstage(element, ion),
      //           calculate_ionrecombcoeff(n, T_e, element, ion, assume_lte, true, printdebug, lower_superlevel_only, per_gmpop) * nne);
      // }
      // fprintf(estimators_file, "\n");

      fprintf(estimators_file, "gamma_R        Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++)
      {
        // const bool printdebug_gammar = (get_element(element) == 26 && get_ionstage(element, ion) == 1);
        const bool printdebug_gammar = false;
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                calculate_iongamma_per_ionpop(n, T_e, element, ion, assume_lte, false, printdebug_gammar));
      }
      fprintf(estimators_file, "\n");

      // fprintf(estimators_file, "gamma_C        Z=%2d", get_element(element));
      // for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //   fprintf(estimators_file, "              ");
      // for (int ion = 0; ion < nions - 1; ion++)
      // {
      //   fprintf(estimators_file, "  %d: %9.3e",
      //           get_ionstage(element, ion),
      //           calculate_iongamma_per_ionpop(n, T_e, element, ion, assume_lte, true, printdebug));
      // }
      // fprintf(estimators_file, "\n");

      if (NT_ON)
      {
        fprintf(estimators_file, "gamma_NT       Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions - 1; ion++)
        {
          const double Y_nt = nt_ionization_ratecoeff(n, element, ion);
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), Y_nt);
        }
        fprintf(estimators_file, "\n");
      }
    }


    const int element = get_elementindex(28);
    const int lowerionstage = 2;
    const int lowerion = lowerionstage - get_ionstage(element, 0);
    const bool printdebug = false;
    const bool lower_superlevel_only = false;
    const bool per_gmpop = false;
    calculate_ionrecombcoeff(n, T_e, element, lowerion + 1, false, false, printdebug, lower_superlevel_only, per_gmpop);

    #ifndef FORCE_LTE
      #if (!NO_LUT_PHOTOION)
        fprintf(estimators_file, "corrphotoionrenorm: ");
        for (int element = 0; element < nelements; element++)
        {
          const int nions = get_nions(element);
          for (int ion = 0; ion < nions; ion++)
          {
            fprintf(estimators_file,"%g ", corrphotoionrenorm[n*nelements*maxion+element*maxion+ion]);
          }
        }
        fprintf(estimators_file, "\n");
        fprintf(estimators_file, "gammaestimator: ");
        for (int element = 0; element < nelements; element++)
        {
          const int nions = get_nions(element);
          for (int ion = 0; ion < nions; ion++)
          {
            fprintf(estimators_file,"%g ", gammaestimator[n*nelements*maxion+element*maxion+ion]);
          }
        }
        fprintf(estimators_file,"\n");
      #endif

      fprintf(estimators_file, "heating: ff %11.5e bf %11.5e coll %11.5e     gamma %11.5e gamma/gamma_dep %.2f\n",
              heatingrates[tid].ff, heatingrates[tid].bf, heatingrates[tid].collisional, heatingrates[tid].gamma, heatingrates[tid].nt_frac_heating);
      fprintf(estimators_file, "cooling: ff %11.5e fb %11.5e coll %11.5e adiabatic %11.5e\n",
              coolingrates[tid].ff, coolingrates[tid].fb, coolingrates[tid].collisional, coolingrates[tid].adiabatic);
    #endif
  }
  else
  {
    // modelgrid cells which are not represented in the simulation grid
    fprintf(estimators_file, "timestep %d modelgridindex %d EMPTYCELL\n", timestep, n);
  }
  fprintf(estimators_file,"\n");

  fflush(estimators_file);
}


void cellhistory_reset(const int modelgridindex, const bool new_timestep)
{
  /// All entries of the cellhistory stack must be flagged as empty at the
  /// onset of the new timestep. Also, boundary crossing?
  /// Calculate the level populations for this cell, and flag the other entries
  /// as empty.
  /// Make known that cellhistory[tid] contains information about the
  /// cell given by cellnumber. (-99 if invalid)
  if ((modelgridindex == cellhistory[tid].cellnumber) && !new_timestep)
  {
    // printout("redundant reset in cell %d\n", modelgridindex);
    return; //TODO: does this break anything?
  }

  cellhistory[tid].cellnumber = modelgridindex;
  //cellhistory[tid].totalcooling = COOLING_UNDEFINED;
  for (int element = 0; element < nelements; element++)
  {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++)
    {
      cellhistory[tid].coolinglist[get_coolinglistoffset(element,ion)].contribution = COOLING_UNDEFINED;
      const int nlevels = get_nlevels(element,ion);
      for (int level = 0; level < nlevels; level++)
      {
        const double population = (new_timestep) ? -99 : calculate_exclevelpop(modelgridindex,element,ion,level);
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].population = population;

        for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
        {
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].sahafact = -99.;
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].spontaneousrecombrate = -99.;
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].bfcooling = -99.;
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].bfheatingcoeff = -99.;
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].corrphotoioncoeff = -99.;
        }
        /// This is the only flag needed for all of the following MA stuff!
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].processrates[MA_ACTION_COLDEEXC] = -99.;
        /*
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_deexc = -99.;
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_recomb = -99.;
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_recomb = -99.;
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_same = -99.;
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_same = -99.;
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_lower = -99.;
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_higher = -99.;

        ndowntrans = get_ndowntrans(element, ion, level);
        nuptrans = get_nuptrans(element, ion, level);
        for (i = 0; i < ndowntrans; i++)
        {
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_rad_deexc[i] = -99.;
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_down_same[i] = -99.;
        }
        for (i = 0; i < nuptrans; i++)
        {
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_up_same[i] = -99.;
        }
        */
      }
    }
  }
  //cellhistory[tid].totalcooling = COOLING_UNDEFINED;
  //cellhistory[tid].phixsflag = PHIXS_UNDEFINED;
}


static void grid_cell_solve_Te_nltepops(const int n, const int nts, const int titer)
{
  const double covergence_tolerance = 0.03;
  for (int nlte_iter = 0; nlte_iter <= NLTEITER; nlte_iter++)
  {
    const time_t sys_time_start_spencerfano = time(NULL);
    if (NT_ON && NT_SOLVE_SPENCERFANO)
    {
      nt_solve_spencerfano(n, nts, nlte_iter);  // depends on the ionization balance, and weakly on nne
    }
    const int duration_solve_spencerfano = time(NULL) - sys_time_start_spencerfano;

    const time_t sys_time_start_partfuncs_or_gamma = time(NULL);
    if (!NLTE_POPS_ON)
    {
      /// These don't depend on T_e, therefore take them out of the T_e iteration
      precalculate_partfuncts(n);
    }
#if (!NO_LUT_PHOTOION)
    else if ((nlte_iter != 0))
    {
      // recalculate the Gammas using the current population estimates
      for (int element = 0; element < nelements; element++)
      {
        const int nions = get_nions(element);
        for (int ion = 0; ion < nions - 1; ion++)
        {
          gammaestimator[n * nelements * maxion + element * maxion + ion] = calculate_iongamma_per_gspop(n, element, ion);
        }
      }
    }
#endif
    const int duration_solve_partfuncs_or_gamma = time(NULL) - sys_time_start_partfuncs_or_gamma;

    /// Find T_e as solution for thermal balance
    const time_t sys_time_start_Te = time(NULL);
    const int nts_for_te = (titer == 0) ? nts - 1 : nts;

    call_T_e_finder(n, nts, time_step[nts_for_te].mid, MINTEMP, MAXTEMP);

    const int duration_solve_T_e = time(NULL) - sys_time_start_Te;

    if (!NLTE_POPS_ON || !NLTE_POPS_ALL_IONS_SIMULTANEOUS) // do this in LTE or NLTE single ion solver mode
    {
      /// Store population values to the grid
      const time_t sys_time_start_pops = time(NULL);
      calculate_populations(n);
      const int duration_solve_pops = time(NULL) - sys_time_start_pops;
      //calculate_cooling_rates(n);
      //calculate_heating_rates(n);
      printout("Grid solver cell %d timestep %d: time spent on: Spencer-Fano %ds, partfuncs/gamma %ds, T_e %ds, populations %ds\n",
               n, nts, duration_solve_spencerfano, duration_solve_partfuncs_or_gamma, duration_solve_T_e, duration_solve_pops);
    }

    if (NLTE_POPS_ON)
    {
      const time_t sys_time_start_nltepops = time(NULL);
      // fractional difference between previous and current iteration's (nne or max(ground state population change))
      double nlte_test;
      for (int element = 0; element < nelements; element++)
      {
        if (NLTE_POPS_ALL_IONS_SIMULTANEOUS)
          solve_nlte_pops_element(element, n, nts, nlte_iter);
        else
        {
          const int nions = get_nions(element);
          for (int ion = 0; ion < nions-1; ion++)
          {
            double trial = fabs(solve_nlte_pops_ion(element, ion, n, nts) - 1);

            if (trial > nlte_test)
              nlte_test = trial;
          }
        }
      }
      const int duration_solve_nltepops = time(NULL) - sys_time_start_nltepops;

      if (NLTE_POPS_ALL_IONS_SIMULTANEOUS)
      {
        const double nne_prev = get_nne(n);
        precalculate_partfuncts(n);
        calculate_electron_densities(n); // sets nne
        nlte_test = fabs((get_nne(n) / nne_prev) - 1);
        printout("NLTE solver cell %d timestep %d iteration %d: time spent on: Spencer-Fano %ds, T_e %ds, NLTE populations %ds\n",
                 n, nts, nlte_iter, duration_solve_spencerfano, duration_solve_T_e, duration_solve_nltepops);
        printout("NLTE (Spencer-Fano/Te/pops) solver cell %d timestep %d iteration %d: previous nne was %g, new nne is %g, fractional difference is %g\n",
                 n, nts, nlte_iter, nne_prev, get_nne(n), nlte_test);
        // damp changes in nne if oscillating to much
        //set_nne(n, (get_nne(n) + nne_prev) / 2.);
      }
      else
      {
        printout("Completed iteration for NLTE population solver in cell %d for timestep %d. Fractional error returned: %g\n",
                 n, nts, nlte_test);
      }

      if (nlte_test <= covergence_tolerance)
      {
        printout("NLTE solver converged to tolerance %g <= %g after %d iterations.\n", nlte_test, covergence_tolerance, nlte_iter + 1);
        break;
      }
      else if (nlte_iter == NLTEITER)
        printout("WARNING: NLTE solver failed to converge after %d iterations. Keeping solution from last iteration\n", nlte_iter + 1);
    }
    else
      break; // no iteration is needed without NLTE_POPS_ON
  }
}


#if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
  static void update_gamma_corrphotoionrenorm_bfheating_estimators(const int n, const double estimator_normfactor)
  {
    #if (!NO_LUT_PHOTOION)
      for (int element = 0; element < nelements; element++)
      {
        const int nions = get_nions(element);
        for (int ion = 0; ion < nions - 1; ion++)
        {
          const int ionestimindex = n * nelements * maxion + element * maxion + ion;
          //printout("mgi %d, element %d, ion %d, gammaest %g\n",n,element,ion,gammaestimator[ionestimindex]);
          gammaestimator[ionestimindex] *= estimator_normfactor / H;
          //printout("mgi %d, element %d, ion %d, gammaest %g\n",n,element,ion,gammaestimator[ionestimindex]);
          #ifdef DO_TITER
            if (gammaestimator_save[ionestimindex] >= 0)
            {
              gammaestimator[ionestimindex] = (gammaestimator[ionestimindex] + gammaestimator_save[ionestimindex]) / 2;
            }
            gammaestimator_save[ionestimindex] = gammaestimator[ionestimindex];
          #endif

          corrphotoionrenorm[ionestimindex] = gammaestimator[ionestimindex] / get_corrphotoioncoeff_ana(element,ion,0,0,n);

          if (!isfinite(corrphotoionrenorm[ionestimindex]))
          {
            printout("[fatal] about to set corrphotoionrenorm = NaN = gammaestimator / get_corrphotoioncoeff_ana(%d,%d,%d,%d,%d)=%g/%g",
                     element,ion,0,0,n,gammaestimator[ionestimindex],get_corrphotoioncoeff_ana(element,ion,0,0,n));
            abort();
          }
        }

      /// 2012-01-11. These loops should terminate here to precalculate *ALL* corrphotoionrenorm values
      /// so that the values are known when required by the call to get_corrphotoioncoeff in the following
      /// loops. Otherwise get_corrphotoioncoeff tries to renormalize by the closest corrphotoionrenorm
      /// in frequency space which can lead to zero contributions to the total photoionsation rate!
      }
    #endif
    #if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
      /// Then reopen the same loops again.
      for (int element = 0; element < nelements; element++)
      {
        const int nions = get_nions(element);
        for (int ion = 0; ion < nions-1; ion++)
        {
          /// Reuse the gammaestimator array as temporary storage of the Gamma values during
          /// the remaining part of the update_grid phase. Afterwards it is reset to record
          /// the next timesteps gamma estimators.
          //nlevels = get_nlevels(element,ion);
          //nlevels = get_ionisinglevels(element,ion);
          const int ionestimindex = n * nelements * maxion + element * maxion + ion;

          #if (!NO_LUT_PHOTOION)
            gammaestimator[ionestimindex] = calculate_iongamma_per_gspop(n, element, ion);
            //printout("mgi %d, element %d, ion %d, Gamma %g\n",n,element,ion,Gamma);
          #endif

          #if (!NO_LUT_BFHEATING)
            bfheatingestimator[ionestimindex] *= estimator_normfactor;
            #ifdef DO_TITER
              if (bfheatingestimator_save[ionestimindex] >= 0)
              {
                bfheatingestimator[ionestimindex] = (bfheatingestimator[ionestimindex]+bfheatingestimator_save[ionestimindex]) / 2;
              }
              bfheatingestimator_save[ionestimindex] = bfheatingestimator[ionestimindex];
            #endif
            /// Now convert bfheatingestimator into the bfheating renormalisation coefficient used in get_bfheating
            /// in the remaining part of update_grid. Later on it's reset and new contributions are added up.

            bfheatingestimator[ionestimindex] = bfheatingestimator[ionestimindex] / get_bfheatingcoeff_ana(element,ion,0,0,n);

            if (!isfinite(bfheatingestimator[ionestimindex]))
            {
              printout("[fatal] about to set bfheatingestimator = NaN = bfheatingestimator / get_bfheatingcoeff_ana(%d,%d,%d,%d,%d)=%g/%g",element,ion,0,0,n,bfheatingestimator[ionestimindex],get_bfheatingcoeff_ana(element,ion,0,0,n));
              abort();
            }

            //printout("cell %d element %d ion %d bfheatingestimator %g\n",n,element,ion,bfheatingestimator[ionestimindex]);
          #endif
        }
      }
    #endif
  }
#endif


#ifdef DO_TITER
#ifndef FORCE_LTE
  static void titer_average_estimators(const int n)
  {
    if (ffheatingestimator_save[n] >= 0)
    {
      ffheatingestimator[n] = (ffheatingestimator[n] + ffheatingestimator_save[n]) / 2;
    }
    ffheatingestimator_save[n] = ffheatingestimator[n];
    if (colheatingestimator_save[n] >= 0)
    {
      colheatingestimator[n] = (colheatingestimator[n] + colheatingestimator_save[n]) / 2;
    }
    colheatingestimator_save[n] = colheatingestimator[n];
  }
#endif
#endif


#if (!NO_LUT_PHOTOION)
  static void zero_gammaestimator(const int modelgridindex)
  {
    for (int element = 0; element < nelements; element++)
    {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++)
        gammaestimator[modelgridindex * nelements * maxion + element * maxion + ion] = 0.;
    }
  }

  static void set_all_corrphotoionrenorm(const int modelgridindex, const double value)
  {
    for (int element = 0; element < nelements; element++)
    {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++)
        corrphotoionrenorm[modelgridindex * nelements * maxion + element * maxion + ion] = value;
    }
  }
#endif


static void update_grid_cell(const int n, const int nts, const int titer, const double tratmid,
                             const double deltaV, const double deltat, double *mps)
// n is the modelgrid index
{
  const int assoc_cells = mg_associated_cells[n];
  if (assoc_cells > 0)
  {
    const time_t sys_time_start_update_cell = time(NULL);
    // const bool log_this_cell = ((n % 50 == 0) || (npts_model < 50));
    const bool log_this_cell = true;

    /// Update current mass density of cell
    //n = nonemptycells[my_rank+ncl*nprocs];
    if (log_this_cell)
      printout("[info] update_grid: working on cell %d before timestep %d ...\n", n, nts);
    //n = nonemptycells[ncl];
    //printout("[debug] update_grid: ncl %d is %d non-empty cell updating grid cell %d ... T_e %g, rho %g\n",ncl,my_rank+ncl*nprocs,n,cell[n].T_e,cell[n].rho);
    modelgrid[n].rho = get_rhoinit(n) / pow(tratmid, 3);
    //cell[n].rho = cell[n].rho_init / pow(tratmid,3);
    //rho = cell[n].rho;
    /// This is done outside update grid now
    //modelgrid[n].totalcooling = COOLING_UNDEFINED;

    if (opacity_case == 4)
    {
      /// Update abundances of radioactive isotopes
      //printout("call update abundances for timestep %d in model cell %d\n",m,n);
      const time_t sys_time_start_update_abundances = time(NULL);
      update_abundances(n, nts, time_step[nts].mid);
      printout("update_abundances for cell %d timestep %d took %d seconds\n", n, nts, time(NULL) - sys_time_start_update_abundances);

      /// For timestep 0 we calculate the level populations straight forward wihout
      /// applying any temperature correction
      if ((nts - itstep) == 0 && titer == 0)
      {
        /// Determine renormalisation factor for corrected photoionization cross-sections
        #ifndef FORCE_LTE
        #if (!NO_LUT_PHOTOION)
          if (!simulation_continued_from_saved)
          {
            set_all_corrphotoionrenorm(n, 1.);
          }
        #endif
        #endif

        /// W == 1 indicates that this modelgrid cell was treated grey in the
        /// last timestep. Therefore it has no valid Gamma estimators and must
        /// be treated in LTE at restart.
        if (modelgrid[n].thick == 0 && get_W(n) == 1)
        {
          if (log_this_cell)
            printout("force modelgrid cell to be grey at restart\n");
          modelgrid[n].thick = 1;
        }
        if (log_this_cell)
        {
          printout("initial_iteration %d\n",initial_iteration);
          printout("modelgrid.thick: %d\n",modelgrid[n].thick);
        }
        precalculate_partfuncts(n);
        //printout("abundance in cell %d is %g\n",n,cell[n].composition[0].abundance);
        if (!simulation_continued_from_saved || !NLTE_POPS_ON)
          calculate_populations(n);  // these were not read from the gridsave file, so calculate them now
        else
        {
          calculate_electron_densities(n);
          // printout("nne: %g\n", get_nne(n));
        }
      }
      else
      /// For all other timesteps temperature corrections have to be applied
      {
        /// we have to calculate the electron density
        /// and all the level populations
        /// Normalise estimators and make sure that they are finite.
        /// Then update T_R and W using the estimators.
        /// (This could in principle also be done for empty cells)

        const double estimator_normfactor = 1 / (deltaV * deltat) / nprocs / assoc_cells;
        const double estimator_normfactor_over4pi = ONEOVER4PI * estimator_normfactor;
        const time_t sys_time_start_temperature_corrections = time(NULL);

        radfield_normalise_J(n, estimator_normfactor_over4pi); // this applies normalisation to the fullspec J
        radfield_set_J_normfactor(n, estimator_normfactor_over4pi); // this stores the factor that will be applied later for the J bins but not fullspec J

        #ifdef DO_TITER
          radfield_titer_J(n);
        #endif

        #ifndef FORCE_LTE
        if (initial_iteration || modelgrid[n].thick == 1)
        #endif
        {
          const double T_R = get_T_R_from_J(n);
          set_TR(n, T_R);
          set_Te(n, T_R);
          set_TJ(n, T_R);
          set_W(n, 1);

          #ifndef FORCE_LTE
            #if (!NO_LUT_PHOTOION)
              set_all_corrphotoionrenorm(n, 1.);
            #endif
          #endif

          precalculate_partfuncts(n);
          calculate_populations(n);
        }
        #ifndef FORCE_LTE
        else // not (initial_iteration || modelgrid[n].thick == 1)
        {
          /// Calculate estimators

          radfield_normalise_nuJ(n, estimator_normfactor_over4pi);

          ffheatingestimator[n] *= estimator_normfactor;
          colheatingestimator[n] *= estimator_normfactor;

          #ifdef DO_TITER
            radfield_titer_nuJ(n);
            titer_average_estimators(n);
          #endif

          // Get radiation field parameters out of the full-spectrum and binned J and nuJ estimators
          radfield_fit_parameters(n, nts);

          #if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
            update_gamma_corrphotoionrenorm_bfheating_estimators(n, estimator_normfactor);
          #endif

          grid_cell_solve_Te_nltepops(n, nts, titer);
        }
        #endif
        printout("Temperature/NLTE solution for cell %d timestep %d took %d seconds\n", n, nts, time(NULL) - sys_time_start_temperature_corrections);
      }

      const float nne = get_nne(n);
      const double compton_optical_depth = SIGMA_T * nne * wid_init(n) * tratmid;

      double radial_pos = modelgrid[n].initial_radial_pos * tratmid / assoc_cells;
      if (grid_type == GRID_SPHERICAL1D)
      {
        const double r_inner = get_cellcoordmin(n, 0) * tratmid / assoc_cells;
        const double r_outer = r_inner + wid_init(n) * tratmid;
        radial_pos = 3./4 * (pow(r_outer, 4.) - pow(r_inner, 4.)) / (pow(r_outer, 3) - pow(r_inner, 3.)); // volume averaged mean radius
        // printout("r_inner %g r_outer %g tratmid %g assoc_cells %d\n", r_inner, r_outer, tratmid, assoc_cells);
      }
      const double grey_optical_deptha = get_kappagrey(n) * get_rho(n) * wid_init(n) * tratmid;
      const double grey_optical_depth = get_kappagrey(n) * get_rho(n) * (rmax * tratmid - radial_pos);
      if (log_this_cell)
      {
        printout("modelgridcell %d, compton optical depth (/propgridcell) %g, grey optical depth (/propgridcell) %g\n", n, compton_optical_depth, grey_optical_deptha);
        printout("radial_pos %g, distance_to_obs %g, tau_dist %g\n", radial_pos, rmax * tratmid - radial_pos, grey_optical_depth);
        //printout("rmax %g, tratmid %g\n",rmax,tratmid);
      }
      modelgrid[n].grey_depth = grey_optical_depth;

//          grey_optical_depth = compton_optical_depth;
      if ((grey_optical_depth > cell_is_optically_thick) && (nts < n_grey_timesteps))
      {
        printout("cell %d is treated in grey approximation (tau %g)\n", n, grey_optical_depth);
        modelgrid[n].thick = 1;
      }
      // else if (grey_optical_depth > cell_is_optically_thick_vpkt)
      //   modelgrid[n].thick = 2;
      else
        modelgrid[n].thick = 0;

      /// Temperatures are are now written to file outside the parallel block!

      /// Reset the volume estimators to zero for the current timestep
      /// Moved next lines to zero_estimators routine (SS)
      /*
      nuJ[n] = 0.;
      J[n] = 0.;
      cell[n].heating_ff = 0.;
      cell[n].heating_bf = 0.;
      */


      //set_kappagrey(n,SIGMA_T);

      /// Cooling rates depend only on cell properties, precalculate total cooling
      /// and ion contributions inside update grid and communicate between MPI tasks
      const time_t sys_time_start_calc_kpkt_rates = time(NULL);

      calculate_kpkt_rates(n);
      printout("calculate_kpkt_rates for cell %d timestep %d took %d seconds\n", n, nts, time(NULL) - sys_time_start_calc_kpkt_rates);
    }
    else if (opacity_case == 3)
    {
      /// MK Begin
      //printout("update_grid: opacity_case 3 ... updating cell[n].kappa_grey"); //MK
      if (get_rho(n) > rho_crit)
      {
        set_kappagrey(n, opcase3_normal * (0.9 * get_ffegrp(n) + 0.1) * rho_crit/get_rho(n));
      }
      else
      {
        set_kappagrey(n, opcase3_normal * (0.9 * get_ffegrp(n) + 0.1));
      }
      /// MK End
    }

    if (do_rlc_est == 2 && get_nne(n) > 0)
    {
      const double cell_len_scale_a = 0.1 / get_nne(n) / SIGMA_T;
      if (cell_len_scale_a < mps[tid])
      {
        mps[tid] = cell_len_scale_a;
      }
      const double cell_len_scale_b = 0.1 / get_rho(n) / GREY_OP;
      if (cell_len_scale_b <  mps[tid])
      {
        mps[tid] = cell_len_scale_b;
      }
    }
    printout("update_grid_cell for cell %d timestep %d took %d seconds\n", n, nts, time(NULL) - sys_time_start_update_cell);
  }
  else
  {
    /// For modelgrid cells that are not represented in the simulation grid,
    /// Set grid properties to zero
    set_TR(n, 0.);
    set_TJ(n, 0.);
    set_Te(n, 0.);
    set_W(n, 0.);

    #ifndef FORCE_LTE
      #if (!NO_LUT_PHOTOION)
        zero_gammaestimator(n);
        set_all_corrphotoionrenorm(n, 0.);
      #endif
    #endif
  }
}


void update_grid(const int nts, const int my_rank, const int nstart, const int ndo, const int titer)
// Subroutine to update the matter quantities in the grid cells at the start
//   of the new timestep.
/// m timestep
{
  /// only needed if all level populations should be printed to the output-file
  //double pop, excitedlevelpops;

  //int samplecell;
  //double deltarho;

  //int first_nonempty_cell = -1000;
  //int lastelement = nelements-1;

  //temprange_paras paras2;
  //int status;
  //double x_0,x_lo,x_hi;
  //int kpkt_cuts_determined = 0;
  //char tempfilename[100],ionfractfilename[100],Gammafilename[100],Alphaspfilename[100];
  //FILE *temperature_file,*ionfract_file,*corrphotoion_file,*thermal_file,*Gamma_file,*Alpha_file,*bfcount_file;
  //FILE *gammaest_file,*gammaana_file;

  /*if (&modelgrid[96].composition[0] == NULL)
  {
   printout("fatal error in ts %d abort\n",m);
   abort();
  }*/

  //printout("[debug] update_grid: starting update for timestep %d...\n",m);
  const double trat = time_step[nts].start / tmin;
  const double tratmid = time_step[nts].mid / tmin;

  double mps[MTHREADS];  /// Thread private substitution of max_path_step. Its minimum is
                         /// assigned to max_path_step after the parallel update_grid finished.
  for (int i = 0; i < nthreads; i++)
  {
    mps[i] = 1.e35;
  }

  /// For debug temperature output
  /*
  sprintf(tempfilename,"d%d_thermal_%.4d.out",m,my_rank);
  //sprintf(tempfilename,"t%d_temperature_%d-%d.out",m,nstart,nstart+nblock-1);
  thermal_file = fopen_required(tempfilename, "w");
  setvbuf(thermal_file, NULL, _IOLBF, 1);
  */

  //printout("[debug] update_grid: time before initialisation of heating file %d\n",time(NULL));
  //#ifndef FORCE_LTE
  //  sprintf(filename,"h%d-%d_heating_%.4d.out",m,titer,my_rank);
  //  heating_file = fopen_required(filename, "w")) == NULL);
  //  setvbuf(heating_file, NULL, _IOLBF, 1);
  //#endif
  //printout("[debug] update_grid: heating file initialised %d\n",time(NULL));

  ///Calculate the critical opacity at which opacity_case 3 switches from a
  ///regime proportional to the density to a regime independent of the density
  ///This is done by solving for tau_sobolev == 1
  ///tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/MNI56 * 3000e-8 * time_step[m].mid;
  rho_crit = ME * CLIGHT * MNI56 / (PI * QE * QE * rho_crit_para * 3000e-8 * time_step[nts].mid);
  printout("update_grid: rho_crit = %g\n", rho_crit);
  //printf("time %ld\n",time(NULL));

  // const double t_current = time_step[nts].start;

  // These values will not be used if nts == 0, but set them anyway
  // nts_prev is the previous timestep, unless this is timestep zero
  const int nts_prev = (titer != 0 || nts == 0) ? nts : nts - 1;
  const double deltat = time_step[nts_prev].width;

  printout("timestep %d, titer %d\n", nts, titer);
  printout("deltat %g\n", deltat);

  /*
  FILE *photoion_file;
  FILE *bf_file;
  char photoion_filename[100],bf_filename[100];
  if (m != itstep)
  {
    sprintf(photoion_filename,"photoion_%.2d.out",m);
    photoion_file = fopen_required(photoion_filename, "w");
    setvbuf(photoion_file, NULL, _IOLBF, 1);
    sprintf(bf_filename,"bf_%.2d.out",m);
    bf_file = fopen_required(bf_filename, "w");
    setvbuf(bf_file, NULL, _IOLBF, 1);

    for (ncl = 0; ncl < nblock; ncl++)
    {
      n = nonemptycells[my_rank+ncl*nprocs];
      fprintf(photoion_file,"%d ",n);
      fprintf(bf_file,"%d ",n);
      for (int i = 0; i < 9; i++)
      {
        fprintf(photoion_file,"%g %g ",cell[n].photoion[i],cell[n].radrecomb[i]);
        fprintf(bf_file,"%g %g ",cell[n].bfabs[i],cell[n].bfem[i]);
      }
      fprintf(photoion_file,"\n");
      fprintf(bf_file,"\n");
    }
    fclose(photoion_file);
    fclose(bf_file);
  }
  */


  /// Compare values of rate coefficients obtained by the tabulated values to
  /// directly integrated values.
  //check_interpolation(T_min,T_max);

//   #ifndef FORCE_LTE
//     //sprintf(Alphaspfilename,"bfcount%.2d-%.2d_%.4d.out",m,titer,my_rank);
//     sprintf(Alphaspfilename,"bfcount%d-%d_%.4d.out",m,titer,my_rank);
//     bfcount_file = fopen_required(Alphaspfilename, "w");
//     setvbuf(bfcount_file, NULL, _IOLBF, 1);
//   #endif

  #ifdef _OPENMP
  #pragma omp parallel
    //copyin(nuJ,J,rhosum,T_Rsum,T_esum,Wsum,associatedcells)
  #endif
  {
    /// Do not use values which are saved in the cellhistory within update_grid
    /// and daughter routines (THREADPRIVATE VARIABLE, THEREFORE HERE!)
    use_cellhist = false;
    cellhistory_reset(-99, true);

    /// Reset histindex to unknown
    //printout("thread%d _ update_grid histindex %d\n",tid,histindex);
    //histindex = -99;
    //printout("thread%d _ update_grid histindex %d\n",tid,histindex);

    /// Updating cell information
    #ifdef _OPENMP
      #pragma omp for schedule(dynamic)
      //T_D,W_D,nne,deltarho_old,deltaT_R_old,deltaT_e_old, deltarho,deltaT_R,deltaT_e,i,rhoindex,T_Rindex,T_eindex,ncl)
    #endif
    //for (n = nstart; n < nstart+nblock; n++)
    //for (ncl = 0; ncl < nblock; ncl++)
    //for (ncl = nstart; ncl < nstart+nblock; ncl++)
    //for (n = nstart; n < nstart+nblock; n++)
    for (int mgi = 0; mgi < npts_model; mgi++)
    {
      /// Check if this task should work on the current model grid cell.
      /// If yes, update the cell and write out the estimators
      if (mgi >= nstart && mgi < nstart + ndo)
      {
        // vol_init should a take a cellindex rather than a modelgridindex
        // but in spherical 1D mode these are the same, and in uniform 3D grid, the parameter is ignored
        const double deltaV = vol_init(mgi) * pow(time_step[nts_prev].mid / tmin, 3);
        update_grid_cell(mgi, nts, titer, tratmid, deltaV, deltat, mps);

        //maybe want to add omp ordered here if the modelgrid cells should be output in order
        #ifdef _OPENMP
        #pragma omp critical(estimators_file)
        #endif
        {
          write_to_estimators_file(mgi, nts);
        }
      }
      else
      {
        /// else, only reset gammaestimator to zero. This allows us to do a global MPI
        /// communication after update_grid to synchronize gammaestimator
        /// and write a contiguous restart file with grid properties
        #if (!NO_LUT_PHOTOION)
          zero_gammaestimator(mgi);
        #endif
      }
    } /// end parallel for loop over all modelgrid cells

    /// Now after all the relevant taks of update_grid have been finished activate
    /// the use of the cellhistory for all OpenMP tasks, in what follows (update_packets)
    use_cellhist = true;

  } /// end OpenMP parallel section

  // alterative way to write out estimators. this keeps the modelgrid cells in order but heatingrates are not valid.
  // #ifdef _OPENMP
  // for (int n = nstart; n < nstart+nblock; n++)
  // {
  //   write_to_estimators_file(n,nts);
  // }
  // #endif

  /// Move that check to outside the update grid routine. These numbers
  /// must be the same as they were during grid_init, but on the
  /// grid as a whole. Therefore densities and kappa must be propagated.
  /*
  //fclose(thermal_file);
  check1=check2=0.;
  int mgi;
  for (n=0; n < ngrid; n++)
  {
    mgi=cell[n].modelgridindex;
    check1 = check1 + (get_kappagrey(mgi)  * get_rho(mgi));
    check2 = check2 + get_rho(mgi);
  }
  printout("Grey normalisation check: %g\n", check1/check2);
  */

  /// Assign the minimum of thread private mps to the global variable max_path_step
  max_path_step = mps[0];
  for (int i = 1; i < nthreads; i++)
  {
    if (mps[i] < max_path_step)
    {
      max_path_step = mps[i];
    }
  }

  if (do_rlc_est == 2)
  {
    if (max_path_step < (wid_init(0) * trat / 10.))
    {
      max_path_step = wid_init(0) / 10. * trat;
    }
  }
  //printout("max_path_step %g\n", max_path_step);

  //#ifndef FORCE_LTE
  //  fclose(heating_file);
  //#endif
}


double calculate_populations(const int modelgridindex)
/// Determines the electron number density for a given cell using one of
/// libgsl's root_solvers and calculates the depending level populations.
{
  /// Initialise the gsl solver
  const gsl_root_fsolver_type *const solvertype = gsl_root_fsolver_brent;
  gsl_root_fsolver *solver;
  solver = gsl_root_fsolver_alloc(solvertype);

  /// and the solution function
  gsl_function f;
  nne_solution_paras paras;
  paras.cellnumber = modelgridindex;
  f.function = &nne_solution_f;
  f.params = &paras;

  neutral_flag = false;

  /// Get temperatures
  const double T_R = get_TR(modelgridindex);
  const float T_e = get_Te(modelgridindex);
  const double W = get_W(modelgridindex);

  double nne_hi = get_rho(modelgridindex) / MH;

  /// The following section of uppermost_ion is (so far) NOT thread safe!!!!!!!!!!!!!!!!!!!!!!!
  int only_neutral_ions = 0;
  int nelements_in_cell = 0;
  for (int element = 0; element < nelements; element++)
  {
    const int nions = get_nions(element);
    //elements[element].uppermost_ion = nions-1;
    elements_uppermost_ion[tid][element] = nions - 1;
    const double abundance = get_abundance(modelgridindex,element);
    if (abundance > 0)
    {
      int uppermost_ion;
#     ifdef FORCE_LTE
      uppermost_ion = get_nions(element) - 1;
#     else
      if (initial_iteration || modelgrid[modelgridindex].thick == 1)
      {
        uppermost_ion = get_nions(element) - 1;
      }
      else
      {
        int ion;
        for (ion = 0; ion < nions-1; ion++)
        {
          //printout("element %d, ion %d, photoionest %g\n",element,ion,photoionestimator[modelgridindex*nelements*maxion+element*maxion+ion]);
          //if (photoionestimator[modelgridindex*nelements*maxion+element*maxion+ion] == 0) break;
          #if NO_LUT_PHOTOION
            const double Gamma = calculate_iongamma_per_gspop(modelgridindex, element, ion);
          #else
            const double Gamma = gammaestimator[modelgridindex * nelements * maxion + element * maxion + ion];
          #endif

          if ((Gamma == 0) &&
             (!NT_ON || ((rpkt_emiss[modelgridindex] == 0.) && (get_modelradioabund(modelgridindex, NUCLIDE_CR48) == 0.) && (get_modelradioabund(modelgridindex, NUCLIDE_NI56) == 0.))))
            break;
        }
        uppermost_ion = ion;
      }
#     endif
      //printout("cell %d, element %d, uppermost_ion by Gamma is %d\n",cellnumber,element,uppermost_ion);

      double factor = 1.;
      int ion;
      for (ion = 0; ion < uppermost_ion; ion++)
      {
        factor *= nne_hi * phi(element, ion, modelgridindex);
        //printout("element %d, ion %d, factor %g\n",element,i,factor);
        if (!isfinite(factor))
        {
          printout("[info] calculate_populations: uppermost_ion limited by phi factors for element Z=%d, ionstage %d in cell %d\n",
                   get_element(element), get_ionstage(element, ion), modelgridindex);
          break;
        }
      }
      uppermost_ion = ion;
      //printout("cell %d, element %d, final uppermost_ion is %d, factor %g\n",modelgridindex,element,uppermost_ion,factor);
      //elements[element].uppermost_ion = uppermost_ion;
      elements_uppermost_ion[tid][element] = uppermost_ion;
      if (uppermost_ion == 0)
        only_neutral_ions++;
      nelements_in_cell++;
    }
  }

  float nne = 0.;
  double nne_tot = 0.;   /// total number of electrons in grid cell which are possible
                         /// targets for compton scattering of gamma rays
  double nntot = 0.;
  if (only_neutral_ions == nelements_in_cell)
  {
    /// Special case of only neutral ions, set nne to some finite value that
    /// packets are not lost in kpkts
    /// Introduce a flag variable which is sent to the T_e solver so that
    /// we get this info only once when T_e is converged and not for each
    /// iteration step.
    neutral_flag = true;
    //printout("[warning] calculate_populations: only neutral ions in cell %d modelgridindex\n",modelgridindex);
    //abort();
    /// Now calculate the ground level populations in nebular approximation and store them to the grid
    for (int element = 0; element < nelements; element++)
    {
      const double abundance = get_abundance(modelgridindex,element);
      /// calculate number density of the current element (abundances are given by mass)
      const double nnelement = abundance / elements[element].mass * get_rho(modelgridindex);
      nne_tot += nnelement * get_element(element);

      const int nions = get_nions(element);
      /// Assign the species population to the neutral ion and set higher ions to MINPOP
      for (int ion = 0; ion < nions; ion++)
      {
        double nnion;
        if (ion == 0)
          nnion = nnelement;
        else if (abundance > 0.)
          nnion = MINPOP;
        else
          nnion = 0.;
        nntot += nnion;
        nne += nnion * (get_ionstage(element,ion)-1);
        modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = (
          nnion * stat_weight(element,ion,0) / modelgrid[modelgridindex].composition[element].partfunct[ion]);

        if (!isfinite(modelgrid[modelgridindex].composition[element].groundlevelpop[ion]))
          printout("[warning] calculate_populations: groundlevelpop infinite in connection with MINPOP\n");
      }
    }
    nntot += nne;
    if (nne < MINPOP)
      nne = MINPOP;
    set_nne(modelgridindex,nne);
  }
  else
  {
    /// Apply solver to get nne
    /// Search solution for nne in [nne_lo,nne_hi]
    //printout("nne@x_lo %g\n", nne_solution_f(nne_lo,f.params));
    //printout("nne@x_hi %g\n", nne_solution_f(nne_hi,f.params));
    //printout("n, x_lo, x_hi, T_R, T_e, W, rho %d, %g, %g, %g, %g, %g, %g\n",modelgridindex,x_lo,x_hi,T_R,T_e,W,cell[modelgridindex].rho);
    double nne_lo = 0.;  //MINPOP;
    if (nne_solution_f(nne_lo,f.params)*nne_solution_f(nne_hi,f.params) > 0)
    {
      printout("n, nne_lo, nne_hi, T_R, T_e, W, rho %d, %g, %g, %g, %g, %g, %g\n",modelgridindex,nne_lo,nne_hi,T_R,T_e,W,get_rho(modelgridindex));
      printout("nne@x_lo %g\n", nne_solution_f(nne_lo,f.params));
      printout("nne@x_hi %g\n", nne_solution_f(nne_hi,f.params));
#     ifndef FORCE_LTE
      for (int element = 0; element < nelements; element++)
      {
        //printout("cell %d, element %d, uppermost_ion is %d\n",modelgridindex,element,elements[element].uppermost_ion);
        printout("cell %d, element %d, uppermost_ion is %d\n",modelgridindex,element,elements_uppermost_ion[tid][element]);
        //for (ion=0; ion <= elements[element].uppermost_ion; ion++)
        #if (!NO_LUT_PHOTOION)
          for (int ion = 0; ion <= elements_uppermost_ion[tid][element]; ion++)
          {
            //printout("element %d, ion %d, photoionest %g\n",element,ion,photoionestimator[modelgridindex*nelements*maxion+element*maxion+ion]);
            printout("element %d, ion %d, photoionest %g\n",element,ion,gammaestimator[modelgridindex*nelements*maxion+element*maxion+ion]);
          }
        #endif
      }
#     endif
    }
    gsl_root_fsolver_set(solver, &f, nne_lo, nne_hi);
    int iter = 0;
    const int maxit = 100;
    const double fractional_accuracy = 1e-3;
    int status;
    do
    {
      iter++;
      gsl_root_fsolver_iterate(solver);
      nne = gsl_root_fsolver_root(solver);
      nne_lo = gsl_root_fsolver_x_lower(solver);
      nne_hi = gsl_root_fsolver_x_upper(solver);
      status = gsl_root_test_interval(nne_lo, nne_hi, 0, fractional_accuracy);
      //if (first_nonempty_cell == -1000) printout("[debug] update_grid:   %d [%g, %g] %g %g\n",iter,x_lo,x_hi,x_0,x_hi-x_lo);
    }
    while (status == GSL_CONTINUE && iter < maxit);

    gsl_root_fsolver_free(solver);

    if (nne < MINPOP)
      nne = MINPOP;

    set_nne(modelgridindex,nne);
    //cell[modelgridindex].nne = nne;
    if (status == GSL_CONTINUE)
      printout("[warning] calculate_populations: nne did not converge within %d iterations\n",maxit);
    //printout("[debug] update_grid:   status = %s\n",gsl_strerror(status));
    //printout("[debug] update_grid:   converged nne %g\n",cell[modelgridindex].nne);

    /// Now calculate the ground level populations in nebular approximation and store them to the grid
    double nne_check = 0.;
    nne_tot = 0.;   /// total number of electrons in grid cell which are possible
                    /// targets for compton scattering of gamma rays

    nntot = nne;
    for (int element = 0; element < nelements; element++)
    {
      const double abundance = get_abundance(modelgridindex,element);
      const int nions = get_nions(element);
      /// calculate number density of the current element (abundances are given by mass)
      const double nnelement = abundance / elements[element].mass * get_rho(modelgridindex);
      nne_tot += nnelement * get_element(element);

      /// Use ionizationfractions to calculate the groundlevel populations
      for (int ion = 0; ion < nions; ion++)
      {
        double nnion;
        //if (ion <= elements[element].uppermost_ion)
        if (ion <= elements_uppermost_ion[tid][element])
        {
          if (abundance > 0)
          {
            nnion = nnelement * ionfract(element,ion,modelgridindex,nne);
            if (nnion < MINPOP)
              nnion = MINPOP;
          }
          else
            nnion = 0.;
        }
        else
          nnion = MINPOP;  /// uppermost_ion is only < nions-1 in cells with nonzero abundance of the given species
        nntot += nnion;
        nne_check += nnion * (get_ionstage(element,ion)-1);
        //if (modelgrid[modelgridindex].composition[element].groundlevelpop[ion] < 0)
        //if (initial_iteration || modelgrid[modelgridindex].thick == 1)
        {
          modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = (nnion *
                stat_weight(element,ion,0) / modelgrid[modelgridindex].composition[element].partfunct[ion]);
          //printout("calculate_populations: setting groundlevelpop of ion %d\n",ion);
        }
        //else
        {
          //printout("calculate_populations: not setting groundlevelpop of ion %d\n",ion);
        }
  /*
  if (element == 21)
    {
      printout("Setting ion %d to have glp %g in cell %d\n", ion, modelgrid[modelgridindex].composition[element].groundlevelpop[ion], modelgridindex);
      printout("ion pop was %g and partfn was %g\n", nnion, modelgrid[modelgridindex].composition[element].partfunct[ion]);
      printout("the ion frac was %g, abundance %g and density %g\n",ionfract(element,ion,modelgridindex,nne), abundance, get_rho(modelgridindex));
    }
  */

        if (!isfinite(modelgrid[modelgridindex].composition[element].groundlevelpop[ion]))
          printout("[warning] calculate_populations: groundlevelpop infinite in connection with MINPOP\n");
      }
    }
  }

  set_nnetot(modelgridindex, nne_tot);
  return nntot;
}


double calculate_electron_densities(const int modelgridindex)
// Determines the free and total electron number densities
// for a given cell and stores them, assuming ion populations (ground_level_pop and partfunc)
// are fixed (determined by NLTE all-ion solver)
{
  double nne_tot = 0.; // total electron density
  float nne = 0.;     // free electron density

  for (int element = 0; element < nelements; element++)
  {
    const double elem_abundance = get_abundance(modelgridindex,element);
    // calculate number density of the current element (abundances are given by mass)
    const double nnelement = elem_abundance / elements[element].mass * get_rho(modelgridindex);
    nne_tot += nnelement * get_element(element);

    // Use ionization fractions to calculate the free electron contributions
    if (elem_abundance > 0)
    {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++)
      {
        //if (ion <= elements[element].uppermost_ion)
        nne += (get_ionstage(element,ion)-1) * ionstagepop(modelgridindex,element,ion);
      }
    }
  }

  set_nne(modelgridindex, nne);
  set_nnetot(modelgridindex, nne_tot);
  return nne_tot;
}

/*
double nuB_nu_integrand(double nu, void *paras)
{
  double T = ((gslintegration_paras *) paras)->T;

  return nu * TWOHOVERCLIGHTSQUARED * pow(nu,3) * 1/(exp(H*nu/KB/T)-1);
}

double B_nu_integrand(double nu, void *paras)
{
  double T = ((gslintegration_paras *) paras)->T;

  return TWOHOVERCLIGHTSQUARED * pow(nu,3) * 1/(exp(H*nu/KB/T)-1);
}
*/


#ifndef FORCE_LTE

/*
double get_ffcooling(int element, int ion, int cellnumber)
{
  double ionstagepop(int cellnumber, int element, int ion);
  int ioncharge;
  double T_e,nne,nncurrention,C_ff;

  ion++;
  nncurrention = ionstagepop(cellnumber,element,ion);
  nne = cell[cellnumber].nne;
  T_e = cell[cellnumber].T_e;

  C_ff = 0.;
  ioncharge = get_ionstage(element,ion) - 1;
  if (ioncharge > 0)
  {
    C_ff += 1.426e-27 * sqrt(T_e) * pow(ioncharge,2) * nncurrention * nne;
  }

  return C_ff;
}

double get_adcool(int cellnumber, int nts)
{

  double t_current,nntot,T_e;
  double p,dV,V;

  t_current = time_step[nts].mid;
  nntot = cell[cellnumber].nn_tot;
  T_e = cell[cellnumber].T_e;
  p = nntot*KB*T_e;
  dV = 3*pow(wid_init/tmin,3)*pow(t_current,2);
  V = pow(wid_init*t_current/tmin,3);

  return p*dV/V;
}

double get_radrecomb(int element, int ion, int level, int cellnumber)
{
  double T_e,nne,nnion,E_threshold,alpha_sp;

  T_e = cell[cellnumber].T_e;
  nne = cell[cellnumber].nne;
  nnion = get_groundlevelpop(cellnumber,element,ion+1);
  E_threshold = epsilon(element,ion+1,0) - epsilon(element,ion,level);
  alpha_sp = interpolate_spontrecombcoeff(element,ion,level,T_e);
  ///// This is ion-wise Alpha_st is the sum of the level-wise alpha_st, so it should be added only once if mulitple levels are present
  //double Alpha_st = stimrecombestimator_save[cellnumber*nelements*maxion+element*maxion+ion];

  //return nnion*nne*(alpha_sp+Alpha_st)*E_threshold;
  return nnion*nne*alpha_sp*E_threshold;
}

*/

/*
double get_thermalratio(int element, int ion, int level, int cellnumber)
{
  double T_e,alpha_sp,alpha_sp_E;

  T_e = cell[cellnumber].T_e;
  alpha_sp = interpolate_spontrecombcoeff(element,ion,level,T_e);
  alpha_sp_E = interpolate_spontrecombcoeff_E(element,ion,level,T_e);

  return alpha_sp/alpha_sp_E;
}*/




/*

double get_Gamma(int cellnumber, int element, int ion)
///photoionization rate: paperII 3.5.2
/// n_1 - occupation number of ground state
{
  int i,level;
  const int nlevels = get_nlevels(element,ion);
  double Gamma = 0.;

  for (i = 0; i < nlevels; i++)
  {
    Gamma += calculate_exclevelpop(cellnumber,element,ion,level) * get_corrphotoioncoeff(element,ion,level,cellnumber);
  }
  Gamma /= calculate_exclevelpop(cellnumber,element,ion,0);

  return Gamma;
}


double get_Gamma_phys(int cellnumber, int element, int ion)
///photoionization rate: paperII 3.5.2
/// n_1 - occupation number of ground state
{
  int i,level;
  const int nlevels = get_nlevels(element,ion);
  double Gamma = 0.;

  //for (i = 0; i < nlevels; i++)
  for (i = 0; i < TAKE_N_BFCONTINUA; i++)
  {
    Gamma += get_corrphotoioncoeff_ana(element,ion,level,cellnumber);
  }

  return Gamma;
}

*/
#endif


void write_grid_restart_data(void)
{
  FILE *restrict gridsave_file = fopen_required("gridsave.dat", "w");

  for (int mgi = 0; mgi < npts_model; mgi++)
  {
    const bool nonemptycell = (mg_associated_cells[mgi] > 0);

    if (nonemptycell)
    {
      fprintf(gridsave_file, "%d %g %g %g %g %hd %lg",
              mgi, get_TR(mgi), get_Te(mgi), get_W(mgi), get_TJ(mgi),
              modelgrid[mgi].thick, rpkt_emiss[mgi]);
    }
    else
    {
      fprintf(gridsave_file, "%d %g %g %g %g %d %lg", mgi, 0., 0., 0., 0., 0, 0.);
    }

    #ifndef FORCE_LTE
      #if (!NO_LUT_PHOTOION)
        for (int element = 0; element < nelements; element++)
        {
          const int nions = get_nions(element);
          for (int ion = 0; ion < nions; ion++)
          {
            const int estimindex = mgi * nelements * maxion + element * maxion + ion;
            fprintf(gridsave_file, " %lg %lg",
                    (nonemptycell ? corrphotoionrenorm[estimindex] : 0.),
                    (nonemptycell ? gammaestimator[estimindex] : 0.));
          }
        }
      #endif
    #endif
    fprintf(gridsave_file,"\n");
  }

  // the order of these calls is very important!
  radfield_write_restart_data(gridsave_file);
  nt_write_restart_data(gridsave_file);
  nltepop_write_restart_data(gridsave_file);
  fclose(gridsave_file);
}


/*static int get_cell(double x, double y, double z, double t)
/// subroutine to identify the cell index from a position and a time.
{
  // Original version of this was general but very slow. Modifying to be
  // faster but only work for regular grid.

  double trat = t / tmin;
  int nx = (x - (cell[0].pos_init[0] * trat))/(wid_init * trat);
  int ny = (y - (cell[0].pos_init[1] * trat))/(wid_init * trat);
  int nz = (z - (cell[0].pos_init[2] * trat))/(wid_init * trat);

  int n = nx + (ncoordgrid[0] * ny) + (ncoordgrid[0] * ncoordgrid[1] * nz);

  // do a check

  if (x < cell[n].pos_init[0] * trat)
  {
    printout("Problem with get_cell (1).\n");
    abort();
  }
  if (x > (cell[n].pos_init[0]+wid_init) * trat)
  {
    printout("Problem with get_cell (2).\n");
    abort();
  }
  if (y < cell[n].pos_init[1] * trat)
  {
    printout("Problem with get_cell (3).\n");
    abort();
  }
  if (y > (cell[n].pos_init[1]+wid_init) * trat)
  {
    printout("Problem with get_cell (4).\n");
    abort();
  }
  if (z < cell[n].pos_init[2] * trat)
  {
    printout("Problem with get_cell (5).\n");
    abort();
  }
  if (z > (cell[n].pos_init[2]+wid_init) * trat)
  {
    printout("Problem with get_cell (6).\n");
    abort();
  }
  return(n);


  // OLD
        //   trat = t / tmin;
        //
        //   for (n = 0; n < ngrid; n++)
        //   {
        //   if (
        //   (x > cell[n].pos_init[0] * trat) &&
        //   (x < (cell[n].pos_init[0] + wid_init) *trat) &&
        //   (y > cell[n].pos_init[1] * trat) &&
        //   (y < (cell[n].pos_init[1] + wid_init) *trat) &&
        //   (z > cell[n].pos_init[2] * trat) &&
        //   (z < (cell[n].pos_init[2] + wid_init) *trat))
        //   {
        //   return(n);
        // }
        // }
  // END OLD

  printout("Failed to find cell (get_cell). \n");
  printout("x %g, y %g, z %g, t %g\n", x, y, z, t);
  printout("xend %g yend %g zend %g\n", cell[ngrid-1].pos_init[0] * trat, cell[ngrid-1].pos_init[1] * trat,cell[ngrid-1].pos_init[2] * trat);
  printout("xend %g yend %g zend %g\n", cell[0].pos_init[0] * trat, cell[0].pos_init[1] * trat,cell[0].pos_init[2] * trat);
  printout("xend %g yend %g zend %g\n", (cell[0].pos_init[0]+wid_init) * trat, (cell[0].pos_init[1]+wid_init) * trat,(cell[0].pos_init[2]+wid_init) * trat);
  printout("xwid %g ywid %g zwid %g\n", (wid_init) * trat, (wid_init) * trat,(wid_init) * trat);

  abort();
}*/

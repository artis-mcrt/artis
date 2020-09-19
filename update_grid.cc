#include <unistd.h>
#include <gsl/gsl_roots.h>
#include "sn3d.h"
#include "atomic.h"
#include "grid_init.h"
#include "kpkt.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "decay.h"
#include "nltepop.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "rpkt.h"
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
      modelgrid[modelgridindex].composition[element].partfunct[ion] = calculate_partfunct(element, ion, modelgridindex);
    }
  }
}

static void write_to_estimators_file(FILE *estimators_file, const int mgi, const int timestep, const int titer, const heatingcoolingrates_t *heatingcoolingrates)
{
  // return; disable for better performance (if estimators files are not needed)
  if (get_numassociatedcells(mgi) > 0)
  {
    const float T_e = get_Te(mgi);
    const float nne = get_nne(mgi);
    //fprintf(estimators_file,"%d %g %g %g %g %d ",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),modelgrid[n].thick);
    //fprintf(estimators_file,"%d %g %g %g %g %g ",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),grey_optical_depth);
    fprintf(estimators_file, "timestep %d modelgridindex %d titeration %d TR %g Te %g W %g TJ %g grey_depth %g nne %g tdays %7.2f\n",
            timestep, mgi, titer, get_TR(mgi), T_e, get_W(mgi), get_TJ(mgi), modelgrid[mgi].grey_depth, nne, time_step[timestep].mid / DAY);
    //fprintf(estimators_file,"%d %g %g %g %g %g %g %g
    //",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),grey_optical_depth,grey_optical_deptha,compton_optical_depth);

    if (NLTE_POPS_ON)  //  && timestep % 2 == 0
    {
      nltepop_write_to_file(mgi, timestep);
    }

    for (int element = 0; element < nelements; element++)
    {
      if (get_abundance(mgi, element) <= 0.) // skip elements with no abundance
        continue;

      fprintf(estimators_file, "populations        Z=%2d", get_element(element));
      const int nions = get_nions(element);
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      double elpop = 0.;
      for (int ion = 0; ion < nions; ion++)
      {
        elpop += ionstagepop(mgi, element, ion);
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ionstagepop(mgi, element, ion));
      }
      fprintf(estimators_file, "  SUM: %9.3e\n", elpop);

      const bool printdebug = false;

      bool assume_lte = true;
      bool per_gmpop = true;
      const bool lower_superlevel_only = false;

      if (timestep % 10 == 0)
      {
        fprintf(estimators_file, "RRC_LTE_Nahar      Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions; ion++)
        {
          fprintf(estimators_file, "  %d: %9.3e",
                  get_ionstage(element, ion),
                  calculate_ionrecombcoeff(-1, T_e, element, ion, assume_lte, false, printdebug, lower_superlevel_only, per_gmpop, false));
        }
        fprintf(estimators_file, "\n");
      }

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

      #if (TRACK_ION_STATS && TRACK_ION_MASTATS)
      fprintf(estimators_file, "MA_IN_RADEXC       Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      double ma_el = 0.;
      for (int ion = 0; ion < nions; ion++)
      {
        const double ma_ion = get_ion_stats(mgi, element, ion, ION_COUNTER_MACROATOM_ENERGYIN_RADEXC);
        ma_el += ma_ion;
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion), ma_ion);
      }
      fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);
      #endif

      #if (TRACK_ION_STATS && TRACK_ION_MASTATS)
      fprintf(estimators_file, "MA_OUT_RADEEXC     Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      ma_el = 0.;
      for (int ion = 0; ion < nions; ion++)
      {
        const double ma_ion = get_ion_stats(mgi, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_RADDEEXC);
        ma_el += ma_ion;
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion), ma_ion);
      }
      fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);
      #endif

      #if (TRACK_ION_STATS && TRACK_ION_MASTATS)
      fprintf(estimators_file, "MA_IN_COLEXC       Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      ma_el = 0.;
      for (int ion = 0; ion < nions; ion++)
      {
        const double ma_ion = get_ion_stats(mgi, element, ion, ION_COUNTER_MACROATOM_ENERGYIN_COLLEXC);
        ma_el += ma_ion;
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion), ma_ion);
      }
      fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);
      #endif

      #if (TRACK_ION_STATS && TRACK_ION_MASTATS)
      fprintf(estimators_file, "MA_OUT_COLDEEXC    Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      ma_el = 0.;
      for (int ion = 0; ion < nions; ion++)
      {
        const double ma_ion = get_ion_stats(mgi, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_COLLDEEXC);
        ma_el += ma_ion;
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion), ma_ion);
      }
      fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);
      #endif

      #if (TRACK_ION_STATS && TRACK_ION_MASTATS)
      fprintf(estimators_file, "MA_IN_PHOTOION     Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      ma_el = 0.;
      for (int ion = 0; ion < nions; ion++)
      {
        const double ma_ion = get_ion_stats(mgi, element, ion, ION_COUNTER_MACROATOM_ENERGYIN_PHOTOION);
        ma_el += ma_ion;
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion), ma_ion);
      }
      fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);
      #endif

      #if (TRACK_ION_STATS && TRACK_ION_MASTATS)
      fprintf(estimators_file, "MA_OUT_RADRECOMB   Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      ma_el = 0.;
      for (int ion = 0; ion < nions; ion++)
      {
        const double ma_ion = get_ion_stats(mgi, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_RADRECOMB);
        ma_el += ma_ion;
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion), ma_ion);
      }
      fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);
      #endif

      #if (TRACK_ION_STATS && TRACK_ION_MASTATS)
      fprintf(estimators_file, "MA_IN_COLION       Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      ma_el = 0.;
      for (int ion = 0; ion < nions; ion++)
      {
        const double ma_ion = get_ion_stats(mgi, element, ion, ION_COUNTER_MACROATOM_ENERGYIN_COLLION);
        ma_el += ma_ion;
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion), ma_ion);
      }
      fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);
      #endif

      #if (TRACK_ION_STATS && TRACK_ION_MASTATS)
      fprintf(estimators_file, "MA_OUT_COLRECOMB   Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      ma_el = 0.;
      for (int ion = 0; ion < nions; ion++)
      {
        const double ma_ion = get_ion_stats(mgi, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_COLLRECOMB);
        ma_el += ma_ion;
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion), ma_ion);
      }
      fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);
      #endif

      #if (TRACK_ION_STATS && TRACK_ION_MASTATS)
      fprintf(estimators_file, "MA_IN_NTCOLION     Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      ma_el = 0.;
      for (int ion = 0; ion < nions; ion++)
      {
        const double ma_ion = get_ion_stats(mgi, element, ion, ION_COUNTER_MACROATOM_ENERGYIN_NTCOLLION);
        ma_el += ma_ion;
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion), ma_ion);
      }
      fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);
      #endif

      #if (TRACK_ION_STATS && TRACK_ION_MASTATS)
      fprintf(estimators_file, "MA_IN_TOTAL        Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      ma_el = 0.;
      for (int ion = 0; ion < nions; ion++)
      {
        const double ma_in_ion = get_ion_stats(mgi, element, ion, ION_COUNTER_MACROATOM_ENERGYIN_TOTAL);
        ma_el += ma_in_ion;
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion), ma_in_ion);
      }
      fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);
      #endif

      #if (TRACK_ION_STATS && TRACK_ION_MASTATS)
      fprintf(estimators_file, "MA_OUT_TOTAL       Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      ma_el = 0.;
      for (int ion = 0; ion < nions; ion++)
      {
        const double ma_ion = get_ion_stats(mgi, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_TOTAL);
        ma_el += ma_ion;
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion), ma_ion);
      }
      fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);
      #endif

      #if (TRACK_ION_STATS && TRACK_ION_MASTATS)
      fprintf(estimators_file, "MA_IN_INTERNAL     Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      ma_el = 0.;
      for (int ion = 0; ion < nions; ion++)
      {
        const double ma_ion = get_ion_stats(mgi, element, ion, ION_COUNTER_MACROATOM_ENERGYIN_INTERNAL);
        ma_el += ma_ion;
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion), ma_ion);
      }
      fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);
      #endif

      #if (TRACK_ION_STATS && TRACK_ION_MASTATS)
      fprintf(estimators_file, "MA_OUT_INTERNAL    Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      ma_el = 0.;
      for (int ion = 0; ion < nions; ion++)
      {
        const double ma_ion = get_ion_stats(mgi, element, ion, ION_COUNTER_MACROATOM_ENERGYOUT_INTERNAL);
        ma_el += ma_ion;
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion), ma_ion);
      }
      fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);
      #endif

      // spontaneous radiative recombination rate coefficient (may or may not include stim. recomb)
      fprintf(estimators_file, "AlphaR*nne         Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions; ion++)
      {
        const bool printdebug = false;
        // const bool printdebug = (get_element(element) >= 26);

        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                calculate_ionrecombcoeff(mgi, T_e, element, ion, assume_lte, false, printdebug, lower_superlevel_only, per_gmpop, false) * nne);
      }
      fprintf(estimators_file, "\n");

      if (timestep % 10 == 0)
      {
        fprintf(estimators_file, "AlphaR_toSL*nne    Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions; ion++)
        {
          const bool printdebug = false;
          // const bool printdebug = (get_element(element) >= 26);

          fprintf(estimators_file, "  %d: %9.3e",
                  get_ionstage(element, ion),
                  calculate_ionrecombcoeff(mgi, T_e, element, ion, assume_lte, false, printdebug, true, per_gmpop, false) * nne);
        }
        fprintf(estimators_file, "\n");
      }

      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "AlphaR_MC_MA*nne   Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions; ion++)
      {
        const double alpha_r_mc = get_ion_stats(mgi, element, ion, ION_COUNTER_RADRECOMB_MACROATOM);
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                alpha_r_mc);
      }
      fprintf(estimators_file, "\n");
      #endif

      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "AlphaR_MC_KPKT*nne Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions; ion++)
      {
        const double alpha_r_mc = get_ion_stats(mgi, element, ion, ION_COUNTER_RADRECOMB_KPKT);
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                alpha_r_mc);
      }
      fprintf(estimators_file, "\n");
      #endif

      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "AlphaR_MC*nne      Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions; ion++)
      {
        const double alpha_r_mc = get_ion_stats(mgi, element, ion, ION_COUNTER_RADRECOMB_MACROATOM) + get_ion_stats(mgi, element, ion, ION_COUNTER_RADRECOMB_KPKT);
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                alpha_r_mc);
      }
      fprintf(estimators_file, "\n");
      #endif

      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "BF_escfrac         Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions; ion++)
      {
        const double alpha_r_mc = get_ion_stats(mgi, element, ion, ION_COUNTER_RADRECOMB_MACROATOM) + get_ion_stats(mgi, element, ion, ION_COUNTER_RADRECOMB_KPKT);
        const double alpha_r_mc_abs = get_ion_stats(mgi, element, ion, ION_COUNTER_RADRECOMB_ABSORBED);
        fprintf(estimators_file, "  %d: %9.3f",
                get_ionstage(element, ion),
                1. - alpha_r_mc_abs / alpha_r_mc);
      }
      fprintf(estimators_file, "\n");
      #endif

      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "BF_escfrac_vpkt    Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions; ion++)
      {
        const double alpha_r_mc = get_ion_stats(mgi, element, ion, ION_COUNTER_RADRECOMB_MACROATOM) + get_ion_stats(mgi, element, ion, ION_COUNTER_RADRECOMB_KPKT);
        const double alpha_r_mc_escaped = get_ion_stats(mgi, element, ion, ION_COUNTER_RADRECOMB_ESCAPED);
        fprintf(estimators_file, "  %d: %9.3f",
                get_ionstage(element, ion),
                alpha_r_mc_escaped / alpha_r_mc);
      }
      fprintf(estimators_file, "\n");
      #endif


      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "BB_escfrac         Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions; ion++)
      {
        const double bb_emitted = get_ion_stats(mgi, element, ion, ION_COUNTER_BOUNDBOUND_MACROATOM);
        const double bb_abs = get_ion_stats(mgi, element, ion, ION_COUNTER_BOUNDBOUND_ABSORBED);
        fprintf(estimators_file, "  %d: %9.3f",
                get_ionstage(element, ion),
                1. - bb_abs / bb_emitted);
      }
      fprintf(estimators_file, "\n");
      #endif

      // stimulated recombination rate coefficient
      // fprintf(estimators_file, "Alpha_stim*nne   Z=%2d", get_element(element));
      // for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //   fprintf(estimators_file, "              ");
      // for (int ion = 0; ion < nions; ion++)
      // {
      //   // const bool printdebug = false;
      //   const bool printdebug = (get_element(element) >= 26);
      //
      //   fprintf(estimators_file, "  %d: %9.3e",
      //           get_ionstage(element, ion),
      //           calculate_ionrecombcoeff(mgi, T_e, element, ion, assume_lte, false, printdebug, lower_superlevel_only, per_gmpop, true) * nne);
      // }
      // fprintf(estimators_file, "\n");

      // thermal collisional recombination rate coefficient
      // fprintf(estimators_file, "Alpha_C*nne      Z=%2d", get_element(element));
      // for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //   fprintf(estimators_file, "              ");
      // for (int ion = 0; ion < nions; ion++)
      // {
      //   fprintf(estimators_file, "  %d: %9.3e",
      //           get_ionstage(element, ion),
      //           calculate_ionrecombcoeff(mgi, T_e, element, ion, assume_lte, true, printdebug, lower_superlevel_only, per_gmpop) * nne);
      // }
      // fprintf(estimators_file, "\n");

      if (timestep % 20 == 0)
      {
        fprintf(estimators_file, "kappa_bf(nuedge)   Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions - 1; ion++)
        {
          double nu_edge = (epsilon(element, ion + 1, 0) - epsilon(element, ion, 0)) / H;
          double kappa_bf = 0.;
          double kappa_fb = 0.;
          calculate_kappa_bf_fb_gammacontr(mgi, nu_edge, &kappa_bf, &kappa_fb);

          fprintf(estimators_file, "  %d: %9.3e",
                  get_ionstage(element, ion),
                  kappa_bf);
        }
        fprintf(estimators_file, "\n");
      }

      // this is expensive
      if (!DETAILED_BF_ESTIMATORS_ON || (timestep % 20) == 0)
      {
        fprintf(estimators_file, "gamma_R_integral   Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions - 1; ion++)
        {
          // const bool printdebug_gammar = (get_element(element) == 26 && get_ionstage(element, ion) == 2);
          const bool printdebug_gammar = false;
          fprintf(estimators_file, "  %d: %9.3e",
                  get_ionstage(element, ion),
                  calculate_iongamma_per_ionpop(mgi, T_e, element, ion, assume_lte, false, printdebug_gammar, false));
        }
        fprintf(estimators_file, "\n");
      }

      if (DETAILED_BF_ESTIMATORS_ON)
      {
        fprintf(estimators_file, "gamma_R_bfest      Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions - 1; ion++)
        {
          // const bool printdebug_gammar = ((get_element(element) == 26 || get_element(element) == 28) && get_ionstage(element, ion) >= 2);
          // const bool printdebug_gammar = (get_element(element) >= 26);
          const bool printdebug_gammar = false;
          fprintf(estimators_file, "  %d: %9.3e",
                  get_ionstage(element, ion),
                  calculate_iongamma_per_ionpop(mgi, T_e, element, ion, assume_lte, false, printdebug_gammar, true));
        }
        fprintf(estimators_file, "\n");
      }

      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC         Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++)
      {
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, ION_COUNTER_PHOTOION));
      }
      fprintf(estimators_file, "\n");
      #endif

      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BF      Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++)
      {
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, ION_COUNTER_PHOTOION_FROMBOUNDFREE));
      }
      fprintf(estimators_file, "\n");
      #endif

      // fprintf(estimators_file, "gamma_R_MC_BFsameZ Z=%2d", get_element(element));
      // for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //   fprintf(estimators_file, "              ");
      // for (int ion = 0; ion < nions; ion++)
      // {
      //   fprintf(estimators_file, "  %d: %9.3e",
      //           get_ionstage(element, ion),
      //           get_ion_stats(mgi, element, ion, ION_COUNTER_PHOTOION_FROMBFSAMEELEMENT]);
      // }
      // fprintf(estimators_file, "\n");

      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BF_i+1  Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++)
      {
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, ION_COUNTER_PHOTOION_FROMBFIONPLUSONE));
      }
      fprintf(estimators_file, "\n");
      #endif

      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BF_i+2  Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++)
      {
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, ION_COUNTER_PHOTOION_FROMBFIONPLUSTWO));
      }
      fprintf(estimators_file, "\n");
      #endif

      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BF_i+3  Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++)
      {
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, ION_COUNTER_PHOTOION_FROMBFIONPLUSTHREE));
      }
      fprintf(estimators_file, "\n");
      #endif

      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BFtoSL  Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++)
      {
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, ION_COUNTER_PHOTOION_FROMBFLOWERSUPERLEVEL));
      }
      fprintf(estimators_file, "\n");
      #endif

      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BB      Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++)
      {
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, ION_COUNTER_PHOTOION_FROMBOUNDBOUND));
      }
      fprintf(estimators_file, "\n");
      #endif

      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BB_i+1  Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++)
      {
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, ION_COUNTER_PHOTOION_FROMBOUNDBOUNDIONPLUSONE));
      }
      fprintf(estimators_file, "\n");
      #endif

      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BB_i+2  Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++)
      {
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, ION_COUNTER_PHOTOION_FROMBOUNDBOUNDIONPLUSTWO));
      }
      fprintf(estimators_file, "\n");
      #endif

      #if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BB_i+3  Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++)
      {
        fprintf(estimators_file, "  %d: %9.3e",
                get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, ION_COUNTER_PHOTOION_FROMBOUNDBOUNDIONPLUSTHREE));
      }
      fprintf(estimators_file, "\n");
      #endif

      // fprintf(estimators_file, "gamma_C          Z=%2d", get_element(element));
      // for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //   fprintf(estimators_file, "              ");
      // for (int ion = 0; ion < nions - 1; ion++)
      // {
      //   fprintf(estimators_file, "  %d: %9.3e",
      //           get_ionstage(element, ion),
      //           calculate_iongamma_per_ionpop(mgi, T_e, element, ion, assume_lte, true, printdebug));
      // }
      // fprintf(estimators_file, "\n");

      if (NT_ON)
      {
        fprintf(estimators_file, "gamma_NT           Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions - 1; ion++)
        {
          const double Y_nt = nt_ionization_ratecoeff(mgi, element, ion);
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), Y_nt);
        }
        fprintf(estimators_file, "\n");

        #if (TRACK_ION_STATS)
        fprintf(estimators_file, "gamma_NT_MC        Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions - 1; ion++)
        {
          fprintf(estimators_file, "  %d: %9.3e",
                  get_ionstage(element, ion),
                  get_ion_stats(mgi, element, ion, ION_COUNTER_NTION));
        }
        fprintf(estimators_file, "\n");
        #endif
      }
    }

    #ifndef FORCE_LTE
      #if (!NO_LUT_PHOTOION)
        fprintf(estimators_file, "corrphotoionrenorm: ");
        for (int element = 0; element < nelements; element++)
        {
          const int nions = get_nions(element);
          for (int ion = 0; ion < nions; ion++)
          {
            fprintf(estimators_file,"%g ", corrphotoionrenorm[mgi*nelements*maxion+element*maxion+ion]);
          }
        }
        fprintf(estimators_file, "\n");
        fprintf(estimators_file, "gammaestimator: ");
        for (int element = 0; element < nelements; element++)
        {
          const int nions = get_nions(element);
          for (int ion = 0; ion < nions; ion++)
          {
            fprintf(estimators_file,"%g ", gammaestimator[mgi*nelements*maxion+element*maxion+ion]);
          }
        }
        fprintf(estimators_file,"\n");
      #endif

      fprintf(estimators_file, "heating: ff %11.5e bf %11.5e coll %11.5e       dep %11.5e heating_dep/total_dep %.2f\n",
              heatingcoolingrates->heating_ff, heatingcoolingrates->heating_bf, heatingcoolingrates->heating_collisional, heatingcoolingrates->heating_gamma, heatingcoolingrates->nt_frac_heating);
      fprintf(estimators_file, "cooling: ff %11.5e fb %11.5e coll %11.5e adiabatic %11.5e\n",
              heatingcoolingrates->cooling_ff, heatingcoolingrates->cooling_fb, heatingcoolingrates->cooling_collisional, heatingcoolingrates->cooling_adiabatic);
    #endif
  }
  else
  {
    // modelgrid cells which are not represented in the simulation grid
    fprintf(estimators_file, "timestep %d modelgridindex %d EMPTYCELL\n", timestep, mgi);
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
    return;
  }

  // force rpkt opacities to be recalculated next time they are accessed
  kappa_rpkt_cont[tid].recalculate_required = true;

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
        for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
        {
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].corrphotoioncoeff = -99.;
#if (SEPARATE_STIMRECOMB)
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].stimrecombcoeff = -99.;
#endif
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


static void solve_Te_nltepops(const int n, const int nts, const int titer, heatingcoolingrates_t *heatingcoolingrates)
{
  // needed for the T_e solver, but only depends on the radiation field, which is fixed during the iterations below
  const time_t sys_time_start_calculate_bfheatingcoeffs = time(NULL);
  calculate_bfheatingcoeffs(n);
  printout("calculate_bfheatingcoeffs for cell %d timestep %d took %ld seconds\n", n, nts, time(NULL) - sys_time_start_calculate_bfheatingcoeffs);

  const double covergence_tolerance = 0.04;
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
    const double prev_T_e = get_Te(n);
    const time_t sys_time_start_Te = time(NULL);
    const int nts_for_te = (titer == 0) ? nts - 1 : nts;

    call_T_e_finder(n, nts, time_step[nts_for_te].mid, MINTEMP, MAXTEMP, heatingcoolingrates);

    const int duration_solve_T_e = time(NULL) - sys_time_start_Te;
    const double fracdiff_T_e = fabs((get_Te(n) / prev_T_e) - 1);

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
      if (NLTE_POPS_ALL_IONS_SIMULTANEOUS)
      {
        for (int element = 0; element < nelements; element++)
        {
          solve_nlte_pops_element(element, n, nts, nlte_iter);
        }
      }
      else
      {
        for (int element = 0; element < nelements; element++)
        {
          const int nions = get_nions(element);
          for (int ion = 0; ion < nions-1; ion++)
          {
            const double trial = fabs(solve_nlte_pops_ion(element, ion, n, nts) - 1);

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
        const double fracdiff_nne = fabs((get_nne(n) / nne_prev) - 1);
        nlte_test = fracdiff_nne;
        printout("NLTE solver cell %d timestep %d iteration %d: time spent on: Spencer-Fano %ds, T_e %ds, NLTE populations %ds\n",
                 n, nts, nlte_iter, duration_solve_spencerfano, duration_solve_T_e, duration_solve_nltepops);
        printout("NLTE (Spencer-Fano/Te/pops) solver cell %d timestep %d iteration %d: prev_iter nne %g, new nne is %g, fracdiff %g, prev T_e %g new T_e %g fracdiff %g\n",
                 n, nts, nlte_iter, nne_prev, get_nne(n), nlte_test, prev_T_e, get_Te(n), fracdiff_T_e);
        // damp changes in nne if oscillating to much
        //set_nne(n, (get_nne(n) + nne_prev) / 2.);
      }
      else
      {
        printout("Completed iteration for NLTE population solver in cell %d for timestep %d. Fractional error returned: %g\n",
                 n, nts, nlte_test);
      }

      if (nlte_test <= covergence_tolerance && fracdiff_T_e <= covergence_tolerance)
      {
        printout("NLTE (Spencer-Fano/Te/pops) solver nne converged to tolerance %g <= %g and T_e to tolerance %g <= %g after %d iterations.\n", nlte_test, covergence_tolerance, fracdiff_T_e, covergence_tolerance, nlte_iter + 1);
        break;
      }
      else if (nlte_iter == NLTEITER)
        printout("WARNING: NLTE solver failed to converge after %d iterations. Keeping solution from last iteration\n", nlte_iter + 1);
    }
    else
      break; // no iteration is needed without NLTE_POPS_ON
  }
}


#if (TRACK_ION_STATS)
static void normalise_ion_estimators(const int mgi, const double deltat, const double deltaV)
{
  for (int element = 0; element < nelements; element++)
  {
    for (int ion = 0; ion < get_nions(element); ion++)
    {
      for (int i = 0; i < ION_COUNTER_COUNT; i++)
      {
        // energy or event count per volume per second
        const double ratedensity = get_ion_stats(mgi, element, ion, (enum ionstatscounters)i) / deltaV / deltat / nprocs;

        if (i < nstatcounters_ratecoeff)
        {
          // convert photon event counters into rate coefficients
          set_ion_stats(mgi, element, ion, (enum ionstatscounters)i, ratedensity / ionstagepop(mgi, element, ion));
        }
        else
        {
          set_ion_stats(mgi, element, ion, (enum ionstatscounters)i, ratedensity);
        }
      }
    }
  }
}
#endif


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

          const double bfheatingcoeff_ana = get_bfheatingcoeff_ana(element, ion, 0, 0, get_TR(n), get_W(n));
          bfheatingestimator[ionestimindex] = bfheatingestimator[ionestimindex] / bfheatingcoeff_ana;

          if (!isfinite(bfheatingestimator[ionestimindex]))
          {
            printout("[fatal] about to set bfheatingestimator = NaN = bfheatingestimator / get_bfheatingcoeff_ana(%d,%d,%d,%d,%d)=%g/%g",element,ion,0,0,n,bfheatingestimator[ionestimindex],bfheatingcoeff_ana);
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


static void update_grid_cell(const int mgi, const int nts, const int nts_prev, const int titer, const double tratmid,
                             const double deltat, double *mps, heatingcoolingrates_t *heatingcoolingrates)
// n is the modelgrid index
{
  const int assoc_cells = get_numassociatedcells(mgi);
  if (assoc_cells > 0)
  {
    // estimators were accumulated in nts_prev, but radiation density, etc should be scaled to the cell volume at nts
    const double deltaV = vol_init_modelcell(mgi) * pow(time_step[nts].mid / tmin, 3);
    const time_t sys_time_start_update_cell = time(NULL);
    // const bool log_this_cell = ((n % 50 == 0) || (npts_model < 50));
    const bool log_this_cell = true;

    /// Update current mass density of cell
    //n = nonemptycells[my_rank+ncl*nprocs];
    if (log_this_cell)
      printout("[info] update_grid: working on cell %d before timestep %d ...\n", mgi, nts);
    //n = nonemptycells[ncl];
    //printout("[debug] update_grid: ncl %d is %d non-empty cell updating grid cell %d ... T_e %g, rho %g\n",ncl,my_rank+ncl*nprocs,n,cell[n].T_e,cell[n].rho);
    modelgrid[mgi].rho = get_rhoinit(mgi) / pow(tratmid, 3);
    //cell[n].rho = cell[n].rho_init / pow(tratmid,3);
    //rho = cell[n].rho;
    /// This is done outside update grid now
    //modelgrid[n].totalcooling = COOLING_UNDEFINED;

    if (opacity_case == 4)
    {
      /// Update abundances of radioactive isotopes
      //printout("call update abundances for timestep %d in model cell %d\n",m,n);
      const time_t sys_time_start_update_abundances = time(NULL);
      update_abundances(mgi, nts, time_step[nts].mid);
      calculate_deposition_rate_density(mgi, nts);
      printout("update_abundances for cell %d timestep %d took %ld seconds\n", mgi, nts, time(NULL) - sys_time_start_update_abundances);

      /// For timestep 0 we calculate the level populations straight forward wihout
      /// applying any temperature correction
      if ((nts - itstep) == 0 && titer == 0)
      {
        /// Determine renormalisation factor for corrected photoionization cross-sections
        #ifndef FORCE_LTE
        #if (!NO_LUT_PHOTOION)
          if (!simulation_continued_from_saved)
          {
            set_all_corrphotoionrenorm(mgi, 1.);
          }
        #endif
        #endif

        /// W == 1 indicates that this modelgrid cell was treated grey in the
        /// last timestep. Therefore it has no valid Gamma estimators and must
        /// be treated in LTE at restart.
        if (modelgrid[mgi].thick == 0 && get_W(mgi) == 1)
        {
          if (log_this_cell)
            printout("force modelgrid cell to be grey at restart\n");
          modelgrid[mgi].thick = 1;
        }
        if (log_this_cell)
        {
          printout("initial_iteration %d\n",initial_iteration);
          printout("modelgrid.thick: %d\n",modelgrid[mgi].thick);
        }

        precalculate_partfuncts(mgi);
        //printout("abundance in cell %d is %g\n",n,cell[n].composition[0].abundance);

        if (!simulation_continued_from_saved || !NLTE_POPS_ON)
        {
          calculate_populations(mgi);  // these were not read from the gridsave file, so calculate them now
        }
        else
        {
          calculate_electron_densities(mgi);
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

        const double estimator_normfactor = 1 / deltaV / deltat / nprocs;
        const double estimator_normfactor_over4pi = ONEOVER4PI * estimator_normfactor;
        const time_t sys_time_start_temperature_corrections = time(NULL);

        radfield_normalise_J(mgi, estimator_normfactor_over4pi); // this applies normalisation to the fullspec J
        radfield_set_J_normfactor(mgi, estimator_normfactor_over4pi); // this stores the factor that will be applied later for the J bins but not fullspec J

        #ifdef DO_TITER
          radfield_titer_J(mgi);
        #endif

        #if (TRACK_ION_STATS)
        normalise_ion_estimators(mgi, deltat, deltaV);
        #endif

        #ifndef FORCE_LTE
        if (initial_iteration || modelgrid[mgi].thick == 1)
        #endif
        {
          const double T_R = get_T_R_from_J(mgi);
          set_TR(mgi, T_R);
          set_Te(mgi, T_R);
          set_TJ(mgi, T_R);
          set_W(mgi, 1);

          #ifndef FORCE_LTE
            #if (!NO_LUT_PHOTOION)
              set_all_corrphotoionrenorm(mgi, 1.);
            #endif
          #endif

          precalculate_partfuncts(mgi);
          calculate_populations(mgi);
        }
        #ifndef FORCE_LTE
        else // not (initial_iteration || modelgrid[n].thick == 1)
        {
          /// Calculate estimators

          radfield_normalise_nuJ(mgi, estimator_normfactor_over4pi);

          ffheatingestimator[mgi] *= estimator_normfactor;
          colheatingestimator[mgi] *= estimator_normfactor;

          #ifdef DO_TITER
            radfield_titer_nuJ(mgi);
            titer_average_estimators(mgi);
          #endif

          // Get radiation field parameters out of the full-spectrum and binned J and nuJ estimators
          radfield_fit_parameters(mgi, nts);

          #if (DETAILED_BF_ESTIMATORS_ON)
          radfield_normalise_bf_estimators(mgi, estimator_normfactor / H);
          #endif

          #if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
          update_gamma_corrphotoionrenorm_bfheating_estimators(mgi, estimator_normfactor);
          #endif

          solve_Te_nltepops(mgi, nts, titer, heatingcoolingrates);
        }
        #endif
        printout("Temperature/NLTE solution for cell %d timestep %d took %ld seconds\n", mgi, nts, time(NULL) - sys_time_start_temperature_corrections);
      }

      const float nne = get_nne(mgi);
      const double compton_optical_depth = SIGMA_T * nne * wid_init(mgi) * tratmid;

      double radial_pos = modelgrid[mgi].initial_radial_pos * tratmid / assoc_cells;
      if (grid_type == GRID_SPHERICAL1D)
      {
        const double r_inner = get_cellcoordmin(mgi, 0) * tratmid / assoc_cells;
        const double r_outer = r_inner + wid_init(mgi) * tratmid;
        radial_pos = 3./4 * (pow(r_outer, 4.) - pow(r_inner, 4.)) / (pow(r_outer, 3) - pow(r_inner, 3.)); // volume averaged mean radius
        // printout("r_inner %g r_outer %g tratmid %g assoc_cells %d\n", r_inner, r_outer, tratmid, assoc_cells);
      }
      const double grey_optical_deptha = get_kappagrey(mgi) * get_rho(mgi) * wid_init(mgi) * tratmid;
      const double grey_optical_depth = get_kappagrey(mgi) * get_rho(mgi) * (rmax * tratmid - radial_pos);
      if (log_this_cell)
      {
        printout("modelgridcell %d, compton optical depth (/propgridcell) %g, grey optical depth (/propgridcell) %g\n", mgi, compton_optical_depth, grey_optical_deptha);
        printout("radial_pos %g, distance_to_obs %g, tau_dist %g\n", radial_pos, rmax * tratmid - radial_pos, grey_optical_depth);
        //printout("rmax %g, tratmid %g\n",rmax,tratmid);
      }
      modelgrid[mgi].grey_depth = grey_optical_depth;

//          grey_optical_depth = compton_optical_depth;
      if ((grey_optical_depth > cell_is_optically_thick) && (nts < n_grey_timesteps))
      {
        printout("cell %d is treated in grey approximation (tau %g)\n", mgi, grey_optical_depth);
        modelgrid[mgi].thick = 1;
      }
      // else if (grey_optical_depth > cell_is_optically_thick_vpkt)
      //   modelgrid[n].thick = 2;
      else
        modelgrid[mgi].thick = 0;

      /// Cooling rates depend only on cell properties, precalculate total cooling
      /// and ion contributions inside update grid and communicate between MPI tasks
      const time_t sys_time_start_calc_kpkt_rates = time(NULL);

      // don't pass pointer to heatingcoolingrates because current populations and rates weren't used to determine T_e
      calculate_cooling_rates(mgi, NULL);

      printout("calculate_kpkt_rates for cell %d timestep %d took %ld seconds\n", mgi, nts, time(NULL) - sys_time_start_calc_kpkt_rates);
    }
    else if (opacity_case == 3)
    {
      /// MK Begin
      //printout("update_grid: opacity_case 3 ... updating cell[n].kappa_grey"); //MK
      if (get_rho(mgi) > rho_crit)
      {
        set_kappagrey(mgi, opcase3_normal * (0.9 * get_ffegrp(mgi) + 0.1) * rho_crit/get_rho(mgi));
      }
      else
      {
        set_kappagrey(mgi, opcase3_normal * (0.9 * get_ffegrp(mgi) + 0.1));
      }
      /// MK End
    }

    if (do_rlc_est == 2 && get_nne(mgi) > 0)
    {
      const double cell_len_scale_a = 0.1 / get_nne(mgi) / SIGMA_T;
      if (cell_len_scale_a < mps[tid])
      {
        mps[tid] = cell_len_scale_a;
      }
      const double cell_len_scale_b = 0.1 / get_rho(mgi) / GREY_OP;
      if (cell_len_scale_b <  mps[tid])
      {
        mps[tid] = cell_len_scale_b;
      }
    }
    printout("update_grid_cell for cell %d timestep %d took %ld seconds\n", mgi, nts, time(NULL) - sys_time_start_update_cell);
  }
  else
  {
    /// For modelgrid cells that are not represented in the simulation grid,
    /// Set grid properties to zero
    set_TR(mgi, 0.);
    set_TJ(mgi, 0.);
    set_Te(mgi, 0.);
    set_W(mgi, 0.);

    #ifndef FORCE_LTE
      #if (!NO_LUT_PHOTOION)
        zero_gammaestimator(mgi);
        set_all_corrphotoionrenorm(mgi, 0.);
      #endif
    #endif
  }
}


void update_grid(FILE *estimators_file, const int nts, const int nts_prev, const int my_rank, const int nstart, const int ndo, const int titer)
// Subroutine to update the matter quantities in the grid cells at the start
//   of the new timestep.
/// m timestep
{
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

  //printout("[debug] update_grid: time before initialisation of heating file %ld\n", time(NULL));
  //#ifndef FORCE_LTE
  //  sprintf(filename,"h%d-%d_heating_%.4d.out",m,titer,my_rank);
  //  heating_file = fopen_required(filename, "w")) == NULL);
  //  setvbuf(heating_file, NULL, _IOLBF, 1);
  //#endif
  //printout("[debug] update_grid: heating file initialised %ld\n", time(NULL));

  ///Calculate the critical opacity at which opacity_case 3 switches from a
  ///regime proportional to the density to a regime independent of the density
  ///This is done by solving for tau_sobolev == 1
  ///tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/nucmass(NUCLIDE_NI56) * 3000e-8 * time_step[m].mid;
  rho_crit = ME * CLIGHT * nucmass(NUCLIDE_NI56) / (PI * QE * QE * rho_crit_para * 3000e-8 * time_step[nts].mid);
  printout("update_grid: rho_crit = %g\n", rho_crit);

  // const double t_current = time_step[nts].start;

  // These values will not be used if nts == 0, but set them anyway
  // nts_prev is the previous timestep, unless this is timestep zero
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
        heatingcoolingrates_t heatingcoolingrates = {};
        update_grid_cell(mgi, nts, nts_prev, titer, tratmid, deltat, mps, &heatingcoolingrates);

        //maybe want to add omp ordered here if the modelgrid cells should be output in order
        const time_t sys_time_start_write_estimators = time(NULL);
        printout("writing to estimators file cell %d timestep %d...\n", mgi, nts);

        #ifdef _OPENMP
        #pragma omp critical(estimators_file)
        #endif
        {
          write_to_estimators_file(estimators_file, mgi, nts, titer, &heatingcoolingrates);
        }

        printout("writing estimators for cell %d timestep %d took %ld seconds\n",
                 mgi, nts, time(NULL) - sys_time_start_write_estimators);
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
  max_path_step = fmin(max_path_step, rmax / 10.);
  printout("max_path_step %g\n", max_path_step);

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
             (!NT_ON || ((rpkt_emiss[modelgridindex] == 0.) && (get_modelinitradioabund(modelgridindex, NUCLIDE_CR48) == 0.) && (get_modelinitradioabund(modelgridindex, NUCLIDE_NI56) == 0.))))
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
    //double nne_check = 0.;
    nne_tot = 0.;   /// total number of electrons in grid cell which are possible
                    /// targets for compton scattering of gamma rays

    nntot = nne;
    for (int element = 0; element < nelements; element++)
    {
      const double abundance = get_abundance(modelgridindex, element);
      const int nions = get_nions(element);
      /// calculate number density of the current element (abundances are given by mass)
      const double nnelement = abundance / elements[element].mass * get_rho(modelgridindex);
      nne_tot += nnelement * get_element(element);

      const int uppermost_ion = elements_uppermost_ion[tid][element];
      double ionfractions[uppermost_ion + 1];
      get_ionfractions(element, modelgridindex, nne, ionfractions, uppermost_ion);

      /// Use ionizationfractions to calculate the groundlevel populations
      for (int ion = 0; ion < nions; ion++)
      {
        double nnion;
        //if (ion <= elements[element].uppermost_ion)
        if (ion <= elements_uppermost_ion[tid][element])
        {
          if (abundance > 0)
          {
            nnion = nnelement * ionfractions[ion];
            if (nnion < MINPOP)
              nnion = MINPOP;
          }
          else
            nnion = 0.;
        }
        else
          nnion = MINPOP;  /// uppermost_ion is only < nions-1 in cells with nonzero abundance of the given species
        nntot += nnion;
        //nne_check += nnion * (get_ionstage(element,ion)-1);
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
      printout("the ion frac was %g, abundance %g and density %g\n",ionfractions[ion], abundance, get_rho(modelgridindex));
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


void write_grid_restart_data(const int timestep)
{
  char filename[100];
  sprintf(filename, "gridsave_ts%d.tmp", timestep);

  const time_t sys_time_start_write_restart = time(NULL);
  printout("Write grid restart data to %s...", filename);

  FILE *gridsave_file = fopen_required(filename, "w");

  fprintf(gridsave_file, "%d ", ntstep);

  for (int i = 0; i < ntstep; i++)
  {
    fprintf(gridsave_file, "%lg %lg ", time_step[i].gamma_dep, time_step[i].positron_dep);
  }

  fprintf(gridsave_file, "%d ", timestep);

  for (int mgi = 0; mgi < npts_model; mgi++)
  {
    const bool nonemptycell = (get_numassociatedcells(mgi) > 0);

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
  printout("done in %ld seconds.\n", time(NULL) - sys_time_start_write_restart);
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

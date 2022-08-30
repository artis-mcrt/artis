#include "update_grid.h"

#include <gsl/gsl_roots.h>
#include <stdlib.h>
#include <cmath>

#include "atomic.h"
#include "decay.h"
#include "grid.h"
#include "kpkt.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "sn3d.h"
#include "stats.h"
#include "thermalbalance.h"
#include "artisoptions.h"
#include "constants.h"
#include "globals.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"

void precalculate_partfuncts(int modelgridindex)
/// The partition functions depend only on T_R and W. This means they don't
/// change during any iteration on T_e. Therefore their precalculation was
/// taken out of calculate_populations to save runtime.
{
  /// Precalculate partition functions for each ion in every cell
  /// this saves a factor 10 in calculation time of Saha-Boltzman populations
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      // printout("precalc element %d, ion %d, mgi %d\n",element,ion,modelgridindex);
      // globals::cell[cellnumber].composition[element].ltepartfunct[ion] = calculate_ltepartfunct(element,ion,T_R);
      grid::modelgrid[modelgridindex].composition[element].partfunct[ion] =
          calculate_partfunct(element, ion, modelgridindex);
    }
  }
}

static void write_to_estimators_file(FILE *estimators_file, const int mgi, const int timestep, const int titer,
                                     const struct heatingcoolingrates *heatingcoolingrates) {
  // return; disable for better performance (if estimators files are not needed)
  if (grid::get_numassociatedcells(mgi) > 0) {
    const float T_e = grid::get_Te(mgi);
    const float nne = grid::get_nne(mgi);
    const double Y_e = grid::get_electronfrac(mgi);
    // fprintf(estimators_file,"%d %g %g %g %g %d
    // ",n,get_TR(n),grid::get_Te(n),get_W(n),get_TJ(n),grid::modelgrid[n].thick); fprintf(estimators_file,"%d %g %g %g
    // %g %g ",n,get_TR(n),grid::get_Te(n),get_W(n),get_TJ(n),grey_optical_depth);
    fprintf(estimators_file,
            "timestep %d modelgridindex %d titeration %d TR %g Te %g W %g TJ %g grey_depth %g thick %d nne %g Ye %g "
            "tdays %7.2f\n",
            timestep, mgi, titer, grid::get_TR(mgi), T_e, grid::get_W(mgi), grid::get_TJ(mgi),
            grid::modelgrid[mgi].grey_depth, grid::modelgrid[mgi].thick, nne, Y_e,
            globals::time_step[timestep].mid / DAY);
    // fprintf(estimators_file,"%d %g %g %g %g %g %g %g
    //",n,get_TR(n),grid::get_Te(n),get_W(n),get_TJ(n),grey_optical_depth,grey_optical_deptha,compton_optical_depth);

    if (NLTE_POPS_ON)  //  && timestep % 2 == 0
    {
      nltepop_write_to_file(mgi, timestep);
    }

    for (int element = 0; element < get_nelements(); element++) {
      if (grid::get_elem_abundance(mgi, element) <= 0.)  // skip elements with no abundance
        continue;

      fprintf(estimators_file, "populations        Z=%2d", get_element(element));
      const int nions = get_nions(element);
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      double elpop = 0.;
      for (int ion = 0; ion < nions; ion++) {
        elpop += ionstagepop(mgi, element, ion);
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ionstagepop(mgi, element, ion));
      }
      if (nions == 0) {
        elpop = grid::get_elem_numberdens(mgi, element);
      }
      fprintf(estimators_file, "  SUM: %9.3e", elpop);

      decay::fprint_nuc_abundances(estimators_file, mgi, globals::time_step[timestep].mid, element);

      if (nions == 0 || elpop <= 0.) {
        // dummy element for nuclear abundances only
        continue;
      }

      const bool printdebug = false;

      bool assume_lte = true;
      bool per_gmpop = true;
      const bool lower_superlevel_only = false;

      // if (timestep % 10 == 0)
      // {
      //   fprintf(estimators_file, "RRC_LTE_Nahar      Z=%2d", get_element(element));
      //   for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //     fprintf(estimators_file, "              ");
      //   for (int ion = 0; ion < nions; ion++)
      //   {
      //     fprintf(estimators_file, "  %d: %9.3e",
      //             get_ionstage(element, ion),
      //             calculate_ionrecombcoeff(-1, T_e, element, ion, assume_lte, false, printdebug,
      //             lower_superlevel_only, per_gmpop, false));
      //   }
      //   fprintf(estimators_file, "\n");
      // }

      per_gmpop = false;

      // fprintf(estimators_file, "AlphaLTE_R*nne Z=%2d", get_element(element));
      // for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //   fprintf(estimators_file, "              ");
      // for (int ion = 0; ion < nions; ion++)
      // {
      //   fprintf(estimators_file, "  %d: %9.3e",
      //           get_ionstage(element, ion),
      //           calculate_ionrecombcoeff(n, T_e, element, ion, assume_lte, false, printdebug, lower_superlevel_only,
      //           per_gmpop) * nne);
      // }
      // fprintf(estimators_file, "\n");

      assume_lte = false;

      if (TRACK_ION_STATS && TRACK_ION_MASTATS) {
        fprintf(estimators_file, "MA_IN_RADEXC       Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        double ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_RADEXC);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_OUT_RADEEXC     Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYOUT_RADDEEXC);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_IN_COLEXC       Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_COLLEXC);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_OUT_COLDEEXC    Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYOUT_COLLDEEXC);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_IN_PHOTOION     Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_PHOTOION);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_OUT_RADRECOMB   Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYOUT_RADRECOMB);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_IN_COLION       Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_COLLION);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_OUT_COLRECOMB   Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYOUT_COLLRECOMB);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_IN_NTCOLION     Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_NTCOLLION);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_IN_TOTAL        Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_in_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_TOTAL);
          ma_el += ma_in_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_in_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_OUT_TOTAL       Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYOUT_TOTAL);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_IN_INTERNAL     Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_INTERNAL);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_OUT_INTERNAL    Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYOUT_INTERNAL);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);
      }

      // // spontaneous radiative recombination rate coefficient (may or may not include stim. recomb)
      // fprintf(estimators_file, "AlphaR*nne         Z=%2d", get_element(element));
      // for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //   fprintf(estimators_file, "              ");
      // for (int ion = 0; ion < nions; ion++)
      // {
      //   const bool printdebug = false;
      //   // const bool printdebug = (get_element(element) >= 26);
      //
      //   fprintf(estimators_file, "  %d: %9.3e",
      //           get_ionstage(element, ion),
      //           calculate_ionrecombcoeff(mgi, T_e, element, ion, assume_lte, false, printdebug,
      //           lower_superlevel_only, per_gmpop, false) * nne);
      // }
      // fprintf(estimators_file, "\n");
      //
      // if (timestep % 10 == 0)
      // {
      //   fprintf(estimators_file, "AlphaR_toSL*nne    Z=%2d", get_element(element));
      //   for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //     fprintf(estimators_file, "              ");
      //   for (int ion = 0; ion < nions; ion++)
      //   {
      //     const bool printdebug = false;
      //     // const bool printdebug = (get_element(element) >= 26);
      //
      //     fprintf(estimators_file, "  %d: %9.3e",
      //             get_ionstage(element, ion),
      //             calculate_ionrecombcoeff(mgi, T_e, element, ion, assume_lte, false, printdebug, true, per_gmpop,
      //             false) * nne);
      //   }
      //   fprintf(estimators_file, "\n");
      // }
      //
      // #if (TRACK_ION_STATS)
      // fprintf(estimators_file, "AlphaR_MC_MA*nne   Z=%2d", get_element(element));
      // for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //   fprintf(estimators_file, "              ");
      // for (int ion = 0; ion < nions; ion++)
      // {
      //   const double alpha_r_mc = get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_MACROATOM);
      //   fprintf(estimators_file, "  %d: %9.3e",
      //           get_ionstage(element, ion),
      //           alpha_r_mc);
      // }
      // fprintf(estimators_file, "\n");
      // #endif
      //
      // #if (TRACK_ION_STATS)
      // fprintf(estimators_file, "AlphaR_MC_KPKT*nne Z=%2d", get_element(element));
      // for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //   fprintf(estimators_file, "              ");
      // for (int ion = 0; ion < nions; ion++)
      // {
      //   const double alpha_r_mc = get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_KPKT);
      //   fprintf(estimators_file, "  %d: %9.3e",
      //           get_ionstage(element, ion),
      //           alpha_r_mc);
      // }
      // fprintf(estimators_file, "\n");
      // #endif

      if (TRACK_ION_STATS) {
        fprintf(estimators_file, "AlphaR_MC*nne      Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions; ion++) {
          const double alpha_r_mc = get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_MACROATOM) +
                                    get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_KPKT);
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), alpha_r_mc);
        }
        fprintf(estimators_file, "\n");

        fprintf(estimators_file, "BF_escfrac         Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions; ion++) {
          const double alpha_r_mc = get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_MACROATOM) +
                                    get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_KPKT);
          const double alpha_r_mc_abs = get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_ABSORBED);
          fprintf(estimators_file, "  %d: %9.3f", get_ionstage(element, ion), 1. - alpha_r_mc_abs / alpha_r_mc);
        }
        fprintf(estimators_file, "\n");

        fprintf(estimators_file, "BF_escfrac_vpkt    Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions; ion++) {
          const double alpha_r_mc = get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_MACROATOM) +
                                    get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_KPKT);
          const double alpha_r_mc_escaped = get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_ESCAPED);
          fprintf(estimators_file, "  %d: %9.3f", get_ionstage(element, ion), alpha_r_mc_escaped / alpha_r_mc);
        }
        fprintf(estimators_file, "\n");

        fprintf(estimators_file, "BB_escfrac         Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions; ion++) {
          const double bb_emitted = get_ion_stats(mgi, element, ion, stats::ION_BOUNDBOUND_MACROATOM);
          const double bb_abs = get_ion_stats(mgi, element, ion, stats::ION_BOUNDBOUND_ABSORBED);
          fprintf(estimators_file, "  %d: %9.3f", get_ionstage(element, ion), 1. - bb_abs / bb_emitted);
        }
        fprintf(estimators_file, "\n");
      }

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
      //           calculate_ionrecombcoeff(mgi, T_e, element, ion, assume_lte, false, printdebug,
      //           lower_superlevel_only, per_gmpop, true) * nne);
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
      //           calculate_ionrecombcoeff(mgi, T_e, element, ion, assume_lte, true, printdebug, lower_superlevel_only,
      //           per_gmpop) * nne);
      // }
      // fprintf(estimators_file, "\n");

      // if (timestep % 20 == 0)
      // {
      //   fprintf(estimators_file, "kappa_bf(nuedge)   Z=%2d", get_element(element));
      //   for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //     fprintf(estimators_file, "              ");
      //   for (int ion = 0; ion < nions - 1; ion++)
      //   {
      //     double nu_edge = (epsilon(element, ion + 1, 0) - epsilon(element, ion, 0)) / H;
      //     double kappa_bf = calculate_kappa_bf_gammacontr(mgi, nu_edge);
      //
      //     fprintf(estimators_file, "  %d: %9.3e",
      //             get_ionstage(element, ion),
      //             kappa_bf);
      //   }
      //   fprintf(estimators_file, "\n");
      // }

      // {
      //   fprintf(estimators_file, "gamma_R            Z=%2d", get_element(element));
      //   for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //     fprintf(estimators_file, "              ");
      //   for (int ion = 0; ion < nions - 1; ion++)
      //   {
      //     // const bool printdebug_gammar = (get_element(element) == 26 && get_ionstage(element, ion) == 2);
      //     const bool printdebug_gammar = false;
      //     fprintf(estimators_file, "  %d: %9.3e",
      //             get_ionstage(element, ion),
      //             calculate_iongamma_per_ionpop(mgi, T_e, element, ion, assume_lte, false, printdebug_gammar, false,
      //             false));
      //   }
      //   fprintf(estimators_file, "\n");
      // }

      if (DETAILED_BF_ESTIMATORS_ON) {
        fprintf(estimators_file, "gamma_R_integral   Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions - 1; ion++) {
          // const bool printdebug_gammar = (get_element(element) == 26 && get_ionstage(element, ion) == 2);
          const bool printdebug_gammar = false;
          fprintf(
              estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
              calculate_iongamma_per_ionpop(mgi, T_e, element, ion, assume_lte, false, printdebug_gammar, false, true));
        }
        fprintf(estimators_file, "\n");
      }

      if (DETAILED_BF_ESTIMATORS_ON) {
        fprintf(estimators_file, "gamma_R_bfest      Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions - 1; ion++) {
          // const bool printdebug_gammar = ((get_element(element) == 26 || get_element(element) == 28) &&
          // get_ionstage(element, ion) >= 2); const bool printdebug_gammar = (get_element(element) >= 26);
          const bool printdebug_gammar = false;
          fprintf(
              estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
              calculate_iongamma_per_ionpop(mgi, T_e, element, ion, assume_lte, false, printdebug_gammar, true, false));
        }
        fprintf(estimators_file, "\n");
      }

#if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC         Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++) {
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, stats::ION_PHOTOION));
      }
      fprintf(estimators_file, "\n");
#endif

#if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BF      Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++) {
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBOUNDFREE));
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
      //           get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBFSAMEELEMENT]);
      // }
      // fprintf(estimators_file, "\n");

#if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BF_i+1  Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++) {
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBFIONPLUSONE));
      }
      fprintf(estimators_file, "\n");
#endif

#if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BF_i+2  Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++) {
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBFIONPLUSTWO));
      }
      fprintf(estimators_file, "\n");
#endif

#if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BF_i+3  Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++) {
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBFIONPLUSTHREE));
      }
      fprintf(estimators_file, "\n");
#endif

#if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BFtoSL  Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++) {
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBFLOWERSUPERLEVEL));
      }
      fprintf(estimators_file, "\n");
#endif

#if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BB      Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++) {
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUND));
      }
      fprintf(estimators_file, "\n");
#endif

#if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BB_i+1  Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++) {
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUNDIONPLUSONE));
      }
      fprintf(estimators_file, "\n");
#endif

#if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BB_i+2  Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++) {
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUNDIONPLUSTWO));
      }
      fprintf(estimators_file, "\n");
#endif

#if (TRACK_ION_STATS)
      fprintf(estimators_file, "gamma_R_MC_BB_i+3  Z=%2d", get_element(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
        fprintf(estimators_file, "              ");
      for (int ion = 0; ion < nions - 1; ion++) {
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUNDIONPLUSTHREE));
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

      if (NT_ON) {
        fprintf(estimators_file, "gamma_NT           Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions - 1; ion++) {
          const double Y_nt = nonthermal::nt_ionization_ratecoeff(mgi, element, ion);
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), Y_nt);
        }
        fprintf(estimators_file, "\n");

#if (TRACK_ION_STATS)
        fprintf(estimators_file, "gamma_NT_MC        Z=%2d", get_element(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
          fprintf(estimators_file, "              ");
        for (int ion = 0; ion < nions - 1; ion++) {
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                  get_ion_stats(mgi, element, ion, stats::ION_NTION));
        }
        fprintf(estimators_file, "\n");
#endif
      }

#if (!defined FORCE_LTE && !NO_LUT_PHOTOION)
      fprintf(estimators_file, "corrphotoionrenorm Z=%2d", get_element(element));
      for (int ion = 0; ion < nions; ion++) {
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                globals::corrphotoionrenorm[mgi * get_nelements() * get_max_nions() + element * get_max_nions() + ion]);
      }
      fprintf(estimators_file, "\n");
      fprintf(estimators_file, "gammaestimator     Z=%2d", get_element(element));
      for (int ion = 0; ion < nions; ion++) {
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                globals::gammaestimator[mgi * get_nelements() * get_max_nions() + element * get_max_nions() + ion]);
      }
      fprintf(estimators_file, "\n");
#endif
    }

#ifndef FORCE_LTE
    fprintf(estimators_file, "heating: ff %11.5e bf %11.5e coll %11.5e       dep %11.5e heating_dep/total_dep %.3f\n",
            heatingcoolingrates->heating_ff, heatingcoolingrates->heating_bf, heatingcoolingrates->heating_collisional,
            heatingcoolingrates->heating_dep, heatingcoolingrates->nt_frac_heating);
    fprintf(estimators_file, "cooling: ff %11.5e fb %11.5e coll %11.5e adiabatic %11.5e\n",
            heatingcoolingrates->cooling_ff, heatingcoolingrates->cooling_fb, heatingcoolingrates->cooling_collisional,
            heatingcoolingrates->cooling_adiabatic);
#endif
  } else {
    // modelgrid cells which are not represented in the simulation grid
    fprintf(estimators_file, "timestep %d modelgridindex %d EMPTYCELL\n", timestep, mgi);
  }
  fprintf(estimators_file, "\n");

  fflush(estimators_file);
}

__host__ __device__ void cellhistory_reset(const int modelgridindex, const bool new_timestep) {
  /// All entries of the cellhistory stack must be flagged as empty at the
  /// onset of the new timestep. Also, boundary crossing?
  /// Calculate the level populations for this cell, and flag the other entries
  /// as empty.
  /// Make known that globals::cellhistory[tid] contains information about the
  /// cell given by cellnumber. (-99 if invalid)
  if ((modelgridindex == globals::cellhistory[tid].cellnumber) && !new_timestep) {
    return;
  }

  const int tid = get_thread_num();

  // force rpkt opacities to be recalculated next time they are accessed
  globals::kappa_rpkt_cont[tid].recalculate_required = true;

  globals::cellhistory[tid].cellnumber = modelgridindex;
  // globals::cellhistory[tid].totalcooling = COOLING_UNDEFINED;
  //  int nlevels_with_processrates = 0;
  //  int nlevels_with_photoioncoeffs = 0;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      globals::cellhistory[tid].cooling_contrib[kpkt::get_coolinglistoffset(element, ion)] = COOLING_UNDEFINED;
      const int nlevels = get_nlevels(element, ion);
      for (int level = 0; level < nlevels; level++) {
        if (modelgridindex >= 0) {
          globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].population =
              calculate_levelpop(modelgridindex, element, ion, level);
        }
        for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++) {
          // bool levelhasvalue = false;
          // if
          // (globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].corrphotoioncoeff
          // >= 0)
          // {
          //   levelhasvalue = true;
          // }
          globals::cellhistory[tid]
              .chelements[element]
              .chions[ion]
              .chlevels[level]
              .chphixstargets[phixstargetindex]
              .corrphotoioncoeff = -99.;
#if (SEPARATE_STIMRECOMB)
          globals::cellhistory[tid]
              .chelements[element]
              .chions[ion]
              .chlevels[level]
              .chphixstargets[phixstargetindex]
              .stimrecombcoeff = -99.;
#endif
        }
        /// This is the only flag needed for all of the following MA stuff!
        // if
        // (globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].processrates[MA_ACTION_COLDEEXC]
        // >= 0)
        //   nlevels_with_processrates++;

        globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].processrates[MA_ACTION_COLDEEXC] =
            -99.;

        // globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_deexc = -99.;
        // globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_recomb = -99.;
        // globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_recomb = -99.;
        // globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_same = -99.;
        // globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_same = -99.;
        // globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_lower = -99.;
        // globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_higher = -99.;
        //
        // ndowntrans = get_ndowntrans(element, ion, level);
        // nuptrans = get_nuptrans(element, ion, level);
        // for (i = 0; i < ndowntrans; i++)
        // {
        //   globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_rad_deexc[i] = -99.;
        //   globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_down_same[i] =
        //   -99.;
        // }
        // for (i = 0; i < nuptrans; i++)
        // {
        //   globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_up_same[i] =
        //   -99.;
        // }
        // if (levelhasvalue)
        //   nlevels_with_photoioncoeffs++;
      }
    }
  }
  // printout("nlevels_with_processrates %d\n", nlevels_with_processrates);
  // printout("nlevels_with_photoioncoeffs %d\n", nlevels_with_photoioncoeffs);

  // globals::cellhistory[tid].totalcooling = COOLING_UNDEFINED;
  // globals::cellhistory[tid].phixsflag = PHIXS_UNDEFINED;
}

static void solve_Te_nltepops(const int n, const int nts, const int titer,
                              struct heatingcoolingrates *heatingcoolingrates)
// n is the modelgridindex (TODO: rename to mgi)
// nts is the timestep number
{
  // bfheating coefficients are needed for the T_e solver, but
  // they only depend on the radiation field, which is fixed during the iterations below
  printout("calculate_bfheatingcoeffs for timestep %d cell %d...", nts, n);
  const time_t sys_time_start_calculate_bfheatingcoeffs = time(NULL);
  calculate_bfheatingcoeffs(n);
  printout("took %ld seconds\n", time(NULL) - sys_time_start_calculate_bfheatingcoeffs);

  const double covergence_tolerance = 0.04;
  for (int nlte_iter = 0; nlte_iter <= NLTEITER; nlte_iter++) {
    const time_t sys_time_start_spencerfano = time(NULL);
    if (NT_ON && NT_SOLVE_SPENCERFANO) {
      nonthermal::solve_spencerfano(n, nts, nlte_iter);  // depends on the ionization balance, and weakly on nne
    }
    const int duration_solve_spencerfano = time(NULL) - sys_time_start_spencerfano;

    const time_t sys_time_start_partfuncs_or_gamma = time(NULL);
    if (!NLTE_POPS_ON) {
      precalculate_partfuncts(n);
    }
#if (!NO_LUT_PHOTOION)
    else if ((nlte_iter != 0)) {
      // recalculate the Gammas using the current population estimates
      for (int element = 0; element < get_nelements(); element++) {
        const int nions = get_nions(element);
        for (int ion = 0; ion < nions - 1; ion++) {
          globals::gammaestimator[n * get_nelements() * get_max_nions() + element * get_max_nions() + ion] =
              calculate_iongamma_per_gspop(n, element, ion);
        }
      }
    }
#endif
    const int duration_solve_partfuncs_or_gamma = time(NULL) - sys_time_start_partfuncs_or_gamma;

    /// Find T_e as solution for thermal balance
    const double prev_T_e = grid::get_Te(n);
    const time_t sys_time_start_Te = time(NULL);
    const int nts_for_te = (titer == 0) ? nts - 1 : nts;

    call_T_e_finder(n, nts, globals::time_step[nts_for_te].mid, MINTEMP, MAXTEMP, heatingcoolingrates);

    const int duration_solve_T_e = time(NULL) - sys_time_start_Te;
    const double fracdiff_T_e = fabs((grid::get_Te(n) / prev_T_e) - 1);

    if (!NLTE_POPS_ON || !NLTE_POPS_ALL_IONS_SIMULTANEOUS)  // do this in LTE or NLTE single ion solver mode
    {
      /// Store population values to the grid
      const time_t sys_time_start_pops = time(NULL);
      calculate_populations(n);
      const int duration_solve_pops = time(NULL) - sys_time_start_pops;
      // calculate_cooling_rates(n);
      // calculate_heating_rates(n);
      printout(
          "Grid solver cell %d timestep %d: time spent on: Spencer-Fano %ds, partfuncs/gamma %ds, T_e %ds, populations "
          "%ds\n",
          n, nts, duration_solve_spencerfano, duration_solve_partfuncs_or_gamma, duration_solve_T_e,
          duration_solve_pops);
    }

    if (NLTE_POPS_ON) {
      const time_t sys_time_start_nltepops = time(NULL);
      // fractional difference between previous and current iteration's (nne or max(ground state population change))
      double nlte_test = 0.;
      if (NLTE_POPS_ALL_IONS_SIMULTANEOUS) {
        for (int element = 0; element < get_nelements(); element++) {
          if (get_nions(element) > 0) {
            solve_nlte_pops_element(element, n, nts, nlte_iter);
          }
        }
      } else {
        for (int element = 0; element < get_nelements(); element++) {
          const int nions = get_nions(element);
          for (int ion = 0; ion < nions - 1; ion++) {
            const double trial = fabs(solve_nlte_pops_ion(element, ion, n, nts) - 1);

            if (trial > nlte_test) {
              nlte_test = trial;
            }
          }
        }
      }
      const int duration_solve_nltepops = time(NULL) - sys_time_start_nltepops;

      if (NLTE_POPS_ALL_IONS_SIMULTANEOUS) {
        const double nne_prev = grid::get_nne(n);
        precalculate_partfuncts(n);
        calculate_electron_densities(n);  // sets nne
        const double fracdiff_nne = fabs((grid::get_nne(n) / nne_prev) - 1);
        nlte_test = fracdiff_nne;
        printout(
            "NLTE solver cell %d timestep %d iteration %d: time spent on: Spencer-Fano %ds, T_e %ds, NLTE populations "
            "%ds\n",
            n, nts, nlte_iter, duration_solve_spencerfano, duration_solve_T_e, duration_solve_nltepops);
        printout(
            "NLTE (Spencer-Fano/Te/pops) solver cell %d timestep %d iteration %d: prev_iter nne %g, new nne is %g, "
            "fracdiff %g, prev T_e %g new T_e %g fracdiff %g\n",
            n, nts, nlte_iter, nne_prev, grid::get_nne(n), nlte_test, prev_T_e, grid::get_Te(n), fracdiff_T_e);
        // damp changes in nne if oscillating to much
        // grid::set_nne(n, (grid::get_nne(n) + nne_prev) / 2.);
      } else {
        printout(
            "Completed iteration for NLTE population solver in cell %d for timestep %d. Fractional error returned: "
            "%g\n",
            n, nts, nlte_test);
      }

      if (nlte_test <= covergence_tolerance && fracdiff_T_e <= covergence_tolerance) {
        printout(
            "NLTE (Spencer-Fano/Te/pops) solver nne converged to tolerance %g <= %g and T_e to tolerance %g <= %g "
            "after %d iterations.\n",
            nlte_test, covergence_tolerance, fracdiff_T_e, covergence_tolerance, nlte_iter + 1);
        break;
      } else if (nlte_iter == NLTEITER) {
        printout("WARNING: NLTE solver failed to converge after %d iterations. Keeping solution from last iteration\n",
                 nlte_iter + 1);
      }
    } else
      break;  // no iteration is needed without NLTE_POPS_ON
  }
}

#if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
static void update_gamma_corrphotoionrenorm_bfheating_estimators(const int n, const double estimator_normfactor) {
#if (!NO_LUT_PHOTOION)
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions - 1; ion++) {
      const int ionestimindex = n * get_nelements() * get_max_nions() + element * get_max_nions() + ion;
      // printout("mgi %d, element %d, ion %d, gammaest %g\n",n,element,ion,globals::gammaestimator[ionestimindex]);
      globals::gammaestimator[ionestimindex] *= estimator_normfactor / H;
// printout("mgi %d, element %d, ion %d, gammaest %g\n",n,element,ion,globals::gammaestimator[ionestimindex]);
#ifdef DO_TITER
      if (globals::gammaestimator_save[ionestimindex] >= 0) {
        globals::gammaestimator[ionestimindex] =
            (globals::gammaestimator[ionestimindex] + globals::gammaestimator_save[ionestimindex]) / 2;
      }
      globals::gammaestimator_save[ionestimindex] = globals::gammaestimator[ionestimindex];
#endif

      globals::corrphotoionrenorm[ionestimindex] =
          globals::gammaestimator[ionestimindex] / get_corrphotoioncoeff_ana(element, ion, 0, 0, n);

      if (!std::isfinite(globals::corrphotoionrenorm[ionestimindex])) {
        printout(
            "[fatal] about to set corrphotoionrenorm = NaN = gammaestimator / "
            "get_corrphotoioncoeff_ana(%d,%d,%d,%d,%d)=%g/%g",
            element, ion, 0, 0, n, globals::gammaestimator[ionestimindex],
            get_corrphotoioncoeff_ana(element, ion, 0, 0, n));
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
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions - 1; ion++) {
      /// Reuse the gammaestimator array as temporary storage of the Gamma values during
      /// the remaining part of the update_grid phase. Afterwards it is reset to record
      /// the next timesteps gamma estimators.
      // nlevels = get_nlevels(element,ion);
      // nlevels = get_ionisinglevels(element,ion);
      const int ionestimindex = n * get_nelements() * get_max_nions() + element * get_max_nions() + ion;

#if (!NO_LUT_PHOTOION)
      globals::gammaestimator[ionestimindex] = calculate_iongamma_per_gspop(n, element, ion);
      // printout("mgi %d, element %d, ion %d, Gamma %g\n",n,element,ion,Gamma);
#endif

#if (!NO_LUT_BFHEATING)
      globals::bfheatingestimator[ionestimindex] *= estimator_normfactor;
#ifdef DO_TITER
      if (bfheatingestimator_save[ionestimindex] >= 0) {
        globals::bfheatingestimator[ionestimindex] =
            (globals::bfheatingestimator[ionestimindex] + bfheatingestimator_save[ionestimindex]) / 2;
      }
      bfheatingestimator_save[ionestimindex] = globals::bfheatingestimator[ionestimindex];
#endif
      /// Now convert bfheatingestimator into the bfheating renormalisation coefficient used in get_bfheating
      /// in the remaining part of update_grid. Later on it's reset and new contributions are added up.

      const double bfheatingcoeff_ana = get_bfheatingcoeff_ana(element, ion, 0, 0, grid::get_TR(n), grid::get_W(n));
      globals::bfheatingestimator[ionestimindex] = globals::bfheatingestimator[ionestimindex] / bfheatingcoeff_ana;

      if (!std::isfinite(globals::bfheatingestimator[ionestimindex])) {
        printout(
            "[fatal] about to set bfheatingestimator = NaN = bfheatingestimator / "
            "get_bfheatingcoeff_ana(%d,%d,%d,%d,%d)=%g/%g",
            element, ion, 0, 0, n, globals::bfheatingestimator[ionestimindex], bfheatingcoeff_ana);
        abort();
      }

      // printout("cell %d element %d ion %d bfheatingestimator %g\n",n,element,ion,bfheatingestimator[ionestimindex]);
#endif
    }
  }
#endif
}
#endif

#ifdef DO_TITER
#ifndef FORCE_LTE
static void titer_average_estimators(const int n) {
  if (globals::ffheatingestimator_save[n] >= 0) {
    globals::ffheatingestimator[n] = (globals::ffheatingestimator[n] + globals::ffheatingestimator_save[n]) / 2;
  }
  globals::ffheatingestimator_save[n] = globals::ffheatingestimator[n];
  if (globals::colheatingestimator_save[n] >= 0) {
    globals::colheatingestimator[n] = (globals::colheatingestimator[n] + globals::colheatingestimator_save[n]) / 2;
  }
  globals::colheatingestimator_save[n] = globals::colheatingestimator[n];
}
#endif
#endif

#if (!NO_LUT_PHOTOION)
static void zero_gammaestimator(const int modelgridindex) {
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++)
      globals::gammaestimator[modelgridindex * get_nelements() * get_max_nions() + element * get_max_nions() + ion] =
          0.;
  }
}

static void set_all_corrphotoionrenorm(const int modelgridindex, const double value) {
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      globals::corrphotoionrenorm[modelgridindex * get_nelements() * get_max_nions() + element * get_max_nions() +
                                  ion] = value;
    }
  }
}
#endif

static void update_grid_cell(const int mgi, const int nts, const int nts_prev, const int titer, const double tratmid,
                             const double deltat, double *mps, struct heatingcoolingrates *heatingcoolingrates)
// n is the modelgrid index
{
  const int assoc_cells = grid::get_numassociatedcells(mgi);
  if (assoc_cells > 0) {
    // estimators were accumulated in nts_prev, but radiation density, etc should be scaled to the cell volume at nts
    const double deltaV = grid::vol_init_modelcell(mgi) * pow(globals::time_step[nts].mid / globals::tmin, 3);
    const time_t sys_time_start_update_cell = time(NULL);
    // const bool log_this_cell = ((n % 50 == 0) || (npts_model < 50));
    const bool log_this_cell = true;

    /// Update current mass density of cell
    // n = nonemptycells[my_rank+ncl*nprocs];
    if (log_this_cell)
      printout("[info] update_grid: working on cell %d before timestep %d titeration %d...\n", mgi, nts, titer);
    // n = nonemptycells[ncl];
    // printout("[debug] update_grid: ncl %d is %d non-empty cell updating grid cell %d ... T_e %g, rho
    // %g\n",ncl,my_rank+ncl*nprocs,n,globals::cell[n].T_e,globals::cell[n].rho);
    grid::modelgrid[mgi].rho = grid::get_rhoinit(mgi) / pow(tratmid, 3);
    // globals::cell[n].rho = globals::cell[n].rho_init / pow(tratmid,3);
    // rho = globals::cell[n].rho;
    /// This is done outside update grid now
    // grid::modelgrid[n].totalcooling = COOLING_UNDEFINED;

    /// Update abundances of radioactive isotopes
    decay::update_abundances(mgi, nts, globals::time_step[nts].mid);
    const double estimator_normfactor = 1 / deltaV / deltat / globals::nprocs;
    const double estimator_normfactor_over4pi = ONEOVER4PI * estimator_normfactor;

    if (globals::opacity_case >= 4) {
      /// For timestep 0 we calculate the level populations straight forward wihout
      /// applying any temperature correction
      if ((nts - globals::itstep) == 0 && titer == 0) {
/// Determine renormalisation factor for corrected photoionization cross-sections
#ifndef FORCE_LTE
#if (!NO_LUT_PHOTOION)
        if (!globals::simulation_continued_from_saved) {
          set_all_corrphotoionrenorm(mgi, 1.);
        }
#endif
#endif

        /// W == 1 indicates that this modelgrid cell was treated grey in the
        /// last timestep. Therefore it has no valid Gamma estimators and must
        /// be treated in LTE at restart.
        if (grid::modelgrid[mgi].thick == 0 && grid::get_W(mgi) == 1) {
          if (log_this_cell) {
            printout("force modelgrid cell %d to be grey at restart\n", mgi);
          }
          grid::modelgrid[mgi].thick = 1;
        }
        if (log_this_cell) {
          printout("initial_iteration %d\n", globals::initial_iteration);
          printout("mgi %d modelgrid.thick: %d\n", mgi, grid::modelgrid[mgi].thick);
        }

        precalculate_partfuncts(mgi);
        // printout("abundance in cell %d is %g\n",n,globals::cell[n].composition[0].abundance);

        if (!globals::simulation_continued_from_saved || !NLTE_POPS_ON) {
          calculate_populations(mgi);  // these were not read from the gridsave file, so calculate them now
        } else {
          calculate_electron_densities(mgi);
          // printout("nne: %g\n", grid::get_nne(n));
        }
      } else
      /// For all other timesteps temperature corrections have to be applied
      {
        /// we have to calculate the electron density
        /// and all the level populations
        /// Normalise estimators and make sure that they are finite.
        /// Then update T_R and W using the estimators.
        /// (This could in principle also be done for empty cells)

        const time_t sys_time_start_temperature_corrections = time(NULL);

        radfield::normalise_J(mgi, estimator_normfactor_over4pi);       // this applies normalisation to the fullspec J
        radfield::set_J_normfactor(mgi, estimator_normfactor_over4pi);  // this stores the factor that will be applied
                                                                        // later for the J bins but not fullspec J

#ifdef DO_TITER
        radfield::titer_J(mgi);
#endif

        if (TRACK_ION_STATS) {
          stats::normalise_ion_estimators(mgi, deltat, deltaV);
        }

#ifndef FORCE_LTE
        // initial_iteration really means either ts 0 or nts < globals::n_lte_timesteps
        if (globals::initial_iteration || grid::modelgrid[mgi].thick == 1)
#endif
        {
          // LTE mode or grey mode (where temperature doesn't matter but is calculated anyway)

          const double T_J = radfield::get_T_J_from_J(mgi);
          grid::set_TR(mgi, T_J);
          grid::set_Te(mgi, T_J);
          grid::set_TJ(mgi, T_J);
          grid::set_W(mgi, 1);

#ifndef FORCE_LTE
#if (!NO_LUT_PHOTOION)
          set_all_corrphotoionrenorm(mgi, 1.);
#endif
#endif

          precalculate_partfuncts(mgi);
          calculate_populations(mgi);
        }
#ifndef FORCE_LTE
        else  // not (initial_iteration || grid::modelgrid[n].thick == 1)
        {
          // non-LTE timesteps with T_e from heating/cooling

          radfield::normalise_nuJ(mgi, estimator_normfactor_over4pi);

          globals::ffheatingestimator[mgi] *= estimator_normfactor;
          globals::colheatingestimator[mgi] *= estimator_normfactor;

#ifdef DO_TITER
          radfield::titer_nuJ(mgi);
          titer_average_estimators(mgi);
#endif

#if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
          update_gamma_corrphotoionrenorm_bfheating_estimators(mgi, estimator_normfactor);
#endif

          // Get radiation field parameters (T_J, T_R, W, and bins if enabled) out of the full-spectrum and binned J and
          // nuJ estimators
          radfield::fit_parameters(mgi, nts);

#if (DETAILED_BF_ESTIMATORS_ON)
          radfield::normalise_bf_estimators(mgi, estimator_normfactor / H);
#endif

          solve_Te_nltepops(mgi, nts, titer, heatingcoolingrates);
        }
#endif
        printout("Temperature/NLTE solution for cell %d timestep %d took %ld seconds\n", mgi, nts,
                 time(NULL) - sys_time_start_temperature_corrections);
      }

      const float nne = grid::get_nne(mgi);
      const double compton_optical_depth = SIGMA_T * nne * grid::wid_init(mgi) * tratmid;

      double radial_pos = grid::modelgrid[mgi].initial_radial_pos_sum * tratmid / assoc_cells;
      if (grid::grid_type == GRID_SPHERICAL1D) {
        const double r_inner = grid::get_cellcoordmin(mgi, 0) * tratmid;
        const double r_outer = r_inner + grid::wid_init(mgi) * tratmid;
        radial_pos = 3. / 4 * (pow(r_outer, 4.) - pow(r_inner, 4.)) /
                     (pow(r_outer, 3) - pow(r_inner, 3.));  // volume averaged mean radius
        // printout("r_inner %g r_outer %g tratmid %g assoc_cells %d\n", r_inner, r_outer, tratmid, assoc_cells);
      }
      const double grey_optical_deptha = grid::get_kappagrey(mgi) * grid::get_rho(mgi) * grid::wid_init(mgi) * tratmid;
      const double grey_optical_depth =
          grid::get_kappagrey(mgi) * grid::get_rho(mgi) * (globals::rmax * tratmid - radial_pos);
      if (log_this_cell) {
        printout("modelgridcell %d, compton optical depth (/propgridcell) %g, grey optical depth (/propgridcell) %g\n",
                 mgi, compton_optical_depth, grey_optical_deptha);
        printout("radial_pos %g, distance_to_obs %g, tau_dist %g\n", radial_pos, globals::rmax * tratmid - radial_pos,
                 grey_optical_depth);
        // printout("rmax %g, tratmid %g\n",rmax,tratmid);
      }
      grid::modelgrid[mgi].grey_depth = grey_optical_depth;

      // grey_optical_depth = compton_optical_depth;

      if ((grey_optical_depth > globals::cell_is_optically_thick) && (nts < globals::n_grey_timesteps)) {
        printout("cell %d is treated in grey approximation (kappa_grey %g [cm2/g], tau %g)\n", mgi,
                 grid::get_kappagrey(mgi), grey_optical_depth);
        grid::modelgrid[mgi].thick = 1;
      }
#ifdef VPKT_ON
      else if (grey_optical_depth > cell_is_optically_thick_vpkt) {
        grid::modelgrid[mgi].thick = 2;
      }
#endif
      else {
        grid::modelgrid[mgi].thick = 0;
      }

      /// Cooling rates depend only on cell properties, precalculate total cooling
      /// and ion contributions inside update grid and communicate between MPI tasks
      const time_t sys_time_start_calc_kpkt_rates = time(NULL);

      // don't pass pointer to heatingcoolingrates because current populations and rates weren't used to determine T_e
      kpkt::calculate_cooling_rates(mgi, NULL);

      printout("calculate_kpkt_rates for cell %d timestep %d took %ld seconds\n", mgi, nts,
               time(NULL) - sys_time_start_calc_kpkt_rates);
    } else {
      // For opacity_case != 4 the opacity treatment is grey. Enforce
      // optically thick treatment in this case (should be equivalent to grey)
      grid::modelgrid[mgi].thick = 1;

      /// Need the total number density of bound and free electrons for Compton scattering
      calculate_electron_densities(mgi);  // if this causes problems, disable the nne calculation (only need nne_tot)

      if (!((nts - globals::itstep) == 0 && titer == 0)) {
        radfield::normalise_J(mgi, estimator_normfactor_over4pi);       // this applies normalisation to the fullspec J
        radfield::set_J_normfactor(mgi, estimator_normfactor_over4pi);  // this stores the factor that will be applied
                                                                        // later for the J bins but not fullspec J

        const double T_J = radfield::get_T_J_from_J(mgi);
        grid::set_TR(mgi, T_J);
        grid::set_Te(mgi, T_J);
        grid::set_TJ(mgi, T_J);
        grid::set_W(mgi, 1);
      }

      if (globals::opacity_case == 3) {
        // printout("update_grid: opacity_case 3 ... updating globals::cell[n].kappa_grey"); //MK
        if (grid::get_rho(mgi) > globals::rho_crit) {
          grid::set_kappagrey(mgi, globals::opcase3_normal * (0.9 * grid::get_ffegrp(mgi) + 0.1) * globals::rho_crit /
                                       grid::get_rho(mgi));
        } else {
          grid::set_kappagrey(mgi, globals::opcase3_normal * (0.9 * grid::get_ffegrp(mgi) + 0.1));
        }
      }
    }

    if (globals::do_rlc_est == 2 && grid::get_nne(mgi) > 0) {
      const double cell_len_scale_a = 0.1 / grid::get_nne(mgi) / SIGMA_T;
      if (cell_len_scale_a < mps[tid]) {
        mps[tid] = cell_len_scale_a;
      }
      const double cell_len_scale_b = 0.1 / grid::get_rho(mgi) / GREY_OP;
      if (cell_len_scale_b < mps[tid]) {
        mps[tid] = cell_len_scale_b;
      }
    }
    printout("update_grid_cell for cell %d timestep %d took %ld seconds\n", mgi, nts,
             time(NULL) - sys_time_start_update_cell);
  } else {
    /// For modelgrid cells that are not represented in the simulation grid,
    /// Set grid properties to zero
    grid::set_TR(mgi, 0.);
    grid::set_TJ(mgi, 0.);
    grid::set_Te(mgi, 0.);
    grid::set_W(mgi, 0.);

#ifndef FORCE_LTE
#if (!NO_LUT_PHOTOION)
    zero_gammaestimator(mgi);
    set_all_corrphotoionrenorm(mgi, 0.);
#endif
#endif
  }
}

void update_grid(FILE *estimators_file, const int nts, const int nts_prev, const int my_rank, const int nstart,
                 const int ndo, const int titer, const time_t real_time_start)
// Subroutine to update the matter quantities in the grid cells at the start
//   of the new timestep.
/// nts timestep
{
  const time_t sys_time_start_update_grid = time(NULL);
  printout("\ntimestep %d: time before update grid %ld (tstart + %ld) simtime ts_mid %g days\n", nts,
           sys_time_start_update_grid, sys_time_start_update_grid - real_time_start, globals::time_step[nts].mid / DAY);

#ifndef FORCE_LTE
#if (!NO_LUT_PHOTOION)
  /// Initialise globals::corrphotoionrenorm[i] to zero before update_grid is called
  /// unless they have been read from file
  if ((!globals::simulation_continued_from_saved) || (nts - globals::itstep != 0) || (titer != 0)) {
    printout("nts %d, titer %d: reset corr photoionrenorm\n", nts, titer);
    for (int i = 0; i < grid::get_npts_model() * get_nelements() * get_max_nions(); i++) {
      globals::corrphotoionrenorm[i] = 0.;
    }
    printout("after nts %d, titer %d: reset corr photoionrenorm\n", nts, titer);
  }
#endif
#endif

  // printout("[debug] update_grid: starting update for timestep %d...\n",m);
  const double tratmid = globals::time_step[nts].mid / globals::tmin;

  double mps[get_max_threads()];  /// Thread private substitution of max_path_step. Its minimum is
                                  /// assigned to max_path_step after the parallel update_grid finished.
  for (int i = 0; i < get_max_threads(); i++) {
    mps[i] = 1.e35;
  }

  /// Calculate the critical opacity at which opacity_case 3 switches from a
  /// regime proportional to the density to a regime independent of the density
  /// This is done by solving for tau_sobolev == 1
  /// tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/nucmass(28, 56) * 3000e-8 * globals::time_step[m].mid;
  globals::rho_crit = ME * CLIGHT * decay::nucmass(28, 56) /
                      (PI * QE * QE * globals::rho_crit_para * 3000e-8 * globals::time_step[nts].mid);
  printout("update_grid: rho_crit = %g\n", globals::rho_crit);

  // These values will not be used if nts == 0, but set them anyway
  // nts_prev is the previous timestep, unless this is timestep zero
  const double deltat = globals::time_step[nts_prev].width;

  // printout("timestep %d, titer %d\n", nts, titer);
  // printout("deltat %g\n", deltat);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    /// Do not use values which are saved in the cellhistory within update_grid
    /// and daughter routines (THREADPRIVATE VARIABLE, THEREFORE HERE!)
    use_cellhist = false;
    cellhistory_reset(-99, true);

/// Updating cell information
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif

    for (int mgi = 0; mgi < grid::get_npts_model(); mgi++) {
      /// Check if this task should work on the current model grid cell.
      /// If yes, update the cell and write out the estimators
      if (mgi >= nstart && mgi < nstart + ndo) {
        struct heatingcoolingrates heatingcoolingrates = {};
        update_grid_cell(mgi, nts, nts_prev, titer, tratmid, deltat, mps, &heatingcoolingrates);

        // maybe want to add omp ordered here if the modelgrid cells should be output in order
        const time_t sys_time_start_write_estimators = time(NULL);
        printout("writing to estimators file cell %d timestep %d...\n", mgi, nts);

#ifdef _OPENMP
#pragma omp critical(estimators_file)
#endif
        { write_to_estimators_file(estimators_file, mgi, nts, titer, &heatingcoolingrates); }

        printout("writing estimators for cell %d timestep %d took %ld seconds\n", mgi, nts,
                 time(NULL) - sys_time_start_write_estimators);
      } else {
/// else, only reset gammaestimator to zero. This allows us to do a global MPI
/// communication after update_grid to synchronize gammaestimator
/// and write a contiguous restart file with grid properties
#if (!NO_LUT_PHOTOION)
        zero_gammaestimator(mgi);
#endif
      }
    }  /// end parallel for loop over all modelgrid cells

    /// Now after all the relevant taks of update_grid have been finished activate
    /// the use of the cellhistory for all OpenMP tasks, in what follows (update_packets)
    use_cellhist = true;
  }  /// end OpenMP parallel section

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
    mgi = grid::get_cell_modelgridindex(n);
    check1 = check1 + (grid::get_kappagrey(mgi)  * grid::get_rho(mgi));
    check2 = check2 + grid::get_rho(mgi);
  }
  printout("Grey normalisation check: %g\n", check1/check2);
  */

  /// Assign the minimum of thread private mps to the global variable max_path_step
  globals::max_path_step = mps[0];
  for (int i = 1; i < get_max_threads(); i++) {
    if (mps[i] < globals::max_path_step) {
      globals::max_path_step = mps[i];
    }
  }

  if (globals::do_rlc_est == 2) {
    const double trat = globals::time_step[nts].start / globals::tmin;
    if (globals::max_path_step < (grid::wid_init(0) * trat / 10.)) {
      globals::max_path_step = grid::wid_init(0) / 10. * trat;
    }
  }
  globals::max_path_step = fmin(globals::max_path_step, globals::rmax / 10.);
  printout("max_path_step %g\n", globals::max_path_step);

  const time_t time_update_grid_end_thisrank = time(NULL);
  printout("finished update grid on this rank at time %ld\n", time_update_grid_end_thisrank);

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  printout("timestep %d: time after update grid for all processes %ld (rank %d took %lds, waited %lds, total %lds)\n",
           nts, time(NULL), my_rank, time_update_grid_end_thisrank - sys_time_start_update_grid,
           time(NULL) - time_update_grid_end_thisrank, time(NULL) - sys_time_start_update_grid);
}

double calculate_populations(const int modelgridindex)
/// Determines the electron number density for a given cell using one of
/// libgsl's root_solvers and calculates the depending level populations.
{
  /// and the solution function
  gsl_function f;
  struct nne_solution_paras paras;
  paras.cellnumber = modelgridindex;
  f.function = &nne_solution_f;
  f.params = &paras;

  neutral_flag = false;

  /// Get temperatures
  const double T_R = grid::get_TR(modelgridindex);
  const float T_e = grid::get_Te(modelgridindex);
  const double W = grid::get_W(modelgridindex);

  double nne_hi = grid::get_rho(modelgridindex) / MH;

  /// The following section of uppermost_ion is (so far) NOT thread safe!!!!!!!!!!!!!!!!!!!!!!!
  int only_neutral_ions = 0;
  int nelements_in_cell = 0;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    // elements[element].uppermost_ion = nions-1;
    grid::set_elements_uppermost_ion(modelgridindex, element, nions - 1);
    const double abundance = grid::get_elem_abundance(modelgridindex, element);
    if (abundance > 0) {
      int uppermost_ion;
#ifdef FORCE_LTE
      uppermost_ion = get_nions(element) - 1;
#else
      if (globals::initial_iteration || grid::modelgrid[modelgridindex].thick == 1) {
        uppermost_ion = get_nions(element) - 1;
      } else {
        int ion;
        for (ion = 0; ion < nions - 1; ion++) {
// printout("element %d, ion %d, photoionest
// %g\n",element,ion,photoionestimator[modelgridindex*get_nelements()*get_max_nions()+element*get_max_nions()+ion]); if
// (photoionestimator[modelgridindex*get_nelements()*get_max_nions()+element*get_max_nions()+ion] == 0) break;
#if NO_LUT_PHOTOION
          const double Gamma = calculate_iongamma_per_gspop(modelgridindex, element, ion);
#else
          const double Gamma = globals::gammaestimator[modelgridindex * get_nelements() * get_max_nions() +
                                                       element * get_max_nions() + ion];
#endif

          if ((Gamma == 0) && (!NT_ON || ((globals::rpkt_emiss[modelgridindex] == 0.) &&
                                          (grid::get_modelinitradioabund(modelgridindex, 24, 48) == 0.) &&
                                          (grid::get_modelinitradioabund(modelgridindex, 28, 56) == 0.))))
            break;
        }
        uppermost_ion = ion;
      }
#endif
      // printout("cell %d, element %d, uppermost_ion by Gamma is %d\n",cellnumber,element,uppermost_ion);

      double factor = 1.;
      int ion;
      for (ion = 0; ion < uppermost_ion; ion++) {
        factor *= nne_hi * phi(element, ion, modelgridindex);
        // printout("element %d, ion %d, factor %g\n",element,i,factor);
        if (!std::isfinite(factor)) {
          printout(
              "[info] calculate_populations: uppermost_ion limited by phi factors for element Z=%d, ionstage %d in "
              "cell %d\n",
              get_element(element), get_ionstage(element, ion), modelgridindex);
          break;
        }
      }
      uppermost_ion = ion;
      // printout("cell %d, element %d, final uppermost_ion is %d, factor
      // %g\n",modelgridindex,element,uppermost_ion,factor); elements[element].uppermost_ion = uppermost_ion;
      grid::set_elements_uppermost_ion(modelgridindex, element, uppermost_ion);
      if (uppermost_ion == 0) only_neutral_ions++;
      nelements_in_cell++;
    }
  }

  float nne = 0.;
  double nne_tot = 0.;  /// total number of electrons in grid cell which are possible
                        /// targets for compton scattering of gamma rays
  double nntot = 0.;
  if (only_neutral_ions == nelements_in_cell) {
    /// Special case of only neutral ions, set nne to some finite value that
    /// packets are not lost in kpkts
    /// Introduce a flag variable which is sent to the T_e solver so that
    /// we get this info only once when T_e is converged and not for each
    /// iteration step.
    neutral_flag = true;
    // printout("[warning] calculate_populations: only neutral ions in cell %d modelgridindex\n",modelgridindex);
    // abort();
    /// Now calculate the ground level populations in nebular approximation and store them to the grid
    for (int element = 0; element < get_nelements(); element++) {
      /// calculate number density of the current element (abundances are given by mass)
      const double nnelement = grid::get_elem_numberdens(modelgridindex, element);
      nne_tot += nnelement * get_element(element);

      const int nions = get_nions(element);
      /// Assign the species population to the neutral ion and set higher ions to MINPOP
      for (int ion = 0; ion < nions; ion++) {
        double nnion;
        if (ion == 0) {
          nnion = nnelement;
        } else if (nnelement > 0.) {
          nnion = MINPOP;
        } else {
          nnion = 0.;
        }
        nntot += nnion;
        nne += nnion * (get_ionstage(element, ion) - 1);
        grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion] =
            (nnion * stat_weight(element, ion, 0) /
             grid::modelgrid[modelgridindex].composition[element].partfunct[ion]);

        if (!std::isfinite(grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion])) {
          printout("[warning] calculate_populations: groundlevelpop infinite in connection with MINPOP\n");
        }
      }
    }
    nntot += nne;
    if (nne < MINPOP) {
      nne = MINPOP;
    }
    grid::set_nne(modelgridindex, nne);
  } else {
    /// Apply solver to get nne
    /// Search solution for nne in [nne_lo,nne_hi]
    // printout("nne@x_lo %g\n", nne_solution_f(nne_lo,f.params));
    // printout("nne@x_hi %g\n", nne_solution_f(nne_hi,f.params));
    // printout("n, x_lo, x_hi, T_R, T_e, W, rho %d, %g, %g, %g, %g, %g,
    // %g\n",modelgridindex,x_lo,x_hi,T_R,T_e,W,globals::cell[modelgridindex].rho);
    double nne_lo = 0.;  // MINPOP;
    if (nne_solution_f(nne_lo, f.params) * nne_solution_f(nne_hi, f.params) > 0) {
      printout("n, nne_lo, nne_hi, T_R, T_e, W, rho %d, %g, %g, %g, %g, %g, %g\n", modelgridindex, nne_lo, nne_hi, T_R,
               T_e, W, grid::get_rho(modelgridindex));
      printout("nne@x_lo %g\n", nne_solution_f(nne_lo, f.params));
      printout("nne@x_hi %g\n", nne_solution_f(nne_hi, f.params));
#ifndef FORCE_LTE
      for (int element = 0; element < get_nelements(); element++) {
        printout("cell %d, element %d, uppermost_ion is %d\n", modelgridindex, element,
                 grid::get_elements_uppermost_ion(modelgridindex, element));

#if (!NO_LUT_PHOTOION)
        for (int ion = 0; ion <= grid::get_elements_uppermost_ion(modelgridindex, element); ion++) {
          // printout("element %d, ion %d, photoionest
          // %g\n",element,ion,photoionestimator[modelgridindex*get_nelements()*get_max_nions()+element*get_max_nions()+ion]);
          printout("element %d, ion %d, gammaionest %g\n", element, ion,
                   globals::gammaestimator[modelgridindex * get_nelements() * get_max_nions() +
                                           element * get_max_nions() + ion]);
        }
#endif
      }
#endif
    }
    gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

    gsl_root_fsolver_set(solver, &f, nne_lo, nne_hi);
    int iter = 0;
    const int maxit = 100;
    const double fractional_accuracy = 1e-3;
    int status;
    do {
      iter++;
      gsl_root_fsolver_iterate(solver);
      nne = gsl_root_fsolver_root(solver);
      nne_lo = gsl_root_fsolver_x_lower(solver);
      nne_hi = gsl_root_fsolver_x_upper(solver);
      status = gsl_root_test_interval(nne_lo, nne_hi, 0, fractional_accuracy);
    } while (status == GSL_CONTINUE && iter < maxit);

    gsl_root_fsolver_free(solver);

    if (nne < MINPOP) {
      nne = MINPOP;
    }

    grid::set_nne(modelgridindex, nne);
    if (status == GSL_CONTINUE)
      printout("[warning] calculate_populations: nne did not converge within %d iterations\n", maxit);
    // printout("[debug] update_grid:   status = %s\n",gsl_strerror(status));
    // printout("[debug] update_grid:   converged nne %g\n",globals::cell[modelgridindex].nne);

    /// Now calculate the ground level populations in nebular approximation and store them to the grid
    // double nne_check = 0.;
    nne_tot = 0.;  /// total number of electrons in grid cell which are possible
                   /// targets for compton scattering of gamma rays

    nntot = nne;
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      /// calculate number density of the current element (abundances are given by mass)
      const double nnelement = grid::get_elem_numberdens(modelgridindex, element);
      nne_tot += nnelement * get_element(element);

      const int uppermost_ion = grid::get_elements_uppermost_ion(modelgridindex, element);
      double ionfractions[uppermost_ion + 1];

      if (nnelement > 0) {
        get_ionfractions(element, modelgridindex, nne, ionfractions, uppermost_ion);
      }

      /// Use ionizationfractions to calculate the groundlevel populations
      for (int ion = 0; ion < nions; ion++) {
        double nnion;
        if (ion <= uppermost_ion) {
          if (nnelement > 0) {
            nnion = nnelement * ionfractions[ion];
            if (nnion < MINPOP) nnion = MINPOP;
          } else {
            nnion = 0.;
          }
        } else {
          nnion = MINPOP;  /// uppermost_ion is only < nions-1 in cells with nonzero abundance of the given species
        }
        nntot += nnion;

        grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion] =
            (nnion * stat_weight(element, ion, 0) /
             grid::modelgrid[modelgridindex].composition[element].partfunct[ion]);

        if (!std::isfinite(grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion])) {
          printout("[warning] calculate_populations: groundlevelpop infinite in connection with MINPOP\n");
        }
      }
    }
  }

  grid::set_nnetot(modelgridindex, nne_tot);
  return nntot;
}

double calculate_electron_densities(const int modelgridindex)
// Determines the free and total electron number densities
// for a given cell and stores them, assuming ion populations (ground_level_pop and partfunc)
// are fixed (determined by NLTE all-ion solver)
{
  double nne_tot = 0.;  // total electron density
  float nne = 0.;       // free electron density

  for (int element = 0; element < get_nelements(); element++) {
    // calculate number density of the current element (abundances are given by mass)
    const double nnelement = grid::get_elem_numberdens(modelgridindex, element);
    nne_tot += nnelement * get_element(element);

    // Use ionization fractions to calculate the free electron contributions
    if (nnelement > 0) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++) {
        // if (ion <= globals::elements[element].uppermost_ion)
        nne += (get_ionstage(element, ion) - 1) * ionstagepop(modelgridindex, element, ion);
      }
    }
  }

  grid::set_nne(modelgridindex, nne);
  grid::set_nnetot(modelgridindex, nne_tot);
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

/*static int get_cell(double x, double y, double z, double t)
/// subroutine to identify the cell index from a position and a time.
{
  // Original version of this was general but very slow. Modifying to be
  // faster but only work for regular grid.

  double trat = t / globals::tmin;
  int nx = (x - (globals::cell[0].pos_min[0] * trat))/(grid::wid_init * trat);
  int ny = (y - (globals::cell[0].pos_min[1] * trat))/(grid::wid_init * trat);
  int nz = (z - (globals::cell[0].pos_min[2] * trat))/(grid::wid_init * trat);

  int n = nx + (grid::ncoordgrid[0] * ny) + (grid::ncoordgrid[0] * grid::ncoordgrid[1] * nz);

  // do a check

  if (x < globals::cell[n].pos_min[0] * trat)
  {
    printout("Problem with get_cell (1).\n");
    abort();
  }
  if (x > (globals::cell[n].pos_min[0]+grid::wid_init) * trat)
  {
    printout("Problem with get_cell (2).\n");
    abort();
  }
  if (y < globals::cell[n].pos_min[1] * trat)
  {
    printout("Problem with get_cell (3).\n");
    abort();
  }
  if (y > (globals::cell[n].pos_min[1]+grid::wid_init) * trat)
  {
    printout("Problem with get_cell (4).\n");
    abort();
  }
  if (z < globals::cell[n].pos_min[2] * trat)
  {
    printout("Problem with get_cell (5).\n");
    abort();
  }
  if (z > (globals::cell[n].pos_min[2]+grid::wid_init) * trat)
  {
    printout("Problem with get_cell (6).\n");
    abort();
  }
  return(n);


  // OLD
        //   trat = t / globals::tmin;
        //
        //   for (n = 0; n < ngrid; n++)
        //   {
        //   if (
        //   (x > globals::cell[n].pos_min[0] * trat) &&
        //   (x < (globals::cell[n].pos_min[0] + grid::wid_init) *trat) &&
        //   (y > globals::cell[n].pos_min[1] * trat) &&
        //   (y < (globals::cell[n].pos_min[1] + grid::wid_init) *trat) &&
        //   (z > globals::cell[n].pos_min[2] * trat) &&
        //   (z < (globals::cell[n].pos_min[2] + grid::wid_init) *trat))
        //   {
        //   return(n);
        // }
        // }
  // END OLD

  printout("Failed to find cell (get_cell). \n");
  printout("x %g, y %g, z %g, t %g\n", x, y, z, t);
  printout("xend %g yend %g zend %g\n", globals::cell[ngrid-1].pos_min[0] * trat, globals::cell[ngrid-1].pos_min[1] *
trat,globals::cell[ngrid-1].pos_min[2] * trat); printout("xend %g yend %g zend %g\n", globals::cell[0].pos_min[0] *
trat, globals::cell[0].pos_min[1] * trat,globals::cell[0].pos_min[2] * trat); printout("xend %g yend %g zend %g\n",
(globals::cell[0].pos_min[0]+grid::wid_init) * trat, (globals::cell[0].pos_min[1]+grid::wid_init) *
trat,(globals::cell[0].pos_min[2]+grid::wid_init) * trat); printout("xwid %g ywid %g zwid %g\n", (grid::wid_init) *
trat, (grid::wid_init) * trat,(grid::wid_init) * trat);

  abort();
}*/

#include "update_grid.h"

#include <cmath>

#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "grid.h"
#include "kpkt.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "rpkt.h"
#include "sn3d.h"
#include "stats.h"
#include "thermalbalance.h"
#include "vpkt.h"

static void write_to_estimators_file(FILE *estimators_file, const int mgi, const int timestep, const int titer,
                                     const struct heatingcoolingrates *heatingcoolingrates) {
  // return; disable for better performance (if estimators files are not needed)
  const time_t sys_time_start_write_estimators = time(nullptr);

  if (grid::get_numassociatedcells(mgi) > 0) {
    printout("writing to estimators file timestep %d cell %d...\n", timestep, mgi);

    const auto T_e = grid::get_Te(mgi);
    const auto nne = grid::get_nne(mgi);
    const auto Y_e = grid::get_electronfrac(mgi);
    // fprintf(estimators_file,"%d %g %g %g %g %d
    // ",n,get_TR(n),grid::get_Te(n),get_W(n),get_TJ(n),grid::modelgrid[n].thick); fprintf(estimators_file,"%d %g %g %g
    // %g %g ",n,get_TR(n),grid::get_Te(n),get_W(n),get_TJ(n),grey_optical_depth);
    fprintf(estimators_file,
            "timestep %d modelgridindex %d titeration %d TR %g Te %g W %g TJ %g grey_depth %g thick %d nne %g Ye %g "
            "tdays %7.2f\n",
            timestep, mgi, titer, grid::get_TR(mgi), T_e, grid::get_W(mgi), grid::get_TJ(mgi),
            grid::modelgrid[mgi].grey_depth, grid::modelgrid[mgi].thick, nne, Y_e,
            globals::timesteps[timestep].mid / DAY);
    // fprintf(estimators_file,"%d %g %g %g %g %g %g %g
    //",n,get_TR(n),grid::get_Te(n),get_W(n),get_TJ(n),grey_optical_depth,grey_optical_deptha,compton_optical_depth);

    if (globals::total_nlte_levels > 0) {
      nltepop_write_to_file(mgi, timestep);
    }

    for (int element = 0; element < get_nelements(); element++) {
      if (grid::get_elem_abundance(mgi, element) <= 0.) {  // skip elements with no abundance
        continue;
      }

      fprintf(estimators_file, "populations        Z=%2d", get_atomicnumber(element));
      const int nions = get_nions(element);
      if (nions > 0) {
        // add spaces for missing lowest ion stages to match other elements
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
      }
      double elpop = 0.;
      for (int ion = 0; ion < nions; ion++) {
        elpop += get_nnion(mgi, element, ion);
        fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), get_nnion(mgi, element, ion));
      }
      if (nions == 0) {
        elpop = grid::get_elem_numberdens(mgi, element);
      }
      fprintf(estimators_file, "  SUM: %9.3e", elpop);

      decay::fprint_nuc_abundances(estimators_file, mgi, globals::timesteps[timestep].mid, element);

      if (nions == 0 || elpop <= 0.) {
        // dummy element for nuclear abundances only
        continue;
      }

      // const bool printdebug = false;

      bool assume_lte = true;
      // bool per_gmpop = true;
      // const bool lower_superlevel_only = false;

      // if (timestep % 10 == 0)
      // {
      //   fprintf(estimators_file, "RRC_LTE_Nahar      Z=%2d", get_atomicnumber(element));
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

      // per_gmpop = false;

      // fprintf(estimators_file, "AlphaLTE_R*nne Z=%2d", get_atomicnumber(element));
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

      if constexpr (TRACK_ION_STATS && TRACK_ION_MASTATS) {
        fprintf(estimators_file, "MA_IN_RADEXC       Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        double ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_RADEXC);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_OUT_RADEEXC     Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYOUT_RADDEEXC);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_IN_COLEXC       Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_COLLEXC);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_OUT_COLDEEXC    Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYOUT_COLLDEEXC);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_IN_PHOTOION     Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_PHOTOION);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_OUT_RADRECOMB   Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYOUT_RADRECOMB);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_IN_COLION       Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_COLLION);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_OUT_COLRECOMB   Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYOUT_COLLRECOMB);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_IN_NTCOLION     Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_NTCOLLION);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_IN_TOTAL        Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_in_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_TOTAL);
          ma_el += ma_in_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_in_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_OUT_TOTAL       Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYOUT_TOTAL);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_IN_INTERNAL     Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_INTERNAL);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);

        fprintf(estimators_file, "MA_OUT_INTERNAL    Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        ma_el = 0.;
        for (int ion = 0; ion < nions; ion++) {
          const double ma_ion = get_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYOUT_INTERNAL);
          ma_el += ma_ion;
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), ma_ion);
        }
        fprintf(estimators_file, "  SUM: %9.3e\n", ma_el);
      }

      // // spontaneous radiative recombination rate coefficient (may or may not include stim. recomb)
      // fprintf(estimators_file, "AlphaR*nne         Z=%2d", get_atomicnumber(element));
      // for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //   fprintf(estimators_file, "              ");
      // for (int ion = 0; ion < nions; ion++)
      // {
      //   const bool printdebug = false;
      //   // const bool printdebug = (get_atomicnumber(element) >= 26);
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
      //   fprintf(estimators_file, "AlphaR_toSL*nne    Z=%2d", get_atomicnumber(element));
      //   for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //     fprintf(estimators_file, "              ");
      //   for (int ion = 0; ion < nions; ion++)
      //   {
      //     const bool printdebug = false;
      //     // const bool printdebug = (get_atomicnumber(element) >= 26);
      //
      //     fprintf(estimators_file, "  %d: %9.3e",
      //             get_ionstage(element, ion),
      //             calculate_ionrecombcoeff(mgi, T_e, element, ion, assume_lte, false, printdebug, true, per_gmpop,
      //             false) * nne);
      //   }
      //   fprintf(estimators_file, "\n");
      // }
      //
      // if constexpr (TRACK_ION_STATS) {
      //   fprintf(estimators_file, "AlphaR_MC_MA*nne   Z=%2d", get_atomicnumber(element));
      //   for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //     fprintf(estimators_file, "              ");
      //   for (int ion = 0; ion < nions; ion++) {
      //     const double alpha_r_mc = get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_MACROATOM);
      //     fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), alpha_r_mc);
      //   }
      //   fprintf(estimators_file, "\n");
      // }
      //
      // if constexpr (TRACK_ION_STATS) {
      //   fprintf(estimators_file, "AlphaR_MC_KPKT*nne Z=%2d", get_atomicnumber(element));
      //   for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //     fprintf(estimators_file, "              ");
      //   for (int ion = 0; ion < nions; ion++) {
      //     const double alpha_r_mc = get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_KPKT);
      //     fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), alpha_r_mc);
      //   }
      //   fprintf(estimators_file, "\n");
      // }

      if constexpr (TRACK_ION_STATS) {
        fprintf(estimators_file, "AlphaR_MC*nne      Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions; ion++) {
          const double alpha_r_mc = get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_MACROATOM) +
                                    get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_KPKT);
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), alpha_r_mc);
        }
        fprintf(estimators_file, "\n");

        fprintf(estimators_file, "BF_escfrac         Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions; ion++) {
          const double alpha_r_mc = get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_MACROATOM) +
                                    get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_KPKT);
          const double alpha_r_mc_abs = get_ion_stats(mgi, element, ion, stats::ION_RADRECOMB_ABSORBED);
          fprintf(estimators_file, "  %d: %9.3f", get_ionstage(element, ion), 1. - alpha_r_mc_abs / alpha_r_mc);
        }
        fprintf(estimators_file, "\n");

        fprintf(estimators_file, "BB_escfrac         Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions; ion++) {
          const double bb_emitted = get_ion_stats(mgi, element, ion, stats::ION_BOUNDBOUND_MACROATOM);
          const double bb_abs = get_ion_stats(mgi, element, ion, stats::ION_BOUNDBOUND_ABSORBED);
          fprintf(estimators_file, "  %d: %9.3f", get_ionstage(element, ion), 1. - bb_abs / bb_emitted);
        }
        fprintf(estimators_file, "\n");
      }

      // stimulated recombination rate coefficient
      // fprintf(estimators_file, "Alpha_stim*nne   Z=%2d", get_atomicnumber(element));
      // for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //   fprintf(estimators_file, "              ");
      // for (int ion = 0; ion < nions; ion++)
      // {
      //   // const bool printdebug = false;
      //   const bool printdebug = (get_atomicnumber(element) >= 26);
      //
      //   fprintf(estimators_file, "  %d: %9.3e",
      //           get_ionstage(element, ion),
      //           calculate_ionrecombcoeff(mgi, T_e, element, ion, assume_lte, false, printdebug,
      //           lower_superlevel_only, per_gmpop, true) * nne);
      // }
      // fprintf(estimators_file, "\n");

      // thermal collisional recombination rate coefficient
      // fprintf(estimators_file, "Alpha_C*nne      Z=%2d", get_atomicnumber(element));
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
      //   fprintf(estimators_file, "chi_bf(nuedge)   Z=%2d", get_atomicnumber(element));
      //   for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //     fprintf(estimators_file, "              ");
      //   for (int ion = 0; ion < nions - 1; ion++)
      //   {
      //     double nu_edge = (epsilon(element, ion + 1, 0) - epsilon(element, ion, 0)) / H;
      //     double chi_bf = calculate_chi_bf_gammacontr(mgi, nu_edge, false);
      //
      //     fprintf(estimators_file, "  %d: %9.3e",
      //             get_ionstage(element, ion),
      //             chi_bf);
      //   }
      //   fprintf(estimators_file, "\n");
      // }

      // {
      //   fprintf(estimators_file, "gamma_R            Z=%2d", get_atomicnumber(element));
      //   for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //     fprintf(estimators_file, "              ");
      //   for (int ion = 0; ion < nions - 1; ion++)
      //   {
      //     // const bool printdebug_gammar = (get_atomicnumber(element) == 26 && get_ionstage(element, ion) == 2);
      //     const bool printdebug_gammar = false;
      //     fprintf(estimators_file, "  %d: %9.3e",
      //             get_ionstage(element, ion),
      //             calculate_iongamma_per_ionpop(mgi, T_e, element, ion, assume_lte, false, printdebug_gammar, false,
      //             false));
      //   }
      //   fprintf(estimators_file, "\n");
      // }

      if (DETAILED_BF_ESTIMATORS_ON) {
        fprintf(estimators_file, "gamma_R_integral   Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions - 1; ion++) {
          // const bool printdebug_gammar = (get_atomicnumber(element) == 26 && get_ionstage(element, ion) == 2);
          const bool printdebug_gammar = false;
          fprintf(
              estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
              calculate_iongamma_per_ionpop(mgi, T_e, element, ion, assume_lte, false, printdebug_gammar, false, true));
        }
        fprintf(estimators_file, "\n");
      }

      if (DETAILED_BF_ESTIMATORS_ON) {
        fprintf(estimators_file, "gamma_R_bfest      Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions - 1; ion++) {
          // const bool printdebug_gammar = ((get_atomicnumber(element) == 26 || get_atomicnumber(element) == 28) &&
          // get_ionstage(element, ion) >= 2); const bool printdebug_gammar = (get_atomicnumber(element) >= 26);
          const bool printdebug_gammar = false;
          fprintf(
              estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
              calculate_iongamma_per_ionpop(mgi, T_e, element, ion, assume_lte, false, printdebug_gammar, true, false));
        }
        fprintf(estimators_file, "\n");
      }

      if constexpr (TRACK_ION_STATS) {
        fprintf(estimators_file, "gamma_R_MC         Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions - 1; ion++) {
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                  get_ion_stats(mgi, element, ion, stats::ION_PHOTOION));
        }
        fprintf(estimators_file, "\n");
      }

      if constexpr (TRACK_ION_STATS) {
        fprintf(estimators_file, "gamma_R_MC_BF      Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions - 1; ion++) {
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                  get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBOUNDFREE));
        }
        fprintf(estimators_file, "\n");
      }

      // fprintf(estimators_file, "gamma_R_MC_BFsameZ Z=%2d", get_atomicnumber(element));
      // for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //   fprintf(estimators_file, "              ");
      // for (int ion = 0; ion < nions; ion++)
      // {
      //   fprintf(estimators_file, "  %d: %9.3e",
      //           get_ionstage(element, ion),
      //           get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBFSAMEELEMENT]);
      // }
      // fprintf(estimators_file, "\n");

      if constexpr (TRACK_ION_STATS) {
        fprintf(estimators_file, "gamma_R_MC_BF_i+1  Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions - 1; ion++) {
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                  get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBFIONPLUSONE));
        }
        fprintf(estimators_file, "\n");
      }

      if constexpr (TRACK_ION_STATS) {
        fprintf(estimators_file, "gamma_R_MC_BF_i+2  Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions - 1; ion++) {
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                  get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBFIONPLUSTWO));
        }
        fprintf(estimators_file, "\n");
      }

      if constexpr (TRACK_ION_STATS) {
        fprintf(estimators_file, "gamma_R_MC_BF_i+3  Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions - 1; ion++) {
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                  get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBFIONPLUSTHREE));
        }
        fprintf(estimators_file, "\n");
      }

      if constexpr (TRACK_ION_STATS) {
        fprintf(estimators_file, "gamma_R_MC_BFtoSL  Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions - 1; ion++) {
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                  get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBFLOWERSUPERLEVEL));
        }
        fprintf(estimators_file, "\n");
      }

      if constexpr (TRACK_ION_STATS) {
        fprintf(estimators_file, "gamma_R_MC_BB      Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions - 1; ion++) {
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                  get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUND));
        }
        fprintf(estimators_file, "\n");
      }

      if constexpr (TRACK_ION_STATS) {
        fprintf(estimators_file, "gamma_R_MC_BB_i+1  Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions - 1; ion++) {
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                  get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUNDIONPLUSONE));
        }
        fprintf(estimators_file, "\n");
      }

      if constexpr (TRACK_ION_STATS) {
        fprintf(estimators_file, "gamma_R_MC_BB_i+2  Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions - 1; ion++) {
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                  get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUNDIONPLUSTWO));
        }
        fprintf(estimators_file, "\n");
      }

      if constexpr (TRACK_ION_STATS) {
        fprintf(estimators_file, "gamma_R_MC_BB_i+3  Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions - 1; ion++) {
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                  get_ion_stats(mgi, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUNDIONPLUSTHREE));
        }
        fprintf(estimators_file, "\n");
      }

      // fprintf(estimators_file, "gamma_C          Z=%2d", get_atomicnumber(element));
      // for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++)
      //   fprintf(estimators_file, "              ");
      // for (int ion = 0; ion < nions - 1; ion++)
      // {
      //   fprintf(estimators_file, "  %d: %9.3e",
      //           get_ionstage(element, ion),
      //           calculate_iongamma_per_ionpop(mgi, T_e, element, ion, assume_lte, true,
      //           printdebug));
      // }
      // fprintf(estimators_file, "\n");

      if (NT_ON) {
        fprintf(estimators_file, "gamma_NT           Z=%2d", get_atomicnumber(element));
        for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
          fprintf(estimators_file, "              ");
        }
        for (int ion = 0; ion < nions - 1; ion++) {
          const double Y_nt = nonthermal::nt_ionization_ratecoeff(mgi, element, ion);
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion), Y_nt);
        }
        fprintf(estimators_file, "\n");

        if constexpr (TRACK_ION_STATS) {
          fprintf(estimators_file, "gamma_NT_MC        Z=%2d", get_atomicnumber(element));
          for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
            fprintf(estimators_file, "              ");
          }
          for (int ion = 0; ion < nions - 1; ion++) {
            fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                    get_ion_stats(mgi, element, ion, stats::ION_NTION));
          }
          fprintf(estimators_file, "\n");
        }
      }

      if (USE_LUT_PHOTOION && globals::nbfcontinua > 0) {
        fprintf(estimators_file, "corrphotoionrenorm Z=%2d", get_atomicnumber(element));
        for (int ion = 0; ion < nions; ion++) {
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                  globals::corrphotoionrenorm[get_ionestimindex(mgi, element, ion)]);
        }
        fprintf(estimators_file, "\n");
        fprintf(estimators_file, "gammaestimator     Z=%2d", get_atomicnumber(element));
        for (int ion = 0; ion < nions; ion++) {
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                  globals::gammaestimator[get_ionestimindex(mgi, element, ion)]);
        }
        fprintf(estimators_file, "\n");
      }
    }

    fprintf(estimators_file, "heating: ff %11.5e bf %11.5e coll %11.5e       dep %11.5e heating_dep/total_dep %.3f\n",
            heatingcoolingrates->heating_ff, heatingcoolingrates->heating_bf, heatingcoolingrates->heating_collisional,
            heatingcoolingrates->heating_dep, heatingcoolingrates->nt_frac_heating);
    fprintf(estimators_file, "cooling: ff %11.5e fb %11.5e coll %11.5e adiabatic %11.5e\n",
            heatingcoolingrates->cooling_ff, heatingcoolingrates->cooling_fb, heatingcoolingrates->cooling_collisional,
            heatingcoolingrates->cooling_adiabatic);

  } else {
    // modelgrid cells which are not represented in the simulation grid
    fprintf(estimators_file, "timestep %d modelgridindex %d EMPTYCELL\n", timestep, mgi);
  }
  fprintf(estimators_file, "\n");

  fflush(estimators_file);

  const int write_estim_duration = time(nullptr) - sys_time_start_write_estimators;
  if (write_estim_duration >= 1) {
    printout("writing estimators for timestep %d cell %d took %d seconds\n", timestep, mgi, write_estim_duration);
  }
}

void cellhistory_reset(const int modelgridindex, const bool new_timestep) {
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
  globals::chi_rpkt_cont[tid].recalculate_required = true;

  globals::cellhistory[tid].cellnumber = modelgridindex;

  //  int nlevels_with_processrates = 0;
  // const double T_e = modelgridindex >= 0 ? grid ::get_Te(modelgridindex) : 0.;
  const int nelements = get_nelements();
  for (int element = 0; element < nelements; element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      globals::cellhistory[tid].cooling_contrib[kpkt::get_coolinglistoffset(element, ion)] = COOLING_UNDEFINED;

      if (modelgridindex >= 0) {
        const int nlevels = get_nlevels(element, ion);
        for (int level = 0; level < nlevels; level++) {
          globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].population =
              calculate_levelpop(modelgridindex, element, ion, level);
        }
      }
    }

    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_nlevels(element, ion);
      for (int level = 0; level < nlevels; level++) {
        for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++) {
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

        // globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_deexc =
        // -99.;
        // globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_recomb =
        // -99.;
        // globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_recomb =
        // -99.;
        // globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_same
        // = -99.;
        // globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_same
        // = -99.;
        // globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_lower
        // = -99.;
        // globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_higher
        // = -99.;
        //
        // ndowntrans = get_ndowntrans(element, ion, level);
        // nuptrans = get_nuptrans(element, ion, level);
        // for (i = 0; i < ndowntrans; i++)
        // {
        //   globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].sum_epstrans_rad_deexc[i]
        //   = -99.;
        //   globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_down_same[i]
        //   = -99.;
        // }
        // for (i = 0; i < nuptrans; i++)
        // {
        //   globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].sum_internal_up_same[i]
        //   = -99.;
        // }
      }
    }
  }

  if (modelgridindex >= 0) {
    const int nbfcont = globals::nbfcontinua;
    std::fill_n(globals::cellhistory[tid].ch_allcont_departureratios, nbfcont, -1);
  }
  // printout("nlevels_with_processrates %d\n", nlevels_with_processrates);

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
  const time_t sys_time_start_calculate_bfheatingcoeffs = time(nullptr);
  calculate_bfheatingcoeffs(n);
  printout("took %ld seconds\n", time(nullptr) - sys_time_start_calculate_bfheatingcoeffs);

  const double covergence_tolerance = 0.04;
  for (int nlte_iter = 0; nlte_iter <= NLTEITER; nlte_iter++) {
    const time_t sys_time_start_spencerfano = time(nullptr);
    if (NT_ON && NT_SOLVE_SPENCERFANO) {
      // SF solution depends on the ionization balance, and weakly on nne
      nonthermal::solve_spencerfano(n, nts, nlte_iter);
    }
    const int duration_solve_spencerfano = time(nullptr) - sys_time_start_spencerfano;

    const time_t sys_time_start_partfuncs_or_gamma = time(nullptr);
    for (int element = 0; element < get_nelements(); element++) {
      if (!elem_has_nlte_levels(element)) {
        calculate_cellpartfuncts(n, element);
      } else if (USE_LUT_PHOTOION && (nlte_iter != 0)) {
        // recalculate the Gammas using the current population estimates
        const int nions = get_nions(element);
        for (int ion = 0; ion < nions - 1; ion++) {
          globals::gammaestimator[get_ionestimindex(n, element, ion)] = calculate_iongamma_per_gspop(n, element, ion);
        }
      }
    }
    const int duration_solve_partfuncs_or_gamma = time(nullptr) - sys_time_start_partfuncs_or_gamma;

    const double prev_T_e = grid::get_Te(n);
    const time_t sys_time_start_Te = time(nullptr);
    const int nts_for_te = (titer == 0) ? nts - 1 : nts;

    /// Find T_e as solution for thermal balance
    call_T_e_finder(n, nts, globals::timesteps[nts_for_te].mid, MINTEMP, MAXTEMP, heatingcoolingrates);

    const int duration_solve_T_e = time(nullptr) - sys_time_start_Te;

    if (globals::total_nlte_levels == 0) {
      const time_t sys_time_start_pops = time(nullptr);
      calculate_ion_balance_nne(n);
      const int duration_solve_pops = time(nullptr) - sys_time_start_pops;

      printout(
          "Grid solver cell %d timestep %d: time spent on: Spencer-Fano %ds, partfuncs/gamma "
          "%ds, T_e %ds, populations %ds\n",
          n, nts, duration_solve_spencerfano, duration_solve_partfuncs_or_gamma, duration_solve_T_e,
          duration_solve_pops);
      break;  // no iteration is needed without nlte pops
    }

    if (globals::total_nlte_levels > 0) {
      const double fracdiff_T_e = fabs((grid::get_Te(n) / prev_T_e) - 1);
      const time_t sys_time_start_nltepops = time(nullptr);
      // fractional difference between previous and current iteration's (nne or max(ground state
      // population change))
      double fracdiff_nne = 0.;
      for (int element = 0; element < get_nelements(); element++) {
        if (get_nions(element) > 0) {
          solve_nlte_pops_element(element, n, nts, nlte_iter);
          calculate_cellpartfuncts(n, element);
        }
      }
      const int duration_solve_nltepops = time(nullptr) - sys_time_start_nltepops;

      const double nne_prev = grid::get_nne(n);
      calculate_ion_balance_nne(n);  // sets nne
      fracdiff_nne = fabs((grid::get_nne(n) / nne_prev) - 1);
      printout(
          "NLTE solver cell %d timestep %d iteration %d: time spent on: Spencer-Fano %ds, T_e "
          "%ds, NLTE populations %ds\n",
          n, nts, nlte_iter, duration_solve_spencerfano, duration_solve_T_e, duration_solve_nltepops);
      printout(
          "NLTE (Spencer-Fano/Te/pops) solver cell %d timestep %d iteration %d: prev_iter nne "
          "%g, new nne is %g, fracdiff %g, prev T_e %g new T_e %g fracdiff %g\n",
          n, nts, nlte_iter, nne_prev, grid::get_nne(n), fracdiff_nne, prev_T_e, grid::get_Te(n), fracdiff_T_e);

      if (fracdiff_nne <= covergence_tolerance && fracdiff_T_e <= covergence_tolerance) {
        printout(
            "NLTE (Spencer-Fano/Te/pops) solver nne converged to tolerance %g <= %g and T_e to "
            "tolerance %g <= %g after %d iterations.\n",
            fracdiff_nne, covergence_tolerance, fracdiff_T_e, covergence_tolerance, nlte_iter + 1);
        break;
      }
      if (nlte_iter == NLTEITER) {
        printout(
            "WARNING: NLTE solver failed to converge after %d iterations. Keeping solution from "
            "last iteration\n",
            nlte_iter + 1);
      }
    }
  }
}

static void update_gamma_corrphotoionrenorm_bfheating_estimators(const int n, const double estimator_normfactor) {
  assert_always(USE_LUT_PHOTOION || USE_LUT_BFHEATING);
  if constexpr (USE_LUT_PHOTOION) {
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions - 1; ion++) {
        const int ionestimindex = get_ionestimindex(n, element, ion);
        globals::gammaestimator[ionestimindex] *= estimator_normfactor / H;
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

      /// 2012-01-11. These loops should terminate here to precalculate *ALL* corrphotoionrenorm
      /// values so that the values are known when required by the call to get_corrphotoioncoeff in
      /// the following loops. Otherwise get_corrphotoioncoeff tries to renormalize by the closest
      /// corrphotoionrenorm in frequency space which can lead to zero contributions to the total
      /// photoionsation rate!
    }
  }
  if constexpr (USE_LUT_PHOTOION || USE_LUT_BFHEATING) {
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions - 1; ion++) {
        /// Reuse the gammaestimator array as temporary storage of the Gamma values during
        /// the remaining part of the update_grid phase. Afterwards it is reset to record
        /// the next timesteps gamma estimators.
        const int ionestimindex = get_ionestimindex(n, element, ion);

        if constexpr (USE_LUT_PHOTOION) {
          globals::gammaestimator[ionestimindex] = calculate_iongamma_per_gspop(n, element, ion);
        }

        if constexpr (USE_LUT_BFHEATING) {
          globals::bfheatingestimator[ionestimindex] *= estimator_normfactor;
#ifdef DO_TITER
          if (bfheatingestimator_save[ionestimindex] >= 0) {
            globals::bfheatingestimator[ionestimindex] =
                (globals::bfheatingestimator[ionestimindex] + bfheatingestimator_save[ionestimindex]) / 2;
          }
          bfheatingestimator_save[ionestimindex] = globals::bfheatingestimator[ionestimindex];
#endif
          /// Now convert bfheatingestimator into the bfheating renormalisation coefficient used in
          /// get_bfheating in the remaining part of update_grid. Later on it's reset and new
          /// contributions are added up.

          const double bfheatingcoeff_ana = get_bfheatingcoeff_ana(element, ion, 0, 0, grid::get_TR(n), grid::get_W(n));
          globals::bfheatingestimator[ionestimindex] = globals::bfheatingestimator[ionestimindex] / bfheatingcoeff_ana;

          if (!std::isfinite(globals::bfheatingestimator[ionestimindex])) {
            printout(
                "[fatal] about to set bfheatingestimator = NaN = bfheatingestimator / "
                "get_bfheatingcoeff_ana(%d,%d,%d,%d,%d)=%g/%g",
                element, ion, 0, 0, n, globals::bfheatingestimator[ionestimindex], bfheatingcoeff_ana);
            abort();
          }
        }
      }
    }
  }
}

#ifdef DO_TITER
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

static void zero_gammaestimator(const int modelgridindex) {
  assert_always(USE_LUT_PHOTOION);
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      globals::gammaestimator[get_ionestimindex(modelgridindex, element, ion)] = 0.;
    }
  }
}

static void set_all_corrphotoionrenorm(const int modelgridindex, const double value) {
  assert_always(USE_LUT_PHOTOION);
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      globals::corrphotoionrenorm[get_ionestimindex(modelgridindex, element, ion)] = value;
    }
  }
}

static void update_grid_cell(const int mgi, const int nts, const int nts_prev, const int titer, const double tratmid,
                             const double deltat, struct heatingcoolingrates *heatingcoolingrates)
// n is the modelgrid index
{
  const int assoc_cells = grid::get_numassociatedcells(mgi);
  if (assoc_cells > 0) {
    const double deltaV =
        grid::get_modelcell_assocvolume_tmin(mgi) * pow(globals::timesteps[nts_prev].mid / globals::tmin, 3);
    const time_t sys_time_start_update_cell = time(nullptr);

    printout("update_grid_cell: working on cell %d before timestep %d titeration %d...\n", mgi, nts, titer);

    /// Update current mass density of cell
    grid::set_rho(mgi, grid::get_rho_tmin(mgi) / pow(tratmid, 3));

    /// Update elemental abundances with radioactive decays
    decay::update_abundances(mgi, nts, globals::timesteps[nts].mid);

    const double estimator_normfactor = 1 / deltaV / deltat / globals::nprocs;
    const double estimator_normfactor_over4pi = ONEOVER4PI * estimator_normfactor;

    if (globals::opacity_case < 4) {
      // various forms of grey opacity
      grid::modelgrid[mgi].thick = 1;

      if (globals::opacity_case == 3) {
        // printout("update_grid: opacity_case 3 ... updating globals::cell[n].chi_grey"); //MK
        if (grid::get_rho(mgi) > globals::rho_crit) {
          grid::set_kappagrey(mgi, globals::opcase3_normal * (0.9 * grid::get_ffegrp(mgi) + 0.1) * globals::rho_crit /
                                       grid::get_rho(mgi));
        } else {
          grid::set_kappagrey(mgi, globals::opcase3_normal * (0.9 * grid::get_ffegrp(mgi) + 0.1));
        }
      }
    }

    if (nts == globals::timestep_initial && titer == 0) {
      // For the initial timestep, temperatures have already been assigned
      // either by trapped energy release calculation, or reading from gridsave file

      if (USE_LUT_PHOTOION && !globals::simulation_continued_from_saved) {
        /// Determine renormalisation factor for corrected photoionization cross-sections
        set_all_corrphotoionrenorm(mgi, 1.);
      }

      /// W == 1 indicates that this modelgrid cell was treated grey in the
      /// last timestep. Therefore it has no valid Gamma estimators and must
      /// be treated in LTE at restart.
      if (grid::modelgrid[mgi].thick == 0 && grid::get_W(mgi) == 1) {
        printout(
            "force modelgrid cell %d to grey/LTE for update grid since existing W == 1. (will not have gamma "
            "estimators)\n",
            mgi);
        grid::modelgrid[mgi].thick = 1;
      }

      printout("lte_iteration %d\n", globals::lte_iteration);
      printout("mgi %d modelgrid.thick: %d (for this grid update only)\n", mgi, grid::modelgrid[mgi].thick);

      for (int element = 0; element < get_nelements(); element++) {
        calculate_cellpartfuncts(mgi, element);
      }
      if (!globals::simulation_continued_from_saved) {
        calculate_ion_balance_nne(mgi);
      }
    } else {
      // For all other timesteps temperature corrections have to be applied

      /// we have to calculate the electron density
      /// and all the level populations
      /// Normalise estimators and make sure that they are finite.
      /// Then update T_R and W using the estimators.
      /// (This could in principle also be done for empty cells)

      const time_t sys_time_start_temperature_corrections = time(nullptr);

      radfield::normalise_J(mgi, estimator_normfactor_over4pi);  // this applies normalisation to the fullspec J
      radfield::set_J_normfactor(mgi,
                                 estimator_normfactor_over4pi);  // this stores the factor that will be applied
                                                                 // later for the J bins but not fullspec J

#ifdef DO_TITER
      radfield::titer_J(mgi);
#endif

      if constexpr (TRACK_ION_STATS) {
        stats::normalise_ion_estimators(mgi, deltat, deltaV);
      }

      // lte_iteration really means either ts 0 or nts < globals::num_lte_timesteps
      if (globals::lte_iteration || grid::modelgrid[mgi].thick == 1) {
        // LTE mode or grey mode (where temperature doesn't matter but is calculated anyway)

        const double T_J = radfield::get_T_J_from_J(mgi);
        grid::set_TR(mgi, T_J);
        grid::set_Te(mgi, T_J);
        grid::set_TJ(mgi, T_J);
        grid::set_W(mgi, 1);

        if constexpr (USE_LUT_PHOTOION) {
          set_all_corrphotoionrenorm(mgi, 1.);
        }

        for (int element = 0; element < get_nelements(); element++) {
          calculate_cellpartfuncts(mgi, element);
        }
        calculate_ion_balance_nne(mgi);
      } else {
        // not lte_iteration and not a thick cell
        // non-LTE timesteps with T_e from heating/cooling

        radfield::normalise_nuJ(mgi, estimator_normfactor_over4pi);

        globals::ffheatingestimator[mgi] *= estimator_normfactor;
        globals::colheatingestimator[mgi] *= estimator_normfactor;

#ifdef DO_TITER
        radfield::titer_nuJ(mgi);
        titer_average_estimators(mgi);
#endif

        if constexpr (USE_LUT_PHOTOION || USE_LUT_BFHEATING) {
          update_gamma_corrphotoionrenorm_bfheating_estimators(mgi, estimator_normfactor);
        }

        // Get radiation field parameters (T_J, T_R, W, and bins if enabled) out of the
        // full-spectrum and binned J and nuJ estimators
        radfield::fit_parameters(mgi, nts);

        if constexpr (DETAILED_BF_ESTIMATORS_ON) {
          radfield::normalise_bf_estimators(mgi, estimator_normfactor / H);
        }

        solve_Te_nltepops(mgi, nts, titer, heatingcoolingrates);
      }
      printout("Temperature/NLTE solution for cell %d timestep %d took %ld seconds\n", mgi, nts,
               time(nullptr) - sys_time_start_temperature_corrections);
    }

    const float nne = grid::get_nne(mgi);
    const double compton_optical_depth = SIGMA_T * nne * grid::wid_init(mgi, 0) * tratmid;

    double const radial_pos = grid::modelgrid[mgi].initial_radial_pos_sum * tratmid / assoc_cells;
    const double grey_optical_deptha = grid::get_kappagrey(mgi) * grid::get_rho(mgi) * grid::wid_init(mgi, 0) * tratmid;
    // cube corners will have radial pos > rmax, so clamp to 0.
    const double dist_to_obs = std::max(0., globals::rmax * tratmid - radial_pos);
    const double grey_optical_depth = grid::get_kappagrey(mgi) * grid::get_rho(mgi) * dist_to_obs;
    printout(
        "modelgridcell %d, compton optical depth (/propgridcell) %g, grey optical depth "
        "(/propgridcell) %g\n",
        mgi, compton_optical_depth, grey_optical_deptha);
    printout("radial_pos %g, distance_to_obs %g, tau_dist %g\n", radial_pos, dist_to_obs, grey_optical_depth);

    grid::modelgrid[mgi].grey_depth = grey_optical_depth;

    // grey_optical_depth = compton_optical_depth;

    if ((grey_optical_depth >= globals::cell_is_optically_thick) && (nts < globals::num_grey_timesteps)) {
      printout("timestep %d cell %d is treated in grey approximation (chi_grey %g [cm2/g], tau %g >= %g)\n", nts, mgi,
               grid::get_kappagrey(mgi), grey_optical_depth, globals::cell_is_optically_thick);
      grid::modelgrid[mgi].thick = 1;
    } else if (VPKT_ON && (grey_optical_depth > cell_is_optically_thick_vpkt)) {
      grid::modelgrid[mgi].thick = 2;
    } else {
      grid::modelgrid[mgi].thick = 0;
    }

    if (grid::modelgrid[mgi].thick == 1) {
      // cooling rates calculation can be skipped for thick cells
      // flag with negative numbers to indicate that the rates are invalid
      grid::modelgrid[mgi].totalcooling = -1.;
      const int element = 0;
      const int ion = 0;
      grid::modelgrid[mgi].cooling_contrib_ion[element][ion] = -1.;
    } else if (globals::simulation_continued_from_saved && nts == globals::timestep_initial) {
      // cooling rates were read from the gridsave file for this timestep
      // make sure they are valid
      assert_always(grid::modelgrid[mgi].totalcooling >= 0.);
      const int element = 0;
      const int ion = 0;
      assert_always(grid::modelgrid[mgi].cooling_contrib_ion[element][ion] >= 0.);
    } else {
      /// Cooling rates depend only on cell properties, precalculate total cooling
      /// and ion contributions inside update grid and communicate between MPI tasks
      const time_t sys_time_start_calc_kpkt_rates = time(nullptr);

      printout("calculate_cooling_rates for timestep %d cell %d...", nts, mgi);

      // don't pass pointer to heatingcoolingrates because current populations and rates weren't
      // used to determine T_e
      kpkt::calculate_cooling_rates(mgi, nullptr);

      printout("took %ld seconds\n", time(nullptr) - sys_time_start_calc_kpkt_rates);
    }

    const int update_grid_cell_seconds = time(nullptr) - sys_time_start_update_cell;
    if (update_grid_cell_seconds > 0) {
      printout("update_grid_cell for cell %d timestep %d took %ld seconds\n", mgi, nts, update_grid_cell_seconds);
    }
  } else {
    /// For modelgrid cells that are not represented in the simulation grid,
    /// Set grid properties to zero
    grid::set_TR(mgi, 0.);
    grid::set_TJ(mgi, 0.);
    grid::set_Te(mgi, 0.);
    grid::set_W(mgi, 0.);

    if constexpr (USE_LUT_PHOTOION) {
      zero_gammaestimator(mgi);
      set_all_corrphotoionrenorm(mgi, 0.);
    }
  }
}

void update_grid(FILE *estimators_file, const int nts, const int nts_prev, const int my_rank, const int nstart,
                 const int ndo, const int titer, const time_t real_time_start)
// Subroutine to update the matter quantities in the grid cells at the start
//   of the new timestep.
/// nts timestep
{
  const time_t sys_time_start_update_grid = time(nullptr);
  printout("\n");
  printout("timestep %d: time before update grid %ld (tstart + %ld) simtime ts_mid %g days\n", nts,
           sys_time_start_update_grid, sys_time_start_update_grid - real_time_start, globals::timesteps[nts].mid / DAY);

  if constexpr (USE_LUT_PHOTOION) {
    /// Initialise globals::corrphotoionrenorm[i] to zero before update_grid is called
    /// unless they have been read from file
    if ((!globals::simulation_continued_from_saved) || (nts - globals::timestep_initial != 0) || (titer != 0)) {
      printout("nts %d, titer %d: reset corr photoionrenorm\n", nts, titer);
      for (int i = 0; i < grid::get_npts_model() * get_nelements() * get_max_nions(); i++) {
        globals::corrphotoionrenorm[i] = 0.;
      }
      printout("after nts %d, titer %d: reset corr photoionrenorm\n", nts, titer);
    }
  }

  // printout("[debug] update_grid: starting update for timestep %d...\n",m);
  const double tratmid = globals::timesteps[nts].mid / globals::tmin;

  /// Thread private substitution of max_path_step. Its minimum is
  /// assigned to max_path_step after the parallel update_grid finished.
  auto mps = std::make_unique<double[]>(get_max_threads());

  for (int i = 0; i < get_max_threads(); i++) {
    mps[i] = 1.e35;
  }

  /// Calculate the critical opacity at which opacity_case 3 switches from a
  /// regime proportional to the density to a regime independent of the density
  /// This is done by solving for tau_sobolev == 1
  /// tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/nucmass(28, 56) * 3000e-8 *
  /// globals::timesteps[m].mid;
  globals::rho_crit = ME * CLIGHT * decay::nucmass(28, 56) /
                      (PI * QE * QE * globals::rho_crit_para * 3000e-8 * globals::timesteps[nts].mid);
  printout("update_grid: rho_crit = %g\n", globals::rho_crit);

  // These values will not be used if nts == 0, but set them anyway
  // nts_prev is the previous timestep, unless this is timestep zero
  const double deltat = globals::timesteps[nts_prev].width;

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
        // use_cellhist = false;
        // cellhistory_reset(-99, true);

        struct heatingcoolingrates heatingcoolingrates = {};
        update_grid_cell(mgi, nts, nts_prev, titer, tratmid, deltat, &heatingcoolingrates);

        // maybe want to add omp ordered here if the modelgrid cells should be output in order
        // use_cellhist = true;
        // cellhistory_reset(mgi, true);
#ifdef _OPENMP
#pragma omp critical(estimators_file)
#endif
        { write_to_estimators_file(estimators_file, mgi, nts, titer, &heatingcoolingrates); }

      } else {
        /// else, only reset gammaestimator to zero. This allows us to do a global MPI
        /// communication after update_grid to synchronize gammaestimator
        /// and write a contiguous restart file with grid properties
        if constexpr (USE_LUT_PHOTOION) {
          zero_gammaestimator(mgi);
        }
      }
    }  /// end parallel for loop over all modelgrid cells

    /// Now after all the relevant taks of update_grid have been finished activate
    /// the use of the cellhistory for all OpenMP tasks, in what follows (update_packets)
    use_cellhist = true;
  }  /// end OpenMP parallel section

  // alterative way to write out estimators. this keeps the modelgrid cells in order but
  // heatingrates are not valid. #ifdef _OPENMP for (int n = nstart; n < nstart+nblock; n++)
  // {
  //   write_to_estimators_file(n,nts);
  // }
  // #endif

  /// Assign the minimum of thread private mps to the global variable max_path_step
  globals::max_path_step = mps[0];
  for (int i = 1; i < get_max_threads(); i++) {
    if (mps[i] < globals::max_path_step) {
      globals::max_path_step = mps[i];
    }
  }

  globals::max_path_step = fmin(globals::max_path_step, globals::rmax / 10.);
  printout("max_path_step %g\n", globals::max_path_step);

  const time_t time_update_grid_end_thisrank = time(nullptr);
  printout("finished update grid on this rank at time %ld\n", time_update_grid_end_thisrank);

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  printout(
      "timestep %d: time after update grid for all processes %ld (rank %d took %lds, waited "
      "%lds, total %lds)\n",
      nts, time(nullptr), my_rank, time_update_grid_end_thisrank - sys_time_start_update_grid,
      time(nullptr) - time_update_grid_end_thisrank, time(nullptr) - sys_time_start_update_grid);
}

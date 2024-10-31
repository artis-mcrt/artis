#include "update_grid.h"

#ifdef MPI_ON
#include <mpi.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "kpkt.h"
#include "ltepop.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "rpkt.h"
#include "sn3d.h"
#include "stats.h"
#include "thermalbalance.h"
#include "vpkt.h"

namespace {

void write_to_estimators_file(FILE *estimators_file, const int mgi, const int timestep, const int titer,
                              const HeatingCoolingRates *heatingcoolingrates) {
  // return; disable for better performance (if estimators files are not needed)

  if (grid::get_numpropcells(mgi) < 1) {
    // modelgrid cells that are not represented in the simulation grid
    fprintf(estimators_file, "timestep %d modelgridindex %d EMPTYCELL\n\n", timestep, mgi);
    fflush(estimators_file);
    return;
  }

  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);

  const auto sys_time_start_write_estimators = std::time(nullptr);

  printout("writing to estimators file timestep %d cell %d...\n", timestep, mgi);

  const auto T_e = grid::get_Te(mgi);
  const auto nne = grid::get_nne(mgi);
  const auto Y_e = grid::get_electronfrac(nonemptymgi);
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

    decay::fprint_nuc_abundances(estimators_file, nonemptymgi, globals::timesteps[timestep].mid, element);

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

    if constexpr (TRACK_ION_STATS) {
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
        fprintf(estimators_file, "  %d: %9.3f", get_ionstage(element, ion), 1. - (alpha_r_mc_abs / alpha_r_mc));
      }
      fprintf(estimators_file, "\n");

      fprintf(estimators_file, "BB_escfrac         Z=%2d", get_atomicnumber(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
        fprintf(estimators_file, "              ");
      }
      for (int ion = 0; ion < nions; ion++) {
        const double bb_emitted = get_ion_stats(mgi, element, ion, stats::ION_BOUNDBOUND_MACROATOM);
        const double bb_abs = get_ion_stats(mgi, element, ion, stats::ION_BOUNDBOUND_ABSORBED);
        fprintf(estimators_file, "  %d: %9.3f", get_ionstage(element, ion), 1. - (bb_abs / bb_emitted));
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

    if (USE_LUT_PHOTOION && globals::nbfcontinua_ground > 0) {
      fprintf(estimators_file, "corrphotoionrenorm Z=%2d", get_atomicnumber(element));
      for (int ion = 0; ion < nions - 1; ion++) {
        const int groundcontindex = globals::elements[element].ions[ion].groundcontindex;
        if (groundcontindex >= 0) {
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                  globals::corrphotoionrenorm[get_ionestimindex_nonemptymgi(nonemptymgi, element, ion)]);
        }
      }
      fprintf(estimators_file, "\n");
      fprintf(estimators_file, "gammaestimator     Z=%2d", get_atomicnumber(element));
      for (int ion = 0; ion < nions - 1; ion++) {
        const int groundcontindex = globals::elements[element].ions[ion].groundcontindex;
        if (groundcontindex >= 0) {
          fprintf(estimators_file, "  %d: %9.3e", get_ionstage(element, ion),
                  globals::gammaestimator[get_ionestimindex_nonemptymgi(nonemptymgi, element, ion)]);
        }
      }
      fprintf(estimators_file, "\n");
    }
  }

  // power densities in erg / s / cm^3
  // ana means analytical at t_mid, i.e. the rates calculated from the nuclear abundances and decay data, not from
  // Monte Carlo
  fprintf(estimators_file, "emission_ana: gamma %11.5e positron %11.5e electron %11.5e alpha %11.5e\n",
          heatingcoolingrates->eps_gamma_ana, heatingcoolingrates->eps_positron_ana,
          heatingcoolingrates->eps_electron_ana, heatingcoolingrates->eps_alpha_ana);
  fprintf(estimators_file, "deposition: gamma %11.5e positron %11.5e electron %11.5e alpha %11.5e\n",
          heatingcoolingrates->dep_gamma, heatingcoolingrates->dep_positron, heatingcoolingrates->dep_electron,
          heatingcoolingrates->dep_alpha);
  fprintf(estimators_file, "heating: ff %11.5e bf %11.5e coll %11.5e       dep %11.5e heating_dep/total_dep %.3f\n",
          heatingcoolingrates->heating_ff, heatingcoolingrates->heating_bf, heatingcoolingrates->heating_collisional,
          heatingcoolingrates->heating_dep, heatingcoolingrates->dep_frac_heating);
  fprintf(estimators_file, "cooling: ff %11.5e fb %11.5e coll %11.5e adiabatic %11.5e\n",
          heatingcoolingrates->cooling_ff, heatingcoolingrates->cooling_fb, heatingcoolingrates->cooling_collisional,
          heatingcoolingrates->cooling_adiabatic);

  fprintf(estimators_file, "\n");

  fflush(estimators_file);

  const auto write_estim_duration = std::time(nullptr) - sys_time_start_write_estimators;
  if (write_estim_duration >= 1) {
    printout("writing estimators for timestep %d cell %d took %ld seconds\n", timestep, mgi, write_estim_duration);
  }
}

void solve_Te_nltepops(const int nonemptymgi, const int nts, const int nts_prev,
                       HeatingCoolingRates *heatingcoolingrates)
// nts is the timestep number
{
  const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  // bfheating coefficients are needed for the T_e solver, but
  // they only depend on the radiation field, which is fixed during the iterations below
  printout("calculate_bfheatingcoeffs for timestep %d cell %d...", nts, mgi);
  const auto sys_time_start_calculate_bfheatingcoeffs = std::time(nullptr);
  thread_local static auto bfheatingcoeffs = std::vector<double>(get_includedlevels());

  calculate_bfheatingcoeffs(nonemptymgi, bfheatingcoeffs);
  printout("took %ld seconds\n", std::time(nullptr) - sys_time_start_calculate_bfheatingcoeffs);

  const double convergence_tolerance = 0.04;
  for (int nlte_iter = 0; nlte_iter <= NLTEITER; nlte_iter++) {
    const auto sys_time_start_spencerfano = std::time(nullptr);
    if (NT_ON && NT_SOLVE_SPENCERFANO) {
      // SF solution depends on the ionization balance, and weakly on nne
      nonthermal::solve_spencerfano(mgi, nts, nlte_iter);
    }
    const int duration_solve_spencerfano = std::time(nullptr) - sys_time_start_spencerfano;

    const auto sys_time_start_partfuncs_or_gamma = std::time(nullptr);
    for (int element = 0; element < get_nelements(); element++) {
      if (!elem_has_nlte_levels(element)) {
        calculate_cellpartfuncts(mgi, element);
      }
    }
    const int duration_solve_partfuncs_or_gamma = std::time(nullptr) - sys_time_start_partfuncs_or_gamma;

    const double prev_T_e = grid::get_Te(mgi);
    const auto sys_time_start_Te = std::time(nullptr);

    // Find T_e as solution for thermal balance
    call_T_e_finder(nonemptymgi, globals::timesteps[nts_prev].mid, MINTEMP, MAXTEMP, heatingcoolingrates,
                    bfheatingcoeffs);

    const int duration_solve_T_e = std::time(nullptr) - sys_time_start_Te;

    if (globals::total_nlte_levels == 0) {
      const auto sys_time_start_pops = std::time(nullptr);
      calculate_ion_balance_nne(mgi);
      const int duration_solve_pops = std::time(nullptr) - sys_time_start_pops;

      printout(
          "Grid solver cell %d timestep %d: time spent on: Spencer-Fano %ds, partfuncs/gamma "
          "%ds, T_e %ds, populations %ds\n",
          mgi, nts, duration_solve_spencerfano, duration_solve_partfuncs_or_gamma, duration_solve_T_e,
          duration_solve_pops);
      break;  // no iteration is needed without nlte pops
    }

    const double fracdiff_T_e = fabs((grid::get_Te(mgi) / prev_T_e) - 1);
    const auto sys_time_start_nltepops = std::time(nullptr);
    // fractional difference between previous and current iteration's (nne or max(ground state
    // population change))
    double fracdiff_nne = 0.;
    for (int element = 0; element < get_nelements(); element++) {
      if (get_nions(element) > 0 && elem_has_nlte_levels(element)) {
        solve_nlte_pops_element(element, mgi, nts, nlte_iter);
        calculate_cellpartfuncts(mgi, element);
      }
    }
    const int duration_solve_nltepops = std::time(nullptr) - sys_time_start_nltepops;

    const double nne_prev = grid::get_nne(mgi);
    calculate_ion_balance_nne(mgi);  // sets nne
    fracdiff_nne = fabs((grid::get_nne(mgi) / nne_prev) - 1);
    printout(
        "NLTE solver cell %d timestep %d iteration %d: time spent on: Spencer-Fano %ds, T_e "
        "%ds, NLTE populations %ds\n",
        mgi, nts, nlte_iter, duration_solve_spencerfano, duration_solve_T_e, duration_solve_nltepops);
    printout(
        "NLTE (Spencer-Fano/Te/pops) solver cell %d timestep %d iteration %d: prev_iter nne "
        "%g, new nne is %g, fracdiff %g, prev T_e %g new T_e %g fracdiff %g\n",
        mgi, nts, nlte_iter, nne_prev, grid::get_nne(mgi), fracdiff_nne, prev_T_e, grid::get_Te(mgi), fracdiff_T_e);

    if (fracdiff_nne <= convergence_tolerance && fracdiff_T_e <= convergence_tolerance) {
      printout(
          "NLTE (Spencer-Fano/Te/pops) solver nne converged to tolerance %g <= %g and T_e to "
          "tolerance %g <= %g after %d iterations.\n",
          fracdiff_nne, convergence_tolerance, fracdiff_T_e, convergence_tolerance, nlte_iter + 1);
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

void update_gamma_corrphotoionrenorm_bfheating_estimators(const int mgi, const double estimator_normfactor) {
  const int nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);
  if constexpr (USE_LUT_PHOTOION) {
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions - 1; ion++) {
        const int groundcontindex = globals::elements[element].ions[ion].groundcontindex;
        if (groundcontindex < 0) {
          continue;
        }
        const int ionestimindex = (nonemptymgi * globals::nbfcontinua_ground) + groundcontindex;

        globals::gammaestimator[ionestimindex] *= estimator_normfactor / H;
#ifdef DO_TITER
        if (globals::gammaestimator_save[ionestimindex] >= 0) {
          globals::gammaestimator[ionestimindex] =
              (globals::gammaestimator[ionestimindex] + globals::gammaestimator_save[ionestimindex]) / 2;
        }
        globals::gammaestimator_save[ionestimindex] = globals::gammaestimator[ionestimindex];
#endif

        globals::corrphotoionrenorm[ionestimindex] =
            globals::gammaestimator[ionestimindex] / get_corrphotoioncoeff_ana(element, ion, 0, 0, mgi);

        if (!std::isfinite(globals::corrphotoionrenorm[ionestimindex])) {
          printout(
              "[fatal] about to set corrphotoionrenorm = NaN = gammaestimator / "
              "get_corrphotoioncoeff_ana(%d,%d,%d,%d,%d)=%g/%g",
              element, ion, 0, 0, mgi, globals::gammaestimator[ionestimindex],
              get_corrphotoioncoeff_ana(element, ion, 0, 0, mgi));
          std::abort();
        }
      }

      // 2012-01-11. These loops should terminate here to precalculate *ALL* corrphotoionrenorm
      // values so that the values are known when required by the call to get_corrphotoioncoeff in
      // the following loops. Otherwise get_corrphotoioncoeff tries to renormalize by the closest
      // corrphotoionrenorm in frequency space which can lead to zero contributions to the total
      // photoionsation rate!
    }
  }
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions - 1; ion++) {
      // Reuse the gammaestimator array as temporary storage of the Gamma values during
      // the remaining part of the update_grid phase. Afterwards it is reset to record
      // the next timesteps gamma estimators.
      const int groundcontindex = globals::elements[element].ions[ion].groundcontindex;
      if (groundcontindex < 0) {
        continue;
      }
      const int ionestimindex = (nonemptymgi * globals::nbfcontinua_ground) + groundcontindex;

      if (!elem_has_nlte_levels(element)) {
        globals::gammaestimator[ionestimindex] = calculate_iongamma_per_gspop(mgi, element, ion);
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
        // Now convert bfheatingestimator into the bfheating renormalisation coefficient used in
        // get_bfheating in the remaining part of update_grid. Later on it's reset and new
        // contributions are added up.

        const double bfheatingcoeff_ana =
            get_bfheatingcoeff_ana(element, ion, 0, 0, grid::get_TR(mgi), grid::get_W(mgi));
        globals::bfheatingestimator[ionestimindex] = globals::bfheatingestimator[ionestimindex] / bfheatingcoeff_ana;

        if (!std::isfinite(globals::bfheatingestimator[ionestimindex])) {
          printout(
              "[fatal] about to set bfheatingestimator = NaN = bfheatingestimator / "
              "get_bfheatingcoeff_ana(%d,%d,%d,%d,%d)=%g/%g",
              element, ion, 0, 0, mgi, globals::bfheatingestimator[ionestimindex], bfheatingcoeff_ana);
          std::abort();
        }
      }
    }
  }
}

#ifdef DO_TITER
static void titer_average_estimators(const int nonemptymgi) {
  if (globals::ffheatingestimator_save[nonemptymgi] >= 0) {
    globals::ffheatingestimator[nonemptymgi] =
        (globals::ffheatingestimator[nonemptymgi] + globals::ffheatingestimator_save[nonemptymgi]) / 2;
  }
  globals::ffheatingestimator_save[nonemptymgi] = globals::ffheatingestimator[nonemptymgi];
  if (globals::colheatingestimator_save[nonemptymgi] >= 0) {
    globals::colheatingestimator[nonemptymgi] =
        (globals::colheatingestimator[nonemptymgi] + globals::colheatingestimator_save[nonemptymgi]) / 2;
  }
  globals::colheatingestimator_save[nonemptymgi] = globals::colheatingestimator[nonemptymgi];
}
#endif

void update_grid_cell(const int mgi, const int nts, const int nts_prev, const int titer, const double tratmid,
                      const double deltat, HeatingCoolingRates *heatingcoolingrates) {
  const int assoc_cells = grid::get_numpropcells(mgi);
  if (assoc_cells < 1) {
    // For modelgrid cells that are not represented in the simulation grid,
    // Set grid properties to zero
    grid::set_TR(mgi, 0.);
    grid::set_TJ(mgi, 0.);
    grid::set_Te(mgi, 0.);
    grid::set_W(mgi, 0.);
    return;
  }

  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);

  const double deltaV =
      grid::get_modelcell_assocvolume_tmin(mgi) * pow(globals::timesteps[nts_prev].mid / globals::tmin, 3);
  const auto sys_time_start_update_cell = std::time(nullptr);

  printout("update_grid_cell: working on cell %d before timestep %d titeration %d...\n", mgi, nts, titer);

  // Update current mass density of cell
  grid::set_rho(mgi, grid::get_rho_tmin(mgi) / pow(tratmid, 3));

  // Update elemental abundances with radioactive decays
  decay::update_abundances(nonemptymgi, nts, globals::timesteps[nts].mid);
  nonthermal::calculate_deposition_rate_density(nonemptymgi, nts, heatingcoolingrates);

  if (globals::opacity_case == 6) {
    grid::calculate_kappagrey();
  }

  const double estimator_normfactor = 1 / deltaV / deltat / globals::nprocs;
  const double estimator_normfactor_over4pi = ONEOVER4PI * estimator_normfactor;

  if (globals::opacity_case < 4) {
    // various forms of grey opacity
    grid::modelgrid[mgi].thick = 1;

    if (globals::opacity_case == 3) {
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
      // Determine renormalisation factor for corrected photoionization cross-sections
      std::fill_n(globals::corrphotoionrenorm + (nonemptymgi * globals::nbfcontinua_ground),
                  globals::nbfcontinua_ground, 1.);
    }

    // W == 1 indicates that this modelgrid cell was treated grey in the
    // last timestep. Therefore it has no valid Gamma estimators and must
    // be treated in LTE at restart.
    if (grid::modelgrid[mgi].thick != 1 && grid::get_W(mgi) == 1) {
      printout(
          "force modelgrid cell %d to grey/LTE thick = 1 for update grid since existing W == 1. (will not have "
          "gamma estimators)\n",
          mgi);
      grid::modelgrid[mgi].thick = 1;
    }

    printout("lte_iteration %d\n", globals::lte_iteration ? 1 : 0);
    printout("mgi %d modelgrid.thick: %d (during grid update)\n", mgi, grid::modelgrid[mgi].thick);

    for (int element = 0; element < get_nelements(); element++) {
      calculate_cellpartfuncts(mgi, element);
    }
    if (!globals::simulation_continued_from_saved) {
      calculate_ion_balance_nne(mgi);
    }
  } else {
    // For all other timesteps temperature corrections have to be applied

    // we have to calculate the electron density
    // and all the level populations
    // Normalise estimators and make sure that they are finite.
    // Then update T_R and W using the estimators.
    // (This could in principle also be done for empty cells)

    const auto sys_time_start_temperature_corrections = std::time(nullptr);

    radfield::normalise_J(mgi, estimator_normfactor_over4pi);  // this applies normalisation to the fullspec J
    // this stores the factor that will be applied later for the J bins but not fullspec J
    radfield::set_J_normfactor(nonemptymgi, estimator_normfactor_over4pi);

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
        std::fill_n(globals::corrphotoionrenorm + (nonemptymgi * globals::nbfcontinua_ground),
                    globals::nbfcontinua_ground, 1.);
      }

      for (int element = 0; element < get_nelements(); element++) {
        calculate_cellpartfuncts(mgi, element);
      }
      calculate_ion_balance_nne(mgi);
    } else {
      // not lte_iteration and not a thick cell
      // non-LTE timesteps with T_e from heating/cooling

      radfield::normalise_nuJ(mgi, estimator_normfactor_over4pi);

      globals::ffheatingestimator[nonemptymgi] *= estimator_normfactor;
      globals::colheatingestimator[nonemptymgi] *= estimator_normfactor;

#ifdef DO_TITER
      radfield::titer_nuJ(mgi);
      titer_average_estimators(mgi);
#endif

      update_gamma_corrphotoionrenorm_bfheating_estimators(mgi, estimator_normfactor);

      // Get radiation field parameters (T_J, T_R, W, and bins if enabled) out of the
      // full-spectrum and binned J and nuJ estimators
      radfield::fit_parameters(mgi, nts);

      solve_Te_nltepops(nonemptymgi, nts, nts_prev, heatingcoolingrates);
    }
    printout("Temperature/NLTE solution for cell %d timestep %d took %ld seconds\n", mgi, nts,
             std::time(nullptr) - sys_time_start_temperature_corrections);
  }

  const float nne = grid::get_nne(mgi);
  const double compton_optical_depth = SIGMA_T * nne * grid::wid_init(mgi, 0) * tratmid;

  const double radial_pos = grid::modelgrid[mgi].initial_radial_pos_sum * tratmid / assoc_cells;
  const double grey_optical_deptha = grid::get_kappagrey(mgi) * grid::get_rho(mgi) * grid::wid_init(mgi, 0) * tratmid;
  // cube corners will have radial pos > rmax, so clamp to 0.
  const double dist_to_obs = std::max(0., (globals::rmax * tratmid) - radial_pos);
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
    grid::ion_cooling_contribs_allcells[(nonemptymgi * get_includedions()) + 0] = -1.;
  } else if (globals::simulation_continued_from_saved && nts == globals::timestep_initial) {
    // cooling rates were read from the gridsave file for this timestep
    // make sure they are valid
    printout("cooling rates read from gridsave file for timestep %d cell %d...", nts, mgi);
    assert_always(grid::modelgrid[mgi].totalcooling >= 0.);
    assert_always(grid::ion_cooling_contribs_allcells[(nonemptymgi * get_includedions()) + 0] >= 0.);
  } else {
    // Cooling rates depend only on cell properties, precalculate total cooling
    // and ion contributions inside update grid and communicate between MPI tasks
    const auto sys_time_start_calc_kpkt_rates = std::time(nullptr);

    printout("calculating cooling_rates for timestep %d cell %d...", nts, mgi);

    // don't pass pointer to heatingcoolingrates because current populations and rates weren't
    // used to determine T_e
    kpkt::calculate_cooling_rates(nonemptymgi, nullptr);

    printout("took %ld seconds\n", std::time(nullptr) - sys_time_start_calc_kpkt_rates);
  }

  if constexpr (EXPANSIONOPACITIES_ON || RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY > 0.) {
    if (grid::modelgrid[mgi].thick != 1) {
      calculate_expansion_opacities(mgi);
    }
  }

  const int update_grid_cell_seconds = std::time(nullptr) - sys_time_start_update_cell;
  if (update_grid_cell_seconds > 0) {
    printout("update_grid_cell for cell %d timestep %d took %d seconds\n", mgi, nts, update_grid_cell_seconds);
  }
}

}  // anonymous namespace

void update_grid(FILE *estimators_file, const int nts, const int nts_prev, const int titer,
                 const std::time_t real_time_start)
// Subroutine to update the matter quantities in the grid cells at the start
//   of the new timestep.
// nts timestep
{
  const auto my_rank = globals::my_rank;
  const auto sys_time_start_update_grid = std::time(nullptr);
  printout("\n");
  printout("timestep %d: time before update grid %ld (tstart + %ld) simtime ts_mid %g days\n", nts,
           sys_time_start_update_grid, sys_time_start_update_grid - real_time_start, globals::timesteps[nts].mid / DAY);

  const double tratmid = globals::timesteps[nts].mid / globals::tmin;

  // Calculate the critical opacity at which opacity_case 3 switches from a
  // regime proportional to the density to a regime independent of the density
  // This is done by solving for tau_sobolev == 1
  // tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/nucmass(28, 56) * 3000e-8 *
  // globals::timesteps[m].mid;
  globals::rho_crit = ME * CLIGHT * decay::nucmass(28, 56) /
                      (PI * QE * QE * globals::rho_crit_para * 3000e-8 * globals::timesteps[nts].mid);
  printout("update_grid: rho_crit = %g\n", globals::rho_crit);

  // These values will not be used if nts == 0, but set them anyway
  // nts_prev is the previous timestep, unless this is timestep zero
  const double deltat = globals::timesteps[nts_prev].width;

  cellcache_change_cell(-99);

  // Do not use values which are saved in the cellcache within update_grid
  use_cellcache = false;

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    radfield::normalise_bf_estimators(nts, nts_prev, titer, deltat);
  }

  const int nstart = grid::get_nstart(my_rank);
  const int ndo = grid::get_ndo(my_rank);
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
// Updating cell information
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif

    for (int mgi = nstart; mgi < nstart + ndo; mgi++) {
      // Check if this task should work on the current model grid cell.
      // If yes, update the cell and write out the estimators
      HeatingCoolingRates heatingcoolingrates{};
      update_grid_cell(mgi, nts, nts_prev, titer, tratmid, deltat, &heatingcoolingrates);

      // maybe want to add omp ordered here if the modelgrid cells should be output in order
#ifdef _OPENMP
#pragma omp critical(estimators_file)
#endif
      {
        write_to_estimators_file(estimators_file, mgi, nts, titer, &heatingcoolingrates);
      }
    }  // end parallel for loop over all modelgrid cells

  }  // end OpenMP parallel section

  // Now after all the relevant tasks of update_grid have been finished activate
  // the use of the cellcache for all OpenMP tasks, in what follows (update_packets)
  use_cellcache = true;

  // alternative way to write out estimators. this keeps the modelgrid cells in order but
  // heatingrates are not valid. #ifdef _OPENMP for (int n = nstart; n < nstart+nblock; n++)
  // {
  //   write_to_estimators_file(n,nts);
  // }
  // #endif

  globals::max_path_step = std::min(1.e35, globals::rmax / 10.);
  printout("max_path_step %g\n", globals::max_path_step);

  const auto time_update_grid_end_thisrank = std::time(nullptr);
  printout("finished update grid on this rank at time %ld\n", time_update_grid_end_thisrank);

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  printout(
      "timestep %d: time after update grid for all processes %ld (rank %d took %lds, waited "
      "%lds, total %lds)\n",
      nts, std::time(nullptr), my_rank, time_update_grid_end_thisrank - sys_time_start_update_grid,
      std::time(nullptr) - time_update_grid_end_thisrank, std::time(nullptr) - sys_time_start_update_grid);
}

void cellcache_change_cell(const int modelgridindex) {
  // All entries of the cellcache stack must be flagged as empty at the
  // onset of the new timestep. Also, boundary crossing?
  // Calculate the level populations for this cell, and flag the other entries
  // as empty.
  auto &cacheslot = globals::cellcache[cellcacheslotid];
  if (modelgridindex == cacheslot.cellnumber) {
    return;
  }

  cacheslot.cellnumber = modelgridindex;
  cacheslot.chi_ff_nnionpart = -1.;

  const int nelements = get_nelements();
  for (int element = 0; element < nelements; element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      cacheslot
          .cooling_contrib[kpkt::get_coolinglistoffset(element, ion) + kpkt::get_ncoolingterms_ion(element, ion) - 1] =
          COOLING_UNDEFINED;

      if (modelgridindex >= 0) {
        const int nlevels = get_nlevels(element, ion);
        auto &chion = cacheslot.chelements[element].chions[ion];
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int level = 0; level < nlevels; level++) {
          chion.chlevels[level].population = calculate_levelpop(modelgridindex, element, ion, level);
        }
      }
    }

    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_nlevels(element, ion);
      auto &chion = cacheslot.chelements[element].chions[ion];
      for (int level = 0; level < nlevels; level++) {
        const auto nphixstargets = get_nphixstargets(element, ion, level);
        auto &chlevel = chion.chlevels[level];
        for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
          chlevel.chphixstargets[phixstargetindex].corrphotoioncoeff = -99.;

#if (SEPARATE_STIMRECOMB)
          chlevel.chphixstargets[phixstargetindex].stimrecombcoeff = -99.;
#endif
        }

        chlevel.processrates[MA_ACTION_INTERNALUPHIGHER] = -99.;
      }
    }
  }

  if (modelgridindex >= 0) {
    std::ranges::fill(cacheslot.ch_allcont_departureratios, -1.);

    const auto nnetot = grid::get_nnetot(modelgridindex);
    for (int i = 0; i < globals::nbfcontinua; i++) {
      const int element = globals::allcont[i].element;
      const int ion = globals::allcont[i].ion;
      const int level = globals::allcont[i].level;
      const auto nnlevel = get_levelpop(modelgridindex, element, ion, level);
      cacheslot.ch_allcont_nnlevel[i] = nnlevel;
      cacheslot.ch_keep_this_cont[i] = nnlevel > 0 && keep_this_cont(element, ion, level, modelgridindex, nnetot);
    }
  }
}

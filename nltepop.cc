#include "nltepop.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <memory>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "gsl/gsl_cblas.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_permutation.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "ratecoeff.h"
#include "sn3d.h"

static FILE *nlte_file = nullptr;

// can save memory by using a combined rate matrix at the cost of diagnostic information
static constexpr bool individual_process_matricies = true;

static inline auto get_nlte_vector_index(const int element, const int ion, const int level) -> int
// this is the index for the NLTE solver that is handling all ions of a single element
// This is NOT an index into grid::modelgrid[modelgridindex].nlte_pops that contains all elements
{
  // have to convert from nlte_pops index to nlte_vector index
  // the difference is that nlte vectors apply to a single element and include ground states
  // The (+ ion) term accounts for the ground state population indicies that are not counted in the NLTE array
  const int gs_index =
      globals::elements[element].ions[ion].first_nlte - globals::elements[element].ions[0].first_nlte + ion;

  // add in level or superlevel number
  const int level_index = gs_index + (is_nlte(element, ion, level) ? level : (get_nlevels_nlte(element, ion) + 1));

  return level_index;
}

static void get_ion_level_of_nlte_vector_index(const int index, const int element, int *ion, int *level) {
  // this could easily be optimized if need be
  for (int dion = 0; dion < get_nions(element); dion++) {
    for (int dlevel = 0; dlevel < get_nlevels(element, dion); dlevel++) {
      if (get_nlte_vector_index(element, dion, dlevel) == index) {
        *ion = dion;
        *level = dlevel;
        return;
      }
    }
  }
}

static void eliminate_nlte_matrix_rowcol(const int index, const int gs_index, gsl_matrix *rate_matrix,
                                         gsl_vector *balance_vector) {
  const gsl_matrix rate_matrix_var = *rate_matrix;

  const int colcount = rate_matrix_var.size2;
  for (int column = 0; column < colcount; column++) {
    gsl_matrix_set(rate_matrix, index, column, 0.0);
  }

  const int rowcount = rate_matrix_var.size1;
  for (int row = 1; row < rowcount; row++) {
    gsl_matrix_set(rate_matrix, row, index, 0.0);
  }

  gsl_matrix_set(rate_matrix, index, gs_index, -1.0);
  gsl_matrix_set(rate_matrix, index, index, 1.0);
  gsl_vector_set(balance_vector, index, 0.0);
}

static void filter_nlte_matrix(const int element, gsl_matrix *rate_matrix, gsl_vector *balance_vector,
                               const gsl_vector * /*pop_norm_factor_vec*/)
// find rows and columns that barely interaction with other levels, and effectively
// removing them by zeroing their interactions and setting their departure
// coeff to 1.0
{
  const gsl_matrix rate_matrix_var = *rate_matrix;
  const int nlte_dimension = rate_matrix_var.size1;
  for (int index = 0; index < nlte_dimension; index++) {
    double row_max = 0.;
    for (int column = 0; column < nlte_dimension; column++) {
      const double element_value = fabs(gsl_matrix_get(rate_matrix, index, column));
      if (element_value > row_max) {
        row_max = element_value;
      }
    }
    double col_max = 0.;
    for (int row = 1; row < nlte_dimension; row++)  // skip the normalisation row 0
    {
      const double element_value = fabs(gsl_matrix_get(rate_matrix, row, index));
      if (element_value > col_max) {
        col_max = element_value;
      }
    }
    int ion = -1;
    int level = -1;
    get_ion_level_of_nlte_vector_index(index, element, &ion, &level);
    // printout("index%4d (ionstage%2d level%4d) row_max %.1e col_max %.1e ",
    //          index,get_ionstage(element,ion),level,row_max,col_max);

    if ((row_max < 1e-100) || (col_max < 1e-100)) {
      if (level == 0) {
        // printout("(Would eliminate but it's a ground state, so keeping it)");
        // printout("(Would eliminate but it's a ground state, so forcing pop=MINPOP=%g)",MINPOP);
        // gsl_vector_set(balance_vector, index, MINPOP / get_vector_get(pop_norm_factor_vec, index));
        // printout("(Eliminating this ground state)");
      } else {
        const int gs_index = get_nlte_vector_index(element, ion, 0);
        eliminate_nlte_matrix_rowcol(index, gs_index, rate_matrix, balance_vector);
        // printout("(forcing LTE population)");
      }
    }
    // printout("\n");
  }
}

static auto get_total_rate(const int index_selected, const gsl_matrix *rate_matrix, const gsl_vector *popvec,
                           const bool into_level, const bool only_levels_below, const bool only_levels_above)
    -> double {
  double total_rate = 0.;
  assert_always(!only_levels_below || !only_levels_above);

  if (into_level) {
    // find rate into selected level
    gsl_vector_const_view row_view = gsl_matrix_const_row(rate_matrix, index_selected);
    const gsl_vector *rates_vec = &row_view.vector;

    // multiply incoming rate coefficients by their corresponding populations to get rates
    if (!only_levels_above)  // add levels below
    {
      for (int index = 0; index < index_selected; index++) {
        total_rate += gsl_vector_get(rates_vec, index) * gsl_vector_get(popvec, index);
      }
    }

    if (!only_levels_below)  // add levels above
    {
      for (unsigned int index = index_selected + 1; index < rate_matrix->size1; index++) {
        total_rate += gsl_vector_get(rates_vec, index) * gsl_vector_get(popvec, index);
      }
    }
  } else {
    // find rate out of selected level
    gsl_vector_const_view col_view = gsl_matrix_const_column(rate_matrix, index_selected);
    const gsl_vector *rates_vec = &col_view.vector;

    // multiply outgoing rate coefficients by the population of the selected level to get rates

    if (!only_levels_above)  // add levels below
    {
      for (int index = 0; index < index_selected; index++) {
        total_rate += gsl_vector_get(rates_vec, index);
      }
    }

    if (!only_levels_below)  // add levels above
    {
      for (unsigned int index = index_selected + 1; index < rate_matrix->size2; index++) {
        total_rate += gsl_vector_get(rates_vec, index);
      }
    }

    total_rate *= gsl_vector_get(popvec, index_selected);
  }

  return total_rate;
}

static auto get_total_rate_in(const int index_to, const gsl_matrix *rate_matrix, const gsl_vector *popvec) -> double {
  return get_total_rate(index_to, rate_matrix, popvec, true, false, false);
}

static auto get_total_rate_out(const int index_from, const gsl_matrix *rate_matrix, const gsl_vector *popvec)
    -> double {
  return get_total_rate(index_from, rate_matrix, popvec, false, false, false);
}

static void print_level_rates_summary(const int element, const int selected_ion, const int selected_level,
                                      const gsl_vector *popvec, const gsl_matrix *rate_matrix_rad_bb,
                                      const gsl_matrix *rate_matrix_coll_bb, const gsl_matrix *rate_matrix_ntcoll_bb,
                                      const gsl_matrix *rate_matrix_rad_bf, const gsl_matrix *rate_matrix_coll_bf,
                                      const gsl_matrix *rate_matrix_ntcoll_bf) {
  const int selected_index = get_nlte_vector_index(element, selected_ion, selected_level);

  for (int i = 0; i <= 3; i++) {
    // rates in from below, in from above, out to below, out to above
    if (i == 0) {
      const int nlevels_nlte = get_nlevels_nlte(element, selected_ion);
      if (ion_has_superlevel(element, selected_ion) && (selected_level == nlevels_nlte + 1)) {
        printout("      superlevel ");
      } else {
        printout("    level%7d ", selected_level);
      }
      printout(" %10.2e ", gsl_vector_get(popvec, selected_index));
    } else {
      printout("                             ");
    }

    const bool into_level = (i <= 1);
    const bool only_levels_below = (i % 2) != 0;
    const bool only_levels_above = !only_levels_below;

    const double rad_bb_total =
        get_total_rate(selected_index, rate_matrix_rad_bb, popvec, into_level, only_levels_below, only_levels_above);
    const double coll_bb_total =
        get_total_rate(selected_index, rate_matrix_coll_bb, popvec, into_level, only_levels_below, only_levels_above);
    const double ntcoll_bb_total =
        get_total_rate(selected_index, rate_matrix_ntcoll_bb, popvec, into_level, only_levels_below, only_levels_above);
    const double rad_bf_total =
        get_total_rate(selected_index, rate_matrix_rad_bf, popvec, into_level, only_levels_below, only_levels_above);
    const double coll_bf_total =
        get_total_rate(selected_index, rate_matrix_coll_bf, popvec, into_level, only_levels_below, only_levels_above);
    const double ntcoll_bf_total =
        get_total_rate(selected_index, rate_matrix_ntcoll_bf, popvec, into_level, only_levels_below, only_levels_above);

    if (into_level) {
      // into this level
      if (only_levels_below) {
        printout(" from below ");
      } else {
        printout(" from above ");
      }
    } else {
      // out of this level
      if (only_levels_below) {
        printout("   to below ");
      } else {
        printout("   to above ");
      }
    }

    printout("%10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n", rad_bb_total, coll_bb_total, ntcoll_bb_total, rad_bf_total,
             coll_bf_total, ntcoll_bf_total);
  }
}

static void print_element_rates_summary(const int element, const int modelgridindex, const int timestep,
                                        const int nlte_iter, const gsl_vector *popvec,
                                        const gsl_matrix *rate_matrix_rad_bb, const gsl_matrix *rate_matrix_coll_bb,
                                        const gsl_matrix *rate_matrix_ntcoll_bb, const gsl_matrix *rate_matrix_rad_bf,
                                        const gsl_matrix *rate_matrix_coll_bf,
                                        const gsl_matrix *rate_matrix_ntcoll_bf) {
  const int nions = get_nions(element);
  for (int ion = 0; ion < nions; ion++) {
    const int nlevels = get_nlevels(element, ion);
    const int nlevels_nlte = get_nlevels_nlte(element, ion);

    const int atomic_number = get_atomicnumber(element);
    const int ionstage = get_ionstage(element, ion);

    const int max_printed_levels = ion_has_superlevel(element, ion) ? nlevels_nlte + 2 : nlevels_nlte + 1;

    for (int level = 0; (level < max_printed_levels) && (level < nlevels) && (level <= nlevels_nlte + 1); level++) {
      if (level == 0) {
        printout("  modelgridindex %d timestep %d NLTE iteration %d Te %g nne %g: NLTE summary for Z=%d ionstage %d:\n",
                 modelgridindex, timestep, nlte_iter, grid::get_Te(modelgridindex), grid::get_nne(modelgridindex),
                 atomic_number, ionstage);
        printout(
            "                         pop       rates     bb_rad     bb_col   bb_ntcol     bf_rad     bf_col   "
            "bf_ntcol\n");
      }

      print_level_rates_summary(element, ion, level, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb,
                                rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    }

    if (ion_has_superlevel(element, ion) && max_printed_levels < (nlevels_nlte + 1)) {
      const int level_superlevel = nlevels_nlte + 1;

      print_level_rates_summary(element, ion, level_superlevel, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb,
                                rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    }
  }
}

static void print_level_rates(const int modelgridindex, const int timestep, const int element, const int selected_ion,
                              const int selected_level, const gsl_vector *popvec, const gsl_matrix *rate_matrix_rad_bb,
                              const gsl_matrix *rate_matrix_coll_bb, const gsl_matrix *rate_matrix_ntcoll_bb,
                              const gsl_matrix *rate_matrix_rad_bf, const gsl_matrix *rate_matrix_coll_bf,
                              const gsl_matrix *rate_matrix_ntcoll_bf) {
  // very detailed output of the NLTE processes for a particular levels

  if (element > get_nelements() - 1 || selected_ion > get_nions(element) - 1 ||
      selected_level >
          (get_nlevels_nlte(element, selected_ion) + (ion_has_superlevel(element, selected_ion) ? 1 : 0))) {
    printout("print_level_rates: invalid element/ion/level arguments\n");
    std::abort();
  }

  if (rate_matrix_rad_bb == rate_matrix_coll_bb) {
    printout(
        "print_level_rates: rate_matrix_rad_bb == rate_matrix_coll_bb. check individual_process_matricies is off\n");
    std::abort();
  }

  const gsl_vector popvector = *popvec;
  const int nlte_dimension = popvector.size;
  const int atomic_number = get_atomicnumber(element);
  const int selected_ionstage = get_ionstage(element, selected_ion);
  const int selected_index = get_nlte_vector_index(element, selected_ion, selected_level);
  const double pop_selectedlevel = gsl_vector_get(popvec, selected_index);
  printout(
      "timestep %d cell %d Te %g nne %g NLTE level diagnostics for Z=%d ionstage %d level %d rates into and out of "
      "this level\n",
      timestep, modelgridindex, grid::get_Te(modelgridindex), grid::get_nne(modelgridindex), atomic_number,
      selected_ionstage, selected_level);

  const double rad_bb_in_total = get_total_rate_in(selected_index, rate_matrix_rad_bb, popvec);
  const double coll_bb_in_total = get_total_rate_in(selected_index, rate_matrix_coll_bb, popvec);
  const double ntcoll_bb_in_total = get_total_rate_in(selected_index, rate_matrix_ntcoll_bb, popvec);
  const double rad_bf_in_total = get_total_rate_in(selected_index, rate_matrix_rad_bf, popvec);
  const double coll_bf_in_total = get_total_rate_in(selected_index, rate_matrix_coll_bf, popvec);
  const double ntcoll_bf_in_total = get_total_rate_in(selected_index, rate_matrix_ntcoll_bf, popvec);
  const double total_rate_in =
      rad_bb_in_total + coll_bb_in_total + rad_bf_in_total + coll_bf_in_total + ntcoll_bf_in_total;
  printout(
      "  TOTAL rates in:             rad_bb_in  %8.2e coll_bb_in  %8.2e ntcoll_bb_in  %8.2e rad_bf_in  %8.2e "
      "coll_bf_in  %8.2e ntcoll_bf_in  %8.2e\n",
      rad_bb_in_total, coll_bb_in_total, ntcoll_bb_in_total, rad_bf_in_total, coll_bf_in_total, ntcoll_bf_in_total);

  const double rad_bb_out_total = get_total_rate_out(selected_index, rate_matrix_rad_bb, popvec);
  const double coll_bb_out_total = get_total_rate_out(selected_index, rate_matrix_coll_bb, popvec);
  const double ntcoll_bb_out_total = get_total_rate_out(selected_index, rate_matrix_ntcoll_bb, popvec);
  const double rad_bf_out_total = get_total_rate_out(selected_index, rate_matrix_rad_bf, popvec);
  const double coll_bf_out_total = get_total_rate_out(selected_index, rate_matrix_coll_bf, popvec);
  const double ntcoll_bf_out_total = get_total_rate_out(selected_index, rate_matrix_ntcoll_bf, popvec);
  const double total_rate_out =
      rad_bb_out_total + coll_bb_out_total + rad_bf_out_total + coll_bf_out_total + ntcoll_bf_out_total;
  printout(
      "  TOTAL rates out:            rad_bb_out %8.2e coll_bb_out %8.2e ntcoll_bb_out %8.2e rad_bf_out %8.2e "
      "coll_bf_out %8.2e ntcoll_bf_out %8.2e\n",
      rad_bb_out_total, coll_bb_out_total, ntcoll_bb_out_total, rad_bf_out_total, coll_bf_out_total,
      ntcoll_bf_out_total);

  for (int index = 0; index < nlte_dimension; index++) {
    if (index == selected_index) {
      continue;
    }
    int ion = -1;
    int level = -1;
    get_ion_level_of_nlte_vector_index(index, element, &ion, &level);
    const int ionstage = get_ionstage(element, ion);
    // in means populating the selected level, out means depopulating the selected level
    const double pop = gsl_vector_get(popvec, index);
    const double rad_bb_in = gsl_matrix_get(rate_matrix_rad_bb, selected_index, index) * pop;
    const double rad_bb_out = gsl_matrix_get(rate_matrix_rad_bb, index, selected_index) * pop_selectedlevel;
    const double coll_bb_in = gsl_matrix_get(rate_matrix_coll_bb, selected_index, index) * pop;
    const double coll_bb_out = gsl_matrix_get(rate_matrix_coll_bb, index, selected_index) * pop_selectedlevel;
    const double ntcoll_bb_in = gsl_matrix_get(rate_matrix_ntcoll_bb, selected_index, index) * pop;
    const double ntcoll_bb_out = gsl_matrix_get(rate_matrix_ntcoll_bb, index, selected_index) * pop_selectedlevel;
    const double rad_bf_in = gsl_matrix_get(rate_matrix_rad_bf, selected_index, index) * pop;
    const double rad_bf_out = gsl_matrix_get(rate_matrix_rad_bf, index, selected_index) * pop_selectedlevel;
    const double coll_bf_in = gsl_matrix_get(rate_matrix_coll_bf, selected_index, index) * pop;
    const double coll_bf_out = gsl_matrix_get(rate_matrix_coll_bf, index, selected_index) * pop_selectedlevel;
    const double ntcoll_bf_in = gsl_matrix_get(rate_matrix_ntcoll_bf, selected_index, index) * pop;
    const double ntcoll_bf_out = gsl_matrix_get(rate_matrix_ntcoll_bf, index, selected_index) * pop_selectedlevel;

    const bool nonzero_rate_in = (fabs(rad_bb_in) > 0. || fabs(coll_bb_in) > 0. || fabs(ntcoll_bb_in) > 0. ||
                                  fabs(rad_bf_in) > 0. || fabs(coll_bf_in) > 0. || fabs(ntcoll_bf_in) > 0.);
    const bool nonzero_rate_out = (fabs(rad_bb_out) > 0. || fabs(coll_bb_out) > 0. || fabs(ntcoll_bb_out) > 0. ||
                                   fabs(rad_bf_out) > 0. || fabs(coll_bf_out) > 0. || fabs(ntcoll_bf_out) > 0.);
    if (nonzero_rate_in || nonzero_rate_out) {
      const double epsilon_trans = fabs(epsilon(element, ion, level) - epsilon(element, selected_ion, selected_level));
      const double nu_trans = epsilon_trans / H;
      const double lambda = 1e8 * CLIGHT / nu_trans;  // should be in Angstroms
      const double level_rate_in = rad_bb_in + coll_bb_in + ntcoll_bb_in + rad_bf_in + coll_bf_in + ntcoll_bf_in;
      const double level_rate_out = rad_bb_out + coll_bb_out + ntcoll_bb_in + rad_bf_out + coll_bf_out + ntcoll_bf_out;
      const double level_percent_in = level_rate_in / total_rate_in * 100.;
      const double level_percent_out = level_rate_out / total_rate_out * 100.;

      printout(
          "  ionstage %d level %4d (%5.1f%% of in)  rad_bb_in  %8.2e coll_bb_in  %8.2e ntcoll_bb_in  %8.2e rad_bf_in  "
          "%8.2e coll_bf_in  %8.2e ntcoll_bf_in  %8.2e lambda %6.0f\n",
          ionstage, level, level_percent_in, rad_bb_in, coll_bb_in, ntcoll_bb_in, rad_bf_in, coll_bf_in, ntcoll_bf_in,
          lambda);
      printout(
          "  ionstage %d level %4d (%5.1f%% of out) rad_bb_out %8.2e coll_bb_out %8.2e ntcoll_bb_out %8.2e rad_bf_out "
          "%8.2e coll_bf_out %8.2e ntcoll_bf_out %8.2e lambda %6.0f\n",
          ionstage, level, level_percent_out, rad_bb_out, coll_bb_out, ntcoll_bb_out, rad_bf_out, coll_bf_out,
          ntcoll_bf_out, lambda);
    }
  }
  printout("\n");
}

static void nltepop_reset_element(const int modelgridindex, const int element) {
  const int nions = get_nions(element);
  for (int ion = 0; ion < nions; ion++) {
    const int nlte_start = globals::elements[element].ions[ion].first_nlte;
    std::fill_n(&grid::modelgrid[modelgridindex].nlte_pops[nlte_start],
                get_nlevels_nlte(element, ion) + (ion_has_superlevel(element, ion) ? 1 : 0), -1.);
  }
}

static auto get_element_superlevelpartfuncs(const int modelgridindex, const int element) -> std::vector<double> {
  const int nions = get_nions(element);
  std::vector<double> superlevel_partfuncs(nions, 0.);
  for (int ion = 0; ion < nions; ion++) {
    if (ion_has_superlevel(element, ion)) {
      const int nlevels_nlte = get_nlevels_nlte(element, ion);
      const int nlevels = get_nlevels(element, ion);
      for (int level = nlevels_nlte + 1; level < nlevels; level++) {
        superlevel_partfuncs[ion] += superlevel_boltzmann(modelgridindex, element, ion, level);
      }
    } else {
      superlevel_partfuncs[ion] = 0.;
    }
  }

  return superlevel_partfuncs;
}

static auto get_element_nlte_dimension(const int element) -> int {
  int nlte_dimension = 0;
  const int nions = get_nions(element);
  for (int ion = 0; ion < nions; ion++) {
    const int nlevels_nlte = get_nlevels_nlte(element, ion);

    nlte_dimension += nlevels_nlte + 1;  // ground state is not counted in nlevels_nlte

    // add super level if it exists
    if (ion_has_superlevel(element, ion)) {
      nlte_dimension++;
    }
  }

  return nlte_dimension;
}

static void nltepop_matrix_add_boundbound(const int modelgridindex, const int element, const int ion,
                                          const double t_mid, const std::vector<double> &s_renorm,
                                          gsl_matrix *rate_matrix_rad_bb, gsl_matrix *rate_matrix_coll_bb,
                                          gsl_matrix *rate_matrix_ntcoll_bb) {
  const auto T_e = grid::get_Te(modelgridindex);
  const float nne = grid::get_nne(modelgridindex);
  const int nlevels = get_nlevels(element, ion);
  for (int level = 0; level < nlevels; level++) {
    const int level_index = get_nlte_vector_index(element, ion, level);
    const double epsilon_level = epsilon(element, ion, level);
    const double statweight = stat_weight(element, ion, level);

    // de-excitation
    const int ndowntrans = get_ndowntrans(element, ion, level);
    for (int i = 0; i < ndowntrans; i++) {
      const auto &downtransition = globals::elements[element].ions[ion].levels[level].downtrans[i];
      const double A_ul = downtransition.einstein_A;
      const int lower = downtransition.targetlevelindex;

      const double epsilon_trans = epsilon_level - epsilon(element, ion, lower);

      const double R = rad_deexcitation_ratecoeff(modelgridindex, element, ion, level, lower, epsilon_trans, A_ul,
                                                  statweight, t_mid) *
                       s_renorm[level];
      const double C =
          col_deexcitation_ratecoeff(T_e, nne, epsilon_trans, element, ion, level, downtransition) * s_renorm[level];

      const int upper_index = level_index;
      const int lower_index = get_nlte_vector_index(element, ion, lower);

      *gsl_matrix_ptr(rate_matrix_rad_bb, upper_index, upper_index) -= R;
      *gsl_matrix_ptr(rate_matrix_rad_bb, lower_index, upper_index) += R;
      *gsl_matrix_ptr(rate_matrix_coll_bb, upper_index, upper_index) -= C;
      *gsl_matrix_ptr(rate_matrix_coll_bb, lower_index, upper_index) += C;
      if ((R < 0) || (C < 0)) {
        printout("  WARNING: Negative de-excitation rate from ionstage %d level %d to level %d\n",
                 get_ionstage(element, ion), level, lower);
      }
    }

    // excitation
    const int nuptrans = get_nuptrans(element, ion, level);
    for (int i = 0; i < nuptrans; i++) {
      const int lineindex = globals::elements[element].ions[ion].levels[level].uptrans[i].lineindex;
      const TransitionLine *line = &globals::linelist[lineindex];
      const int upper = line->upperlevelindex;
      const double epsilon_trans = epsilon(element, ion, upper) - epsilon_level;

      const double R =
          rad_excitation_ratecoeff(modelgridindex, element, ion, level, i, epsilon_trans, lineindex, t_mid) *
          s_renorm[level];
      assert_always(R >= 0);
      assert_always(std::isfinite(R));

      const double C =
          col_excitation_ratecoeff(T_e, nne, element, ion, level, i, epsilon_trans, statweight) * s_renorm[level];
      assert_always(C >= 0);
      assert_always(std::isfinite(C));

      const double NTC =
          nonthermal::nt_excitation_ratecoeff(modelgridindex, element, ion, level, i, epsilon_trans, lineindex) *
          s_renorm[level];

      // if ((Z == 26) && (ionstage == 1) && (level == 0) && (upper <= 5))
      // {
      //   const double tau_sobolev = get_tau_sobolev(modelgridindex, lineindex, t_mid);
      //   const double nu_trans = epsilon_trans / H;
      //   const double lambda = 1e8 * CLIGHT / nu_trans; // should be in Angstroms
      //   printout("Z=%d ionstage %d lower %d upper %d lambda %6.1fÅ tau_sobolev=%g einstein_A %g osc_strength %g
      //   coll_str %g\n",
      //            Z, ionstage, level, upper, lambda, tau_sobolev,
      //            linelist[lineindex].einstein_A, linelist[lineindex].osc_strength, linelist[lineindex].coll_str);
      // }

      const int lower_index = level_index;
      const int upper_index = get_nlte_vector_index(element, ion, upper);

      *gsl_matrix_ptr(rate_matrix_rad_bb, lower_index, lower_index) -= R;
      *gsl_matrix_ptr(rate_matrix_rad_bb, upper_index, lower_index) += R;
      *gsl_matrix_ptr(rate_matrix_coll_bb, lower_index, lower_index) -= C;
      *gsl_matrix_ptr(rate_matrix_coll_bb, upper_index, lower_index) += C;
      *gsl_matrix_ptr(rate_matrix_ntcoll_bb, lower_index, lower_index) -= NTC;
      *gsl_matrix_ptr(rate_matrix_ntcoll_bb, upper_index, lower_index) += NTC;
      if ((R < 0) || (C < 0)) {
        printout("  WARNING: Negative excitation rate from ion %d level %d to level %d\n", get_ionstage(element, ion),
                 level, upper);
      }
    }
  }
}

static void nltepop_matrix_add_ionisation(const int modelgridindex, const int element, const int ion,
                                          const std::vector<double> &s_renorm, gsl_matrix *rate_matrix_rad_bf,
                                          gsl_matrix *rate_matrix_coll_bf) {
  assert_always(ion + 1 < get_nions(element));  // can't ionise the top ion
  const auto T_e = grid::get_Te(modelgridindex);
  const float nne = grid::get_nne(modelgridindex);
  const int nionisinglevels = get_ionisinglevels(element, ion);
  const int maxrecombininglevel = get_maxrecombininglevel(element, ion + 1);

  for (int level = 0; level < nionisinglevels; level++) {
    const int lower_index = get_nlte_vector_index(element, ion, level);

    // thermal collisional ionization, photoionisation and recombination processes
    const double epsilon_current = epsilon(element, ion, level);

    for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++) {
      const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
      const int upper_index = get_nlte_vector_index(element, ion + 1, upper);
      const double epsilon_trans = epsilon(element, ion + 1, upper) - epsilon_current;

      // ionization

      // the R part is slow!
      const double R_ionisation = get_corrphotoioncoeff(element, ion, level, phixstargetindex, modelgridindex);
      const double C_ionisation =
          col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);

      *gsl_matrix_ptr(rate_matrix_rad_bf, lower_index, lower_index) -= R_ionisation * s_renorm[level];
      *gsl_matrix_ptr(rate_matrix_rad_bf, upper_index, lower_index) += R_ionisation * s_renorm[level];
      *gsl_matrix_ptr(rate_matrix_coll_bf, lower_index, lower_index) -= C_ionisation * s_renorm[level];
      *gsl_matrix_ptr(rate_matrix_coll_bf, upper_index, lower_index) += C_ionisation * s_renorm[level];

      if ((R_ionisation < 0) || (C_ionisation < 0)) {
        printout("  WARNING: Negative ionization rate from ionstage %d level %d phixstargetindex %d\n",
                 get_ionstage(element, ion), level, phixstargetindex);
      }

      // recombination
      if (upper <= maxrecombininglevel)  // we can skip this part if the functions below will return zero anyway
      {
        const double R_recomb = rad_recombination_ratecoeff(T_e, nne, element, ion + 1, upper, level, modelgridindex);
        const double C_recomb =
            col_recombination_ratecoeff(modelgridindex, element, ion + 1, upper, level, epsilon_trans);

        *gsl_matrix_ptr(rate_matrix_rad_bf, upper_index, upper_index) -= R_recomb * s_renorm[upper];
        *gsl_matrix_ptr(rate_matrix_rad_bf, lower_index, upper_index) += R_recomb * s_renorm[upper];
        *gsl_matrix_ptr(rate_matrix_coll_bf, upper_index, upper_index) -= C_recomb * s_renorm[upper];
        *gsl_matrix_ptr(rate_matrix_coll_bf, lower_index, upper_index) += C_recomb * s_renorm[upper];

        if ((R_recomb < 0) || (C_recomb < 0)) {
          printout("  WARNING: Negative recombination rate to ionstage %d level %d phixstargetindex %d\n",
                   get_ionstage(element, ion), level, phixstargetindex);
        }
      }
    }
  }
}

static void nltepop_matrix_add_nt_ionisation(const int modelgridindex, const int element, const int ion,
                                             const std::vector<double> &s_renorm, gsl_matrix *rate_matrix_ntcoll_bf) {
  // collisional ionization by non-thermal electrons

  assert_always(ion + 1 < get_nions(element));  // can't ionise the top ion
  const double Y_nt = nonthermal::nt_ionization_ratecoeff(modelgridindex, element, ion);
  if (Y_nt < 0.) {
    printout("  WARNING: Negative NT_ionization rate from ionstage %d\n", get_ionstage(element, ion));
  }

  const int nlevels = get_nlevels(element, ion);

  for (int upperion = ion + 1; upperion <= nonthermal::nt_ionisation_maxupperion(element, ion); upperion++) {
    const double Y_nt_thisupperion =
        Y_nt * nonthermal::nt_ionization_upperion_probability(modelgridindex, element, ion, upperion, false);

    if (Y_nt_thisupperion > 0.) {
      const int upper_groundstate_index = get_nlte_vector_index(element, upperion, 0);
      for (int level = 0; level < nlevels; level++) {
        const int lower_index = get_nlte_vector_index(element, ion, level);

        *gsl_matrix_ptr(rate_matrix_ntcoll_bf, lower_index, lower_index) -= Y_nt_thisupperion * s_renorm[level];
        *gsl_matrix_ptr(rate_matrix_ntcoll_bf, upper_groundstate_index, lower_index) +=
            Y_nt_thisupperion * s_renorm[level];
      }
    }
  }
}

static void nltepop_matrix_normalise(const int modelgridindex, const int element, gsl_matrix *rate_matrix,
                                     gsl_vector *pop_norm_factor_vec) {
  const size_t nlte_dimension = pop_norm_factor_vec->size;
  assert_always(pop_norm_factor_vec->size == nlte_dimension);
  assert_always(rate_matrix->size1 == nlte_dimension);
  assert_always(rate_matrix->size2 == nlte_dimension);

  // TODO: consider replacing normalisation by LTE populations with
  // GSL's gsl_linalg_balance_matrix(gsl_matrix * A, gsl_vector * D) function instead
  for (size_t column = 0; column < nlte_dimension; column++) {
    int ion = -1;
    int level = -1;
    get_ion_level_of_nlte_vector_index(column, element, &ion, &level);

    gsl_vector_set(pop_norm_factor_vec, column, calculate_levelpop_lte(modelgridindex, element, ion, level));

    if ((level != 0) && (!is_nlte(element, ion, level))) {
      // level is a superlevel, so add populations of higher levels to the norm factor
      for (int dummylevel = level + 1; dummylevel < get_nlevels(element, ion); dummylevel++) {
        if (!is_nlte(element, ion, dummylevel)) {
          *gsl_vector_ptr(pop_norm_factor_vec, column) +=
              calculate_levelpop_lte(modelgridindex, element, ion, dummylevel);
        }
      }
      // NOTE: above calculation is not always equal to the sum of LTE populations
      // since calculate_levelpop_lte imposes MINPOP minimum
      // printout("superlevel norm factor index %d is %g, partfunc is %g, partfunc*levelpop(SL)/g(SL) %g\n",
      //          column, gsl_vector_get(pop_norm_factor_vec, column), superlevel_partfunc[ion],
      //          superlevel_partfunc[ion] * calculate_levelpop_lte(modelgridindex,element,ion,level) /
      //          stat_weight(element,ion,level));
    }

    // apply the normalisation factor to this column in the rate_matrix
    gsl_vector_view column_view = gsl_matrix_column(rate_matrix, column);
    gsl_vector_scale(&column_view.vector, gsl_vector_get(pop_norm_factor_vec, column));
  }
}

static void set_element_pops_lte(const int modelgridindex, const int element) {
  nltepop_reset_element(modelgridindex, element);  // set NLTE pops as invalid so that LTE pops will be used instead

  calculate_cellpartfuncts(modelgridindex, element);
  set_groundlevelpops(modelgridindex, element, grid::get_nne(modelgridindex), true);
}

static auto lumatrix_is_singular(const gsl_matrix *LU, const int element) -> bool {
  size_t const n = LU->size1;
  bool is_singular = false;

  for (size_t i = 0; i < n; i++) {
    const double u = gsl_matrix_get(LU, i, i);
    if (u == 0) {
      int ion = -1;
      int level = -1;
      get_ion_level_of_nlte_vector_index(i, element, &ion, &level);
      if (is_nlte(element, ion, level)) {
        printout("NLTE disconnected level: Z=%d ionstage %d level %d\n", get_atomicnumber(element),
                 get_ionstage(element, ion), level);
      } else {
        printout("NLTE disconnected superlevel: Z=%d ionstage %d\n", get_atomicnumber(element),
                 get_ionstage(element, ion));
      }
      is_singular = true;
    }
  }

  return is_singular;
}

static auto nltepop_matrix_solve(const int element, const gsl_matrix *rate_matrix, const gsl_vector *balance_vector,
                                 gsl_vector *popvec, const gsl_vector *pop_normfactor_vec) -> bool
// solve rate_matrix * x = balance_vector,
// then popvec[i] = x[i] / pop_norm_factor_vec[i]
{
  bool completed_solution = false;
  const size_t nlte_dimension = balance_vector->size;
  assert_always(pop_normfactor_vec->size == nlte_dimension);
  assert_always(rate_matrix->size1 == nlte_dimension);
  assert_always(rate_matrix->size2 == nlte_dimension);

  gsl_vector *x = gsl_vector_alloc(nlte_dimension);  // population solution vector (normalised)

  // make a copy of the rate matrix for the LU decomp
  gsl_matrix *rate_matrix_LU_decomp = gsl_matrix_alloc(nlte_dimension, nlte_dimension);
  gsl_matrix_memcpy(rate_matrix_LU_decomp, rate_matrix);

  gsl_permutation *p = gsl_permutation_alloc(nlte_dimension);

  int s = 0;  // sign of the transformation
  gsl_linalg_LU_decomp(rate_matrix_LU_decomp, p, &s);

  if (lumatrix_is_singular(rate_matrix_LU_decomp, element)) {
    printout("ERROR: NLTE matrix is singular for element Z=%d!\n", get_atomicnumber(element));
    // abort();
    completed_solution = false;
  } else {
    gsl_error_handler_t *previous_handler = gsl_set_error_handler(gsl_error_handler_printout);

    // solve matrix equation: rate_matrix * x = balance_vector for x (population vector)
    gsl_linalg_LU_solve(rate_matrix_LU_decomp, p, balance_vector, x);

    gsl_set_error_handler(previous_handler);

    // gsl_linalg_HH_solve (&m.matrix, &b.vector, x);

    const double TOLERANCE = 1e-40;

    gsl_vector *gsl_work_vector = gsl_vector_alloc(nlte_dimension);
    double error_best = -1.;
    gsl_vector *x_best = gsl_vector_alloc(nlte_dimension);  // population solution vector with lowest error
    gsl_vector *residual_vector = gsl_vector_alloc(nlte_dimension);
    int iteration = 0;
    for (iteration = 0; iteration < 10; iteration++) {
      if (iteration > 0) {
        gsl_linalg_LU_refine(rate_matrix, rate_matrix_LU_decomp, p, balance_vector, x, gsl_work_vector);
      }

      gsl_vector_memcpy(residual_vector, balance_vector);
      gsl_blas_dgemv(CblasNoTrans, 1.0, rate_matrix, x, -1.0, residual_vector);  // calculate Ax - b = residual
      const double error = fabs(
          gsl_vector_get(residual_vector, gsl_blas_idamax(residual_vector)));  // value of the largest absolute residual

      if (error < error_best || error_best < 0.) {
        gsl_vector_memcpy(x_best, x);
        error_best = error;
      }
      // printout("Linear algebra solver iteration %d has a maximum residual of %g\n",iteration,error);
      if (error < TOLERANCE) {
        break;
      }
    }
    if (error_best >= 0.) {
      // printout("  NLTE solver matrix LU_refine: After %d iterations, keeping solution vector with a max residual of
      // %g\n",iteration,error_best);
      if (error_best > 1e-10) {
        printout(
            "  NLTE solver matrix LU_refine: After %d iterations, best solution vector has a max residual of %g "
            "(WARNING!)\n",
            iteration, error_best);
      }

      gsl_vector_memcpy(x, x_best);
    }

    gsl_vector_free(x_best);
    gsl_vector_free(gsl_work_vector);

    // get the real populations using the x vector and the normalisation factors
    gsl_vector_memcpy(popvec, x);
    gsl_vector_mul(popvec, pop_normfactor_vec);
    // popvec will be used contains the real population densities

    for (size_t row = 0; row < nlte_dimension; row++) {
      double recovered_balance_vector_elem = 0.;
      gsl_vector_const_view row_view = gsl_matrix_const_row(rate_matrix, row);
      gsl_blas_ddot(&row_view.vector, x, &recovered_balance_vector_elem);

      int ion = 0;
      int level = 0;
      get_ion_level_of_nlte_vector_index(row, element, &ion, &level);

      // printout("index %4d (ionstage %d level%4d): residual %+.2e recovered balance: %+.2e normed pop %.2e pop %.2e
      // departure ratio %.4f\n",
      //          row,get_ionstage(element,ion),level, gsl_vector_get(residual_vector,row),
      //          recovered_balance_vector_elem, gsl_vector_get(x,row),
      //          gsl_vector_get(popvec, row),
      //          gsl_vector_get(x, row) / gsl_vector_get(x,get_nlte_vector_index(element,ion,0)));

      if (gsl_vector_get(popvec, row) < 0.0) {
        printout(
            "  WARNING: NLTE solver gave negative population to index %zud (Z=%d ionstage %d level %d), pop = %g. "
            "Replacing with LTE pop of %g\n",
            row, get_atomicnumber(element), get_ionstage(element, ion), level,
            gsl_vector_get(x, row) * gsl_vector_get(pop_normfactor_vec, row), gsl_vector_get(pop_normfactor_vec, row));
        gsl_vector_set(popvec, row, gsl_vector_get(pop_normfactor_vec, row));
      }
    }

    gsl_vector_free(residual_vector);
    completed_solution = true;
  }

  gsl_vector_free(x);
  gsl_matrix_free(rate_matrix_LU_decomp);
  gsl_permutation_free(p);

  return completed_solution;
}

void solve_nlte_pops_element(const int element, const int modelgridindex, const int timestep, const int nlte_iter)
// solves the statistical balance equations to find NLTE level populations for all ions of an element
// (ionisation balance follows from this too)
{
  const int atomic_number = get_atomicnumber(element);

  if (grid::get_elem_abundance(modelgridindex, element) <= 0.) {
    // abundance of this element is zero, so do not store any NLTE populations
    printout("Not solving for NLTE populations in cell %d at timestep %d for element Z=%d due to zero abundance\n",
             modelgridindex, timestep, atomic_number);

    nltepop_reset_element(modelgridindex, element);
    return;
  }

  const auto sys_time_start_nltesolver = std::time(nullptr);

  const double t_mid = globals::timesteps[timestep].mid;
  const int nions = get_nions(element);
  const double nnelement = grid::get_elem_numberdens(modelgridindex, element);

  printout(
      "Solving for NLTE populations in cell %d at timestep %d NLTE iteration %d for element Z=%d (mass fraction %.2e, "
      "nnelement %.2e cm^-3)\n",
      modelgridindex, timestep, nlte_iter, atomic_number, grid::get_elem_abundance(modelgridindex, element), nnelement);

  auto superlevel_partfunc = std::vector<double>(nions) = get_element_superlevelpartfuncs(modelgridindex, element);
  const int nlte_dimension = get_element_nlte_dimension(element);

  // printout("NLTE: the vector dimension is %d", nlte_dimension);

  gsl_matrix *rate_matrix = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
  gsl_matrix *rate_matrix_rad_bb = nullptr;
  gsl_matrix *rate_matrix_coll_bb = nullptr;
  gsl_matrix *rate_matrix_ntcoll_bb = nullptr;
  gsl_matrix *rate_matrix_rad_bf = nullptr;
  gsl_matrix *rate_matrix_coll_bf = nullptr;
  gsl_matrix *rate_matrix_ntcoll_bf = nullptr;
  if (individual_process_matricies) {
    rate_matrix_rad_bb = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    rate_matrix_coll_bb = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    rate_matrix_ntcoll_bb = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    rate_matrix_rad_bf = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    rate_matrix_coll_bf = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    rate_matrix_ntcoll_bf = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
  } else {
    // alias the single matrix accounting for all processes
    rate_matrix_rad_bb = rate_matrix;
    rate_matrix_coll_bb = rate_matrix;
    rate_matrix_ntcoll_bb = rate_matrix;
    rate_matrix_rad_bf = rate_matrix;
    rate_matrix_coll_bf = rate_matrix;
    rate_matrix_ntcoll_bf = rate_matrix;
  }

  gsl_vector *const balance_vector = gsl_vector_calloc(nlte_dimension);

  assert_always(balance_vector != nullptr);

  // printout("  Adding rates for ion stages:");
  for (int ion = 0; ion < nions; ion++) {
    // const int ionstage = get_ionstage(element, ion);
    // printout(" %d", ionstage);

    const int nlevels = get_nlevels(element, ion);
    const int nlevels_nlte = get_nlevels_nlte(element, ion);

    auto s_renorm = std::vector<double>(nlevels);
    for (int level = 0; level <= nlevels_nlte; level++) {
      s_renorm[level] = 1.;
    }

    for (int level = (nlevels_nlte + 1); level < nlevels; level++) {
      s_renorm[level] = superlevel_boltzmann(modelgridindex, element, ion, level) / superlevel_partfunc[ion];
    }

    nltepop_matrix_add_boundbound(modelgridindex, element, ion, t_mid, s_renorm, rate_matrix_rad_bb,
                                  rate_matrix_coll_bb, rate_matrix_ntcoll_bb);

    if (ion < nions - 1) {
      // this is the slowest component
      nltepop_matrix_add_ionisation(modelgridindex, element, ion, s_renorm, rate_matrix_rad_bf, rate_matrix_coll_bf);
      if (NT_ON) {
        nltepop_matrix_add_nt_ionisation(modelgridindex, element, ion, s_renorm, rate_matrix_ntcoll_bf);
      }
    }
  }
  // printout("\n");

  if (individual_process_matricies) {
    // sum the matricies for each transition process to get a total rate matrix
    gsl_matrix_add(rate_matrix, rate_matrix_rad_bb);
    gsl_matrix_add(rate_matrix, rate_matrix_coll_bb);
    gsl_matrix_add(rate_matrix, rate_matrix_ntcoll_bb);
    gsl_matrix_add(rate_matrix, rate_matrix_rad_bf);
    gsl_matrix_add(rate_matrix, rate_matrix_coll_bf);
    gsl_matrix_add(rate_matrix, rate_matrix_ntcoll_bf);
  }

  // replace the first row of the matrix and balance vector with the normalisation
  // constraint on the total element population
  gsl_vector_view first_row_view = gsl_matrix_row(rate_matrix, 0);
  gsl_vector_set_all(&first_row_view.vector, 1.0);
  // set first balance vector entry to the element population (all other entries will be zero)
  gsl_vector_set(balance_vector, 0, nnelement);

  if (FORCE_SAHA_ION_BALANCE(atomic_number)) {
    const auto ionfractions = calculate_ionfractions(element, modelgridindex, grid::get_nne(modelgridindex), true);
    const int uppermost_ion = static_cast<int>(ionfractions.size() - 1);
    for (int ion = 1; ion <= uppermost_ion; ion++) {
      // replace matrix row for ion's ground state with sum of this ion's level populations is equal to the ion
      // population
      const double nnion = nnelement * ionfractions[ion];
      const int index_ion_ground = get_nlte_vector_index(element, ion, 0);
      const int index_ion_toplevel = get_nlte_vector_index(element, ion, get_nlevels(element, ion));
      gsl_vector_view ion_ground_row_view = gsl_matrix_row(rate_matrix, index_ion_ground);
      gsl_vector_set_all(&ion_ground_row_view.vector, 0.);
      for (int index = index_ion_ground; index <= index_ion_toplevel; index++) {
        gsl_vector_set(&ion_ground_row_view.vector, index, 1.);
      }

      gsl_vector_set(balance_vector, get_nlte_vector_index(element, ion, index_ion_ground), nnion);
    }
  }

  // calculate the normalisation factors and apply them to the matrix
  // columns and balance vector elements
  gsl_vector *pop_norm_factor_vec = gsl_vector_calloc(nlte_dimension);
  nltepop_matrix_normalise(modelgridindex, element, rate_matrix, pop_norm_factor_vec);

  // printout("Rate matrix | balance vector:\n");
  // for (int row = 0; row < nlte_dimension; row++)
  // {
  //   for (int column = 0; column < nlte_dimension; column++)
  //   {
  //     char str[15];
  //     snprintf(str, 15, "%+.1e ", gsl_matrix_get(rate_matrix, row, column));
  //     printout(str);
  //   }
  //   printout("| ");
  //   char str[15];
  //   snprintf(str, 15, "%+.1e\n", gsl_vector_get(balance_vector, row));
  //   printout(str);
  // }
  // printout("\n");

  // eliminate barely-interacting levels from the NLTE matrix by removing
  // their interactions and setting their normalised populations (probably departure coeff) to 1.0
  // filter_nlte_matrix(element, rate_matrix, balance_vector, pop_norm_factor_vec);

  gsl_vector *popvec = gsl_vector_alloc(nlte_dimension);  // the true population densities

  const bool matrix_solve_success =
      nltepop_matrix_solve(element, rate_matrix, balance_vector, popvec, pop_norm_factor_vec);

  if (!matrix_solve_success) {
    printout(
        "WARNING: Can't solve for NLTE populations in cell %d at timestep %d for element Z=%d due to singular matrix. "
        "Attempting to use LTE solution instead\n",
        modelgridindex, timestep, atomic_number);
    set_element_pops_lte(modelgridindex, element);
  } else {
    // check calculated NLTE populations are valid
    for (int index = 0; index < nlte_dimension; index++) {
      assert_always(std::isfinite(gsl_vector_get(popvec, index)));
      assert_always(gsl_vector_get(popvec, index) >= 0.);
    }

    for (int ion = 0; ion < nions; ion++) {
      const int nlevels_nlte = get_nlevels_nlte(element, ion);
      const int index_gs = get_nlte_vector_index(element, ion, 0);
      // const int ionstage = get_ionstage(element, ion);
      // printout("  [ionstage %d]\n", ionstage);
      //
      // printout("    For ionstage %d, the ground state populations are %g (function) and %g (matrix result with
      // normed pop %g, ltepopnormfactor %g)\n",get_ionstage(element,ion),
      //          get_groundlevelpop(modelgridindex, element, ion), gsl_vector_get(popvec, index_gs),
      //          gsl_vector_get(x, index_gs), gsl_vector_get(pop_norm_factor_vec, index_gs));

      // store the NLTE level populations
      const int nlte_start = globals::elements[element].ions[ion].first_nlte;
      // double solution_ion_pop = 0.;
      for (int level = 1; level <= nlevels_nlte; level++) {
        const int index = get_nlte_vector_index(element, ion, level);
        grid::modelgrid[modelgridindex].nlte_pops[nlte_start + level - 1] =
            gsl_vector_get(popvec, index) / grid::get_rho(modelgridindex);
        // solution_ion_pop += gsl_vector_get(popvec, index);
      }

      // store the superlevel population if there is one
      if (ion_has_superlevel(element, ion))  // a superlevel exists
      {
        const int index_sl = get_nlte_vector_index(element, ion, nlevels_nlte + 1);
        grid::modelgrid[modelgridindex].nlte_pops[nlte_start + nlevels_nlte] =
            (gsl_vector_get(popvec, index_sl) / grid::modelgrid[modelgridindex].rho / superlevel_partfunc[ion]);
      }

      // store the ground level population
      grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = gsl_vector_get(popvec, index_gs);
      // solution_ion_pop += gsl_vector_get(popvec, index_gs);

      calculate_cellpartfuncts(modelgridindex, element);
    }

    const double elem_pop_matrix = gsl_blas_dasum(popvec);
    const double elem_pop_error_percent = fabs((nnelement / elem_pop_matrix) - 1) * 100;
    if (elem_pop_error_percent > 1.0) {
      printout(
          "  WARNING: The Z=%d element population is: %g (from abundance) and %g (from matrix solution sum of level "
          "pops), error: %.1f%%. Forcing element pops to LTE.\n",
          atomic_number, nnelement, elem_pop_matrix, elem_pop_error_percent);
      set_element_pops_lte(modelgridindex, element);
    }

    if (individual_process_matricies && (timestep % 5 == 0) &&
        (nlte_iter == 0))  // output NLTE stats every nth timestep for the first NLTE iteration only
    {
      print_element_rates_summary(element, modelgridindex, timestep, nlte_iter, popvec, rate_matrix_rad_bb,
                                  rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf,
                                  rate_matrix_ntcoll_bf);
    }

    const bool print_detailed_level_stats = false;

    // if ((atomic_number == 26) && ((timestep % 5) == 0) && (nlte_iter == 0))
    // {
    //   print_detailed_level_stats = true;
    // }

    if (individual_process_matricies && print_detailed_level_stats) {
      const int ionstage = 2;
      const int ion = ionstage - get_ionstage(element, 0);

      for (int level = 0; level < get_nlevels_nlte(element, ion); level++) {
        print_level_rates(modelgridindex, timestep, element, ion, level, popvec, rate_matrix_rad_bb,
                          rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf,
                          rate_matrix_ntcoll_bf);
      }

      if (ion_has_superlevel(element, ion)) {
        const int slindex = get_nlevels_nlte(element, ion) + 1;
        print_level_rates(modelgridindex, timestep, element, ion, slindex, popvec, rate_matrix_rad_bb,
                          rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf,
                          rate_matrix_ntcoll_bf);
      }
    }
  }

  if (individual_process_matricies) {
    gsl_matrix_free(rate_matrix_rad_bb);
    gsl_matrix_free(rate_matrix_coll_bb);
    gsl_matrix_free(rate_matrix_ntcoll_bb);
    gsl_matrix_free(rate_matrix_rad_bf);
    gsl_matrix_free(rate_matrix_coll_bf);
    gsl_matrix_free(rate_matrix_ntcoll_bf);
  }

  gsl_vector_free(popvec);

  gsl_matrix_free(rate_matrix);
  gsl_vector_free(balance_vector);
  gsl_vector_free(pop_norm_factor_vec);
  const int duration_nltesolver = std::time(nullptr) - sys_time_start_nltesolver;
  if (duration_nltesolver > 2) {
    printout("NLTE population solver call for Z=%d took %d seconds\n", get_atomicnumber(element), duration_nltesolver);
  }
}

auto superlevel_boltzmann(const int modelgridindex, const int element, const int ion, const int level) -> double {
  const int superlevel_index = get_nlevels_nlte(element, ion) + 1;
  const double T_exc = LTEPOP_EXCITATION_USE_TJ ? grid::get_TJ(modelgridindex) : grid::get_Te(modelgridindex);
  const double E_level = epsilon(element, ion, level);
  const double E_superlevel = epsilon(element, ion, superlevel_index);

  return stat_weight(element, ion, level) / stat_weight(element, ion, superlevel_index) *
         exp(-(E_level - E_superlevel) / KB / T_exc);
}

void nltepop_open_file(const int my_rank) {
  char filename[MAXFILENAMELENGTH];
  snprintf(filename, MAXFILENAMELENGTH, "nlte_%.4d.out", my_rank);
  assert_always(nlte_file == nullptr);
  nlte_file = fopen_required(filename, "w");
  fprintf(nlte_file, "timestep modelgridindex Z ionstage level n_LTE n_NLTE ion_popfrac\n");
}

void nltepop_close_file() {
  if (nlte_file != nullptr) {
    fclose(nlte_file);
    nlte_file = nullptr;
  }
}

void nltepop_write_to_file(const int modelgridindex, const int timestep) {
  if (globals::lte_iteration || grid::modelgrid[modelgridindex].thick == 1) {  // NLTE solver hasn't been run yet
    return;
  }

  assert_always(nlte_file != nullptr);
  // fprintf(nlte_file,"#timestep %d modelgridindex %d T_R %g T_e %g W %g T_J %g nne %g\n",
  //         timestep, n, grid::get_TR(n), grid::get_Te(n), grid::get_W(n), grid::get_TJ(n), grid::get_nne(n));

  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    const int atomic_number = get_atomicnumber(element);

    for (int ion = 0; ion < nions; ion++) {
      const int nlevels_nlte = get_nlevels_nlte(element, ion);
      const int ion_first_nlte = globals::elements[element].ions[ion].first_nlte;
      const int ionstage = get_ionstage(element, ion);

      const int nsuperlevels = ion_has_superlevel(element, ion) ? 1 : 0;

      for (int level = 0; level <= nlevels_nlte + nsuperlevels; level++) {
        double nnlevellte = calculate_levelpop_lte(modelgridindex, element, ion, level);
        double nnlevelnlte{NAN};

        fprintf(nlte_file, "%d %d %d %d ", timestep, modelgridindex, atomic_number, ionstage);
        if (level <= nlevels_nlte) {
          fprintf(nlte_file, "%d ", level);

          if (level == 0) {
            nnlevelnlte = get_groundlevelpop(modelgridindex, element, ion);
          } else {
            nnlevelnlte = (grid::modelgrid[modelgridindex].nlte_pops[ion_first_nlte + level - 1] *
                           grid::modelgrid[modelgridindex].rho);
          }
        } else {
          // superlevel, so add the populations of all other levels in the superlevel
          const double slpopfactor = (grid::modelgrid[modelgridindex].nlte_pops[ion_first_nlte + nlevels_nlte] *
                                      grid::modelgrid[modelgridindex].rho);

          nnlevellte = 0;
          double superlevel_partfunc = 0;
          fprintf(nlte_file, "%d ", -1);
          for (int level_sl = nlevels_nlte + 1; level_sl < get_nlevels(element, ion); level_sl++) {
            nnlevellte += calculate_levelpop_lte(modelgridindex, element, ion, level_sl);
            superlevel_partfunc += superlevel_boltzmann(modelgridindex, element, ion, level_sl);
          }

          nnlevelnlte = slpopfactor * superlevel_partfunc;

          // printout("nltepop_write_to_file: The Z=%d ionstage %d superlevel population is %g with rho %g and
          // superlevel_partfunc %g Te %g scaled pop stored as %g\n", get_atomicnumber(element), get_ionstage(element,
          // ion), nnlevelnlte, grid::modelgrid[modelgridindex].rho, superlevel_partfunc, grid::get_Te(modelgridindex),
          // grid::modelgrid[modelgridindex].nlte_pops[ion_first_nlte + nlevels_nlte]);
        }

        const double ion_popfrac = nnlevelnlte / get_nnion(modelgridindex, element, ion);
        fprintf(nlte_file, "%.5e %.5e %.5e\n", nnlevellte, nnlevelnlte, ion_popfrac);
      }
    }
  }

  fflush(nlte_file);
}

void nltepop_write_restart_data(FILE *restart_file) {
  printout("populations, ");

  fprintf(restart_file, "%d\n", 75618527);  // special number marking the beginning of nlte data

  fprintf(restart_file, "%d\n", globals::total_nlte_levels);

  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    const int modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    fprintf(restart_file, "%d %la\n", modelgridindex, grid::modelgrid[modelgridindex].totalcooling);
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++) {
        fprintf(restart_file, "%d %a %a %la\n", ion,
                grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion],
                grid::modelgrid[modelgridindex].composition[element].partfunct[ion],
                grid::modelgrid[modelgridindex].cooling_contrib_ion[element][ion]);
      }
    }
    for (int nlteindex = 0; nlteindex < globals::total_nlte_levels; nlteindex++) {
      fprintf(restart_file, "%la ", grid::modelgrid[modelgridindex].nlte_pops[nlteindex]);
    }
  }
}

void nltepop_read_restart_data(FILE *restart_file) {
  printout("Reading restart data for populations\n");

  int code_check = 0;
  assert_always(fscanf(restart_file, "%d\n", &code_check) == 1);
  if (code_check != 75618527) {
    printout("ERROR: Beginning of NLTE restart data not found!\n");
    std::abort();
  }

  int total_nlte_levels_in = 0;
  assert_always(fscanf(restart_file, "%d\n", &total_nlte_levels_in) == 1);
  if (total_nlte_levels_in != globals::total_nlte_levels) {
    printout("ERROR: Expected %d NLTE levels but found %d in restart file\n", globals::total_nlte_levels,
             total_nlte_levels_in);
    std::abort();
  }

  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    const int modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    int mgi_in = 0;
    assert_always(fscanf(restart_file, "%d %la\n", &mgi_in, &grid::modelgrid[modelgridindex].totalcooling) == 2);
    if (mgi_in != modelgridindex) {
      printout("ERROR: expected data for cell %d but found cell %d\n", modelgridindex, mgi_in);
      std::abort();
    }

    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++) {
        int ion_in = 0;
        assert_always(fscanf(restart_file, "%d %a %a %la\n", &ion_in,
                             &grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion],
                             &grid::modelgrid[modelgridindex].composition[element].partfunct[ion],
                             &grid::modelgrid[modelgridindex].cooling_contrib_ion[element][ion]) == 4);
        if (ion_in != ion) {
          printout("ERROR: expected data for ion %d but found ion %d\n", ion, ion_in);
          std::abort();
        }
      }
    }
    for (int nlteindex = 0; nlteindex < globals::total_nlte_levels; nlteindex++) {
#ifdef MPI_ON
      if (globals::rank_in_node != 0) {
        assert_always(fscanf(restart_file, "%*a ") == 0);  // discard value (master rank of this node will set it)
      } else
#endif
        assert_always(fscanf(restart_file, "%la ", &grid::modelgrid[modelgridindex].nlte_pops[nlteindex]) == 1);
    }
  }
}
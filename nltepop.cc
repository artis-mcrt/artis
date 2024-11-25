#include "nltepop.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_double.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <ranges>
#include <span>
#include <tuple>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "ratecoeff.h"
#include "sn3d.h"

namespace {
FILE *nlte_file{};

// can save memory by using a combined rate matrix at the cost of diagnostic information
constexpr bool individual_process_matrices = true;

// this is the index for the NLTE solver that is handling all ions of a single element
// This is NOT an index into grid::modelgrid[nonemptymgi].nlte_pops that contains all elements
auto get_nlte_vector_index(const int element, const int ion, const int level) -> int {
  // have to convert from nlte_pops index to nlte_vector index
  // the difference is that nlte vectors apply to a single element and include ground states
  // The (+ ion) term accounts for the ground state population indices that are not counted in the NLTE array
  const int gs_index =
      globals::elements[element].ions[ion].first_nlte - globals::elements[element].ions[0].first_nlte + ion;

  // add in level or superlevel number
  const int level_index = gs_index + (is_nlte(element, ion, level) ? level : (get_nlevels_nlte(element, ion) + 1));

  return level_index;
}

[[nodiscard]] auto get_ion_level_of_nlte_vector_index(const int index, const int element) -> std::tuple<int, int> {
  // this could easily be optimized if need be
  for (int dion = 0; dion < get_nions(element); dion++) {
    for (int dlevel = 0; dlevel < get_nlevels(element, dion); dlevel++) {
      if (get_nlte_vector_index(element, dion, dlevel) == index) {
        return {dion, dlevel};
      }
    }
  }
  assert_always(false);
  return {-1, -1};
}

void eliminate_nlte_matrix_rowcol(const int index, const int gs_index, gsl_matrix *rate_matrix,
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

// find rows and columns that barely interaction with other levels, and effectively
// removing them by zeroing their interactions and setting their departure
// coeff to 1.0
void filter_nlte_matrix(const int element, gsl_matrix *rate_matrix, gsl_vector *balance_vector,
                        const gsl_vector * /*pop_norm_factor_vec*/) {
  const gsl_matrix rate_matrix_var = *rate_matrix;
  const int nlte_dimension = rate_matrix_var.size1;
  for (int index = 0; index < nlte_dimension; index++) {
    double row_max = 0.;
    for (int column = 0; column < nlte_dimension; column++) {
      const double element_value = fabs(gsl_matrix_get(rate_matrix, index, column));
      row_max = std::max(element_value, row_max);
    }
    double col_max = 0.;
    for (int row = 1; row < nlte_dimension; row++)  // skip the normalisation row 0
    {
      const double element_value = fabs(gsl_matrix_get(rate_matrix, row, index));
      col_max = std::max(element_value, col_max);
    }
    const auto [ion, level] = get_ion_level_of_nlte_vector_index(index, element);
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

[[nodiscard]] auto get_total_rate(const int index_selected, const gsl_matrix *rate_matrix, const gsl_vector *popvec,
                                  const bool into_level, const bool only_levels_below, const bool only_levels_above)
    -> double {
  double total_rate = 0.;
  assert_always(!only_levels_below || !only_levels_above);

  if (into_level) {
    // find rate into selected level
    auto row_vec = gsl_matrix_const_row(rate_matrix, index_selected).vector;

    // multiply incoming rate coefficients by their corresponding populations to get rates
    if (!only_levels_above)  // add levels below
    {
      for (int index = 0; index < index_selected; index++) {
        total_rate += gsl_vector_get(&row_vec, index) * gsl_vector_get(popvec, index);
      }
    }

    if (!only_levels_below)  // add levels above
    {
      for (unsigned int index = index_selected + 1; index < rate_matrix->size1; index++) {
        total_rate += gsl_vector_get(&row_vec, index) * gsl_vector_get(popvec, index);
      }
    }
  } else {
    // find rate out of selected level
    auto col_vec = gsl_matrix_const_column(rate_matrix, index_selected).vector;

    // multiply outgoing rate coefficients by the population of the selected level to get rates

    if (!only_levels_above)  // add levels below
    {
      for (int index = 0; index < index_selected; index++) {
        total_rate += gsl_vector_get(&col_vec, index);
      }
    }

    if (!only_levels_below)  // add levels above
    {
      for (unsigned int index = index_selected + 1; index < rate_matrix->size2; index++) {
        total_rate += gsl_vector_get(&col_vec, index);
      }
    }

    total_rate *= gsl_vector_get(popvec, index_selected);
  }

  return total_rate;
}

auto get_total_rate_in(const int index_to, const gsl_matrix *rate_matrix, const gsl_vector *popvec) -> double {
  return get_total_rate(index_to, rate_matrix, popvec, true, false, false);
}

auto get_total_rate_out(const int index_from, const gsl_matrix *rate_matrix, const gsl_vector *popvec) -> double {
  return get_total_rate(index_from, rate_matrix, popvec, false, false, false);
}

void print_level_rates_summary(const int element, const int selected_ion, const int selected_level,
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

void print_element_rates_summary(const int element, const int modelgridindex, const int timestep, const int nlte_iter,
                                 const gsl_vector *popvec, const gsl_matrix *rate_matrix_rad_bb,
                                 const gsl_matrix *rate_matrix_coll_bb, const gsl_matrix *rate_matrix_ntcoll_bb,
                                 const gsl_matrix *rate_matrix_rad_bf, const gsl_matrix *rate_matrix_coll_bf,
                                 const gsl_matrix *rate_matrix_ntcoll_bf) {
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
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
                 modelgridindex, timestep, nlte_iter, grid::get_Te(nonemptymgi), grid::get_nne(nonemptymgi),
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

void print_level_rates(const int modelgridindex, const int timestep, const int element, const int selected_ion,
                       const int selected_level, const gsl_vector *popvec, const gsl_matrix *rate_matrix_rad_bb,
                       const gsl_matrix *rate_matrix_coll_bb, const gsl_matrix *rate_matrix_ntcoll_bb,
                       const gsl_matrix *rate_matrix_rad_bf, const gsl_matrix *rate_matrix_coll_bf,
                       const gsl_matrix *rate_matrix_ntcoll_bf) {
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  // very detailed output of the NLTE processes for a particular levels

  if (element > get_nelements() - 1 || selected_ion > get_nions(element) - 1 ||
      selected_level >
          (get_nlevels_nlte(element, selected_ion) + (ion_has_superlevel(element, selected_ion) ? 1 : 0))) {
    printout("print_level_rates: invalid element/ion/level arguments\n");
    std::abort();
  }

  if (rate_matrix_rad_bb == rate_matrix_coll_bb) {
    printout(
        "print_level_rates: rate_matrix_rad_bb == rate_matrix_coll_bb. check individual_process_matrices is off\n");
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
      timestep, modelgridindex, grid::get_Te(nonemptymgi), grid::get_nne(nonemptymgi), atomic_number, selected_ionstage,
      selected_level);

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
    const auto [ion, level] = get_ion_level_of_nlte_vector_index(index, element);
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

void nltepop_reset_element(const int nonemptymgi, const int element) {
  const int nions = get_nions(element);
  for (int ion = 0; ion < nions; ion++) {
    std::fill_n(&grid::nltepops_allcells[(nonemptymgi * globals::total_nlte_levels)],
                get_nlevels_nlte(element, ion) + (ion_has_superlevel(element, ion) ? 1 : 0), -1.);
  }
}

auto get_element_superlevelpartfuncs(const int nonemptymgi, const int element) -> std::vector<double> {
  const int nions = get_nions(element);
  auto superlevel_partfuncs = std::vector<double>(nions, 0.);
  for (int ion = 0; ion < nions; ion++) {
    if (ion_has_superlevel(element, ion)) {
      const int nlevels_nlte = get_nlevels_nlte(element, ion);
      const int nlevels = get_nlevels(element, ion);
      for (int level = nlevels_nlte + 1; level < nlevels; level++) {
        superlevel_partfuncs[ion] += superlevel_boltzmann(nonemptymgi, element, ion, level);
      }
    }
  }

  return superlevel_partfuncs;
}

[[nodiscard]] auto get_element_nlte_dimension(const int element) -> int {
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

// get the maximum NLTE dimension for any of the included elements
[[nodiscard]] auto get_max_nlte_dimension() {
  int max_nlte_dimension = 0;
  for (int element = 0; element < get_nelements(); element++) {
    max_nlte_dimension = std::max(max_nlte_dimension, get_element_nlte_dimension(element));
  }
  return max_nlte_dimension;
}

void nltepop_matrix_add_boundbound(const int nonemptymgi, const int element, const int ion, const double t_mid,
                                   const std::vector<double> &s_renorm, gsl_matrix *rate_matrix_rad_bb,
                                   gsl_matrix *rate_matrix_coll_bb, gsl_matrix *rate_matrix_ntcoll_bb) {
  const auto T_e = grid::get_Te(nonemptymgi);
  const float nne = grid::get_nne(nonemptymgi);
  const int nlevels = get_nlevels(element, ion);
  const auto levels = std::ranges::iota_view{0, nlevels};
  std::for_each(levels.begin(), levels.end(), [&](const auto level) {
    const int level_index = get_nlte_vector_index(element, ion, level);
    const double epsilon_level = epsilon(element, ion, level);
    const double statweight = stat_weight(element, ion, level);
    const auto nnlevel = get_levelpop(nonemptymgi, element, ion, level);

    // de-excitation
    const int ndowntrans = get_ndowntrans(element, ion, level);
    const auto leveldowntrans = std::span(get_downtranslist(element, ion, level), ndowntrans);
    std::for_each(leveldowntrans.begin(), leveldowntrans.end(), [&](const auto &downtransition) {
      const double A_ul = downtransition.einstein_A;
      const int lower = downtransition.targetlevelindex;

      const double epsilon_trans = epsilon_level - epsilon(element, ion, lower);

      const double R = rad_deexcitation_ratecoeff(nonemptymgi, element, ion, lower, epsilon_trans, A_ul, statweight,
                                                  nnlevel, t_mid) *
                       s_renorm[level];
      const double C =
          col_deexcitation_ratecoeff(T_e, nne, epsilon_trans, element, ion, level, downtransition) * s_renorm[level];

      const int upper_index = level_index;
      const int lower_index = get_nlte_vector_index(element, ion, lower);

      atomicadd(*gsl_matrix_ptr(rate_matrix_rad_bb, upper_index, upper_index), -R);
      atomicadd(*gsl_matrix_ptr(rate_matrix_rad_bb, lower_index, upper_index), R);
      atomicadd(*gsl_matrix_ptr(rate_matrix_coll_bb, upper_index, upper_index), -C);
      atomicadd(*gsl_matrix_ptr(rate_matrix_coll_bb, lower_index, upper_index), C);
      if ((R < 0) || (C < 0)) {
        printout("  WARNING: Negative de-excitation rate from ionstage %d level %d to level %d\n",
                 get_ionstage(element, ion), level, lower);
      }
    });

    // excitation
    const int nuptrans = get_nuptrans(element, ion, level);
    const auto *const leveluptranslist = get_uptranslist(element, ion, level);
    const auto nuptransindices = std::ranges::iota_view{0, nuptrans};
    std::for_each(nuptransindices.begin(), nuptransindices.end(), [&](const auto i) {
      const int lineindex = leveluptranslist[i].lineindex;
      const int upper = leveluptranslist[i].targetlevelindex;
      const double epsilon_trans = epsilon(element, ion, upper) - epsilon_level;

      const double R =
          rad_excitation_ratecoeff(nonemptymgi, element, ion, level, i, epsilon_trans, nnlevel, lineindex, t_mid) *
          s_renorm[level];
      assert_always(R >= 0);
      assert_always(std::isfinite(R));

      const double C =
          col_excitation_ratecoeff(T_e, nne, element, ion, level, i, epsilon_trans, statweight) * s_renorm[level];
      assert_always(C >= 0);
      assert_always(std::isfinite(C));

      const double NTC =
          nonthermal::nt_excitation_ratecoeff(nonemptymgi, element, ion, level, i, lineindex) * s_renorm[level];

      const int lower_index = level_index;
      const int upper_index = get_nlte_vector_index(element, ion, upper);

      atomicadd(*gsl_matrix_ptr(rate_matrix_rad_bb, lower_index, lower_index), -R);
      atomicadd(*gsl_matrix_ptr(rate_matrix_rad_bb, upper_index, lower_index), R);
      atomicadd(*gsl_matrix_ptr(rate_matrix_coll_bb, lower_index, lower_index), -C);
      atomicadd(*gsl_matrix_ptr(rate_matrix_coll_bb, upper_index, lower_index), C);
      atomicadd(*gsl_matrix_ptr(rate_matrix_ntcoll_bb, lower_index, lower_index), -NTC);
      atomicadd(*gsl_matrix_ptr(rate_matrix_ntcoll_bb, upper_index, lower_index), NTC);

      if ((R < 0) || (C < 0)) {
        printout("  WARNING: Negative excitation rate from ion %d level %d to level %d\n", get_ionstage(element, ion),
                 level, upper);
      }
    });
  });
}

void nltepop_matrix_add_ionisation(const int modelgridindex, const int element, const int ion,
                                   const std::vector<double> &s_renorm, gsl_matrix *rate_matrix_rad_bf,
                                   gsl_matrix *rate_matrix_coll_bf) {
  assert_always((ion + 1) < get_nions(element));  // can't ionise the top ion
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  const auto T_e = grid::get_Te(nonemptymgi);
  const float nne = grid::get_nne(nonemptymgi);
  const int nionisinglevels = get_nlevels_ionising(element, ion);
  const int maxrecombininglevel = get_maxrecombininglevel(element, ion + 1);

  const auto levels = std::ranges::iota_view{0, nionisinglevels};
  std::for_each(EXEC_PAR levels.begin(), levels.end(), [&](const auto level) {
    const int lower_index = get_nlte_vector_index(element, ion, level);

    // thermal collisional ionization, photoionisation and recombination processes
    const double epsilon_current = epsilon(element, ion, level);

    const auto nphixstargets = get_nphixstargets(element, ion, level);
    for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
      const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
      const int upper_index = get_nlte_vector_index(element, ion + 1, upper);
      const double epsilon_trans = epsilon(element, ion + 1, upper) - epsilon_current;

      // ionization

      // the R part is slow!
      const double R_ionisation = get_corrphotoioncoeff(element, ion, level, phixstargetindex, nonemptymgi);
      const double C_ionisation =
          col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);

      atomicadd(*gsl_matrix_ptr(rate_matrix_rad_bf, lower_index, lower_index), -R_ionisation * s_renorm[level]);
      atomicadd(*gsl_matrix_ptr(rate_matrix_rad_bf, upper_index, lower_index), R_ionisation * s_renorm[level]);
      atomicadd(*gsl_matrix_ptr(rate_matrix_coll_bf, lower_index, lower_index), -C_ionisation * s_renorm[level]);
      atomicadd(*gsl_matrix_ptr(rate_matrix_coll_bf, upper_index, lower_index), C_ionisation * s_renorm[level]);

      if ((R_ionisation < 0) || (C_ionisation < 0)) {
        printout("  WARNING: Negative ionization rate from ionstage %d level %d phixstargetindex %d\n",
                 get_ionstage(element, ion), level, phixstargetindex);
      }

      // recombination
      if (upper <= maxrecombininglevel)  // we can skip this part if the functions below will return zero anyway
      {
        const double R_recomb = rad_recombination_ratecoeff(T_e, nne, element, ion + 1, upper, level, nonemptymgi);
        const double C_recomb = col_recombination_ratecoeff(T_e, nne, element, ion + 1, upper, level, epsilon_trans);

        atomicadd(*gsl_matrix_ptr(rate_matrix_rad_bf, upper_index, upper_index), -R_recomb * s_renorm[upper]);
        atomicadd(*gsl_matrix_ptr(rate_matrix_rad_bf, lower_index, upper_index), R_recomb * s_renorm[upper]);
        atomicadd(*gsl_matrix_ptr(rate_matrix_coll_bf, upper_index, upper_index), -C_recomb * s_renorm[upper]);
        atomicadd(*gsl_matrix_ptr(rate_matrix_coll_bf, lower_index, upper_index), C_recomb * s_renorm[upper]);

        if ((R_recomb < 0) || (C_recomb < 0)) {
          printout("  WARNING: Negative recombination rate to ionstage %d level %d phixstargetindex %d\n",
                   get_ionstage(element, ion), level, phixstargetindex);
        }
      }
    }
  });
}

void nltepop_matrix_add_nt_ionisation(const int nonemptymgi, const int element, const int ion,
                                      const std::vector<double> &s_renorm, gsl_matrix *rate_matrix_ntcoll_bf) {
  // collisional ionization by non-thermal electrons

  assert_always(ion + 1 < get_nions(element));  // can't ionise the top ion
  const double Y_nt = nonthermal::nt_ionization_ratecoeff(nonemptymgi, element, ion);
  if (Y_nt < 0.) {
    printout("  WARNING: Negative NT_ionization rate from ionstage %d\n", get_ionstage(element, ion));
  }

  const int nlevels = get_nlevels(element, ion);

  for (int upperion = ion + 1; upperion <= nonthermal::nt_ionisation_maxupperion(element, ion); upperion++) {
    const double Y_nt_thisupperion =
        Y_nt * nonthermal::nt_ionization_upperion_probability(nonemptymgi, element, ion, upperion, false);

    if (Y_nt_thisupperion > 0.) {
      const int upper_groundstate_index = get_nlte_vector_index(element, upperion, 0);
      for (int level = 0; level < nlevels; level++) {
        const int lower_index = get_nlte_vector_index(element, ion, level);

        atomicadd(*gsl_matrix_ptr(rate_matrix_ntcoll_bf, lower_index, lower_index),
                  -Y_nt_thisupperion * s_renorm[level]);
        atomicadd(*gsl_matrix_ptr(rate_matrix_ntcoll_bf, upper_groundstate_index, lower_index),
                  Y_nt_thisupperion * s_renorm[level]);
      }
    }
  }
}

void nltepop_matrix_normalise(const int nonemptymgi, const int element, gsl_matrix *rate_matrix,
                              gsl_vector *pop_norm_factor_vec) {
  const size_t nlte_dimension = pop_norm_factor_vec->size;
  assert_always(pop_norm_factor_vec->size == nlte_dimension);
  assert_always(rate_matrix->size1 == nlte_dimension);
  assert_always(rate_matrix->size2 == nlte_dimension);

  // TODO: consider replacing normalisation by LTE populations with
  // GSL's gsl_linalg_balance_matrix(gsl_matrix * A, gsl_vector * D) function instead
  for (size_t column = 0; column < nlte_dimension; column++) {
    const auto [ion, level] = get_ion_level_of_nlte_vector_index(column, element);

    gsl_vector_set(pop_norm_factor_vec, column, calculate_levelpop_lte(nonemptymgi, element, ion, level));

    if ((level != 0) && (!is_nlte(element, ion, level))) {
      // level is a superlevel, so add populations of higher levels to the norm factor
      for (int dummylevel = level + 1; dummylevel < get_nlevels(element, ion); dummylevel++) {
        if (!is_nlte(element, ion, dummylevel)) {
          *gsl_vector_ptr(pop_norm_factor_vec, column) += calculate_levelpop_lte(nonemptymgi, element, ion, dummylevel);
        }
      }
      // NOTE: above calculation is not always equal to the sum of LTE populations
      // since calculate_levelpop_lte imposes MINPOP minimum
      // printout("superlevel norm factor index %d is %g, partfunc is %g, partfunc*levelpop(SL)/g(SL) %g\n",
      //          column, gsl_vector_get(pop_norm_factor_vec, column), superlevel_partfunc[ion],
      //          superlevel_partfunc[ion] * calculate_levelpop_lte(nonemptymgi,element,ion,level) /
      //          stat_weight(element,ion,level));
    }

    // apply the normalisation factor to this column in the rate_matrix
    gsl_vector_view column_view = gsl_matrix_column(rate_matrix, column);
    gsl_vector_scale(&column_view.vector, gsl_vector_get(pop_norm_factor_vec, column));
  }
}

void set_element_pops_lte(const int nonemptymgi, const int element) {
  nltepop_reset_element(nonemptymgi, element);  // set NLTE pops as invalid so that LTE pops will be used instead
  calculate_cellpartfuncts(nonemptymgi, element);
  // Recall find_uppermost_ion with force_lte = true so uppermost ion used in set_groundlevelpops
  // is reset based on LTE phi factors instead of coming from NLTE phi factors
  const double nne_hi = grid::get_rho(nonemptymgi) / MH;
  const bool force_lte = true;
  const int uppermost_ion = find_uppermost_ion(nonemptymgi, element, nne_hi, force_lte);
  grid::set_elements_uppermost_ion(nonemptymgi, element, uppermost_ion);
  set_groundlevelpops(nonemptymgi, element, grid::get_nne(nonemptymgi), true);
}

[[nodiscard]] auto lumatrix_is_singular(const gsl_matrix *LU, const int element) -> bool {
  size_t const n = LU->size1;

  for (size_t i = 0; i < n; i++) {
    const double u = gsl_matrix_get(LU, i, i);
    if (u == 0) {
      const auto [ion, level] = get_ion_level_of_nlte_vector_index(i, element);
      if (is_nlte(element, ion, level)) {
        printout("NLTE disconnected level: Z=%d ionstage %d level %d\n", get_atomicnumber(element),
                 get_ionstage(element, ion), level);
      } else {
        printout("NLTE disconnected superlevel: Z=%d ionstage %d\n", get_atomicnumber(element),
                 get_ionstage(element, ion));
      }
      return true;
    }
  }

  return false;
}

// solve rate_matrix * x = balance_vector,
// then popvec[i] = x[i] / pop_norm_factor_vec[i]
// return true if the solution is successful, or false if the matrix is singular
[[nodiscard]] auto nltepop_matrix_solve(const int element, const gsl_matrix *rate_matrix,
                                        const gsl_vector *balance_vector, gsl_vector *popvec,
                                        const gsl_vector *pop_normfactor_vec, const int max_nlte_dimension) -> bool {
  const size_t nlte_dimension = balance_vector->size;
  assert_always(pop_normfactor_vec->size == nlte_dimension);
  assert_always(rate_matrix->size1 == nlte_dimension);
  assert_always(rate_matrix->size2 == nlte_dimension);

  // backing storage for gsl vectors
  THREADLOCALONHOST std::vector<double> vec_x;
  vec_x.resize(max_nlte_dimension, 0.);
  gsl_vector x = gsl_vector_view_array(vec_x.data(), nlte_dimension).vector;

  THREADLOCALONHOST std::vector<double> vec_rate_matrix_LU_decomp;
  vec_rate_matrix_LU_decomp.resize(max_nlte_dimension * max_nlte_dimension, 0.);
  // make a copy of the rate matrix for the LU decomp
  gsl_matrix rate_matrix_LU_decomp =
      gsl_matrix_view_array(vec_rate_matrix_LU_decomp.data(), nlte_dimension, nlte_dimension).matrix;
  gsl_matrix_memcpy(&rate_matrix_LU_decomp, rate_matrix);

  THREADLOCALONHOST std::vector<size_t> vec_permutation;
  vec_permutation.resize(max_nlte_dimension, 0);
  gsl_permutation p{.size = nlte_dimension, .data = vec_permutation.data()};
  gsl_permutation_init(&p);

  int s = 0;  // sign of the transformation
  gsl_linalg_LU_decomp(&rate_matrix_LU_decomp, &p, &s);

  if (lumatrix_is_singular(&rate_matrix_LU_decomp, element)) {
    printout("ERROR: NLTE matrix is singular for element Z=%d!\n", get_atomicnumber(element));
    // abort();
    return false;
  }

  gsl_error_handler_t *previous_handler = gsl_set_error_handler(gsl_error_handler_printout);

  // solve matrix equation: rate_matrix * x = balance_vector for x (population vector)
  gsl_linalg_LU_solve(&rate_matrix_LU_decomp, &p, balance_vector, &x);

  gsl_set_error_handler(previous_handler);

  // gsl_linalg_HH_solve (&m.matrix, &b.vector, x);

  const double TOLERANCE = 1e-40;
  THREADLOCALONHOST std::vector<double> vec_work;
  vec_work.resize(max_nlte_dimension, 0.);
  gsl_vector gsl_work_vector = gsl_vector_view_array(vec_work.data(), nlte_dimension).vector;

  double error_best = -1.;

  // population solution vector with lowest error
  THREADLOCALONHOST std::vector<double> vec_x_best;
  vec_x_best.resize(max_nlte_dimension, 0.);
  gsl_vector gsl_x_best = gsl_vector_view_array(vec_x_best.data(), nlte_dimension).vector;

  THREADLOCALONHOST std::vector<double> vec_residual;
  vec_residual.resize(max_nlte_dimension, 0.);
  gsl_vector gsl_vec_residual = gsl_vector_view_array(vec_residual.data(), nlte_dimension).vector;

  int iteration = 0;
  for (iteration = 0; iteration < 10; iteration++) {
    if (iteration > 0) {
      gsl_linalg_LU_refine(rate_matrix, &rate_matrix_LU_decomp, &p, balance_vector, &x, &gsl_work_vector);
    }

    gsl_vector_memcpy(&gsl_vec_residual, balance_vector);
    gsl_blas_dgemv(CblasNoTrans, 1.0, rate_matrix, &x, -1.0, &gsl_vec_residual);  // calculate Ax - b = residual
    const double error = fabs(gsl_vector_get(
        &gsl_vec_residual, gsl_blas_idamax(&gsl_vec_residual)));  // value of the largest absolute residual

    if (error < error_best || error_best < 0.) {
      gsl_vector_memcpy(&gsl_x_best, &x);
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

    gsl_vector_memcpy(&x, &gsl_x_best);
  }

  // get the real populations using the x vector and the normalisation factors
  gsl_vector_memcpy(popvec, &x);
  gsl_vector_mul(popvec, pop_normfactor_vec);
  // popvec will be used contains the real population densities

  size_t row_ground_state = 0;
  for (size_t row = 0; row < nlte_dimension; row++) {
    double recovered_balance_vector_elem = 0.;
    gsl_vector_const_view row_view = gsl_matrix_const_row(rate_matrix, row);
    gsl_blas_ddot(&row_view.vector, &x, &recovered_balance_vector_elem);

    const auto [ion, level] = get_ion_level_of_nlte_vector_index(row, element);
    if (level == 0) {
      row_ground_state = row;
    }

    // printout("index %4d (ionstage %d level%4d): residual %+.2e recovered balance: %+.2e normed pop %.2e pop %.2e
    // departure ratio %.4f\n",
    //          row,get_ionstage(element,ion),level, gsl_vector_get(residual_vector,row),
    //          recovered_balance_vector_elem, gsl_vector_get(x,row),
    //          gsl_vector_get(popvec, row),
    //          gsl_vector_get(x, row) / gsl_vector_get(x,get_nlte_vector_index(element,ion,0)));

    // Checking that groundpop is greater than MINPOP here - if it is then use LTE populations for element instead
    if (gsl_vector_get(popvec, row) < MINPOP && row == row_ground_state) {
      printout(
          "  WARNING: NLTE solver gave ground population less than MINPOP for index %zud (Z=%d ionstage %d level %d), "
          "pop = %g. "
          "Returning nltepop_matrix_solve fail (using LTE pops instead)\n",
          row, get_atomicnumber(element), get_ionstage(element, ion), level,
          gsl_vector_get(&x, row) * gsl_vector_get(pop_normfactor_vec, row));

        return false;
    }
    if (gsl_vector_get(popvec, row) < 0.0) {
      printout("  WARNING: NLTE solver gave negative population for index %zud (Z=%d ionstage %d level %d), pop = %g",
               row, get_atomicnumber(element), get_ionstage(element, ion), level,
               gsl_vector_get(&x, row) * gsl_vector_get(pop_normfactor_vec, row));
      if (gsl_vector_get(popvec, row) < -1*MINPOP) {
          printout(
              "  WARNING: negative pop = %g less than -1*MINPOP (-%g) unlikely a rounding error to zero so returning "
              "nltepop_matrix_solve fail (using LTE pops instead)\n", gsl_vector_get(popvec, row), MINPOP);

      return false;
      }
      printout(
              "  WARNING: negative pop = %g greater than -1*MINPOP (-%g) likely a rounding error to zero so continue "
              "with NLTE pops \n", gsl_vector_get(popvec, row), MINPOP);
    }
    if (row != row_ground_state &&
        gsl_vector_get(popvec, row_ground_state) <
        (stat_weight(element, ion, 0) / stat_weight(element, ion, level)) * gsl_vector_get(popvec, row)) {
            printout("[debug] WARNING: pop inversion: (g_pop %g)/(e_pop %g) = %g is less than (g_sw %g)/(e_sw %g) = %g "
            "for index %zud Z=%d ionstage %d level %d (factor %g inversion) - ",
            gsl_vector_get(popvec, row_ground_state), gsl_vector_get(popvec, row),
            gsl_vector_get(popvec, row_ground_state) / gsl_vector_get(popvec, row),
            stat_weight(element, ion, 0), stat_weight(element, ion, level), stat_weight(element, ion, 0) / stat_weight(element, ion, level),
            row, get_atomicnumber(element), get_ionstage(element, ion), level,
            (stat_weight(element, ion, 0) / stat_weight(element, ion, level)) / (gsl_vector_get(popvec, row_ground_state) / gsl_vector_get(popvec, row)));

            if (gsl_vector_get(popvec, row_ground_state) * 10000. <
                (stat_weight(element, ion, 0) / stat_weight(element, ion, level)) * gsl_vector_get(popvec, row)) {
              printout(
                  "large pop inversion (ground_pop * 10000 < ([g_gs / g_es] * excited_pop) - return matrix solve "
                  "fail and use LTE pops for element \n");
              return false;
            }
        if (gsl_vector_get(popvec, row_ground_state) * 10. < (stat_weight(element, ion, 0) / stat_weight(element, ion, level)) * gsl_vector_get(popvec, row)) {
          printout("more substantial pop inversion (ground_pop * 10 < ([g_gs / g_es] * excited_pop) - "
              "but continue with NLTE solution\n");
        }
        else {
          printout("relatively small pop inversion (ground_pop * 10 > ([g_gs / g_es] * excited_pop) - "
              "continue with NLTE solution\n");
        }
      }
  }

  return true;
}

} // anonymous namespace

void solve_nlte_pops_element(const int element, const int nonemptymgi, const int timestep, const int nlte_iter)
// solves the statistical balance equations to find NLTE level populations for all ions of an element
// (ionisation balance follows from this too)
{
  const int atomic_number = get_atomicnumber(element);
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);

  if (grid::get_elem_abundance(nonemptymgi, element) <= 0.) {
    // abundance of this element is zero, so do not store any NLTE populations
    printout("Not solving for NLTE populations in cell %d at timestep %d for element Z=%d due to zero abundance\n",
             modelgridindex, timestep, atomic_number);

    nltepop_reset_element(nonemptymgi, element);
    return;
  }

  const double cell_Te = grid::get_Te(nonemptymgi);

  if (cell_Te == MINTEMP) {
    printout(
        "Not solving for NLTE populations in cell %d at timestep %d for element Z=%d due to low temperature Te=%g\n",
        modelgridindex, timestep, atomic_number, cell_Te);
    set_element_pops_lte(nonemptymgi, element);
    return;
  }

  const auto sys_time_start_nltesolver = std::time(nullptr);

  const double t_mid = globals::timesteps[timestep].mid;
  const int nions = get_nions(element);
  const double nnelement = grid::get_elem_numberdens(nonemptymgi, element);

  printout(
      "Solving for NLTE populations in cell %d at timestep %d NLTE iteration %d for element Z=%d (mass fraction %.2e, "
      "nnelement %.2e cm^-3)\n",
      modelgridindex, timestep, nlte_iter, atomic_number, grid::get_elem_abundance(nonemptymgi, element), nnelement);

  const auto superlevel_partfunc = get_element_superlevelpartfuncs(nonemptymgi, element);
  const int nlte_dimension = get_element_nlte_dimension(element);

  // printout("NLTE: the vector dimension is %d", nlte_dimension);

  const auto max_nlte_dimension = get_max_nlte_dimension();

  THREADLOCALONHOST std::vector<double> vec_rate_matrix;
  vec_rate_matrix.resize(max_nlte_dimension * max_nlte_dimension, 0.);
  auto rate_matrix = gsl_matrix_view_array(vec_rate_matrix.data(), nlte_dimension, nlte_dimension).matrix;
  gsl_matrix_set_all(&rate_matrix, 0.);

  gsl_matrix rate_matrix_rad_bb;
  gsl_matrix rate_matrix_coll_bb;
  gsl_matrix rate_matrix_ntcoll_bb;
  gsl_matrix rate_matrix_rad_bf;
  gsl_matrix rate_matrix_coll_bf;
  gsl_matrix rate_matrix_ntcoll_bf;

  THREADLOCALONHOST std::vector<double> vec_rate_matrix_rad_bb;
  THREADLOCALONHOST std::vector<double> vec_rate_matrix_coll_bb;
  THREADLOCALONHOST std::vector<double> vec_rate_matrix_ntcoll_bb;
  THREADLOCALONHOST std::vector<double> vec_rate_matrix_rad_bf;
  THREADLOCALONHOST std::vector<double> vec_rate_matrix_coll_bf;
  THREADLOCALONHOST std::vector<double> vec_rate_matrix_ntcoll_bf;
  if constexpr (individual_process_matrices) {
    vec_rate_matrix_rad_bb.resize(max_nlte_dimension * max_nlte_dimension, 0.);
    rate_matrix_rad_bb = gsl_matrix_view_array(vec_rate_matrix_rad_bb.data(), nlte_dimension, nlte_dimension).matrix;
    gsl_matrix_set_all(&rate_matrix_rad_bb, 0.);

    vec_rate_matrix_coll_bb.resize(max_nlte_dimension * max_nlte_dimension, 0.);
    rate_matrix_coll_bb = gsl_matrix_view_array(vec_rate_matrix_coll_bb.data(), nlte_dimension, nlte_dimension).matrix;
    gsl_matrix_set_all(&rate_matrix_coll_bb, 0.);

    vec_rate_matrix_ntcoll_bb.resize(max_nlte_dimension * max_nlte_dimension, 0.);
    rate_matrix_ntcoll_bb =
        gsl_matrix_view_array(vec_rate_matrix_ntcoll_bb.data(), nlte_dimension, nlte_dimension).matrix;
    gsl_matrix_set_all(&rate_matrix_ntcoll_bb, 0.);

    vec_rate_matrix_rad_bf.resize(max_nlte_dimension * max_nlte_dimension, 0.);
    rate_matrix_rad_bf = gsl_matrix_view_array(vec_rate_matrix_rad_bf.data(), nlte_dimension, nlte_dimension).matrix;
    gsl_matrix_set_all(&rate_matrix_rad_bf, 0.);

    vec_rate_matrix_coll_bf.resize(max_nlte_dimension * max_nlte_dimension, 0.);
    rate_matrix_coll_bf = gsl_matrix_view_array(vec_rate_matrix_coll_bf.data(), nlte_dimension, nlte_dimension).matrix;
    gsl_matrix_set_all(&rate_matrix_coll_bf, 0.);

    vec_rate_matrix_ntcoll_bf.resize(max_nlte_dimension * max_nlte_dimension, 0.);
    rate_matrix_ntcoll_bf =
        gsl_matrix_view_array(vec_rate_matrix_ntcoll_bf.data(), nlte_dimension, nlte_dimension).matrix;
    gsl_matrix_set_all(&rate_matrix_ntcoll_bf, 0.);
  } else {
    // if not individual_process_matrices, alias a single matrix for all transition types
    // the "gsl_matrix" structs are independent, but the data is shared
    rate_matrix_rad_bb = rate_matrix;
    rate_matrix_coll_bb = rate_matrix;
    rate_matrix_ntcoll_bb = rate_matrix;
    rate_matrix_rad_bf = rate_matrix;
    rate_matrix_coll_bf = rate_matrix;
    rate_matrix_ntcoll_bf = rate_matrix;
  }

  // printout("  Adding rates for ion stages:");
  const auto ions = std::ranges::iota_view{0, nions};
  std::for_each(EXEC_PAR ions.begin(), ions.end(), [&](const auto ion) {
    // const int ionstage = get_ionstage(element, ion);
    // printout(" %d", ionstage);

    const int nlevels = get_nlevels(element, ion);
    const int nlevels_nlte = get_nlevels_nlte(element, ion);  // does not count the ground state!

    auto s_renorm = std::vector<double>(nlevels);
    std::fill_n(s_renorm.begin(), nlevels_nlte + 1, 1.);

    for (int level = (nlevels_nlte + 1); level < nlevels; level++) {
      s_renorm[level] = superlevel_boltzmann(nonemptymgi, element, ion, level) / superlevel_partfunc[ion];
    }

    nltepop_matrix_add_boundbound(nonemptymgi, element, ion, t_mid, s_renorm, &rate_matrix_rad_bb, &rate_matrix_coll_bb,
                                  &rate_matrix_ntcoll_bb);

    if (ion < (nions - 1)) {
      // this is the slowest component
      nltepop_matrix_add_ionisation(modelgridindex, element, ion, s_renorm, &rate_matrix_rad_bf, &rate_matrix_coll_bf);
      if (NT_ON) {
        nltepop_matrix_add_nt_ionisation(nonemptymgi, element, ion, s_renorm, &rate_matrix_ntcoll_bf);
      }
    }
  });
  // printout("\n");

  if (individual_process_matrices) {
    // sum the matrices for each transition type to get a total rate matrix
    gsl_matrix_add(&rate_matrix, &rate_matrix_rad_bb);
    gsl_matrix_add(&rate_matrix, &rate_matrix_coll_bb);
    gsl_matrix_add(&rate_matrix, &rate_matrix_ntcoll_bb);
    gsl_matrix_add(&rate_matrix, &rate_matrix_rad_bf);
    gsl_matrix_add(&rate_matrix, &rate_matrix_coll_bf);
    gsl_matrix_add(&rate_matrix, &rate_matrix_ntcoll_bf);
  }

  // replace the first row of the matrix and balance vector with the normalisation
  // constraint on the total element population
  gsl_vector_view first_row_view = gsl_matrix_row(&rate_matrix, 0);
  gsl_vector_set_all(&first_row_view.vector, 1.0);

  THREADLOCALONHOST std::vector<double> vec_balance_vector;
  vec_balance_vector.resize(max_nlte_dimension, 0.);
  auto balance_vector = gsl_vector_view_array(vec_balance_vector.data(), nlte_dimension).vector;
  gsl_vector_set_all(&balance_vector, 0.);
  // set first balance vector entry to the element population (all other entries will be zero)
  gsl_vector_set(&balance_vector, 0, nnelement);

  if (FORCE_SAHA_ION_BALANCE(atomic_number)) {
    const auto ionfractions = calculate_ionfractions(element, nonemptymgi, grid::get_nne(nonemptymgi), true);
    const int uppermost_ion = static_cast<int>(ionfractions.size() - 1);
    for (int ion = 1; ion <= uppermost_ion; ion++) {
      // replace matrix row for ion's ground state with sum of this ion's level populations is equal to the ion
      // population
      const double nnion = nnelement * ionfractions[ion];
      const int index_ion_ground = get_nlte_vector_index(element, ion, 0);
      const int index_ion_toplevel = get_nlte_vector_index(element, ion, get_nlevels(element, ion));
      gsl_vector_view ion_ground_row_view = gsl_matrix_row(&rate_matrix, index_ion_ground);
      gsl_vector_set_all(&ion_ground_row_view.vector, 0.);
      for (int index = index_ion_ground; index <= index_ion_toplevel; index++) {
        gsl_vector_set(&ion_ground_row_view.vector, index, 1.);
      }

      gsl_vector_set(&balance_vector, get_nlte_vector_index(element, ion, index_ion_ground), nnion);
    }
  }

  // calculate the normalisation factors and apply them to the matrix
  // columns and balance vector elements
  THREADLOCALONHOST std::vector<double> vec_pop_norm_factor_vec;
  vec_pop_norm_factor_vec.resize(max_nlte_dimension, 0.);
  auto pop_norm_factor_vec = gsl_vector_view_array(vec_pop_norm_factor_vec.data(), nlte_dimension).vector;
  gsl_vector_set_all(&pop_norm_factor_vec, 1.0);

  nltepop_matrix_normalise(nonemptymgi, element, &rate_matrix, &pop_norm_factor_vec);

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

  // the true population densities
  THREADLOCALONHOST std::vector<double> vec_pop;
  vec_pop.resize(max_nlte_dimension, 0.);
  auto popvec = gsl_vector_view_array(vec_pop.data(), nlte_dimension).vector;

  const bool matrix_solve_success =
      nltepop_matrix_solve(element, &rate_matrix, &balance_vector, &popvec, &pop_norm_factor_vec, max_nlte_dimension);

  if (!matrix_solve_success) {
    printout(
        "WARNING: Can't solve for NLTE populations in cell %d at timestep %d for element Z=%d due to singular matrix. "
        ", negative population or large population inversion. Attempting to use LTE solution instead\n",
        modelgridindex, timestep, atomic_number);
    set_element_pops_lte(nonemptymgi, element);
  } else {
    // check calculated NLTE populations are valid
    for (int index = 0; index < nlte_dimension; index++) {
      assert_always(std::isfinite(gsl_vector_get(&popvec, index)));
      assert_always(gsl_vector_get(&popvec, index) >= 0.);
    }

    for (int ion = 0; ion < nions; ion++) {
      const int nlevels_nlte = get_nlevels_nlte(element, ion);
      const int index_gs = get_nlte_vector_index(element, ion, 0);
      // const int ionstage = get_ionstage(element, ion);
      // printout("  [ionstage %d]\n", ionstage);
      //
      // printout("    For ionstage %d, the ground state populations are %g (function) and %g (matrix result with
      // normed pop %g, ltepopnormfactor %g)\n",get_ionstage(element,ion),
      //          get_groundlevelpop(nonemptymgi, element, ion), gsl_vector_get(popvec, index_gs),
      //          gsl_vector_get(x, index_gs), gsl_vector_get(pop_norm_factor_vec, index_gs));

      // store the NLTE level populations
      // double solution_ion_pop = 0.;
      for (int level = 1; level <= nlevels_nlte; level++) {
        const int index = get_nlte_vector_index(element, ion, level);
        set_nlte_levelpop_over_rho(nonemptymgi, element, ion, level,
                                   gsl_vector_get(&popvec, index) / grid::get_rho(nonemptymgi));
        // solution_ion_pop += gsl_vector_get(popvec, index);
      }

      // store the superlevel population if there is one
      if (ion_has_superlevel(element, ion))  // a superlevel exists
      {
        const int index_sl = get_nlte_vector_index(element, ion, nlevels_nlte + 1);
        set_nlte_superlevelpop_over_rho(
            nonemptymgi, element, ion,
            gsl_vector_get(&popvec, index_sl) / grid::get_rho(nonemptymgi) / superlevel_partfunc[ion]);
      }

      // store the ground level population
      grid::ion_groundlevelpops_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) +
                                         get_uniqueionindex(element, ion)] = gsl_vector_get(&popvec, index_gs);
      // solution_ion_pop += gsl_vector_get(popvec, index_gs);

      calculate_cellpartfuncts(nonemptymgi, element);
    }

    const double elem_pop_matrix = gsl_blas_dasum(&popvec);
    const double elem_pop_error_percent = fabs((nnelement / elem_pop_matrix) - 1) * 100;
    if (elem_pop_error_percent > 1.0) {
      printout(
          "  WARNING: The Z=%d element population is: %g (from abundance) and %g (from matrix solution sum of level "
          "pops), error: %.1f%%. Forcing element pops to LTE.\n",
          atomic_number, nnelement, elem_pop_matrix, elem_pop_error_percent);
      set_element_pops_lte(nonemptymgi, element);
    }

    if (individual_process_matrices && (timestep % 5 == 0) &&
        (nlte_iter == 0))  // output NLTE stats every nth timestep for the first NLTE iteration only
    {
      print_element_rates_summary(element, modelgridindex, timestep, nlte_iter, &popvec, &rate_matrix_rad_bb,
                                  &rate_matrix_coll_bb, &rate_matrix_ntcoll_bb, &rate_matrix_rad_bf,
                                  &rate_matrix_coll_bf, &rate_matrix_ntcoll_bf);
    }

    const bool print_detailed_level_stats = false;

    // if ((atomic_number == 26) && ((timestep % 5) == 0) && (nlte_iter == 0))
    // {
    //   print_detailed_level_stats = true;
    // }

    if (individual_process_matrices && print_detailed_level_stats) {
      const int ionstage = 2;
      const int ion = ionstage - get_ionstage(element, 0);

      for (int level = 0; level < get_nlevels_nlte(element, ion); level++) {
        print_level_rates(modelgridindex, timestep, element, ion, level, &popvec, &rate_matrix_rad_bb,
                          &rate_matrix_coll_bb, &rate_matrix_ntcoll_bb, &rate_matrix_rad_bf, &rate_matrix_coll_bf,
                          &rate_matrix_ntcoll_bf);
      }

      if (ion_has_superlevel(element, ion)) {
        const int slindex = get_nlevels_nlte(element, ion) + 1;
        print_level_rates(modelgridindex, timestep, element, ion, slindex, &popvec, &rate_matrix_rad_bb,
                          &rate_matrix_coll_bb, &rate_matrix_ntcoll_bb, &rate_matrix_rad_bf, &rate_matrix_coll_bf,
                          &rate_matrix_ntcoll_bf);
      }
    }
  }

  const int duration_nltesolver = std::time(nullptr) - sys_time_start_nltesolver;
  if (duration_nltesolver > 2) {
    printout("NLTE population solver call for Z=%d took %d seconds\n", get_atomicnumber(element), duration_nltesolver);
  }
}

// Get a Boltzman factor for a level within the super level (combined Non-LTE level)
__host__ __device__ auto superlevel_boltzmann(const int nonemptymgi, const int element, const int ion, const int level)
    -> double {
  assert_testmodeonly(level_isinsuperlevel(element, ion, level));
  const int superlevel_index = get_nlevels_nlte(element, ion) + 1;
  const double T_exc = LTEPOP_EXCITATION_USE_TJ ? grid::get_TJ(nonemptymgi) : grid::get_Te(nonemptymgi);
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

void nltepop_write_to_file(const int nonemptymgi, const int timestep) {
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  if (globals::lte_iteration || grid::modelgrid[nonemptymgi].thick == 1) {  // NLTE solver hasn't been run yet
    return;
  }

  assert_always(nlte_file != nullptr);
  // fprintf(nlte_file,"#timestep %d modelgridindex %d T_R %g T_e %g W %g T_J %g nne %g\n",
  //         timestep, n, grid::get_TR(n), grid::get_Te(n), grid::get_W(n), grid::get_TJ(n), grid::get_nne(n));

  for (int element = 0; element < get_nelements(); element++) {
    if (!elem_has_nlte_levels(element)) {
      continue;
    }
    const int nions = get_nions(element);
    const int atomic_number = get_atomicnumber(element);

    for (int ion = 0; ion < nions; ion++) {
      const int nlevels_nlte = get_nlevels_nlte(element, ion);
      const int ionstage = get_ionstage(element, ion);

      const int nsuperlevels = ion_has_superlevel(element, ion) ? 1 : 0;

      for (int level = 0; level <= nlevels_nlte + nsuperlevels; level++) {
        double nnlevellte = calculate_levelpop_lte(nonemptymgi, element, ion, level);
        double nnlevelnlte{NAN};

        fprintf(nlte_file, "%d %d %d %d ", timestep, modelgridindex, atomic_number, ionstage);
        if (level <= nlevels_nlte) {
          fprintf(nlte_file, "%d ", level);

          if (level == 0) {
            nnlevelnlte = get_groundlevelpop(nonemptymgi, element, ion);
          } else {
            nnlevelnlte =
                get_nlte_levelpop_over_rho(nonemptymgi, element, ion, level) * grid::modelgrid[nonemptymgi].rho;
          }
        } else {
          // superlevel, so add the populations of all other levels in the superlevel
          const double slpopfactor =
              get_nlte_superlevelpop_over_rho(nonemptymgi, element, ion) * grid::modelgrid[nonemptymgi].rho;

          nnlevellte = 0;
          double superlevel_partfunc = 0;
          fprintf(nlte_file, "%d ", -1);
          for (int level_sl = nlevels_nlte + 1; level_sl < get_nlevels(element, ion); level_sl++) {
            nnlevellte += calculate_levelpop_lte(nonemptymgi, element, ion, level_sl);
            superlevel_partfunc += superlevel_boltzmann(nonemptymgi, element, ion, level_sl);
          }

          nnlevelnlte = slpopfactor * superlevel_partfunc;
        }

        const double ion_popfrac = nnlevelnlte / get_nnion(nonemptymgi, element, ion);
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
  const auto nincludedions = get_includedions();

  for (ptrdiff_t nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    const int modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    fprintf(restart_file, "%d %la\n", modelgridindex, grid::modelgrid[nonemptymgi].totalcooling);
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++) {
        const int uniqueionindex = get_uniqueionindex(element, ion);
        fprintf(restart_file, "%d %a %a %la\n", ion,
                grid::ion_groundlevelpops_allcells[(nonemptymgi * nincludedions) + uniqueionindex],
                grid::ion_partfuncts_allcells[(nonemptymgi * nincludedions) + uniqueionindex],
                grid::ion_cooling_contribs_allcells[(nonemptymgi * nincludedions) + uniqueionindex]);
      }
    }
    for (int nlteindex = 0; nlteindex < globals::total_nlte_levels; nlteindex++) {
      fprintf(restart_file, "%la ", grid::nltepops_allcells[(nonemptymgi * globals::total_nlte_levels) + nlteindex]);
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
  const auto nincludedions = get_includedions();

  for (ptrdiff_t nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    int mgi_in = 0;
    assert_always(fscanf(restart_file, "%d %la\n", &mgi_in, &grid::modelgrid[nonemptymgi].totalcooling) == 2);
    assert_always(mgi_in == grid::get_mgi_of_nonemptymgi(nonemptymgi));

    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++) {
        int ion_in = 0;
        const int uniqueionindex = get_uniqueionindex(element, ion);
        assert_always(fscanf(restart_file, "%d %a %a %la\n", &ion_in,
                             &grid::ion_groundlevelpops_allcells[(nonemptymgi * nincludedions) + uniqueionindex],
                             &grid::ion_partfuncts_allcells[(nonemptymgi * nincludedions) + uniqueionindex],
                             &grid::ion_cooling_contribs_allcells[(nonemptymgi * nincludedions) + uniqueionindex]) ==
                      4);
        assert_always(ion_in == ion);
      }
    }
    for (int nlteindex = 0; nlteindex < globals::total_nlte_levels; nlteindex++) {
      assert_always(fscanf(restart_file, "%la ",
                           &grid::nltepops_allcells[(nonemptymgi * globals::total_nlte_levels) + nlteindex]) == 1);
    }
  }
}

auto get_nlte_levelpop_over_rho(const int nonemptymgi, const int element, const int ion, const int level) -> double {
  assert_testmodeonly(level <= get_nlevels_nlte(element, ion));
  return grid::nltepops_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * globals::total_nlte_levels) +
                                 globals::elements[element].ions[ion].first_nlte + level - 1];
}

auto get_nlte_superlevelpop_over_rho(const int nonemptymgi, const int element, const int ion) -> double {
  assert_testmodeonly(ion_has_superlevel(element, ion));
  const int sl_nlte_index = globals::elements[element].ions[ion].first_nlte + get_nlevels_nlte(element, ion);
  return grid::nltepops_allcells[(nonemptymgi * globals::total_nlte_levels) + sl_nlte_index];
}

void set_nlte_levelpop_over_rho(const int nonemptymgi, const int element, const int ion, const int level,
                                const double value) {
  assert_testmodeonly(level <= get_nlevels_nlte(element, ion));
  grid::nltepops_allcells[(nonemptymgi * globals::total_nlte_levels) + globals::elements[element].ions[ion].first_nlte +
                          level - 1] = value;
}

void set_nlte_superlevelpop_over_rho(const int nonemptymgi, const int element, const int ion, const double value) {
  assert_testmodeonly(ion_has_superlevel(element, ion));
  const int sl_nlte_index = globals::elements[element].ions[ion].first_nlte + get_nlevels_nlte(element, ion);
  grid::nltepops_allcells[(nonemptymgi * globals::total_nlte_levels) + sl_nlte_index] = value;
}

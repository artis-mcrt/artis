#include "nltepop.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <format>
#include <fstream>
#include <ios>
#include <numeric>
#include <ranges>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_double.h>
#pragma clang unsafe_buffer_usage end

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

static_assert(STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING >= 1);
static_assert(STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING <
              STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL);
namespace {
std::fstream nlte_file;

struct RateMatrices {
  int max_nlte_dimension{0};
  int used_nlte_dimension{0};

  std::vector<double> vec_rate_matrix;
  std::vector<double> vec_rate_matrix_rad_bb;
  std::vector<double> vec_rate_matrix_coll_bb;
  std::vector<double> vec_rate_matrix_ntcoll_bb;
  std::vector<double> vec_rate_matrix_rad_bf;
  std::vector<double> vec_rate_matrix_coll_bf;
  std::vector<double> vec_rate_matrix_ntcoll_bf;
  std::vector<double> vec_rate_matrix_autoion;

  gsl_matrix rad_bb{};
  gsl_matrix coll_bb{};
  gsl_matrix ntcoll_bb{};
  gsl_matrix rad_bf{};
  gsl_matrix coll_bf{};
  gsl_matrix ntcoll_bf{};
  gsl_matrix autoion{};

  explicit RateMatrices(int max_nlte_dimension_in) : max_nlte_dimension(max_nlte_dimension_in) {
    // allocation of the maximum required size is done once,
    // while the used_nlte_dimension is set later
    const auto max_dim_square = max_nlte_dimension * max_nlte_dimension;
    resize_exactly(vec_rate_matrix, max_dim_square);
    resize_exactly(vec_rate_matrix_rad_bb, max_dim_square);
    resize_exactly(vec_rate_matrix_coll_bb, max_dim_square);
    resize_exactly(vec_rate_matrix_ntcoll_bb, max_dim_square);
    resize_exactly(vec_rate_matrix_rad_bf, max_dim_square);
    resize_exactly(vec_rate_matrix_coll_bf, max_dim_square);
    resize_exactly(vec_rate_matrix_ntcoll_bf, max_dim_square);
    resize_exactly(vec_rate_matrix_autoion, max_dim_square);
  }

  void set_used_dimension(int nlte_dimension) {
    assert_always(std::cmp_less_equal(nlte_dimension, max_nlte_dimension));
    used_nlte_dimension = nlte_dimension;
    std::fill_n(vec_rate_matrix_rad_bb.data(), used_nlte_dimension * used_nlte_dimension, 0.);
    std::fill_n(vec_rate_matrix_coll_bb.data(), used_nlte_dimension * used_nlte_dimension, 0.);
    std::fill_n(vec_rate_matrix_ntcoll_bb.data(), used_nlte_dimension * used_nlte_dimension, 0.);
    std::fill_n(vec_rate_matrix_rad_bf.data(), used_nlte_dimension * used_nlte_dimension, 0.);
    std::fill_n(vec_rate_matrix_coll_bf.data(), used_nlte_dimension * used_nlte_dimension, 0.);
    std::fill_n(vec_rate_matrix_ntcoll_bf.data(), used_nlte_dimension * used_nlte_dimension, 0.);
    std::fill_n(vec_rate_matrix_autoion.data(), used_nlte_dimension * used_nlte_dimension, 0.);

    rad_bb = gsl_matrix_view_array(vec_rate_matrix_rad_bb.data(), used_nlte_dimension, used_nlte_dimension).matrix;
    coll_bb = gsl_matrix_view_array(vec_rate_matrix_coll_bb.data(), used_nlte_dimension, used_nlte_dimension).matrix;
    ntcoll_bb =
        gsl_matrix_view_array(vec_rate_matrix_ntcoll_bb.data(), used_nlte_dimension, used_nlte_dimension).matrix;
    rad_bf = gsl_matrix_view_array(vec_rate_matrix_rad_bf.data(), used_nlte_dimension, used_nlte_dimension).matrix;
    coll_bf = gsl_matrix_view_array(vec_rate_matrix_coll_bf.data(), used_nlte_dimension, used_nlte_dimension).matrix;
    ntcoll_bf =
        gsl_matrix_view_array(vec_rate_matrix_ntcoll_bf.data(), used_nlte_dimension, used_nlte_dimension).matrix;
    autoion = gsl_matrix_view_array(vec_rate_matrix_autoion.data(), used_nlte_dimension, used_nlte_dimension).matrix;
  }

  [[nodiscard]] auto get_summed_rate_matrix() -> gsl_matrix {
    // sum the matrices for each transition type to get a total rate matrix
    for (int i = 0; i < used_nlte_dimension * used_nlte_dimension; i++) {
      vec_rate_matrix[i] = vec_rate_matrix_rad_bb[i] + vec_rate_matrix_coll_bb[i] + vec_rate_matrix_ntcoll_bb[i] +
                           vec_rate_matrix_rad_bf[i] + vec_rate_matrix_coll_bf[i] + vec_rate_matrix_ntcoll_bf[i] +
                           vec_rate_matrix_autoion[i];
    }

    return gsl_matrix_view_array(vec_rate_matrix.data(), used_nlte_dimension, used_nlte_dimension).matrix;
  }
};

// this is the matrix/vector index for the NLTE solver that is handling all ions of a single element
auto get_nlte_vector_index(const int element, const int ion, const int level, const int first_ion_used) -> int {
  // have to convert from nlte_pops index to nlte_vector index
  // the difference is that nlte vectors apply to a single element and include ground states
  // The (+ ion) term accounts for the ground state population indices that are not counted in the NLTE array
  int offset_autoion = 0;
  for (int dion = first_ion_used; dion < ion; dion++) {
    offset_autoion += get_nlevels_autoion(element, dion);
  }
  assert_testmodeonly(first_ion_used >= 0);
  assert_testmodeonly(first_ion_used < get_nions(element));
  const int gs_index = globals::elements[element].ions[ion].first_nlte -
                       globals::elements[element].ions[first_ion_used].first_nlte + (ion - first_ion_used) +
                       offset_autoion;
  const int first_autoion = get_nlevels(element, ion) - get_nlevels_autoion(element, ion);
  if (ion_has_superlevel(element, ion) && level > (first_autoion - 1)) {
    return (gs_index + get_nlevels_excited_nlte(element, ion) + level - first_autoion + 2);
  }
  // add in level or superlevel number
  const int level_index =
      gs_index + (is_nlte(element, ion, level) ? level : (get_nlevels_excited_nlte(element, ion) + 1));
  return level_index;
}

[[nodiscard]] auto get_ion_level_of_nlte_vector_index(const int index, const int element, const int first_ion_used,
                                                      const int nions_used) -> std::tuple<int, int> {
  // this could easily be optimized if need be
  for (int dion = first_ion_used; dion < first_ion_used + nions_used; dion++) {
    for (int dlevel = 0; dlevel < get_nlevels(element, dion); dlevel++) {
      if (get_nlte_vector_index(element, dion, dlevel, first_ion_used) == index) {
        return {dion, dlevel};
      }
    }
  }
  assert_always(false);
  return {-1, -1};
}

[[nodiscard]] auto get_total_rate(const size_t index_selected, const gsl_matrix &rate_matrix,
                                  const std::span<const double> popvec, const bool into_level,
                                  const bool only_levels_below, const bool only_levels_above) -> double {
  double total_rate = 0.;
  assert_always(!only_levels_below || !only_levels_above);

  if (into_level) {
    // find rate into selected level
    const auto row_vec = gsl_matrix_const_row(&rate_matrix, index_selected).vector;

    // multiply incoming rate coefficients by their corresponding populations to get rates
    if (!only_levels_above) {  // add levels below
      for (size_t index = 0; index < index_selected; index++) {
        total_rate += gsl_vector_get(&row_vec, index) * popvec[index];
      }
    }

    if (!only_levels_below) {  // add levels above
      for (size_t index = index_selected + 1; index < rate_matrix.size1; index++) {
        total_rate += gsl_vector_get(&row_vec, index) * popvec[index];
      }
    }
  } else {
    // find rate out of selected level
    const auto col_vec = gsl_matrix_const_column(&rate_matrix, index_selected).vector;

    // multiply outgoing rate coefficients by the population of the selected level to get rates

    if (!only_levels_above) {  // add levels below
      for (size_t index = 0; index < index_selected; index++) {
        total_rate += gsl_vector_get(&col_vec, index);
      }
    }

    if (!only_levels_below) {  // add levels above
      for (size_t index = index_selected + 1; index < rate_matrix.size2; index++) {
        total_rate += gsl_vector_get(&col_vec, index);
      }
    }

    total_rate *= popvec[index_selected];
  }

  return total_rate;
}

auto get_total_rate_in(const int index_to, const gsl_matrix &rate_matrix, const std::span<const double> popvec)
    -> double {
  return get_total_rate(index_to, rate_matrix, popvec, true, false, false);
}

auto get_total_rate_out(const int index_from, const gsl_matrix &rate_matrix, const std::span<const double> popvec)
    -> double {
  return get_total_rate(index_from, rate_matrix, popvec, false, false, false);
}

void print_level_rates_summary(const int element, const int selected_ion, const int selected_level,
                               std::span<const double> popvec, const RateMatrices &rate_matrices,
                               const int first_ion_used) {
  const int selected_index = get_nlte_vector_index(element, selected_ion, selected_level, first_ion_used);

  for (int i = 0; i <= 3; i++) {
    // rates in from below, in from above, out to below, out to above
    if (i == 0) {
      const int nlevels_nlte = get_nlevels_excited_nlte(element, selected_ion);
      if (ion_has_superlevel(element, selected_ion) && (selected_level == nlevels_nlte + 1)) {
        printout("      superlevel ");
      } else {
        printout("    level%7d ", selected_level);
      }
      printout(" %10.2e ", popvec[selected_index]);
    } else {
      printout("                             ");
    }

    const bool into_level = (i <= 1);
    const bool only_levels_below = (i % 2) != 0;
    const bool only_levels_above = !only_levels_below;

    const double rad_bb_total =
        get_total_rate(selected_index, rate_matrices.rad_bb, popvec, into_level, only_levels_below, only_levels_above);
    const double coll_bb_total =
        get_total_rate(selected_index, rate_matrices.coll_bb, popvec, into_level, only_levels_below, only_levels_above);
    const double ntcoll_bb_total = get_total_rate(selected_index, rate_matrices.ntcoll_bb, popvec, into_level,
                                                  only_levels_below, only_levels_above);
    const double rad_bf_total =
        get_total_rate(selected_index, rate_matrices.rad_bf, popvec, into_level, only_levels_below, only_levels_above);
    const double coll_bf_total =
        get_total_rate(selected_index, rate_matrices.coll_bf, popvec, into_level, only_levels_below, only_levels_above);
    const double ntcoll_bf_total = get_total_rate(selected_index, rate_matrices.ntcoll_bf, popvec, into_level,
                                                  only_levels_below, only_levels_above);

    if (into_level) {
      // into this level
      printout(" from ");
    } else {
      // out of this level
      printout("   to ");
    }

    if (only_levels_below) {
      printout("below ");
    } else {
      printout("above ");
    }

    printout("%10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n", rad_bb_total, coll_bb_total, ntcoll_bb_total, rad_bf_total,
             coll_bf_total, ntcoll_bf_total);
  }
}

void print_element_rates_summary(const int element, const int modelgridindex, const int timestep, const int nlte_iter,
                                 std::span<const double> popvec, const RateMatrices &rate_matrices,
                                 const int first_ion_used, const int nions_used) {
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  for (int ion = first_ion_used; ion < first_ion_used + nions_used - 1; ion++) {
    const int nlevels = get_nlevels(element, ion);
    const int nlevels_nlte = get_nlevels_excited_nlte(element, ion);

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

      print_level_rates_summary(element, ion, level, popvec, rate_matrices, first_ion_used);
    }

    if (ion_has_superlevel(element, ion) && max_printed_levels < (nlevels_nlte + 1)) {
      const int level_superlevel = nlevels_nlte + 1;

      print_level_rates_summary(element, ion, level_superlevel, popvec, rate_matrices, first_ion_used);
    }
  }
}

void print_level_rates(const int nonemptymgi, const int timestep, const int element, const int selected_ion,
                       const int selected_level, const std::span<const double> popvec,
                       const RateMatrices &rate_matrices, const int first_ion_used, const int nions_used) {
  // very detailed output of each transition for a particular level

  assert_always(selected_ion >= first_ion_used);
  assert_always(selected_ion < first_ion_used + nions_used - 1);
  assert_always(selected_level <= (get_nlevels_excited_nlte(element, selected_ion) +
                                   (ion_has_superlevel(element, selected_ion) ? 1 : 0)));

  const int nlte_dimension = popvec.size();
  const int atomic_number = get_atomicnumber(element);
  const int selected_ionstage = get_ionstage(element, selected_ion);
  const int selected_index = get_nlte_vector_index(element, selected_ion, selected_level, first_ion_used);
  const double pop_selectedlevel = popvec[selected_index];
  printout(
      "timestep %d cell %d Te %g nne %g NLTE level diagnostics for Z=%d ionstage %d level %d rates into and out of "
      "this level\n",
      timestep, grid::get_mgi_of_nonemptymgi(nonemptymgi), grid::get_Te(nonemptymgi), grid::get_nne(nonemptymgi),
      atomic_number, selected_ionstage, selected_level);

  const double rad_bb_in_total = get_total_rate_in(selected_index, rate_matrices.rad_bb, popvec);
  const double coll_bb_in_total = get_total_rate_in(selected_index, rate_matrices.coll_bb, popvec);
  const double ntcoll_bb_in_total = get_total_rate_in(selected_index, rate_matrices.ntcoll_bb, popvec);
  const double rad_bf_in_total = get_total_rate_in(selected_index, rate_matrices.rad_bf, popvec);
  const double coll_bf_in_total = get_total_rate_in(selected_index, rate_matrices.coll_bf, popvec);
  const double ntcoll_bf_in_total = get_total_rate_in(selected_index, rate_matrices.ntcoll_bf, popvec);
  const double total_rate_in =
      rad_bb_in_total + coll_bb_in_total + rad_bf_in_total + coll_bf_in_total + ntcoll_bf_in_total;
  printout(
      "  TOTAL rates in:             rad_bb_in  %8.2e coll_bb_in  %8.2e ntcoll_bb_in  %8.2e rad_bf_in  %8.2e "
      "coll_bf_in  %8.2e ntcoll_bf_in  %8.2e\n",
      rad_bb_in_total, coll_bb_in_total, ntcoll_bb_in_total, rad_bf_in_total, coll_bf_in_total, ntcoll_bf_in_total);

  const double rad_bb_out_total = get_total_rate_out(selected_index, rate_matrices.rad_bb, popvec);
  const double coll_bb_out_total = get_total_rate_out(selected_index, rate_matrices.coll_bb, popvec);
  const double ntcoll_bb_out_total = get_total_rate_out(selected_index, rate_matrices.ntcoll_bb, popvec);
  const double rad_bf_out_total = get_total_rate_out(selected_index, rate_matrices.rad_bf, popvec);
  const double coll_bf_out_total = get_total_rate_out(selected_index, rate_matrices.coll_bf, popvec);
  const double ntcoll_bf_out_total = get_total_rate_out(selected_index, rate_matrices.ntcoll_bf, popvec);
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
    const auto [ion, level] = get_ion_level_of_nlte_vector_index(index, element, first_ion_used, nions_used);
    const int ionstage = get_ionstage(element, ion);
    // in means populating the selected level, out means depopulating the selected level
    const double pop = popvec[index];
    const double rad_bb_in = gsl_matrix_get(&rate_matrices.rad_bb, selected_index, index) * pop;
    const double rad_bb_out = gsl_matrix_get(&rate_matrices.rad_bb, index, selected_index) * pop_selectedlevel;
    const double coll_bb_in = gsl_matrix_get(&rate_matrices.coll_bb, selected_index, index) * pop;
    const double coll_bb_out = gsl_matrix_get(&rate_matrices.coll_bb, index, selected_index) * pop_selectedlevel;
    const double ntcoll_bb_in = gsl_matrix_get(&rate_matrices.ntcoll_bb, selected_index, index) * pop;
    const double ntcoll_bb_out = gsl_matrix_get(&rate_matrices.ntcoll_bb, index, selected_index) * pop_selectedlevel;
    const double rad_bf_in = gsl_matrix_get(&rate_matrices.rad_bf, selected_index, index) * pop;
    const double rad_bf_out = gsl_matrix_get(&rate_matrices.rad_bf, index, selected_index) * pop_selectedlevel;
    const double coll_bf_in = gsl_matrix_get(&rate_matrices.coll_bf, selected_index, index) * pop;
    const double coll_bf_out = gsl_matrix_get(&rate_matrices.coll_bf, index, selected_index) * pop_selectedlevel;
    const double ntcoll_bf_in = gsl_matrix_get(&rate_matrices.ntcoll_bf, selected_index, index) * pop;
    const double ntcoll_bf_out = gsl_matrix_get(&rate_matrices.ntcoll_bf, index, selected_index) * pop_selectedlevel;

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
  for (int ion = 0; ion < get_nions(element); ion++) {
    const int nlte_start = globals::elements[element].ions[ion].first_nlte;
    std::fill_n(&grid::nltepops_allcells[(nonemptymgi * globals::total_nlte_levels) + nlte_start],
                get_nlevels_excited_nlte(element, ion) + (ion_has_superlevel(element, ion) ? 1 : 0), -1.);
  }
}

auto get_element_superlevelpartfuncs(const int nonemptymgi, const int element) -> std::vector<double> {
  std::vector<double> superlevel_partfuncs;
  resize_exactly(superlevel_partfuncs, get_nions(element));
  for (int ion = 0; ion < get_nions(element); ion++) {
    superlevel_partfuncs[ion] = std::ranges::fold_left(
        std::views::iota(get_nlevels_excited_nlte(element, ion) + 1,
                         get_nlevels(element, ion) - get_nlevels_autoion(element, ion)),
        0.0, [&](double sum, int level) { return sum + superlevel_boltzmann(nonemptymgi, element, ion, level); });
  }

  return superlevel_partfuncs;
}

[[nodiscard]] auto get_element_nlte_dimension(const int element, const int first_ion_used, const int nions_used)
    -> int {
  assert_testmodeonly(nions_used >= 0);
  assert_testmodeonly(nions_used <= get_nions(element));
  assert_testmodeonly(first_ion_used >= 0);
  assert_testmodeonly(first_ion_used < get_nions(element));
  assert_testmodeonly((first_ion_used + nions_used - 1) < get_nions(element));
  int nlte_dimension = 0;
  for (int ion = first_ion_used; ion < (first_ion_used + nions_used); ion++) {
    // add super level if it exists
    if (ion_has_superlevel(element, ion)) {
      nlte_dimension += get_nlevels_excited_nlte(element, ion) + get_nlevels_autoion(element, ion) + 2;
      // printout("Here 1: For element %d ion %d adding %d to nlte_dimension. \n", element, ion, nlevels_nlte +
      // get_nlevels_autoion(element, ion) + 2); printout("checks: %d %d\n", nlevels_nlte, get_nlevels_autoion(element,
      // ion)); if it has a superlevel then need + 2 for ground state and super and to add autionising levels
    } else {  // if it doesn't have a superlevel
      nlte_dimension += get_nlevels(element, ion);
      // printout("Here 2: For element %d ion %d adding %d to nlte_dimension. \n", element, ion,
      // get_nlevels(element,ion));
    }
  }

  return nlte_dimension;
}

// get the maximum NLTE dimension for any of the included elements
[[nodiscard]] auto get_max_nlte_dimension() -> int {
  int max_nlte_dimension = 0;
  for (int element = 0; element < get_nelements(); element++) {
    const int first_ion_used = 0;
    const int nions_used = get_nions(element);
    max_nlte_dimension = std::max(max_nlte_dimension, get_element_nlte_dimension(element, first_ion_used, nions_used));
  }
  return max_nlte_dimension;
}

void nltepop_matrix_add_boundbound(const int nonemptymgi, const int element, const int ion, const double t_mid,
                                   const std::span<double> s_renorm, RateMatrices &rate_matrices,
                                   const int first_ion_used) {
  const auto T_e = grid::get_Te(nonemptymgi);
  const auto nne = grid::get_nne(nonemptymgi);
  const int nlevels = get_nlevels(element, ion);
  const auto ionuniquelevelindexstart = globals::elements[element].ions[ion].uniquelevelindexstart;
  const auto nlte_dimension = rate_matrices.used_nlte_dimension;

  const auto levels = std::views::iota(0, nlevels);
  std::for_each(levels.begin(), levels.end(), [&](const auto level) {
    const int level_index = get_nlte_vector_index(element, ion, level, first_ion_used);
    const auto matrix_index_level_level = (level_index * nlte_dimension) + level_index;
    const auto uniquelevelindex = ionuniquelevelindexstart + level;
    const double epsilon_level = epsilon(uniquelevelindex);
    const double statweight = stat_weight(uniquelevelindex);
    const auto nnlevel = get_levelpop(nonemptymgi, uniquelevelindex);

    // de-excitation
    const auto alltrans_startdown = get_alltrans_startdown(uniquelevelindex);
    const int ndowntrans = get_ndowntrans(uniquelevelindex);
    const auto ndowntransindices = std::views::iota(0, ndowntrans);
    std::for_each(ndowntransindices.begin(), ndowntransindices.end(), [&](const auto &i) {
      const auto alltransindex = alltrans_startdown + i;
      const int lower = globals::alltrans.targetlevelindex[alltransindex];
      const auto lower_uniquelevelindex = ionuniquelevelindexstart + lower;
      const auto lower_statweight = stat_weight(lower_uniquelevelindex);

      const double epsilon_trans = epsilon_level - epsilon(lower_uniquelevelindex);
      const double R = rad_deexcitation_ratecoeff(nonemptymgi, lower_uniquelevelindex, epsilon_trans,
                                                  globals::alltrans.einstein_A[alltransindex], statweight,
                                                  lower_statweight, nnlevel, t_mid) *
                       s_renorm[level];
      const double C =
          col_deexcitation_ratecoeff(T_e, nne, epsilon_trans, statweight, lower_statweight, alltransindex) *
          s_renorm[level];

      const int lower_index = get_nlte_vector_index(element, ion, lower, first_ion_used);
      const auto matrix_index_upper_upper = matrix_index_level_level;
      const auto matrix_index_lower_upper = (lower_index * nlte_dimension) + level_index;

      atomicadd(rate_matrices.vec_rate_matrix_rad_bb[matrix_index_upper_upper], -R);
      atomicadd(rate_matrices.vec_rate_matrix_rad_bb[matrix_index_lower_upper], R);
      atomicadd(rate_matrices.vec_rate_matrix_coll_bb[matrix_index_upper_upper], -C);
      atomicadd(rate_matrices.vec_rate_matrix_coll_bb[matrix_index_lower_upper], C);
    });

    // excitation
    const int nuptrans = get_nuptrans(uniquelevelindex);
    const auto alltrans_startup = get_alltrans_startup(uniquelevelindex);
    const auto nuptransindices = std::views::iota(0, nuptrans);
    std::for_each(nuptransindices.begin(), nuptransindices.end(), [&](const auto i) {
      const auto alltransindex = alltrans_startup + i;
      const int upper = globals::alltrans.targetlevelindex[alltransindex];
      const auto upper_uniquelevelindex = ionuniquelevelindexstart + upper;
      const double epsilon_trans = epsilon(upper_uniquelevelindex) - epsilon_level;
      const auto upper_statweight = stat_weight(upper_uniquelevelindex);

      const double R = rad_excitation_ratecoeff(nonemptymgi, upper_uniquelevelindex, upper_statweight,
                                                globals::alltrans.einstein_A[alltransindex], epsilon_trans, nnlevel,
                                                statweight, alltransindex, t_mid) *
                       s_renorm[level];

      const double C = col_excitation_ratecoeff(T_e, nne, upper_statweight, alltransindex, epsilon_trans, statweight) *
                       s_renorm[level];

      const double NTC =
          nonthermal::nt_excitation_ratecoeff(nonemptymgi, level, upper, alltransindex) * s_renorm[level];

      const int upper_index = get_nlte_vector_index(element, ion, upper, first_ion_used);
      const auto matrix_index_lower_lower = matrix_index_level_level;
      const auto matrix_index_upper_lower = (upper_index * nlte_dimension) + level_index;

      atomicadd(rate_matrices.vec_rate_matrix_rad_bb[matrix_index_lower_lower], -R);
      atomicadd(rate_matrices.vec_rate_matrix_rad_bb[matrix_index_upper_lower], R);
      atomicadd(rate_matrices.vec_rate_matrix_coll_bb[matrix_index_lower_lower], -C);
      atomicadd(rate_matrices.vec_rate_matrix_coll_bb[matrix_index_upper_lower], C);
      atomicadd(rate_matrices.vec_rate_matrix_ntcoll_bb[matrix_index_lower_lower], -NTC);
      atomicadd(rate_matrices.vec_rate_matrix_ntcoll_bb[matrix_index_upper_lower], NTC);
    });
  });
}

void nltepop_matrix_add_ionisation(const int nonemptymgi, const int element, const int ion,
                                   const std::span<double> s_renorm, RateMatrices &rate_matrices,
                                   const int first_ion_used, const int nions_used) {
  assert_always((ion + 1) < (nions_used + first_ion_used));  // can't ionise top ion stage
  const auto T_e = grid::get_Te(nonemptymgi);
  const float nne = grid::get_nne(nonemptymgi);
  const int nionisinglevels = get_nlevels_ionising(element, ion);
  const int maxrecombininglevel = get_maxrecombininglevel(element, ion + 1);
  const auto nlte_dimension = rate_matrices.used_nlte_dimension;

  const auto levels = std::views::iota(0, nionisinglevels);
  std::for_each(EXEC_PAR levels.begin(), levels.end(), [&](const auto level) {
    const int lower_index = get_nlte_vector_index(element, ion, level, first_ion_used);
    const auto matrix_index_lower_lower = (lower_index * nlte_dimension) + lower_index;

    const double epsilon_current = epsilon(element, ion, level);

    const auto nphixstargets = get_nphixstargets(element, ion, level);
    for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
      const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
      const int upper_index = get_nlte_vector_index(element, ion + 1, upper, first_ion_used);
      const double epsilon_trans = epsilon(element, ion + 1, upper) - epsilon_current;

      // photoionization and collisional ionization
      const double R_ionisation = get_corrphotoioncoeff(element, ion, level, phixstargetindex, nonemptymgi);
      const double C_ionisation =
          col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);

      const auto matrix_index_upper_lower = (upper_index * nlte_dimension) + lower_index;

      atomicadd(rate_matrices.vec_rate_matrix_rad_bf[matrix_index_lower_lower], -R_ionisation * s_renorm[level]);
      atomicadd(rate_matrices.vec_rate_matrix_rad_bf[matrix_index_upper_lower], R_ionisation * s_renorm[level]);
      atomicadd(rate_matrices.vec_rate_matrix_coll_bf[matrix_index_lower_lower], -C_ionisation * s_renorm[level]);
      atomicadd(rate_matrices.vec_rate_matrix_coll_bf[matrix_index_upper_lower], C_ionisation * s_renorm[level]);

      if ((R_ionisation < 0) || (C_ionisation < 0)) {
        printout("  WARNING: Negative ionization rate from ionstage %d level %d phixstargetindex %d\n",
                 get_ionstage(element, ion), level, phixstargetindex);
      }

      // recombination
      if (upper <= maxrecombininglevel) {
        const double R_recomb = rad_recombination_ratecoeff(T_e, nne, element, ion + 1, upper, level, nonemptymgi);
        const double C_recomb = col_recombination_ratecoeff(T_e, nne, element, ion + 1, upper, level, epsilon_trans);

        const auto matrix_index_upper_upper = (upper_index * nlte_dimension) + upper_index;
        const auto matrix_index_lower_upper = (lower_index * nlte_dimension) + upper_index;

        atomicadd(rate_matrices.vec_rate_matrix_rad_bf[matrix_index_upper_upper], -R_recomb * s_renorm[upper]);
        atomicadd(rate_matrices.vec_rate_matrix_rad_bf[matrix_index_lower_upper], R_recomb * s_renorm[upper]);
        atomicadd(rate_matrices.vec_rate_matrix_coll_bf[matrix_index_upper_upper], -C_recomb * s_renorm[upper]);
        atomicadd(rate_matrices.vec_rate_matrix_coll_bf[matrix_index_lower_upper], C_recomb * s_renorm[upper]);

        if ((R_recomb < 0) || (C_recomb < 0)) {
          printout("  WARNING: Negative recombination rate to ionstage %d level %d phixstargetindex %d\n",
                   get_ionstage(element, ion), level, phixstargetindex);
        }
      }
    }
  });
}

void nltepop_matrix_add_nt_ionisation(const int nonemptymgi, const int element, const int ion,
                                      const std::span<double> s_renorm, RateMatrices &rate_matrices,
                                      const int first_ion_used, const int nions_used) {
  // collisional ionization by non-thermal electrons
  const int max_ion_used = first_ion_used + nions_used - 1;
  assert_always(ion < max_ion_used);  // can't ionise top ion stage
  const double Y_nt = nonthermal::nt_ionization_ratecoeff(nonemptymgi, element, ion);

  const int nlevels = get_nlevels(element, ion);
  const auto nlte_dimension = rate_matrices.used_nlte_dimension;

  for (int upperion = ion + 1; upperion <= nonthermal::nt_ionisation_maxupperion(element, ion); upperion++) {
    const double Y_nt_thisupperion =
        Y_nt * nonthermal::nt_ionization_upperion_probability(nonemptymgi, element, ion, upperion, false);

    if (Y_nt_thisupperion > 0.) {
      // ensure if upperion here is past the upper ion used in the NLTE matrix the rates
      // to this ion are instead added to the ground state of the uppermost ion used in the NLTE matrix
      const int upperion_clamped = std::min(upperion, max_ion_used);
      const int upper_groundstate_index = get_nlte_vector_index(element, upperion_clamped, 0, first_ion_used);
      for (int level = 0; level < nlevels; level++) {
        const int lower_index = get_nlte_vector_index(element, ion, level, first_ion_used);

        atomicadd(rate_matrices.vec_rate_matrix_ntcoll_bf[(lower_index * nlte_dimension) + lower_index],
                  -Y_nt_thisupperion * s_renorm[level]);
        atomicadd(rate_matrices.vec_rate_matrix_ntcoll_bf[(upper_groundstate_index * nlte_dimension) + lower_index],
                  Y_nt_thisupperion * s_renorm[level]);
      }
    }
  }
}

void nltepop_matrix_add_autoionisation(const int nonemptymgi, const int element, const int ion,
                                       const std::vector<double> &s_renorm, RateMatrices &rate_matrices,
                                       const int first_ion_used, const int nions_used) {
  // Autoionization and inverse (i.e. collisional capture part of di-el)
  const auto nlte_dimension = rate_matrices.used_nlte_dimension;
  const int max_ion_used = first_ion_used + nions_used - 1;
  assert_always(ion < max_ion_used);  // can't ionise top ion stage
  const auto T_e = grid::get_Te(nonemptymgi);
  const float nne = grid::get_nne(nonemptymgi);
  const int nlevels = get_nlevels(element, ion);
  for (int level = 0; level < nlevels; level++) {
    const int level_index = get_nlte_vector_index(element, ion, level, first_ion_used);
    const double epsilon_level = epsilon(element, ion, level);
    const double statweight = stat_weight(element, ion, level);

    const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
    const int nautoiondowntrans = get_nautoiondowntrans(uniquelevelindex);
    for (int i = 0; i < nautoiondowntrans; i++) {
      // autoionization (which is a de-excitation propcess)
      const auto &autoiontransition = globals::allautoion[globals::alllevels.allautoion_start[uniquelevelindex] + i];
      const double A_a = autoiontransition.autoion_A;
      const int target_ion = autoiontransition.upperionindex;
      const int target_level = autoiontransition.upperlevelindex;

      const double epsilon_trans = epsilon_level - epsilon(element, target_ion, target_level);

      double R = A_a * s_renorm[level];

      const int upper_index = level_index;
      const int lower_index = get_nlte_vector_index(element, target_ion, target_level, first_ion_used);

      rate_matrices.vec_rate_matrix_autoion[(upper_index * nlte_dimension) + upper_index] -= R;
      rate_matrices.vec_rate_matrix_autoion[(lower_index * nlte_dimension) + upper_index] += R;
      if ((R < 0)) {
        printout("  WARNING: Negative autoionization rate from ionstage %d level %d to level %d\n",
                 get_ionstage(element, ion), level, target_level);
      }

      // capture (which is an excitation process, and the first part of di-electronic recomb)
      R = nne * A_a * statweight / stat_weight(element, target_ion, target_level) * SAHACONST * pow(T_e, -1.5) *
          exp(-1. * epsilon_trans / KB / T_e);
      // renorm??

      rate_matrices.vec_rate_matrix_autoion[(lower_index * nlte_dimension) + lower_index] -= R;
      rate_matrices.vec_rate_matrix_autoion[(upper_index * nlte_dimension) + lower_index] += R;
      if ((R < 0)) {
        printout("  WARNING: Negative autoionization rate from ionstage %d level %d to level %d\n",
                 get_ionstage(element, ion), level, target_level);
      }
    }
  }
}

void nltepop_matrix_normalise(const int nonemptymgi, const int element, gsl_matrix *rate_matrix,
                              std::span<double> pop_norm_factors, const int first_ion_used, const int nions_used) {
  const size_t nlte_dimension = rate_matrix->size1;
  assert_always(pop_norm_factors.size() == nlte_dimension);
  assert_always(rate_matrix->size2 == nlte_dimension);

  for (size_t column = 0; column < nlte_dimension; column++) {
    const auto [ion, level] = get_ion_level_of_nlte_vector_index(column, element, first_ion_used, nions_used);

    pop_norm_factors[column] = calculate_levelpop_boltzmann(nonemptymgi, element, ion, level);

    if (level_isinsuperlevel(element, ion, level)) {
      // levels in the superlevel get combined together
      for (int dummylevel = level + 1; dummylevel < get_nlevels(element, ion); dummylevel++) {
        if (level_isinsuperlevel(element, ion, dummylevel)) {
          pop_norm_factors[column] += calculate_levelpop_boltzmann(nonemptymgi, element, ion, dummylevel);
        }
      }
    }

    // apply the normalisation factor to this column in the rate_matrix
    gsl_vector_view column_view = gsl_matrix_column(rate_matrix, column);
    gsl_vector_scale(&column_view.vector, pop_norm_factors[column]);
  }
}

void set_nlte_levelpop_over_rho(const int nonemptymgi, const int element, const int ion, const int level,
                                const double value) {
  assert_testmodeonly(level > 0);  // ground state is stored separately
  assert_testmodeonly(level <= get_nlevels_excited_nlte(element, ion));
  grid::nltepops_allcells[(nonemptymgi * globals::total_nlte_levels) + globals::elements[element].ions[ion].first_nlte +
                          level - 1] = value;
}

void set_nlte_superlevelpop_over_rho_over_slpartfunc(const int nonemptymgi, const int element, const int ion,
                                                     const double value) {
  assert_testmodeonly(ion_has_superlevel(element, ion));
  const int sl_nlte_index = globals::elements[element].ions[ion].first_nlte + get_nlevels_excited_nlte(element, ion);
  grid::nltepops_allcells[(nonemptymgi * globals::total_nlte_levels) + sl_nlte_index] = value;
}

void set_element_pops_lte(const int nonemptymgi, const int element) {
  // set NLTE level pops as invalid so that Boltzmann pops will be used instead
  nltepop_reset_element(nonemptymgi, element);
  calculate_cellpartfuncts(nonemptymgi, element);
  // set_groundlevelpops uses uppermost_ion. Previously, this was set based on the NLTE phi factors. Therefore we need
  // to call find_uppermost_ion with force_saha = true so the uppermost ion used in set_groundlevelpops is changed to
  // the one based on the correct LTE phi factors instead
  if (NLTEPOP_FAILURE_USE_FIND_UPPERMOST_ION_FOR_LTE_RESET) {
    printout("NLTEPOP_FAILURE_USE_FIND_UPPERMOST_ION_FOR_LTE_RESET for element %d\n", element);
    const double nne_hi = grid::get_rho(nonemptymgi) / MH;
    const bool force_saha = true;
    const int uppermost_ion = find_uppermost_ion(nonemptymgi, element, nne_hi, force_saha);
    grid::set_elements_uppermost_ion(nonemptymgi, element, uppermost_ion);
  }
  set_groundlevelpops(nonemptymgi, element, grid::get_nne(nonemptymgi), true);
}

[[nodiscard]] auto lumatrix_is_singular(const gsl_matrix *LU, const int element, const int first_ion_used,
                                        const int nions_used) -> bool {
  for (size_t i = 0; i < LU->size1; i++) {
    // diagonal elements of LU matrix should not be zero
    // if they are, then the matrix is singular and the NLTE solution will fail
    if (gsl_matrix_get(LU, i, i) == 0) {
      const auto [ion, level] = get_ion_level_of_nlte_vector_index(i, element, first_ion_used, nions_used);
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

[[nodiscard]] auto solution_pops_are_valid(const int nonemptymgi, const int element, std::span<double> popvec,
                                           std::span<const double> pop_normfactors, const int first_ion_used,
                                           const int nions_used) -> bool {
  const size_t nlte_dimension = popvec.size();
  const auto superlevel_partfuncs = get_element_superlevelpartfuncs(nonemptymgi, element);
  for (size_t index = 0; index < nlte_dimension; index++) {
    const auto [ion, level] = get_ion_level_of_nlte_vector_index(index, element, first_ion_used, nions_used);
    const auto ionstage = get_ionstage(element, ion);
    const auto population = popvec[index];

    if constexpr (!STRICT_POPULATION_CHECKING) {
      if (population < 0.0) {
        printout(
            "  WARNING: NLTE solver gave negative population to index %zu (Z=%d ionstage %d level %d), pop = %g. "
            "Replacing with LTE pop of %g\n",
            index, get_atomicnumber(element), ionstage, level, population, pop_normfactors[index]);
        popvec[index] = pop_normfactors[index];
      }
    } else {
      // if groundpop is below MINPOP then the NLTE solution fails
      const size_t index_ion_ground = get_nlte_vector_index(element, ion, 0, first_ion_used);
      if (index == index_ion_ground && population < MINPOP) {
        printout(
            "  WARNING: NLTE solver gave ground pop less than MINPOP for index %zu (Z=%d ionstage %d level %d), "
            "pop = %g. Returning nltepop_matrix_solve fail\n",
            index, get_atomicnumber(element), ionstage, level, population);
        return false;
      }

      if (population < 0.0) {
        printout("  WARNING: NLTE solver gave negative population for index %zu (Z=%d ionstage %d level %d), pop = %g",
                 index, get_atomicnumber(element), ionstage, level, population);
        if (population < -MINPOP) {
          printout(
              "  WARNING: negative pop = %g less than -MINPOP (-%g) unlikely a rounding error to zero so "
              "returning nltepop_matrix_solve fail\n",
              population, MINPOP);
          return false;
        }
        printout(
            "  WARNING: negative pop = %g greater than -MINPOP (-%g) likely a rounding error to zero so continue "
            "with NLTE pops but set this level to MINPOP\n",
            population, MINPOP);
        popvec[index] = MINPOP;
      }

      if (index != index_ion_ground) {
        // Check for population inversions
        const double ground_pop = popvec[index_ion_ground];
        const auto statweight_ground = stat_weight(element, ion, 0);
        const auto statweight = stat_weight(element, ion, level);
        if (!level_isinsuperlevel(element, ion, level)) {
          const double inversion_factor = (population / statweight) / (ground_pop / statweight_ground);

          if (inversion_factor > STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING) {
            printout(
                "[debug] WARNING: pop inversion greater than factor %g: (g_pop %g)/(e_pop %g) = %g is less than (g_sw "
                "%g)/(e_sw %g) = %g for index %zu Z=%d ionstage %d level %d (factor %g inversion) - ",
                STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING, ground_pop, population,
                ground_pop / population, statweight_ground, statweight, statweight_ground / statweight, index,
                get_atomicnumber(element), ionstage, level, inversion_factor);

            if (inversion_factor > STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL) {
              printout("large pop inversion - matrix solve failed\n");
              return false;
            }
            printout("relatively small pop inversion so continue with NLTE solution\n");
          }
        } else {
          // check if the first sublevel in the superlevel is inverted relative to the ground state. Only need to check
          // the first level as the superlevel levels are related by Boltzmann factors. Therefore if the first level in
          // the superlevel isn't inverted relative to the ground, none of the superlevel levels will be.
          const double sublevel_pop =
              population / superlevel_partfuncs[ion] * superlevel_boltzmann(nonemptymgi, element, ion, level);
          const double inversion_factor = (sublevel_pop / statweight) / (ground_pop / statweight_ground);
          if (inversion_factor > STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING) {
            assert_testmodeonly(ion_has_superlevel(element, ion));
            printout(
                "[debug] WARNING: superlevel pop inversion greater than factor %g: (g_pop %g)/(SL_first_level_pop %g) "
                "= %g is less than (g_sw %g)/(SL_first_level_sw %g) = %g for index %zu Z=%d ionstage %d level %d "
                "(factor %g inversion) - ",
                STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING, ground_pop, sublevel_pop,
                ground_pop / sublevel_pop, statweight_ground, statweight, statweight_ground / statweight, index,
                get_atomicnumber(element), ionstage, level, inversion_factor);

            if (inversion_factor > STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL) {
              printout("large pop inversion for superlevel - matrix solve failed\n");
              return false;
            }
            printout("relatively small pop inversion for superlevel so continue with NLTE solution\n");
          }
        }
      }
    }
  }
  return true;
}

// solve rate_matrix * x = balance_vector, so that popvec[i] = x[i] * pop_norm_factors[i]
// return true if the solution is successful, or false if the matrix is singular, contains negative or inverted
// populations.
[[nodiscard]] auto nltepop_matrix_solve(const int element, const int nonemptymgi, const gsl_matrix *rate_matrix,
                                        std::span<double> balance_vector, std::span<double> popvec,
                                        std::span<const double> pop_normfactors, const int max_nlte_dimension,
                                        const int first_ion_used, const int nions_used) -> bool {
  const size_t nlte_dimension = balance_vector.size();
  assert_always(pop_normfactors.size() == nlte_dimension);
  assert_always(rate_matrix->size1 == nlte_dimension);
  assert_always(rate_matrix->size2 == nlte_dimension);
  assert_always(std::cmp_greater_equal(max_nlte_dimension, nlte_dimension));
  if (lumatrix_is_singular(rate_matrix, element, first_ion_used, nions_used)) {
    printout("ERROR: NLTE matrix is singular for element Z=%d!\n", get_atomicnumber(element));
    return false;
  }

  // backing storage for gsl vectors
  THREADLOCALONHOST std::vector<double> vec_x;
  vec_x.reserve(max_nlte_dimension);
  vec_x.resize(nlte_dimension);
  gsl_vector x = gsl_vector_view_array(vec_x.data(), nlte_dimension).vector;

  THREADLOCALONHOST std::vector<double> vec_rate_matrix_LU_decomp;
  resize_exactly(vec_rate_matrix_LU_decomp, max_nlte_dimension * max_nlte_dimension);
  // make a copy of the rate matrix for the LU decomp
  gsl_matrix rate_matrix_LU_decomp =
      gsl_matrix_view_array(vec_rate_matrix_LU_decomp.data(), nlte_dimension, nlte_dimension).matrix;
  gsl_matrix_memcpy(&rate_matrix_LU_decomp, rate_matrix);

  THREADLOCALONHOST std::vector<size_t> vec_permutation;
  resize_exactly(vec_permutation, max_nlte_dimension * max_nlte_dimension);
  gsl_permutation_struct p{.size = nlte_dimension, .data = vec_permutation.data()};
  gsl_permutation_init(&p);

  auto gsl_balance_vector = gsl_vector_view_array(balance_vector.data(), nlte_dimension).vector;

  int s = 0;  // sign of the transformation
  assert_always(gsl_linalg_LU_decomp(&rate_matrix_LU_decomp, &p, &s) == GSL_SUCCESS);

#if !USE_SIMPSON_INTEGRATOR
  gsl_error_handler_t *previous_handler = gsl_set_error_handler(gsl_error_handler_printout);
#endif

  // solve matrix equation: rate_matrix * x = balance_vector for x (population vector)
  gsl_linalg_LU_solve(&rate_matrix_LU_decomp, &p, &gsl_balance_vector, &x);

#if !USE_SIMPSON_INTEGRATOR
  gsl_set_error_handler(previous_handler);
#endif

  constexpr double TOLERANCE = 1e-40;
  THREADLOCALONHOST std::vector<double> vec_work;
  resize_exactly(vec_work, max_nlte_dimension);
  gsl_vector gsl_work_vector = gsl_vector_view_array(vec_work.data(), nlte_dimension).vector;

  double error_best = -1.;

  // population solution vector with lowest error
  THREADLOCALONHOST std::vector<double> vec_x_best;
  vec_x_best.reserve(max_nlte_dimension);
  vec_x_best.resize(nlte_dimension);

  THREADLOCALONHOST std::vector<double> vec_residual;
  resize_exactly(vec_residual, max_nlte_dimension);
  gsl_vector gsl_vec_residual = gsl_vector_view_array(vec_residual.data(), nlte_dimension).vector;

  int iteration = 0;
  for (iteration = 0; iteration < 10; iteration++) {
    if (iteration > 0) {
      gsl_linalg_LU_refine(rate_matrix, &rate_matrix_LU_decomp, &p, &gsl_balance_vector, &x, &gsl_work_vector);
    }

    gsl_vector_memcpy(&gsl_vec_residual, &gsl_balance_vector);
    gsl_blas_dgemv(CblasNoTrans, 1.0, rate_matrix, &x, -1.0, &gsl_vec_residual);  // calculate Ax - b = residual
    const double error = fabs(gsl_vector_get(
        &gsl_vec_residual, gsl_blas_idamax(&gsl_vec_residual)));  // value of the largest absolute residual

    if (error < error_best || error_best < 0.) {
      std::ranges::copy(vec_x, vec_x_best.begin());
      error_best = error;
    }
    if (error < TOLERANCE) {
      break;
    }
  }
  if (error_best >= 0.) {
    if (error_best > 1e-8) {
      printout(
          "  NLTE solver matrix LU_refine: After %d iterations, best solution vector has a max residual of %g "
          "(WARNING!)\n",
          iteration, error_best);
    }

    std::ranges::copy(vec_x_best, vec_x.begin());
  }

  // get the unnormalised populations from the x solution vector and the normalisation factors
  for (size_t i = 0; i < nlte_dimension; i++) {
    popvec[i] = gsl_vector_get(&x, i) * pop_normfactors[i];
  }

  return solution_pops_are_valid(nonemptymgi, element, popvec, pop_normfactors, first_ion_used, nions_used);
}

auto can_remove_ion(const int element, const int ion, const int first_ion_used, const int nions_used,
                    const double nnelement, std::span<const double> popvec) -> bool {
  if (nions_used <= 1) {
    // can't remove an ion if there is only one ion stage used in the NLTE matrix
    printout("  WARNING: can't remove ion stage %d from NLTE matrix for element Z=%d as only one ion stage used\n",
             get_ionstage(element, ion), get_atomicnumber(element));
    return false;
  }

  const int max_ion_used = first_ion_used + nions_used - 1;
  assert_always(ion == first_ion_used || ion == max_ion_used);
  const std::string ionname = (ion == first_ion_used) ? "bottom" : "top";

  const int index_gs = get_nlte_vector_index(element, ion, 0, first_ion_used);
  const double ground_pop = popvec[index_gs];
  if (ground_pop / nnelement > NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION) {
    // ground pop relative to the element population must not exceed the limit
    printout(
        "  WARNING: %s ion ground state population too large to remove ion (ground_pop / "
        "nnelement (%g/%g) > (%g) NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION)\n",
        ionname.c_str(), ground_pop, nnelement, NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION);

    return false;
  }

  double nlte_excited_pop_sum = 0.0;
  // excited populations relative to the total element population must not exceed the limit
  const int nlevels_nlte_excited = get_nlevels_excited_nlte(element, ion);
  for (int level = 1; level <= nlevels_nlte_excited; level++) {
    const int index_level = get_nlte_vector_index(element, ion, level, first_ion_used);
    const double levelpop = popvec[index_level];
    nlte_excited_pop_sum += fabs(levelpop);
    if (levelpop / nnelement > NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION) {
      printout(
          "  WARNING: %s ion excited state (level %d) population too large to remove ion "
          "(nlte_excited_pop_bottom_ion "
          "/ nnelement (%g/%g) > (%g) NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION)\n",
          ionname.c_str(), level, levelpop, nnelement, NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION);
      return false;
    }
  }

  double superlevel_pop = 0.0;
  if (ion_has_superlevel(element, ion)) {
    // superlevel populations relative to the total element population must not exceed the limit
    const int index_superlevel = get_nlte_vector_index(element, ion, nlevels_nlte_excited + 1, first_ion_used);
    superlevel_pop = popvec[index_superlevel];
    if (superlevel_pop / nnelement > NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION) {
      printout(
          "  WARNING: %s ion superlevel population too large to remove ion (superlevel_pop_bottom_ion / "
          "nnelement (%g/%g) > (%g) NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION)\n",
          ionname.c_str(), superlevel_pop, nnelement, NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION);
      return false;
    }
  }

  printout(
      "  Passed checks: removing %s ion from NLTE matrix and attempting re-solve (ionstage %d, "
      "ground pop %g, fabs(excited pops) sum %g, superlevel pop %g, nnelement %g) matrix\n",
      ionname.c_str(), get_ionstage(element, ion), ground_pop, nlte_excited_pop_sum, superlevel_pop, nnelement);

  return true;
}

}  // anonymous namespace

void solve_nlte_pops_element(const int element, const int nonemptymgi, const int timestep, const int nlte_iter) {
  // solve the statistical balance equations to find NLTE level populations for all ions of an element
  // (ionisation balance follows from this too). Failure modes can lead to quasi-LTE populations being set instead.

  const int atomic_number = get_atomicnumber(element);
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);

  if (grid::get_elem_abundance(nonemptymgi, element) <= 0.) {
    // abundance of this element is zero, so do not store any NLTE populations
    printout(
        "Not solving for NLTE populations in cell %d at timestep %d NLTE iteration %d for element Z=%d due to zero "
        "abundance\n",
        modelgridindex, timestep, nlte_iter, atomic_number);

    nltepop_reset_element(nonemptymgi, element);
    return;
  }

  const double cell_Te = grid::get_Te(nonemptymgi);

  if (cell_Te == MINTEMP) {
    printout(
        "Not solving for NLTE populations in cell %d at timestep %d NLTE iteration %d for element Z=%d due to low "
        "temperature Te=MINTEMP=%g K\n",
        modelgridindex, timestep, atomic_number, nlte_iter, cell_Te);
    set_element_pops_lte(nonemptymgi, element);
    return;
  }

  const auto sys_time_start_nltesolver = std::time(nullptr);

  const double t_mid = globals::timesteps[timestep].mid;
  const int nions = get_nions(element);
  const double nnelement = grid::get_elem_numberdens(nonemptymgi, element);

  printout(
      "Solving for NLTE populations in cell %d at timestep %d NLTE iteration %d for element Z=%d (mass fraction "
      "%.2e, nnelement %.2e cm^-3)\n",
      modelgridindex, timestep, nlte_iter, atomic_number, grid::get_elem_abundance(nonemptymgi, element), nnelement);
  const auto superlevel_partfuncs = get_element_superlevelpartfuncs(nonemptymgi, element);
  int nions_used = nions;
  int first_ion_used = 0;
  bool matrix_solve_success = false;
  const auto max_nlte_dimension = get_max_nlte_dimension();

  // will hold the un-normalised population densities [cm^-3]
  THREADLOCALONHOST std::vector<double> popvec;
  popvec.reserve(max_nlte_dimension);

  THREADLOCALONHOST RateMatrices rate_matrices{max_nlte_dimension};

  bool matrix_solve_required = true;
  while (matrix_solve_required) {
    const int nlte_dimension = get_element_nlte_dimension(element, first_ion_used, nions_used);
    rate_matrices.set_used_dimension(nlte_dimension);
    popvec.resize(nlte_dimension);

    const int max_ion_used = first_ion_used + nions_used - 1;
    const auto ions = std::views::iota(first_ion_used, max_ion_used + 1);
    std::for_each(ions.begin(), ions.end(), [&](const auto ion) {
      const int nlevels = get_nlevels(element, ion);
      const int level_superlevel_start = get_nlevels_excited_nlte(element, ion) + 1;

      auto s_renorm = std::vector<double>(nlevels);
      std::fill_n(s_renorm.begin(), level_superlevel_start, 1.);

      // nlevels_nlte is the lowest superlevel index
      for (int level = level_superlevel_start; level < nlevels; level++) {
        s_renorm[level] = superlevel_boltzmann(nonemptymgi, element, ion, level) / superlevel_partfuncs[ion];
      }

      nltepop_matrix_add_boundbound(nonemptymgi, element, ion, t_mid, s_renorm, rate_matrices, first_ion_used);

      if (ion < max_ion_used) {
        nltepop_matrix_add_ionisation(nonemptymgi, element, ion, s_renorm, rate_matrices, first_ion_used, nions_used);
        if (NT_ON) {
          nltepop_matrix_add_nt_ionisation(nonemptymgi, element, ion, s_renorm, rate_matrices, first_ion_used,
                                           nions_used);
        }
        nltepop_matrix_add_autoionisation(nonemptymgi, element, ion, s_renorm, rate_matrices, first_ion_used,
                                          nions_used);
      }
    });

    // replace the zeroth row of the matrix and balance vector with the normalisation
    // constraint (sum of levelpops = total element population)

    auto rate_matrix = rate_matrices.get_summed_rate_matrix();
    gsl_vector_view first_row_view = gsl_matrix_row(&rate_matrix, 0);
    gsl_vector_set_all(&first_row_view.vector, 1.0);

    THREADLOCALONHOST std::vector<double> balance_vector;
    balance_vector.reserve(max_nlte_dimension);
    balance_vector.resize(nlte_dimension);
    std::ranges::fill(balance_vector, 0.0);  // statistical equilibrium means balance vector is zero
    // except for zeroth element used to normalise the total element population
    balance_vector[0] = nnelement;

    if (FORCE_SAHA_ION_BALANCE(atomic_number)) {
      const auto ionfractions = calculate_ionfractions(element, nonemptymgi, grid::get_nne(nonemptymgi), true);
      const int uppermost_ion = static_cast<int>(ionfractions.size() - 1);
      for (int ion = 1; ion <= uppermost_ion; ion++) {
        // replace matrix row for ion's ground state with sum of this ion's level populations is equal to the ion
        // population
        const double nnion = nnelement * ionfractions[ion];
        const int index_ion_ground = get_nlte_vector_index(element, ion, 0, first_ion_used);
        const int index_ion_toplevel = get_nlte_vector_index(element, ion, get_nlevels(element, ion), first_ion_used);
        gsl_vector_view ion_ground_row_view = gsl_matrix_row(&rate_matrix, index_ion_ground);
        gsl_vector_set_all(&ion_ground_row_view.vector, 0.);
        for (int index = index_ion_ground; index <= index_ion_toplevel; index++) {
          gsl_vector_set(&ion_ground_row_view.vector, index, 1.);
        }

        balance_vector[get_nlte_vector_index(element, ion, index_ion_ground, first_ion_used)] = nnion;
      }
    }

    // calculate the normalisation factors and apply them to the matrix
    // columns and balance vector elements
    THREADLOCALONHOST std::vector<double> pop_norm_factors;
    pop_norm_factors.reserve(max_nlte_dimension);
    pop_norm_factors.resize(nlte_dimension);
    std::ranges::fill(pop_norm_factors, 1.0);
    nltepop_matrix_normalise(nonemptymgi, element, &rate_matrix, pop_norm_factors, first_ion_used, nions_used);

    matrix_solve_success = nltepop_matrix_solve(element, nonemptymgi, &rate_matrix, balance_vector, popvec,
                                                pop_norm_factors, max_nlte_dimension, first_ion_used, nions_used);

    matrix_solve_required = false;  // will be set to true if we need to retry with a different ion range
    if (matrix_solve_success) {
      if (nions_used < nions) {
        printout(
            "successfully solved NLTE matrix when reducing ions used for element to Z=%d ionstage=%d to ionstage=%d\n",
            atomic_number, get_ionstage(element, first_ion_used), get_ionstage(element, max_ion_used));
      }
    } else if (NLTE_LIMIT_ION_STAGES_AFTER_FAILURE) {
      printout("  WARNING: NLTE matrix solution failed for element Z=%d using ionstage %d to %d\n", atomic_number,
               get_ionstage(element, first_ion_used), get_ionstage(element, max_ion_used));

      if (can_remove_ion(element, max_ion_used, first_ion_used, nions_used, nnelement, popvec)) {
        // we can remove the top ion from the solution
        nions_used--;
        matrix_solve_required = true;
      } else if (can_remove_ion(element, first_ion_used, first_ion_used, nions_used, nnelement, popvec)) {
        // we can remove the bottom ion from the solution
        first_ion_used++;
        nions_used--;
        matrix_solve_required = true;
      } else {
        printout(
            "  WARNING: can't remove top or bottom ion stage from NLTE equations, therefore unable to find an NLTE "
            "solution for this element\n");
      }
    }
  }

  if (!matrix_solve_success) {
    printout(
        "WARNING: Can't solve for NLTE populations in cell %d at timestep %d for element Z=%d due to singular "
        "matrix, negative pop or large pop inversion and unable to recover solution by reducing ion range. "
        "Attempting to use LTE solution instead\n",
        modelgridindex, timestep, atomic_number);
    set_element_pops_lte(nonemptymgi, element);
  } else {
    // check calculated NLTE populations are valid
    for (const auto pop : popvec) {
      assert_always(std::isfinite(pop));
      assert_always(pop >= 0.);
    }

    // set the ground level, excited level and possible superlevel populations for this element
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels_excited_nlte = get_nlevels_excited_nlte(element, ion);
      const int index_gs = get_nlte_vector_index(element, ion, 0, first_ion_used);

      if (ion < first_ion_used || ion >= (first_ion_used + nions_used)) {
        printout("  WARNING: Z=%d ionstage %d removed from NLTE rate matrix. Setting all levelpops for ion to zero \n",
                 get_atomicnumber(element), get_ionstage(element, ion));

        grid::ion_groundlevelpops_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) +
                                           get_uniqueionindex(element, ion)] = 0.;

        for (int level = 1; level <= nlevels_excited_nlte; level++) {
          set_nlte_levelpop_over_rho(nonemptymgi, element, ion, level, 0.);
        }

        if (ion_has_superlevel(element, ion)) {
          set_nlte_superlevelpop_over_rho_over_slpartfunc(nonemptymgi, element, ion, 0.);
        }
      } else {
        grid::ion_groundlevelpops_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) +
                                           get_uniqueionindex(element, ion)] = static_cast<float>(popvec[index_gs]);

        for (int level = 1; level <= nlevels_excited_nlte; level++) {
          const int index = get_nlte_vector_index(element, ion, level, first_ion_used);
          set_nlte_levelpop_over_rho(nonemptymgi, element, ion, level, popvec[index] / grid::get_rho(nonemptymgi));
        }

        if (ion_has_superlevel(element, ion)) {
          const int index_sl = get_nlte_vector_index(element, ion, nlevels_excited_nlte + 1, first_ion_used);
          set_nlte_superlevelpop_over_rho_over_slpartfunc(
              nonemptymgi, element, ion, popvec[index_sl] / grid::get_rho(nonemptymgi) / superlevel_partfuncs[ion]);
        }
      }
    }
    calculate_cellpartfuncts(nonemptymgi, element);

    const double elem_pop_matrix = std::accumulate(popvec.begin(), popvec.end(), 0.0,
                                                   [](const double sum, const double pop) { return sum + fabs(pop); });
    const double elem_pop_error_percent = fabs((nnelement / elem_pop_matrix) - 1) * 100;
    if (elem_pop_error_percent > 1.0) {
      printout(
          "  WARNING: timestep %d nlteiter %d Z=%d element population is: %g (from abundance) and %g (from matrix "
          "solution), error: %.2f%%. Forcing element pops to LTE.\n",
          timestep, nlte_iter, atomic_number, nnelement, elem_pop_matrix, elem_pop_error_percent);
      set_element_pops_lte(nonemptymgi, element);
    }

    // output NLTE stats every nth timestep for the first NLTE iteration only
    if ((timestep % 5 == 0) && (nlte_iter == 0)) {
      print_element_rates_summary(element, modelgridindex, timestep, nlte_iter, popvec, rate_matrices, first_ion_used,
                                  nions_used);
    }

    // detailed levels stats (very verbose, only for debugging)
    // const bool print_detailed_stats = ((atomic_number == 26) && ((timestep % 5) == 0) && (nlte_iter == 0))
    const bool print_detailed_level_stats = false;
    if (print_detailed_level_stats) {
      const int ionstage = 2;
      const int ion = ionstage - get_ionstage(element, 0);

      const int nlevels = get_nlevels_excited_nlte(element, ion) + (ion_has_superlevel(element, ion) ? 1 : 0);
      for (int level = 0; level <= nlevels; level++) {
        print_level_rates(nonemptymgi, timestep, element, ion, level, popvec, rate_matrices, first_ion_used,
                          nions_used);
      }
    }
  }

  if (const auto duration_nltesolver = std::time(nullptr) - sys_time_start_nltesolver; duration_nltesolver > 2) {
    printout("NLTE population solver call for Z=%d took %ld seconds\n", get_atomicnumber(element), duration_nltesolver);
  }
}

// Get a Boltzman factor for a level within the super level (combined Non-LTE level)
__host__ __device__ auto superlevel_boltzmann(const int nonemptymgi, const int element, const int ion, const int level)
    -> double {
  assert_testmodeonly(level_isinsuperlevel(element, ion, level));
  const int level_superlevel_start = get_nlevels_excited_nlte(element, ion) + 1;
  const double T_exc = LTEPOP_EXCITATION_USE_TJ ? grid::get_TJ(nonemptymgi) : grid::get_Te(nonemptymgi);
  const double E_level = epsilon(element, ion, level);
  const double E_superlevel = epsilon(element, ion, level_superlevel_start);

  return stat_weight(element, ion, level) / stat_weight(element, ion, level_superlevel_start) *
         exp(-(E_level - E_superlevel) / KB / T_exc);
}

void nltepop_open_file(const int my_rank) {
  nlte_file = fstream_required(std::format("nlte_{:04d}.out", my_rank), std::ios::out | std::ios::trunc);
  nlte_file << "timestep modelgridindex Z ionstage level n_LTE n_NLTE ion_popfrac\n";
}

void nltepop_write_to_file(const int nonemptymgi, const int timestep) {
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  if (globals::lte_iteration || grid::modelgrid[nonemptymgi].thick == 1) {  // NLTE solver hasn't been run yet
    return;
  }

  for (int element = 0; element < get_nelements(); element++) {
    if (!elem_has_nlte_levels(element)) {
      continue;
    }

    const auto superlevel_partfuncs = get_element_superlevelpartfuncs(nonemptymgi, element);

    for (int ion = 0; ion < get_nions(element); ion++) {
      const int nlevels_excited_nlte = get_nlevels_excited_nlte(element, ion);
      const int nsuperlevels = ion_has_superlevel(element, ion) ? 1 : 0;
      const auto nnion = get_nnion(nonemptymgi, element, ion);

      for (int level = 0; level <= nlevels_excited_nlte + nsuperlevels; level++) {
        double nnlevellte = calculate_levelpop_boltzmann(nonemptymgi, element, ion, level);
        double nnlevelnlte{NAN};

        nlte_file << timestep << ' ' << modelgridindex << ' ' << get_atomicnumber(element) << ' '
                  << get_ionstage(element, ion) << ' ';
        if (level <= nlevels_excited_nlte) {
          nlte_file << level << ' ';

          if (level == 0) {
            nnlevelnlte = get_groundlevelpop(nonemptymgi, element, ion);
          } else {
            nnlevelnlte =
                get_nlte_levelpop_over_rho(nonemptymgi, element, ion, level) * grid::modelgrid[nonemptymgi].rho;
          }
        } else {
          // superlevel, so add the populations of all other levels in the superlevel
          const double slpopfactor = get_nlte_superlevelpop_over_rho_over_slpartfunc(nonemptymgi, element, ion) *
                                     grid::modelgrid[nonemptymgi].rho;

          nnlevellte = 0;
          nlte_file << -1 << ' ';
          for (int level_sl = nlevels_excited_nlte + 1; level_sl < get_nlevels(element, ion); level_sl++) {
            if (level_isinsuperlevel(element, ion, level_sl)) {
              nnlevellte += calculate_levelpop_boltzmann(nonemptymgi, element, ion, level_sl);
            }
          }

          nnlevelnlte = slpopfactor * superlevel_partfuncs[ion];
        }

        nlte_file << std::format("{:.5e} {:.5e} {:.5e}\n", nnlevellte, nnlevelnlte, nnlevelnlte / nnion);
      }
    }
  }

  nlte_file.flush();
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
      for (int ion = 0; ion < get_nions(element); ion++) {
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
  assert_always(code_check == 75618527);

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

__host__ __device__ auto get_nlte_levelpop_over_rho(const int nonemptymgi, const int element, const int ion,
                                                    const int level) -> double {
  assert_testmodeonly(level > 0);  // ground state is stored separately
  assert_testmodeonly(level <= get_nlevels_excited_nlte(element, ion));
  return grid::nltepops_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * globals::total_nlte_levels) +
                                 globals::elements[element].ions[ion].first_nlte + level - 1];
}

[[nodiscard]] __host__ __device__ auto get_nlte_superlevelpop_over_rho_over_slpartfunc(const int nonemptymgi,
                                                                                       const int element, const int ion)
    -> double {
  assert_testmodeonly(ion_has_superlevel(element, ion));
  const int sl_nlte_index = globals::elements[element].ions[ion].first_nlte + get_nlevels_excited_nlte(element, ion);
  return grid::nltepops_allcells[(nonemptymgi * globals::total_nlte_levels) + sl_nlte_index];
}

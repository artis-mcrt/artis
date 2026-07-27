// The full NLTE level population solver: assembles the statistical equilibrium rate matrix
// (radiative and collisional bound-bound and bound-free rates, plus non-thermal ionisation)
// for all levels of all ions of an element and solves for the level populations.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <fstream>
#include <ios>
#include <numeric>
#include <print>
#include <ranges>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#ifdef EIGEN_OFF
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_double.h>
#else
#include <Eigen/Core>
#include <Eigen/Dense>
#endif
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "kpkt.h"
#include "ltepop.h"
#include "macroatom.h"
#include "mpi_logging.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "ratecoeff.h"
#include "sn3d.h"

static_assert(STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING >= 1);
static_assert(STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING <
              STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL);
namespace {
std::fstream nlte_file;

struct RateMatrices {
  int used_nlte_dimension{0};

  std::vector<double> rad_bb;
  std::vector<double> coll_bb;
  std::vector<double> ntcoll_bb;
  std::vector<double> rad_bf;
  std::vector<double> coll_bf;
  std::vector<double> ntcoll_bf;
  std::vector<double> autoion;

  explicit RateMatrices(int max_nlte_dimension) {
    // allocation of the maximum required size is done once,
    // while the used_nlte_dimension is set later
    const auto max_dim_squared = static_cast<ptrdiff_t>(max_nlte_dimension) * max_nlte_dimension;
    summed_rates.reserve(max_dim_squared);
    rad_bb.reserve(max_dim_squared);
    coll_bb.reserve(max_dim_squared);
    ntcoll_bb.reserve(max_dim_squared);
    rad_bf.reserve(max_dim_squared);
    coll_bf.reserve(max_dim_squared);
    ntcoll_bf.reserve(max_dim_squared);
    autoion.reserve(max_dim_squared);
  }

  void set_used_dimension(int used_nlte_dimension_in) {
    used_nlte_dimension = used_nlte_dimension_in;
    const auto used_dim_squared = used_nlte_dimension * used_nlte_dimension;

    assert_always(std::cmp_less_equal(used_dim_squared, summed_rates.capacity()));
    summed_rates.resize(used_dim_squared);
    assert_always(std::cmp_less_equal(used_dim_squared, rad_bb.capacity()));
    rad_bb.resize(used_dim_squared);
    assert_always(std::cmp_less_equal(used_dim_squared, coll_bb.capacity()));
    coll_bb.resize(used_dim_squared);
    assert_always(std::cmp_less_equal(used_dim_squared, ntcoll_bb.capacity()));
    ntcoll_bb.resize(used_dim_squared);
    assert_always(std::cmp_less_equal(used_dim_squared, rad_bf.capacity()));
    rad_bf.resize(used_dim_squared);
    assert_always(std::cmp_less_equal(used_dim_squared, coll_bf.capacity()));
    coll_bf.resize(used_dim_squared);
    assert_always(std::cmp_less_equal(used_dim_squared, ntcoll_bf.capacity()));
    ntcoll_bf.resize(used_dim_squared);
    assert_always(std::cmp_less_equal(used_dim_squared, autoion.capacity()));
    autoion.resize(used_dim_squared);

    std::ranges::fill(rad_bb, 0.);
    std::ranges::fill(coll_bb, 0.);
    std::ranges::fill(ntcoll_bb, 0.);
    std::ranges::fill(rad_bf, 0.);
    std::ranges::fill(coll_bf, 0.);
    std::ranges::fill(ntcoll_bf, 0.);
    std::ranges::fill(autoion, 0.);
  }

  [[nodiscard]] auto get_summed_rate_matrix() -> std::span<double> {
    // sum the matrices for each transition type to get a total rate matrix
    assert_always(summed_rates.size() == rad_bb.size());
    assert_always(summed_rates.size() == coll_bb.size());
    assert_always(summed_rates.size() == ntcoll_bb.size());
    assert_always(summed_rates.size() == rad_bf.size());
    assert_always(summed_rates.size() == coll_bf.size());
    assert_always(summed_rates.size() == ntcoll_bf.size());
    assert_always(summed_rates.size() == autoion.size());
    for (auto i = 0U; i < summed_rates.size(); i++) {
      summed_rates[i] = rad_bb[i] + coll_bb[i] + ntcoll_bb[i] + rad_bf[i] + coll_bf[i] + ntcoll_bf[i] + autoion[i];
    }

    return summed_rates;
  }

 private:
  // keep this private to ensure that the sum has been calculated before passing a reference
  std::vector<double> summed_rates;
};

// this is the matrix/vector index for the NLTE solver that is handling all ions of a single element
auto get_nlte_vector_index(const int element, const int ion, const int level, const int first_ion_used) -> int {
  // have to convert from nlte_pops index to nlte_vector index
  // the difference is that nlte vectors apply to a single element and include ground and autoionising states
  int offset_autoion = 0;
  for (int dion = first_ion_used; dion < ion; dion++) {
    if (ion_has_superlevel(element, dion)) {
      offset_autoion += get_nlevels_autoion(element, dion);
    }
  }
  assert_testmodeonly(first_ion_used >= 0);
  assert_testmodeonly(first_ion_used < get_nions(element));
  // The (ion - first_ion_used) term accounts for the ground state population indices that are not counted in the NLTE
  // array
  const int gs_index = get_allnltelevelsindexstart(element, ion) -
                       get_allnltelevelsindexstart(element, first_ion_used) + (ion - first_ion_used) + offset_autoion;
  const int first_autoion = get_nlevels_nonautoion(element, ion);
  if (level >= first_autoion) {
    if (ion_has_superlevel(element, ion)) {
      // autoionising levels are stored after the ground state, the NLTE levels, and the superlevel
      return (gs_index + get_nlevels_excited_nlte(element, ion) + level - first_autoion + 2);
    }
    // without a superlevel every level of the ion has its own slot (see get_element_nlte_dimension),
    // so autoionising levels must not collapse onto the (non-existent) superlevel index
    return gs_index + level;
  }
  // add in level or superlevel number
  const int level_index =
      gs_index + (is_nlte(element, ion, level) ? level : (get_nlevels_excited_nlte(element, ion) + 1));
  return level_index;
}

[[nodiscard]] auto get_ion_level_of_nlte_vector_index(const std::ptrdiff_t index, const int element,
                                                      const int first_ion_used, const int nions_used)
    -> std::tuple<int, int> {
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

[[nodiscard]] auto get_total_rate(const ptrdiff_t index_selected, const std::span<const double> rate_matrix,
                                  const std::span<const double> popvec, const bool into_level,
                                  const bool only_levels_below, const bool only_levels_above) -> double {
  double total_rate = 0.;
  assert_always(!only_levels_below || !only_levels_above);
  const auto nlte_dimension = std::ssize(popvec);
  assert_testmodeonly(std::ssize(rate_matrix) == (nlte_dimension * nlte_dimension));
  assert_testmodeonly(index_selected >= 0);
  assert_testmodeonly(index_selected < nlte_dimension);

  if (into_level) {
    // find rate into selected level

    // multiply incoming rate coefficients by their corresponding populations to get rates
    if (!only_levels_above) {  // add levels below
      for (auto index = 0Z; index < index_selected; index++) {
        total_rate += rate_matrix[(index_selected * nlte_dimension) + index] * popvec[index];
      }
    }

    if (!only_levels_below) {  // add levels above
      for (auto index = index_selected + 1; index < nlte_dimension; index++) {
        total_rate += rate_matrix[(index_selected * nlte_dimension) + index] * popvec[index];
      }
    }
  } else {
    // find rate out of selected level

    // multiply outgoing rate coefficients by the population of the selected level to get rates

    if (!only_levels_above) {  // add levels below
      for (auto index = 0; index < index_selected; index++) {
        total_rate += rate_matrix[(index * nlte_dimension) + index_selected];
      }
    }

    if (!only_levels_below) {  // add levels above
      for (auto index = index_selected + 1; index < nlte_dimension; index++) {
        total_rate += rate_matrix[(index * nlte_dimension) + index_selected];
      }
    }

    total_rate *= popvec[index_selected];
  }

  return total_rate;
}

auto get_total_rate_in(const int index_to, const std::span<const double> rate_matrix,
                       const std::span<const double> popvec) -> double {
  return get_total_rate(index_to, rate_matrix, popvec, true, false, false);
}

auto get_total_rate_out(const int index_from, const std::span<const double> rate_matrix,
                        const std::span<const double> popvec) -> double {
  return get_total_rate(index_from, rate_matrix, popvec, false, false, false);
}

void print_level_rates_summary(const int element, const int selected_ion, const int selected_level,
                               std::span<const double> popvec, const RateMatrices& rate_matrices,
                               const int first_ion_used) {
  const int selected_index = get_nlte_vector_index(element, selected_ion, selected_level, first_ion_used);

  for (int i = 0; i <= 3; i++) {
    // rates in from above, in from below, out to above, out to below
    if (i == 0) {
      const int nlevels_nlte = get_nlevels_excited_nlte(element, selected_ion);
      if (ion_has_superlevel(element, selected_ion) && (selected_level == nlevels_nlte + 1)) {
        printlog("      superlevel ");
      } else {
        printlog("    level{:7} ", selected_level);
      }
      printlog(" {:10.2e} ", popvec[selected_index]);
    } else {
      printlog("                             ");
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
    const double autoion_total =
        get_total_rate(selected_index, rate_matrices.autoion, popvec, into_level, only_levels_below, only_levels_above);

    if (into_level) {
      // into this level
      printlog(" from ");
    } else {
      // out of this level
      printlog("   to ");
    }

    if (only_levels_below) {
      printlog("below ");
    } else {
      printlog("above ");
    }

    printlnlog("{:10.2e} {:10.2e} {:10.2e} {:10.2e} {:10.2e} {:10.2e} {:10.2e}", rad_bb_total, coll_bb_total,
               ntcoll_bb_total, rad_bf_total, coll_bf_total, ntcoll_bf_total, autoion_total);
  }
}

void print_element_rates_summary(const int element, const int modelgridindex, const int timestep, const int nlte_iter,
                                 std::span<const double> popvec, const RateMatrices& rate_matrices,
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
        printlnlog(
            "  modelgridindex {} timestep {} NLTE iteration {} Te {:g} nne {:g}: NLTE summary for Z={} ionstage {}:",
            modelgridindex, timestep, nlte_iter, grid::get_Te(nonemptymgi), grid::get_nne(nonemptymgi), atomic_number,
            ionstage);
        printlnlog(
            "                         pop       rates     bb_rad     bb_col   bb_ntcol     bf_rad     bf_col   "
            "bf_ntcol autoion");
      }

      print_level_rates_summary(element, ion, level, popvec, rate_matrices, first_ion_used);
    }
  }
}

void print_level_rates(const int nonemptymgi, const int timestep, const int element, const int selected_ion,
                       const int selected_level, const std::span<const double> popvec,
                       const RateMatrices& rate_matrices, const int first_ion_used, const int nions_used) {
  // very detailed output of each transition for a particular level

  assert_always(selected_ion >= first_ion_used);
  assert_always(selected_ion < first_ion_used + nions_used - 1);
  assert_always(selected_level <= (get_nlevels_excited_nlte(element, selected_ion) +
                                   (ion_has_superlevel(element, selected_ion) ? 1 : 0)));

  const auto nltedim = std::ssize(popvec);
  const int atomic_number = get_atomicnumber(element);
  const int selected_ionstage = get_ionstage(element, selected_ion);
  const auto selected_index = get_nlte_vector_index(element, selected_ion, selected_level, first_ion_used);
  const double pop_selectedlevel = popvec[selected_index];
  printlnlog(
      "timestep {} cell {} Te {:g} nne {:g} NLTE level diagnostics for Z={} ionstage {} level {} rates into and out of "
      "this level",
      timestep, grid::get_mgi_of_nonemptymgi(nonemptymgi), grid::get_Te(nonemptymgi), grid::get_nne(nonemptymgi),
      atomic_number, selected_ionstage, selected_level);

  const double rad_bb_in_total = get_total_rate_in(selected_index, rate_matrices.rad_bb, popvec);
  const double coll_bb_in_total = get_total_rate_in(selected_index, rate_matrices.coll_bb, popvec);
  const double ntcoll_bb_in_total = get_total_rate_in(selected_index, rate_matrices.ntcoll_bb, popvec);
  const double rad_bf_in_total = get_total_rate_in(selected_index, rate_matrices.rad_bf, popvec);
  const double coll_bf_in_total = get_total_rate_in(selected_index, rate_matrices.coll_bf, popvec);
  const double ntcoll_bf_in_total = get_total_rate_in(selected_index, rate_matrices.ntcoll_bf, popvec);
  const double total_rate_in =
      rad_bb_in_total + coll_bb_in_total + ntcoll_bb_in_total + rad_bf_in_total + coll_bf_in_total + ntcoll_bf_in_total;
  printlnlog(
      "  TOTAL rates in:             rad_bb_in  {:8.2e} coll_bb_in  {:8.2e} ntcoll_bb_in  {:8.2e} rad_bf_in  {:8.2e} "
      "coll_bf_in  {:8.2e} ntcoll_bf_in  {:8.2e}",
      rad_bb_in_total, coll_bb_in_total, ntcoll_bb_in_total, rad_bf_in_total, coll_bf_in_total, ntcoll_bf_in_total);

  const double rad_bb_out_total = get_total_rate_out(selected_index, rate_matrices.rad_bb, popvec);
  const double coll_bb_out_total = get_total_rate_out(selected_index, rate_matrices.coll_bb, popvec);
  const double ntcoll_bb_out_total = get_total_rate_out(selected_index, rate_matrices.ntcoll_bb, popvec);
  const double rad_bf_out_total = get_total_rate_out(selected_index, rate_matrices.rad_bf, popvec);
  const double coll_bf_out_total = get_total_rate_out(selected_index, rate_matrices.coll_bf, popvec);
  const double ntcoll_bf_out_total = get_total_rate_out(selected_index, rate_matrices.ntcoll_bf, popvec);
  const double total_rate_out = rad_bb_out_total + coll_bb_out_total + ntcoll_bb_out_total + rad_bf_out_total +
                                coll_bf_out_total + ntcoll_bf_out_total;
  printlnlog(
      "  TOTAL rates out:            rad_bb_out {:8.2e} coll_bb_out {:8.2e} ntcoll_bb_out {:8.2e} rad_bf_out {:8.2e} "
      "coll_bf_out {:8.2e} ntcoll_bf_out {:8.2e}",
      rad_bb_out_total, coll_bb_out_total, ntcoll_bb_out_total, rad_bf_out_total, coll_bf_out_total,
      ntcoll_bf_out_total);

  for (auto index = 0Z; index < nltedim; index++) {
    if (index == selected_index) {
      continue;
    }
    const auto [ion, level] = get_ion_level_of_nlte_vector_index(index, element, first_ion_used, nions_used);
    const int ionstage = get_ionstage(element, ion);
    // in means populating the selected level, out means depopulating the selected level
    const double pop = popvec[index];
    const auto ij_selectedindex_index = (selected_index * nltedim) + index;
    const auto ij_index_selectedindex = (index * nltedim) + selected_index;
    const double rad_bb_in = rate_matrices.rad_bb[ij_selectedindex_index] * pop;
    const double rad_bb_out = rate_matrices.rad_bb[ij_index_selectedindex] * pop_selectedlevel;
    const double coll_bb_in = rate_matrices.coll_bb[ij_selectedindex_index] * pop;
    const double coll_bb_out = rate_matrices.coll_bb[ij_index_selectedindex] * pop_selectedlevel;
    const double ntcoll_bb_in = rate_matrices.ntcoll_bb[ij_selectedindex_index] * pop;
    const double ntcoll_bb_out = rate_matrices.ntcoll_bb[ij_index_selectedindex] * pop_selectedlevel;
    const double rad_bf_in = rate_matrices.rad_bf[ij_selectedindex_index] * pop;
    const double rad_bf_out = rate_matrices.rad_bf[ij_index_selectedindex] * pop_selectedlevel;
    const double coll_bf_in = rate_matrices.coll_bf[ij_selectedindex_index] * pop;
    const double coll_bf_out = rate_matrices.coll_bf[ij_index_selectedindex] * pop_selectedlevel;
    const double ntcoll_bf_in = rate_matrices.ntcoll_bf[ij_selectedindex_index] * pop;
    const double ntcoll_bf_out = rate_matrices.ntcoll_bf[ij_index_selectedindex] * pop_selectedlevel;

    const bool nonzero_rate_in = (fabs(rad_bb_in) > 0. || fabs(coll_bb_in) > 0. || fabs(ntcoll_bb_in) > 0. ||
                                  fabs(rad_bf_in) > 0. || fabs(coll_bf_in) > 0. || fabs(ntcoll_bf_in) > 0.);
    const bool nonzero_rate_out = (fabs(rad_bb_out) > 0. || fabs(coll_bb_out) > 0. || fabs(ntcoll_bb_out) > 0. ||
                                   fabs(rad_bf_out) > 0. || fabs(coll_bf_out) > 0. || fabs(ntcoll_bf_out) > 0.);
    if (nonzero_rate_in || nonzero_rate_out) {
      const double epsilon_trans = fabs(epsilon(element, ion, level) - epsilon(element, selected_ion, selected_level));
      const double nu_trans = epsilon_trans / H;
      const double lambda = 1e8 * CLIGHT / nu_trans;  // [Angstroms]
      const double level_rate_in = rad_bb_in + coll_bb_in + ntcoll_bb_in + rad_bf_in + coll_bf_in + ntcoll_bf_in;
      const double level_rate_out = rad_bb_out + coll_bb_out + ntcoll_bb_out + rad_bf_out + coll_bf_out + ntcoll_bf_out;
      const double level_percent_in = level_rate_in / total_rate_in * 100.;
      const double level_percent_out = level_rate_out / total_rate_out * 100.;

      printlnlog(
          "  ionstage {} level {:4} ({:5.1f}% of in)  rad_bb_in  {:8.2e} coll_bb_in  {:8.2e} ntcoll_bb_in  {:8.2e} "
          "rad_bf_in  {:8.2e} coll_bf_in  {:8.2e} ntcoll_bf_in  {:8.2e} lambda {:6.0f}",
          ionstage, level, level_percent_in, rad_bb_in, coll_bb_in, ntcoll_bb_in, rad_bf_in, coll_bf_in, ntcoll_bf_in,
          lambda);
      printlnlog(
          "  ionstage {} level {:4} ({:5.1f}% of out) rad_bb_out {:8.2e} coll_bb_out {:8.2e} ntcoll_bb_out {:8.2e} "
          "rad_bf_out {:8.2e} coll_bf_out {:8.2e} ntcoll_bf_out {:8.2e} lambda {:6.0f}",
          ionstage, level, level_percent_out, rad_bb_out, coll_bb_out, ntcoll_bb_out, rad_bf_out, coll_bf_out,
          ntcoll_bf_out, lambda);
    }
  }
  printlnlog("");
}

void nltepop_reset_element(const int nonemptymgi, const int element) {
  for (int ion = 0; ion < get_nions(element); ion++) {
    const int nlte_start = get_allnltelevelsindexstart(element, ion);
    std::fill_n(&nltepops_allcells[(nonemptymgi * globals::total_nlte_levels) + nlte_start],
                get_nlevels_excited_nlte(element, ion) + (ion_has_superlevel(element, ion) ? 1 : 0), -1.);
  }
}

auto get_element_superlevelpartfuncs(const int nonemptymgi, const int element) -> std::vector<double> {
  std::vector<double> superlevel_partfuncs;
  reserve_resize(superlevel_partfuncs, get_nions(element));
  for (int ion = 0; ion < get_nions(element); ion++) {
    if (ion_has_superlevel(element, ion)) {
      superlevel_partfuncs[ion] = 0.;
      for (int level = get_nlevels_excited_nlte(element, ion) + 1; level < get_nlevels_nonautoion(element, ion);
           level++) {
        superlevel_partfuncs[ion] += superlevel_boltzmann(nonemptymgi, element, ion, level);
      }
    } else {
      superlevel_partfuncs[ion] = -1.;
    }
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
      // if it has a superlevel then need +1 for ground state, +1 for superlevel and add autoionising levels
      nlte_dimension += get_nlevels_excited_nlte(element, ion) + get_nlevels_autoion(element, ion) + 2;
    } else {  // if it doesn't have a superlevel
      nlte_dimension += get_nlevels(element, ion);
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
                                   const std::span<const double> s_renorm, RateMatrices& rate_matrices,
                                   const int first_ion_used) {
  const auto T_e = grid::get_Te(nonemptymgi);
  const auto clumpednne = grid::get_clumpfactor(nonemptymgi) * grid::get_nne(nonemptymgi);
  const int nlevels = get_nlevels(element, ion);
  const auto ionuniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);
  const auto nlte_dimension = rate_matrices.used_nlte_dimension;

  const auto levels = std::views::iota(0, nlevels);
  std::for_each(levels.begin(), levels.end(), [&](const auto level) {
    const int level_index = get_nlte_vector_index(element, ion, level, first_ion_used);
    const auto matrix_index_level_level = (level_index * nlte_dimension) + level_index;
    const auto uniquelevelindex = ionuniquelevelindexstart + level;
    const double epsilon_level = epsilon(uniquelevelindex);
    const double statweight = stat_weight(uniquelevelindex);
    const auto nnlevel = calculate_levelpop(nonemptymgi, element, ion, level);

    // de-excitation
    const auto alltrans_startdown = get_alltrans_startdown(uniquelevelindex);
    const int ndowntrans = get_ndowntrans(uniquelevelindex);
    const auto ndowntransindices = std::views::iota(0, ndowntrans);
    std::for_each(ndowntransindices.begin(), ndowntransindices.end(), [&](const auto& i) {
      const auto alltransindex = alltrans_startdown + i;
      const int lower = globals::alltrans.targetlevelindex[alltransindex];
      const auto lower_uniquelevelindex = ionuniquelevelindexstart + lower;
      const auto lower_statweight = stat_weight(lower_uniquelevelindex);
      const auto nnlevel_lower = calculate_levelpop(nonemptymgi, element, ion, lower);

      const double epsilon_trans = epsilon_level - epsilon(lower_uniquelevelindex);
      const double R = rad_deexcitation_ratecoeff(epsilon_trans, globals::alltrans.einstein_A[alltransindex],
                                                  statweight, lower_statweight, nnlevel, nnlevel_lower, t_mid) *
                       s_renorm[level];
      const double C =
          col_deexcitation_ratecoeff(T_e, clumpednne, epsilon_trans, statweight, lower_statweight, alltransindex) *
          s_renorm[level];

      const int lower_index = get_nlte_vector_index(element, ion, lower, first_ion_used);
      const auto matrix_index_upper_upper = matrix_index_level_level;
      const auto matrix_index_lower_upper = (lower_index * nlte_dimension) + level_index;

      atomicadd(rate_matrices.rad_bb[matrix_index_upper_upper], -R);
      atomicadd(rate_matrices.rad_bb[matrix_index_lower_upper], R);
      atomicadd(rate_matrices.coll_bb[matrix_index_upper_upper], -C);
      atomicadd(rate_matrices.coll_bb[matrix_index_lower_upper], C);
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
      const auto nnlevel_upper = calculate_levelpop(nonemptymgi, element, ion, upper);

      const double R =
          rad_excitation_ratecoeff(nonemptymgi, upper_statweight, globals::alltrans.einstein_A[alltransindex],
                                   epsilon_trans, nnlevel, nnlevel_upper, statweight, alltransindex, t_mid) *
          s_renorm[level];

      const double C =
          col_excitation_ratecoeff(T_e, clumpednne, upper_statweight, alltransindex, epsilon_trans, statweight) *
          s_renorm[level];

      const double NTC =
          nonthermal::nt_excitation_ratecoeff(nonemptymgi, level, upper, alltransindex) * s_renorm[level];

      const int upper_index = get_nlte_vector_index(element, ion, upper, first_ion_used);
      const auto matrix_index_lower_lower = matrix_index_level_level;
      const auto matrix_index_upper_lower = (upper_index * nlte_dimension) + level_index;

      atomicadd(rate_matrices.rad_bb[matrix_index_lower_lower], -R);
      atomicadd(rate_matrices.rad_bb[matrix_index_upper_lower], R);
      atomicadd(rate_matrices.coll_bb[matrix_index_lower_lower], -C);
      atomicadd(rate_matrices.coll_bb[matrix_index_upper_lower], C);
      atomicadd(rate_matrices.ntcoll_bb[matrix_index_lower_lower], -NTC);
      atomicadd(rate_matrices.ntcoll_bb[matrix_index_upper_lower], NTC);
    });
  });
}

// add photoionisation and thermal collisional ionisation to NLTE matrix
void nltepop_matrix_add_ionisation(const int nonemptymgi, const int element, const int ion,
                                   const std::span<const double> s_renorm,
                                   const std::span<const double> s_renorm_upperion, RateMatrices& rate_matrices,
                                   const int first_ion_used, const int nions_used) {
  assert_always((ion + 1) < (nions_used + first_ion_used));  // can't ionise top ion stage
  const auto T_e = grid::get_Te(nonemptymgi);
  const float clumpednne = grid::get_clumpfactor(nonemptymgi) * grid::get_nne(nonemptymgi);
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

      // photoionisation and collisional ionisation
      const double R_ionisation = get_corrphotoioncoeff(element, ion, level, phixstargetindex, nonemptymgi, false);
      const double C_ionisation =
          col_ionisation_ratecoeff(T_e, clumpednne, element, ion, level, phixstargetindex, epsilon_trans);

      const auto matrix_index_upper_lower = (upper_index * nlte_dimension) + lower_index;

      atomicadd(rate_matrices.rad_bf[matrix_index_lower_lower], -R_ionisation * s_renorm[level]);
      atomicadd(rate_matrices.rad_bf[matrix_index_upper_lower], R_ionisation * s_renorm[level]);
      atomicadd(rate_matrices.coll_bf[matrix_index_lower_lower], -C_ionisation * s_renorm[level]);
      atomicadd(rate_matrices.coll_bf[matrix_index_upper_lower], C_ionisation * s_renorm[level]);

      assert_always((R_ionisation >= 0) && (C_ionisation >= 0));

      // recombination
      if (upper <= maxrecombininglevel) {
        const double R_recomb = rad_recombination_ratecoeff(T_e, clumpednne, element, ion + 1, upper, level);
        const double C_recomb =
            col_recombination_ratecoeff(T_e, clumpednne, element, ion + 1, upper, level, epsilon_trans);

        const auto matrix_index_upper_upper = (upper_index * nlte_dimension) + upper_index;
        const auto matrix_index_lower_upper = (lower_index * nlte_dimension) + upper_index;

        atomicadd(rate_matrices.rad_bf[matrix_index_upper_upper], -R_recomb * s_renorm_upperion[upper]);
        atomicadd(rate_matrices.rad_bf[matrix_index_lower_upper], R_recomb * s_renorm_upperion[upper]);
        atomicadd(rate_matrices.coll_bf[matrix_index_upper_upper], -C_recomb * s_renorm_upperion[upper]);
        atomicadd(rate_matrices.coll_bf[matrix_index_lower_upper], C_recomb * s_renorm_upperion[upper]);

        assert_always((R_recomb >= 0) && (C_recomb >= 0));
      }
    }
  });
}

// add collisional ionisation by non-thermal electrons to NLTE matrix
void nltepop_matrix_add_nt_ionisation(const int nonemptymgi, const int element, const int ion,
                                      const std::span<const double> s_renorm, RateMatrices& rate_matrices,
                                      const int first_ion_used, const int nions_used) {
  if (!NT_ON) {
    return;
  }
  const int max_ion_used = first_ion_used + nions_used - 1;
  assert_always(ion < max_ion_used);  // can't ionise top ion stage
  const double Y_nt = nonthermal::nt_ionisation_ratecoeff(nonemptymgi, element, ion);

  const int nlevels = get_nlevels(element, ion);
  const auto nlte_dimension = rate_matrices.used_nlte_dimension;

  for (int upperion = ion + 1; upperion <= nonthermal::nt_ionisation_maxupperion(element, ion); upperion++) {
    const double Y_nt_thisupperion =
        Y_nt * nonthermal::nt_ionisation_upperion_probability(nonemptymgi, element, ion, upperion, false);

    if (Y_nt_thisupperion > 0.) {
      // ensure if upperion here is past the upper ion used in the NLTE matrix the rates
      // to this ion are instead added to the ground state of the uppermost ion used in the NLTE matrix
      const int upperion_clamped = std::min(upperion, max_ion_used);
      const int upper_groundstate_index = get_nlte_vector_index(element, upperion_clamped, 0, first_ion_used);
      for (int level = 0; level < nlevels; level++) {
        const int lower_index = get_nlte_vector_index(element, ion, level, first_ion_used);

        atomicadd(rate_matrices.ntcoll_bf[(lower_index * nlte_dimension) + lower_index],
                  -Y_nt_thisupperion * s_renorm[level]);
        atomicadd(rate_matrices.ntcoll_bf[(upper_groundstate_index * nlte_dimension) + lower_index],
                  Y_nt_thisupperion * s_renorm[level]);
      }
    }
  }
}

// Add autoionisation and inverse (i.e. collisional capture part of di-el)
void nltepop_matrix_add_autoionisation(const int nonemptymgi, const int element, const int ion,
                                       const std::span<const std::vector<double>> s_renorm_allions,
                                       RateMatrices& rate_matrices, const int first_ion_used, const int nions_used) {
  if (get_nlevels_autoion(element, ion) == 0) {
    return;  // no autoionising levels for this ion
  }
  const auto& s_renorm = s_renorm_allions[ion];

  const auto nlte_dimension = rate_matrices.used_nlte_dimension;
  const int max_ion_used = first_ion_used + nions_used - 1;
  assert_always(ion < max_ion_used);  // can't ionise top ion stage
  const auto T_e = grid::get_Te(nonemptymgi);
  const float clumpednne = grid::get_clumpfactor(nonemptymgi) * grid::get_nne(nonemptymgi);
  const int nlevels = get_nlevels(element, ion);
  for (int level = 0; level < nlevels; level++) {
    const int level_index = get_nlte_vector_index(element, ion, level, first_ion_used);
    const double epsilon_level = epsilon(element, ion, level);
    const double statweight = stat_weight(element, ion, level);

    const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
    const int nautoiondowntrans = get_nautoiondowntrans(uniquelevelindex);
    for (int i = 0; i < nautoiondowntrans; i++) {
      // autoionisation (which is a de-excitation propcess)
      const auto& autoiontransition = globals::allautoion[globals::alllevels.allautoion_start[uniquelevelindex] + i];
      const double A_a = autoiontransition.autoion_A;
      const int target_ion = autoiontransition.upperionindex;
      const int target_level = autoiontransition.upperlevelindex;
      assert_always(target_ion <= max_ion_used);  // the target ion must be included in the NLTE matrix

      const double epsilon_trans = epsilon_level - epsilon(element, target_ion, target_level);

      double R = A_a * s_renorm[level];

      const int upper_index = level_index;
      const int lower_index = get_nlte_vector_index(element, target_ion, target_level, first_ion_used);

      rate_matrices.autoion[(upper_index * nlte_dimension) + upper_index] -= R;
      rate_matrices.autoion[(lower_index * nlte_dimension) + upper_index] += R;
      if ((R < 0)) {
        printlnlog("  WARNING: Negative autoionisation rate from ionstage {} level {} to level {}",
                   get_ionstage(element, ion), level, target_level);
      }

      // capture (which is an excitation process, and the first part of di-electronic recomb).
      // the rate is out of the target ion's level, so if that level is folded into the target ion's
      // superlevel, weight the rate by its Boltzmann fraction of the superlevel population
      R = clumpednne * A_a * statweight / stat_weight(element, target_ion, target_level) * SAHACONST * pow(T_e, -1.5) *
          exp(-epsilon_trans / KB / T_e) * s_renorm_allions[target_ion][target_level];

      rate_matrices.autoion[(lower_index * nlte_dimension) + lower_index] -= R;
      rate_matrices.autoion[(upper_index * nlte_dimension) + lower_index] += R;
      if ((R < 0)) {
        printlnlog("  WARNING: Negative autoionisation rate from ionstage {} level {} to level {}",
                   get_ionstage(element, ion), level, target_level);
      }
    }
  }
}

// Calculate the l2 norm of row or column `i` of `n*n` matrix
auto matrix_row_or_col_norm(std::span<const double> matrix, const int i, const bool row, const size_t n) -> double {
  double result = 0.;

  for (size_t j = 0; j < n; j++) {
    const double val = row ? matrix[(i * n) + j] : matrix[(j * n) + i];
    result += pow2(val);
  }

  return std::sqrt(result);
}

// Multiply row or column `i` by `val`
void matrix_scale_row_or_col(std::span<double> matrix, const int i, const bool row, const size_t n, const double val) {
  for (size_t j = 0; j < n; j++) {
    const size_t idx = row ? (i * n) + j : (j * n) + i;
    matrix[idx] *= val;
  }
}

void nltepop_matrix_normalise(std::span<double> rate_matrix, std::span<double> balance_vector,
                              std::span<double> pop_normfactors) {
  const auto nlte_dimension = std::ssize(pop_normfactors);
  assert_always(std::ssize(rate_matrix) == (nlte_dimension * nlte_dimension));

  for (int iter = 0; iter < 10; iter++) {
    bool changed = false;

    for (auto i = 0; i < nlte_dimension; i++) {
      const double row_norm = matrix_row_or_col_norm(rate_matrix, i, true, nlte_dimension);
      const double col_norm = matrix_row_or_col_norm(rate_matrix, i, false, nlte_dimension);
      if (row_norm == 0 || col_norm == 0) {
        continue;
      }

      const double f = std::sqrt(col_norm / row_norm);
      if (std::abs(f - 1.) > 1e-3) {
        matrix_scale_row_or_col(rate_matrix, i, true, nlte_dimension, f);
        matrix_scale_row_or_col(rate_matrix, i, false, nlte_dimension, 1. / f);
        pop_normfactors[i] /= f;
        changed = true;
      }
    }

    if (!changed) {
      break;
    }
  }

  for (int i = 0; i < nlte_dimension; i++) {
    balance_vector[i] /= pop_normfactors[i];
  }
}

void set_nlte_levelpop_over_rho(const int nonemptymgi, const int element, const int ion, const int level,
                                const double value) {
  assert_testmodeonly(level > 0);  // ground state is stored separately
  assert_testmodeonly(level <= get_nlevels_excited_nlte(element, ion));
  nltepops_allcells[(nonemptymgi * globals::total_nlte_levels) + get_allnltelevelsindexstart(element, ion) + level -
                    1] = value;
}

void set_nlte_superlevelpop_over_rho_over_slpartfunc(const int nonemptymgi, const int element, const int ion,
                                                     const double value) {
  assert_testmodeonly(ion_has_superlevel(element, ion));
  const int sl_nlte_index = get_allnltelevelsindexstart(element, ion) + get_nlevels_excited_nlte(element, ion);
  nltepops_allcells[(nonemptymgi * globals::total_nlte_levels) + sl_nlte_index] = value;
}

void set_element_pops_lte(const int nonemptymgi, const int element) {
  // set NLTE level pops as invalid so that Boltzmann pops will be used instead
  nltepop_reset_element(nonemptymgi, element);
  calculate_cellpartfuncts(nonemptymgi, element);
  // set_groundlevelpops uses uppermost_ion. Previously, this was set based on the NLTE phi factors. Therefore we need
  // to call find_uppermost_ion with force_saha = true so the uppermost ion used in set_groundlevelpops is changed to
  // the one based on the correct LTE phi factors instead
  printlnlog("  WARNING: NLTEPOP setting element Z={} level pops to LTE", get_atomicnumber(element));
  const double nne_hi = grid::get_rho(nonemptymgi) / MH;
  constexpr bool force_saha = true;
  grid::set_elements_uppermost_ion(nonemptymgi, element, find_uppermost_ion(nonemptymgi, element, nne_hi, force_saha));
  set_groundlevelpops(nonemptymgi, element, grid::get_nne(nonemptymgi), force_saha);
}

[[nodiscard]] auto solution_pops_are_valid(const int nonemptymgi, const int element, std::span<double> popvec,
                                           const int first_ion_used, const int nions_used) -> bool {
  const size_t nlte_dimension = popvec.size();
  const auto superlevel_partfuncs = get_element_superlevelpartfuncs(nonemptymgi, element);
  for (auto index = 0ZU; index < nlte_dimension; index++) {
    const auto [ion, level] = get_ion_level_of_nlte_vector_index(index, element, first_ion_used, nions_used);
    const auto ionstage = get_ionstage(element, ion);
    const auto population = popvec[index];

    if constexpr (!STRICT_POPULATION_CHECKING) {
      if (population < 0.0 || !std::isfinite(population)) {
        double ltepop = calculate_levelpop_boltzmann(nonemptymgi, element, ion, level);

        if (level_isinsuperlevel(element, ion, level)) {
          for (int dummylevel = level + 1; dummylevel < get_nlevels(element, ion); dummylevel++) {
            if (level_isinsuperlevel(element, ion, dummylevel)) {
              ltepop += calculate_levelpop_boltzmann(nonemptymgi, element, ion, dummylevel);
            }
          }
        }

        if (!std::isfinite(ltepop) || ltepop <= 0.) {
          printlnlog(
              "  WARNING: Boltzmann population for Z={} ionstage {} level {} is non-finite or non-positive ({}). "
              "Setting to MINPOP",
              get_atomicnumber(element), get_ionstage(element, ion), level, ltepop);
          ltepop = MINPOP;
        }

        printlnlog(
            "  WARNING: NLTE solver gave negative/non-finite population to index {} (Z={} ionstage {} level {}), pop = "
            "{:g}. Replacing with LTE pop of {:g}",
            index, get_atomicnumber(element), ionstage, level, population, ltepop);
        popvec[index] = ltepop;
      }
    } else {
      // if groundpop is below MINPOP then the NLTE solution fails
      const size_t index_ion_ground = get_nlte_vector_index(element, ion, 0, first_ion_used);
      if (index == index_ion_ground && population < MINPOP) {
        printlnlog(
            "  WARNING: NLTE solver gave ground pop less than MINPOP for index {} (Z={} ionstage {} level {}), pop = "
            "{:g}. Returning nltepop_matrix_solve fail",
            index, get_atomicnumber(element), ionstage, level, population);
        return false;
      }

      if (population < 0.0 || !std::isfinite(population)) {
        printlnlog(
            "  WARNING: NLTE solver gave negative/non-finite population for index {} (Z={} ionstage {} level {}), pop "
            "= {:g}",
            index, get_atomicnumber(element), ionstage, level, population);
        if (population < -MINPOP) {
          printlnlog(
              "  WARNING: negative pop = {:g} less than -MINPOP (-{:g}) unlikely a rounding error to zero so returning "
              "nltepop_matrix_solve fail",
              population, MINPOP);
          return false;
        }
        printlnlog(
            "  WARNING: negative pop = {:g} greater than -MINPOP (-{:g}) likely a rounding error to zero so continue "
            "with NLTE pops but set this level to MINPOP",
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
            printlog(
                " WARNING: pop inversion greater than factor {:g}: (g_pop {:g})/(e_pop {:g}) = {:g} is less "
                "than (g_sw {:g})/(e_sw {:g}) = {:g} for index {} Z={} ionstage {} level {} (factor {:g} inversion) - ",
                STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING, ground_pop, population,
                ground_pop / population, statweight_ground, statweight, statweight_ground / statweight, index,
                get_atomicnumber(element), ionstage, level, inversion_factor);

            if (inversion_factor > STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL) {
              printlnlog("large pop inversion - matrix solve failed");
              return false;
            }
            printlnlog("relatively small pop inversion so continue with NLTE solution");
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
            printlog(
                " WARNING: superlevel pop inversion greater than factor {:g}: (g_pop {:g})/(SL_first_level_pop "
                "{:g}) = {:g} is less than (g_sw {:g})/(SL_first_level_sw {:g}) = {:g} for index {} Z={} ionstage {} "
                "level {} (factor {:g} inversion) - ",
                STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING, ground_pop, sublevel_pop,
                ground_pop / sublevel_pop, statweight_ground, statweight, statweight_ground / statweight, index,
                get_atomicnumber(element), ionstage, level, inversion_factor);

            if (inversion_factor > STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL) {
              printlnlog("large pop inversion for superlevel - matrix solve failed");
              return false;
            }
            printlnlog("relatively small pop inversion for superlevel so continue with NLTE solution");
          }
        }
      }
    }
  }
  return true;
}

// solve rate_matrix * x = balance_vector, so that popvec[i] = x[i] * pop_normfactors[i]
// return true if the solution is successful, or false if the matrix is singular, contains negative or inverted
// populations.
[[nodiscard]] auto nltepop_matrix_solve(const int element, const int nonemptymgi,
                                        const std::span<const double> rate_matrix,
                                        std::span<const double> balance_vector, std::span<double> popvec,
                                        std::span<const double> pop_normfactors, const int max_nlte_dimension,
                                        const int first_ion_used, const int nions_used) -> bool {
  const size_t nlte_dimension = balance_vector.size();
  assert_always(pop_normfactors.size() == nlte_dimension);
  assert_always(rate_matrix.size() == (nlte_dimension * nlte_dimension));
  assert_always(std::cmp_greater_equal(max_nlte_dimension, nlte_dimension));

  // solution vector for the matrix equation
  THREADLOCALONHOST std::vector<double> vec_x;
  vec_x.reserve(max_nlte_dimension);
  vec_x.resize(nlte_dimension);

#ifdef EIGEN_OFF

  THREADLOCALONHOST std::vector<double> rate_matrix_LU_decomp;
  rate_matrix_LU_decomp.reserve(static_cast<ptrdiff_t>(max_nlte_dimension) * max_nlte_dimension);
  rate_matrix_LU_decomp.resize(rate_matrix.size());

  // make a copy of the rate matrix for the LU decomp call as gsl_linalg_LU_decomp modifies the input matrix
  std::ranges::copy(rate_matrix, rate_matrix_LU_decomp.begin());
  auto gsl_rate_matrix_LU_decomp =
      gsl_matrix_view_array(rate_matrix_LU_decomp.data(), nlte_dimension, nlte_dimension).matrix;
  auto gsl_x = gsl_vector_view_array(vec_x.data(), nlte_dimension).vector;
  auto gsl_rate_matrix = gsl_matrix_const_view_array(rate_matrix.data(), nlte_dimension, nlte_dimension).matrix;

  THREADLOCALONHOST std::vector<size_t> vec_permutation;
  vec_permutation.reserve(max_nlte_dimension);
  vec_permutation.resize(nlte_dimension);
  gsl_permutation_struct p{.size = nlte_dimension, .data = vec_permutation.data()};
  gsl_permutation_init(&p);

  auto gsl_balance_vector = gsl_vector_const_view_array(balance_vector.data(), nlte_dimension).vector;

  int s = 0;  // sign of the transformation
  assert_always(gsl_linalg_LU_decomp(&gsl_rate_matrix_LU_decomp, &p, &s) == GSL_SUCCESS);

  bool lumatrix_is_singular = false;
  for (auto i = 0ZU; i < nlte_dimension; i++) {
    // diagonal elements of LU matrix must be non-zero and finite
    const double& matrixelement = rate_matrix_LU_decomp[(i * nlte_dimension) + i];
    if (matrixelement == 0. || !std::isfinite(matrixelement)) {
      const auto [ion, level] = get_ion_level_of_nlte_vector_index(i, element, first_ion_used, nions_used);
      if (!lumatrix_is_singular) {
        printlog("ERROR: NLTE disconnected: Z={}", get_atomicnumber(element));
      }
      if (is_nlte(element, ion, level)) {
        printlog(" ionstage {} level {}", get_ionstage(element, ion), level);
      } else {
        printlog(" ionstage {} superlevel ", get_ionstage(element, ion));
      }

      lumatrix_is_singular = true;
    }
  }
  if (lumatrix_is_singular) {
    printlnlog("");
    printlnlog("  ERROR: NLTE matrix is singular for element Z={}!", get_atomicnumber(element));
    return false;
  }

  // solve matrix equation: rate_matrix * x = balance_vector for x (population vector)
  gsl_linalg_LU_solve(&gsl_rate_matrix_LU_decomp, &p, &gsl_balance_vector, &gsl_x);

#else

  auto eigen_vec_x = Eigen::Map<Eigen::VectorX<double>>(vec_x.data(), nlte_dimension);
  const auto eigen_balance_vector = Eigen::Map<const Eigen::VectorX<double>>(balance_vector.data(), nlte_dimension);
  assert_always(eigen_balance_vector.allFinite());

  const auto eigen_rate_matrix =
      Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
          rate_matrix.data(), nlte_dimension, nlte_dimension);
  assert_always(eigen_rate_matrix.allFinite());

  THREADLOCALONHOST Eigen::PartialPivLU<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
      eigen_rate_matrix_lu;

  eigen_rate_matrix_lu.compute(eigen_rate_matrix);
  bool lumatrix_is_singular = false;
  for (auto i = 0ZU; i < nlte_dimension; i++) {
    // diagonal elements of LU matrix must be non-zero and finite
    const double& matrixelement = eigen_rate_matrix_lu.matrixLU().diagonal()[i];
    if (matrixelement == 0. || !std::isfinite(matrixelement)) {
      const auto [ion, level] = get_ion_level_of_nlte_vector_index(i, element, first_ion_used, nions_used);
      if (!lumatrix_is_singular) {
        printlog("ERROR: NLTE disconnected: Z={}", get_atomicnumber(element));
      }
      if (is_nlte(element, ion, level)) {
        printlog(" ionstage {} level {}", get_ionstage(element, ion), level);
      } else {
        printlog(" ionstage {} superlevel ", get_ionstage(element, ion));
      }

      lumatrix_is_singular = true;
    }
  }
  if (lumatrix_is_singular) {
    printlnlog("");
    printlnlog("ERROR: NLTE matrix is singular for element Z={}!", get_atomicnumber(element));
    return false;
  }

  eigen_vec_x = eigen_rate_matrix_lu.solve(eigen_balance_vector);
  assert_always(eigen_vec_x.allFinite());

#endif

  double error_best = -1.;

  // population solution vector with lowest error
  THREADLOCALONHOST std::vector<double> vec_x_best;
  vec_x_best.reserve(max_nlte_dimension);
  vec_x_best.resize(nlte_dimension);

  THREADLOCALONHOST std::vector<double> vec_residual;
  vec_residual.reserve(max_nlte_dimension);
  vec_residual.resize(nlte_dimension);

#ifdef EIGEN_OFF
  auto gsl_vec_residual = gsl_vector_view_array(vec_residual.data(), nlte_dimension).vector;
#else
  auto eigen_vec_residual = Eigen::Map<Eigen::VectorX<double>>(vec_residual.data(), nlte_dimension);
#endif

  int iteration = 0;
  for (iteration = 0; iteration < 10; iteration++) {
#ifdef EIGEN_OFF
    if (iteration > 0) {
      // use gsl_vec_residual as the temporary work vector
      gsl_linalg_LU_refine(&gsl_rate_matrix, &gsl_rate_matrix_LU_decomp, &p, &gsl_balance_vector, &gsl_x,
                           &gsl_vec_residual);
    }

    std::ranges::copy(balance_vector, vec_residual.begin());
    gsl_blas_dgemv(CblasNoTrans, 1.0, &gsl_rate_matrix, &gsl_x, -1.0,
                   &gsl_vec_residual);  // calculate Ax - b = residual

    const double error = fabs(gsl_vector_get(
        &gsl_vec_residual, gsl_blas_idamax(&gsl_vec_residual)));  // value of the largest absolute residual
#else
    if (iteration > 0) {
      eigen_vec_x += eigen_rate_matrix_lu.solve(eigen_vec_residual);
    }
    eigen_vec_residual = eigen_balance_vector - eigen_rate_matrix * eigen_vec_x;
    const double error = eigen_vec_residual.cwiseAbs().maxCoeff();
#endif

    if ((std::isfinite(error) && (error < error_best)) || error_best < 0.) {
      std::ranges::copy(vec_x, vec_x_best.begin());
      error_best = error;
    }
    if (error < 1e-40) {
      break;
    }
  }

  if (!(error_best >= 0.)) {
    printlnlog("  ERROR: NLTE solver failed to find a solution with finite error for element Z={}",
               get_atomicnumber(element));
    return false;
  }
  std::ranges::copy(vec_x_best, vec_x.begin());

  if (error_best > 1e-8) {
    printlnlog(
        "  NLTE solver matrix LU_refine: After {} iterations, best solution vector has a max residual of {:g} "
        "(WARNING!)",
        iteration, error_best);
  }

  // get the unnormalised populations from the x solution vector and the normalisation factors
  for (auto i = 0ZU; i < nlte_dimension; i++) {
    popvec[i] = vec_x[i] * pop_normfactors[i];
  }

  return solution_pops_are_valid(nonemptymgi, element, popvec, first_ion_used, nions_used);
}

auto can_remove_ion(const int element, const int ion, const int first_ion_used, const int nions_used,
                    const double nnelement, std::span<const double> popvec) -> bool {
  if (nions_used <= 1) {
    // can't remove an ion if there is only one ion stage used in the NLTE matrix
    printlnlog("  WARNING: can't remove ion stage {} from NLTE matrix for element Z={} as only one ion stage used",
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
    printlnlog(
        "  WARNING: {} ion ground state population too large to remove ion (ground_pop / nnelement) = ({}/{}) = {} "
        "> {}",
        ionname, ground_pop, nnelement, ground_pop / nnelement,
        NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION);

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
      printlnlog(
          "  WARNING: {} ion excited state (level {}) population too large to remove ion (nlte_excited_pop_{}_ion "
          "/ nnelement) = ({}/{}) > ({}))",
          ionname, ionname, level, levelpop, nnelement, NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION);
      return false;
    }
  }

  double superlevel_pop = 0.0;
  if (ion_has_superlevel(element, ion)) {
    // superlevel populations relative to the total element population must not exceed the limit
    const int index_superlevel = get_nlte_vector_index(element, ion, nlevels_nlte_excited + 1, first_ion_used);
    superlevel_pop = popvec[index_superlevel];
    if (superlevel_pop / nnelement > NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION) {
      printlnlog(
          "  WARNING: {} ion superlevel population too large to remove ion (superlevel_pop_{}_ion / nnelement "
          "({:g}/{:g}) > ({:g}) NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION)",
          ionname, ionname, superlevel_pop, nnelement, NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION);
      return false;
    }
  }

  printlnlog(
      "  Passed checks: removing {} ion from NLTE matrix and attempting re-solve (ionstage {}, ground pop {:g}, "
      "fabs(excited pops) sum {:g}, superlevel pop {:g}, nnelement {:g}) matrix",
      ionname, get_ionstage(element, ion), ground_pop, nlte_excited_pop_sum, superlevel_pop, nnelement);

  return true;
}

// Assemble and solve the NLTE rate matrix, retrying with a reduced range of ion stages when the
// solve fails and NLTE_LIMIT_ION_STAGES_AFTER_FAILURE allows an edge ion to be removed.
// Returns whether a solution was found and the range of ions included in it.
auto nltepop_solve_matrix_with_ion_reduction(const int element, const int nonemptymgi, const double t_mid,
                                             const double nnelement,
                                             const std::vector<std::vector<double>>& s_renorm_allions,
                                             RateMatrices& rate_matrices, std::vector<double>& popvec)
    -> std::tuple<bool, int, int> {
  const int atomic_number = get_atomicnumber(element);
  const int nions = get_nions(element);
  const auto max_nlte_dimension = get_max_nlte_dimension();
  int nions_used = nions;
  int first_ion_used = 0;
  bool matrix_solve_success = false;

  bool matrix_solve_required = true;
  while (matrix_solve_required) {
    const auto nlte_dimension = get_element_nlte_dimension(element, first_ion_used, nions_used);
    rate_matrices.set_used_dimension(nlte_dimension);
    popvec.resize(nlte_dimension);
    // popvec is a reused (thread-local) buffer. Clear it so that if the solve fails before writing a
    // solution (e.g. singular matrix), can_remove_ion() reads zeros for this element rather than stale
    // populations left over from a previous element or cell.
    std::ranges::fill(popvec, 0.);

    const int max_ion_used = first_ion_used + nions_used - 1;

    const auto ions = std::views::iota(first_ion_used, max_ion_used + 1);
    std::for_each(ions.begin(), ions.end(), [&](const auto ion) {
      const auto& s_renorm = s_renorm_allions[ion];

      nltepop_matrix_add_boundbound(nonemptymgi, element, ion, t_mid, s_renorm, rate_matrices, first_ion_used);

      if (ion < max_ion_used) {
        nltepop_matrix_add_ionisation(nonemptymgi, element, ion, s_renorm, s_renorm_allions[ion + 1], rate_matrices,
                                      first_ion_used, nions_used);
        nltepop_matrix_add_nt_ionisation(nonemptymgi, element, ion, s_renorm, rate_matrices, first_ion_used,
                                         nions_used);
        nltepop_matrix_add_autoionisation(nonemptymgi, element, ion, s_renorm_allions, rate_matrices, first_ion_used,
                                          nions_used);
      }
    });

    // replace the zeroth row of the matrix and balance vector with the normalisation
    // constraint (sum of levelpops = total element population)

    const auto rate_matrix = rate_matrices.get_summed_rate_matrix();
    std::ranges::fill(std::span{rate_matrix}.first(nlte_dimension), 1.0);

    THREADLOCALONHOST std::vector<double> balance_vector;
    balance_vector.reserve(max_nlte_dimension);
    balance_vector.resize(nlte_dimension);
    std::ranges::fill(balance_vector, 0.0);  // statistical equilibrium means balance vector is zero
    // except for zeroth element used to normalise the total element population
    balance_vector[0] = nnelement;

    if (FORCE_SAHA_ION_BALANCE(atomic_number)) {
      const auto ionfractions = calculate_ionfractions(element, nonemptymgi, grid::get_nne(nonemptymgi), true);
      // ssize() to avoid unsigned wraparound to a huge positive value when the vector is empty
      const int uppermost_ion = static_cast<int>(std::ssize(ionfractions)) - 1;
      for (int ion = first_ion_used + 1; ion <= std::min(uppermost_ion, max_ion_used); ion++) {
        // replace matrix row for ion's ground state with:
        // sum of this ion's level populations is equal to the ion population
        const double nnion = nnelement * ionfractions[ion];
        const int index_ion_ground = get_nlte_vector_index(element, ion, 0, first_ion_used);
        const int index_ion_toplevel =
            get_nlte_vector_index(element, ion, get_nlevels(element, ion) - 1, first_ion_used);
        for (int index = 0; index < nlte_dimension; index++) {
          rate_matrix[(index_ion_ground * nlte_dimension) + index] =
              (index >= index_ion_ground && index <= index_ion_toplevel) ? 1.0 : 0.0;
        }

        balance_vector[index_ion_ground] = nnion;
      }
    }

    // calculate the normalisation factors and apply them to the matrix
    // columns and balance vector elements
    THREADLOCALONHOST std::vector<double> pop_normfactors;
    pop_normfactors.reserve(max_nlte_dimension);
    pop_normfactors.resize(nlte_dimension);
    std::ranges::fill(pop_normfactors, 1.0);
    nltepop_matrix_normalise(rate_matrix, balance_vector, pop_normfactors);

    matrix_solve_success = nltepop_matrix_solve(element, nonemptymgi, rate_matrix, balance_vector, popvec,
                                                pop_normfactors, max_nlte_dimension, first_ion_used, nions_used);

    matrix_solve_required = false;  // will be set to true if we need to retry with a different ion range
    if (matrix_solve_success) {
      if (nions_used < nions) {
        printlnlog(
            "successfully solved NLTE matrix when reducing ions used for element to Z={} ionstage={} to ionstage={}",
            atomic_number, get_ionstage(element, first_ion_used), get_ionstage(element, max_ion_used));
      }
    } else if (NLTE_LIMIT_ION_STAGES_AFTER_FAILURE) {
      printlnlog("  WARNING: NLTE matrix solution failed for element Z={} using ionstage {} to {}", atomic_number,
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
        printlnlog(
            "  WARNING: can't remove top or bottom ion stage from NLTE equations, therefore unable to find an NLTE "
            "solution for this element");
      }
    }
  }

  return {matrix_solve_success, first_ion_used, nions_used};
}

// Store the solved populations (ground level, NLTE excited levels, and superlevel) for every
// ion of the element, then check that the total element population matches the element
// abundance and print the periodic rates summaries.
void nltepop_apply_solution(const int element, const int nonemptymgi, const int timestep, const int nlte_iter,
                            const double nnelement, const std::vector<double>& superlevel_partfuncs,
                            const std::vector<double>& popvec, const RateMatrices& rate_matrices,
                            const int first_ion_used, const int nions_used) {
  const int atomic_number = get_atomicnumber(element);
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const int nions = get_nions(element);

  // check calculated NLTE populations are valid
  for (const auto pop : popvec) {
    assert_always(std::isfinite(pop));
    assert_always(pop >= 0.);
  }

  // set the ground level, excited level and possible superlevel populations for this element
  for (int ion = 0; ion < nions; ion++) {
    const int nlevels_excited_nlte = get_nlevels_excited_nlte(element, ion);

    if (ion < first_ion_used || ion >= (first_ion_used + nions_used)) {
      printlnlog("  WARNING: Z={} ionstage {} removed from NLTE rate matrix. Setting all levelpops for ion to zero ",
                 get_atomicnumber(element), get_ionstage(element, ion));

      set_groundlevelpop(nonemptymgi, element, ion, 0.);

      for (int level = 1; level <= nlevels_excited_nlte; level++) {
        set_nlte_levelpop_over_rho(nonemptymgi, element, ion, level, 0.);
      }

      if (ion_has_superlevel(element, ion)) {
        set_nlte_superlevelpop_over_rho_over_slpartfunc(nonemptymgi, element, ion, 0.);
      }
    } else {
      // only valid for ions inside the solved range [first_ion_used, first_ion_used + nions_used)
      const int index_gs = get_nlte_vector_index(element, ion, 0, first_ion_used);
      set_groundlevelpop(nonemptymgi, element, ion, static_cast<float>(popvec[index_gs]));

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
    printlnlog(
        "  WARNING: timestep {} nlteiter {} Z={} element population is: {:g} (from abundance) and {:g} (from matrix "
        "solution), error: {:.2f}%. Forcing element pops to LTE.",
        timestep, nlte_iter, atomic_number, nnelement, elem_pop_matrix, elem_pop_error_percent);
    set_element_pops_lte(nonemptymgi, element);
  }

  // output NLTE stats every nth timestep for the first NLTE iteration only
  if ((timestep % 5 == 0) && (nlte_iter == 0)) {
    print_element_rates_summary(element, modelgridindex, timestep, nlte_iter, popvec, rate_matrices, first_ion_used,
                                nions_used);
  }

  // detailed levels stats (very verbose, only for debugging)
  const bool print_detailed_level_stats = false;
  if (print_detailed_level_stats) {
    const int ionstage = 2;
    const int ion = ionstage - get_ionstage(element, 0);

    const int nlevels = get_nlevels_excited_nlte(element, ion) + (ion_has_superlevel(element, ion) ? 1 : 0);
    for (int level = 0; level <= nlevels; level++) {
      print_level_rates(nonemptymgi, timestep, element, ion, level, popvec, rate_matrices, first_ion_used, nions_used);
    }
  }
}

}  // anonymous namespace

void solve_nlte_pops_element(const int element, const int nonemptymgi, const int timestep, const int nlte_iter) {
  // solve the statistical balance equations to find NLTE level populations for all ions of an element
  // (ionisation balance follows from this too). Failure modes can lead to quasi-LTE populations being set instead.

  const int atomic_number = get_atomicnumber(element);
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);

  if (grid::get_elem_massfrac(nonemptymgi, element) <= 0.) {
    // abundance of this element is zero, so do not store any NLTE populations
    printlnlog(
        "Not solving for NLTE populations in cell {} at timestep {} NLTE iteration {} for element Z={} due to zero "
        "abundance",
        modelgridindex, timestep, nlte_iter, atomic_number);

    nltepop_reset_element(nonemptymgi, element);
    return;
  }

  const double cell_Te = grid::get_Te(nonemptymgi);

  if (cell_Te == MINTEMP) {
    printlnlog(
        "Not solving for NLTE populations in cell {} at timestep {} NLTE iteration {} for element Z={} due to low "
        "temperature Te=MINTEMP={:g} K",
        modelgridindex, timestep, nlte_iter, atomic_number, cell_Te);
    set_element_pops_lte(nonemptymgi, element);
    return;
  }

  const auto sys_time_start_nltesolver = std::chrono::steady_clock::now();

  const double t_mid = globals::timesteps[timestep].mid;
  const int nions = get_nions(element);
  const double nnelement = grid::get_elem_numberdens(nonemptymgi, element);

  printlnlog(
      "Solving for NLTE populations in cell {} at timestep {} NLTE iteration {} for element Z={} (mass fraction "
      "{:.2e}, nnelement {:.2e} cm^-3)",
      modelgridindex, timestep, nlte_iter, atomic_number, grid::get_elem_massfrac(nonemptymgi, element), nnelement);

  const auto superlevel_partfuncs = get_element_superlevelpartfuncs(nonemptymgi, element);

  const auto max_nlte_dimension = get_max_nlte_dimension();

  // will hold the un-normalised population densities [cm^-3]
  THREADLOCALONHOST std::vector<double> popvec;
  popvec.reserve(max_nlte_dimension);

  // superlevel population renormalisation factors for every ion of the element:
  // superlevel members get their Boltzmann fraction of the superlevel population, while levels
  // with their own matrix entries (ground, NLTE excited, and autoionising levels) get a weight of 1.
  THREADLOCALONHOST std::vector<std::vector<double>> s_renorm_allions;
  s_renorm_allions.reserve(get_max_nions());
  s_renorm_allions.resize(nions);
  for (int ion = 0; ion < nions; ion++) {
    const int nlevels = get_nlevels(element, ion);
    s_renorm_allions[ion].resize(nlevels);
    std::ranges::fill(s_renorm_allions[ion], 1.0);  // default to 1. for all levels
    const int level_superlevel_start = get_nlevels_excited_nlte(element, ion) + 1;
    const int first_autoion = get_nlevels_nonautoion(element, ion);

    // set the superlevel renormalisation factors for all levels in the superlevel
    for (int level = level_superlevel_start; level < first_autoion; level++) {
      s_renorm_allions[ion][level] = superlevel_boltzmann(nonemptymgi, element, ion, level) / superlevel_partfuncs[ion];
    }
  }

  THREADLOCALONHOST RateMatrices rate_matrices{max_nlte_dimension};

  const auto [matrix_solve_success, first_ion_used, nions_used] = nltepop_solve_matrix_with_ion_reduction(
      element, nonemptymgi, t_mid, nnelement, s_renorm_allions, rate_matrices, popvec);

  if (!matrix_solve_success) {
    printlnlog(
        "WARNING: Can't solve for NLTE populations in cell {} at timestep {} for element Z={} due to singular matrix, "
        "negative pop or large pop inversion and unable to recover solution by reducing ion range. Attempting to use "
        "LTE solution instead",
        modelgridindex, timestep, atomic_number);
    set_element_pops_lte(nonemptymgi, element);
  } else {
    nltepop_apply_solution(element, nonemptymgi, timestep, nlte_iter, nnelement, superlevel_partfuncs, popvec,
                           rate_matrices, first_ion_used, nions_used);
  }

  if (const auto duration_nltesolver =
          std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_nltesolver).count();
      duration_nltesolver > 2.) {
    printlnlog("NLTE population solver call for Z={} took {:.1f} seconds", get_atomicnumber(element),
               duration_nltesolver);
  }
}

// Get a Boltzman factor for a level within the super level (combined Non-LTE level)
DEVICE_FUNC auto superlevel_boltzmann(const int nonemptymgi, const int element, const int ion, const int level)
    -> double {
  // autoionising levels are also allowed here: outside the NLTE solver their populations are
  // attached to the superlevel (see the comment in read_autoion_data())
  assert_testmodeonly(level_isinsuperlevel(element, ion, level) || level_isautoionising(element, ion, level));
  const int level_superlevel_start = get_nlevels_excited_nlte(element, ion) + 1;
  const double T_exc = LTEPOP_EXCITATION_USE_TJ ? grid::get_TJ(nonemptymgi) : grid::get_Te(nonemptymgi);
  const double E_level = epsilon(element, ion, level);
  const double E_superlevel = epsilon(element, ion, level_superlevel_start);

  return stat_weight(element, ion, level) / stat_weight(element, ion, level_superlevel_start) *
         exp(-(E_level - E_superlevel) / KB / T_exc);
}

void nltepop_open_file() {
  nlte_file = fstream_required(std::format("nlte_{:04d}.out", globals::my_rank), std::ios::out | std::ios::trunc);
  std::println(nlte_file, "timestep modelgridindex Z ionstage level n_LTE n_NLTE ion_popfrac");
}

void nltepop_write_to_file(const int nonemptymgi, const int timestep) {
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  if (globals::lte_iteration || grid::thick_allcells[nonemptymgi] == 1) {  // NLTE solver hasn't been run yet
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

        std::print(nlte_file, "{} {} {} {} ", timestep, modelgridindex, get_atomicnumber(element),
                   get_ionstage(element, ion));
        if (level <= nlevels_excited_nlte) {
          std::print(nlte_file, "{} ", level);

          if (level == 0) {
            nnlevelnlte = get_groundlevelpop(nonemptymgi, element, ion);
          } else {
            nnlevelnlte = get_nlte_levelpop_over_rho(nonemptymgi, element, ion, level) * grid::get_rho(nonemptymgi);
          }
        } else {
          // superlevel, so add the populations of all other levels in the superlevel
          const double slpopfactor =
              get_nlte_superlevelpop_over_rho_over_slpartfunc(nonemptymgi, element, ion) * grid::get_rho(nonemptymgi);

          nnlevellte = 0;
          std::print(nlte_file, "-1 ");
          for (int level_sl = nlevels_excited_nlte + 1; level_sl < get_nlevels(element, ion); level_sl++) {
            if (level_isinsuperlevel(element, ion, level_sl)) {
              nnlevellte += calculate_levelpop_boltzmann(nonemptymgi, element, ion, level_sl);
            }
          }

          nnlevelnlte = slpopfactor * superlevel_partfuncs[ion];
        }

        std::println(nlte_file, "{:.5e} {:.5e} {:.5e}", nnlevellte, nnlevelnlte, nnlevelnlte / nnion);
      }
    }
  }

  nlte_file.flush();
}

void nltepop_write_restart_data(FILE* restart_file) {
  printlog("populations, ");

  fprintf(restart_file, "%d\n", 75618527);  // special number marking the beginning of nlte data

  fprintf(restart_file, "%d\n", globals::total_nlte_levels);
  const auto nincludedions = get_includedions();

  for (auto nonemptymgi = 0Z; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    const int modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    fprintf(restart_file, "%d\n", modelgridindex);
    for (int element = 0; element < get_nelements(); element++) {
      for (int ion = 0; ion < get_nions(element); ion++) {
        const int uniqueionindex = get_uniqueionindex(element, ion);
        fprintf(restart_file, "%d %a %a %la\n", ion,
                grid::ion_groundlevelpops_allcells[(nonemptymgi * nincludedions) + uniqueionindex],
                grid::ion_partfuncts_allcells[(nonemptymgi * nincludedions) + uniqueionindex],
                kpkt::ion_cooling_contribs_allcells[(nonemptymgi * nincludedions) + uniqueionindex]);
      }
    }
    for (int nlteindex = 0; nlteindex < globals::total_nlte_levels; nlteindex++) {
      fprintf(restart_file, "%la ", nltepops_allcells[(nonemptymgi * globals::total_nlte_levels) + nlteindex]);
    }
  }
}

void nltepop_read_restart_data(FILE* restart_file) {
  printlnlog("Reading restart data for populations");

  int code_check = 0;
  assert_always(fscanf(restart_file, "%d\n", &code_check) == 1);
  assert_always(code_check == 75618527);

  int total_nlte_levels_in = 0;
  assert_always(fscanf(restart_file, "%d\n", &total_nlte_levels_in) == 1);
  if (total_nlte_levels_in != globals::total_nlte_levels) {
    printlnlog("ERROR: Expected {} NLTE levels but found {} in restart file", globals::total_nlte_levels,
               total_nlte_levels_in);
    std::abort();
  }
  const auto nincludedions = get_includedions();

  for (auto nonemptymgi = 0Z; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    int mgi_in = 0;
    assert_always(fscanf(restart_file, "%d\n", &mgi_in) == 1);
    assert_always(mgi_in == grid::get_mgi_of_nonemptymgi(nonemptymgi));

    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++) {
        int ion_in = 0;
        const int uniqueionindex = get_uniqueionindex(element, ion);
        assert_always(fscanf(restart_file, "%d %a %a %la\n", &ion_in,
                             &grid::ion_groundlevelpops_allcells[(nonemptymgi * nincludedions) + uniqueionindex],
                             &grid::ion_partfuncts_allcells[(nonemptymgi * nincludedions) + uniqueionindex],
                             &kpkt::ion_cooling_contribs_allcells[(nonemptymgi * nincludedions) + uniqueionindex]) ==
                      4);
        assert_always(ion_in == ion);
      }
    }
    for (int nlteindex = 0; nlteindex < globals::total_nlte_levels; nlteindex++) {
      assert_always(fscanf(restart_file, "%la ",
                           &nltepops_allcells[(nonemptymgi * globals::total_nlte_levels) + nlteindex]) == 1);
    }
  }
}

DEVICE_FUNC auto get_nlte_levelpop_over_rho(const int nonemptymgi, const int element, const int ion, const int level)
    -> double {
  assert_testmodeonly(level > 0);  // ground state is stored separately
  assert_testmodeonly(level <= get_nlevels_excited_nlte(element, ion));
  return nltepops_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * globals::total_nlte_levels) +
                           get_allnltelevelsindexstart(element, ion) + level - 1];
}

[[nodiscard]] DEVICE_FUNC auto get_nlte_superlevelpop_over_rho_over_slpartfunc(const int nonemptymgi, const int element,
                                                                               const int ion) -> double {
  assert_testmodeonly(ion_has_superlevel(element, ion));
  const int sl_nlte_index = get_allnltelevelsindexstart(element, ion) + get_nlevels_excited_nlte(element, ion);
  return nltepops_allcells[(nonemptymgi * globals::total_nlte_levels) + sl_nlte_index];
}

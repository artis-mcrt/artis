// The full NLTE level population solver: assembles the statistical equilibrium rate matrix
// (radiative and collisional bound-bound and bound-free rates, plus non-thermal ionisation)
// for all levels of all ions of an element and solves for the level populations.
//
// The NLTE ionisation and population solver is described by Shingles et al. (2020), MNRAS, 492,
// 2029, section 2.3, doi:10.1093/mnras/stz3412.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <fstream>
#include <numeric>
#include <optional>
#include <print>
#include <ranges>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <Eigen/Core>
#include <Eigen/LU>
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

// context so that log messages from the NLTE solver helpers can identify the cell, timestep, and
// iteration of the solve they were emitted from (set by solve_nlte_pops_element)
struct NLTELogContext {
  int modelgridindex = -1;
  int timestep = -1;
  int nlte_iter = -1;
};
THREADLOCALONHOST NLTELogContext nltelog{};

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
  // inverse of get_nlte_vector_index() by linear search. Called once per vector entry when validating or
  // reporting a solution, never from the rate matrix assembly, so the O(nlevels) scan is not on a hot path.
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

// log " ionstage {} level {}" or " ionstage {} superlevel" (no newline) identifying an NLTE vector index in a message
void printlog_nlte_vector_index(const std::ptrdiff_t index, const int element, const int first_ion_used,
                                const int nions_used) {
  const auto [ion, level] = get_ion_level_of_nlte_vector_index(index, element, first_ion_used, nions_used);
  if (is_nlte(element, ion, level)) {
    printlog(" ionstage {} level {}", get_ionstage(element, ion), level);
  } else {
    printlog(" ionstage {} superlevel", get_ionstage(element, ion));
  }
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

    printlnlog("{}{}{:10.2e} {:10.2e} {:10.2e} {:10.2e} {:10.2e} {:10.2e} {:10.2e}", into_level ? " from " : "   to ",
               only_levels_below ? "below " : "above ", rad_bb_total, coll_bb_total, ntcoll_bb_total, rad_bf_total,
               coll_bf_total, ntcoll_bf_total, autoion_total);
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
            modelgridindex, timestep, nlte_iter, grid::Te_allcells[nonemptymgi], grid::get_nne(nonemptymgi),
            atomic_number, ionstage);
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
      timestep, grid::get_mgi_of_nonemptymgi(nonemptymgi), grid::Te_allcells[nonemptymgi], grid::get_nne(nonemptymgi),
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
  const auto T_e = grid::Te_allcells[nonemptymgi];
  const auto clumpednne = grid::get_clumpfactor(nonemptymgi) * grid::get_nne(nonemptymgi);
  const int nlevels = get_nlevels(element, ion);
  const auto ionuniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);
  const auto nlte_dimension = rate_matrices.used_nlte_dimension;

  // the level populations and matrix indices are constants of the fill but are needed for both
  // endpoints of every transition, so look them up once per level instead of twice per transition
  std::vector<double> levelpops(nlevels);
  std::vector<int> levelindices(nlevels);
  for (int level = 0; level < nlevels; level++) {
    levelpops[level] = calculate_levelpop(nonemptymgi, element, ion, level);
    levelindices[level] = get_nlte_vector_index(element, ion, level, first_ion_used);
  }

  const auto levels = std::views::iota(0, nlevels);
  std::for_each(levels.begin(), levels.end(), [&](const auto level) {
    const int level_index = levelindices[level];
    const auto matrix_index_level_level = (level_index * nlte_dimension) + level_index;
    const auto uniquelevelindex = ionuniquelevelindexstart + level;
    const double epsilon_level = epsilon(uniquelevelindex);
    const double statweight = stat_weight(uniquelevelindex);
    const auto nnlevel = levelpops[level];

    // de-excitation
    const auto alltrans_startdown = get_alltrans_startdown(uniquelevelindex);
    const int ndowntrans = get_ndowntrans(uniquelevelindex);
    const auto ndowntransindices = std::views::iota(0, ndowntrans);
    std::for_each(ndowntransindices.begin(), ndowntransindices.end(), [&](const auto& i) {
      const auto alltransindex = alltrans_startdown + i;
      const int lower = globals::alltrans.targetlevelindex[alltransindex];
      assert_testmodeonly(lower >= 0 && lower < nlevels);
      const int lower_index = levelindices[lower];
      if (lower_index == level_index) {
        // both endpoints are folded into the same superlevel, so every rate below would be subtracted from
        // and added back to the same matrix element and cancel. Such a transition only redistributes
        // population within the superlevel, whose internal level ratios are held fixed by construction.
        return;
      }
      const auto lower_uniquelevelindex = ionuniquelevelindexstart + lower;
      const auto lower_statweight = stat_weight(lower_uniquelevelindex);
      const auto nnlevel_lower = levelpops[lower];

      const double epsilon_trans = epsilon_level - epsilon(lower_uniquelevelindex);
      const double R = rad_deexcitation_ratecoeff(epsilon_trans, globals::alltrans.einstein_A[alltransindex],
                                                  statweight, lower_statweight, nnlevel, nnlevel_lower, t_mid) *
                       s_renorm[level];
      const double C =
          col_deexcitation_ratecoeff(T_e, clumpednne, epsilon_trans, statweight, lower_statweight, alltransindex) *
          s_renorm[level];

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
      assert_testmodeonly(upper >= 0 && upper < nlevels);
      const int upper_index = levelindices[upper];
      if (upper_index == level_index) {
        // superlevel-internal transition: see the matching check in the de-excitation loop above
        return;
      }
      const auto upper_uniquelevelindex = ionuniquelevelindexstart + upper;
      const double epsilon_trans = epsilon(upper_uniquelevelindex) - epsilon_level;
      const auto upper_statweight = stat_weight(upper_uniquelevelindex);
      const auto nnlevel_upper = levelpops[upper];

      const double R =
          rad_excitation_ratecoeff(nonemptymgi, upper_statweight, globals::alltrans.einstein_A[alltransindex],
                                   epsilon_trans, nnlevel, nnlevel_upper, statweight, alltransindex, t_mid) *
          s_renorm[level];

      const double C =
          col_excitation_ratecoeff(T_e, clumpednne, upper_statweight, alltransindex, epsilon_trans, statweight) *
          s_renorm[level];

      const double NTC =
          nonthermal::nt_excitation_ratecoeff(nonemptymgi, level, upper, alltransindex) * s_renorm[level];

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
  const auto T_e = grid::Te_allcells[nonemptymgi];
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
  const auto T_e = grid::Te_allcells[nonemptymgi];
  const float clumpednne = grid::get_clumpfactor(nonemptymgi) * grid::get_nne(nonemptymgi);
  const int nlevels = get_nlevels(element, ion);
  for (int level = 0; level < nlevels; level++) {
    const int level_index = get_nlte_vector_index(element, ion, level, first_ion_used);
    const double epsilon_level = epsilon(element, ion, level);
    const double statweight = stat_weight(element, ion, level);

    const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
    const int nautoiondowntrans = get_nautoiondowntrans(uniquelevelindex);
    for (int i = 0; i < nautoiondowntrans; i++) {
      // autoionisation (which is a de-excitation process)
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
        printlnlog("  [warning] cell {} ts {}: negative autoionisation rate from Z={} ionstage {} level {} to level {}",
                   nltelog.modelgridindex, nltelog.timestep, get_atomicnumber(element), get_ionstage(element, ion),
                   level, target_level);
      }

      // capture (which is an excitation process, and the first part of di-electronic recomb).
      // the rate is out of the target ion's level, so if that level is folded into the target ion's
      // superlevel, weight the rate by its Boltzmann fraction of the superlevel population
      R = clumpednne * A_a * statweight / stat_weight(element, target_ion, target_level) * SAHACONST * pow(T_e, -1.5) *
          exp(-epsilon_trans / KB / T_e) * s_renorm_allions[target_ion][target_level];

      rate_matrices.autoion[(lower_index * nlte_dimension) + lower_index] -= R;
      rate_matrices.autoion[(upper_index * nlte_dimension) + lower_index] += R;
      if ((R < 0)) {
        printlnlog("  [warning] cell {} ts {}: negative autoionisation rate from Z={} ionstage {} level {} to level {}",
                   nltelog.modelgridindex, nltelog.timestep, get_atomicnumber(element), get_ionstage(element, ion),
                   level, target_level);
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
  printlnlog("  [warning] cell {} ts {}: NLTEPOP setting element Z={} level pops to LTE", nltelog.modelgridindex,
             nltelog.timestep, get_atomicnumber(element));
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
              "  [warning] cell {} ts {}: Boltzmann population for Z={} ionstage {} level {} is non-finite or "
              "non-positive ({}). Setting to MINPOP",
              nltelog.modelgridindex, nltelog.timestep, get_atomicnumber(element), get_ionstage(element, ion), level,
              ltepop);
          ltepop = MINPOP;
        }

        printlnlog(
            "  [warning] cell {} ts {}: NLTE solver gave negative/non-finite population to index {} (Z={} ionstage {} "
            "level {}), pop = {:g}. Replacing with LTE pop of {:g}",
            nltelog.modelgridindex, nltelog.timestep, index, get_atomicnumber(element), ionstage, level, population,
            ltepop);
        popvec[index] = ltepop;
      }
    } else {
      // A non-finite population (reachable via the popvec[i] = vec_x[i] * pop_normfactors[i] denormalisation
      // product overflowing - the solver itself rejects a non-finite solution vector through the residual score)
      // fails every comparison in the triage below, so without this rejection it would pass through unhandled and
      // trip the assert_always(std::isfinite(pop)) in nltepop_apply_solution, aborting the run instead of letting
      // the caller drop an ion stage and retry.
      if (!std::isfinite(population)) {
        printlnlog(
            "  [warning] cell {} ts {}: NLTE solver gave non-finite population for index {} (Z={} ionstage {} level "
            "{}), pop = {:g}. Returning nltepop_matrix_solve failure",
            nltelog.modelgridindex, nltelog.timestep, index, get_atomicnumber(element), ionstage, level, population);
        // clear the failed attempt so that can_remove_ion's triage reads zeros, as it does after a singular-matrix
        // failure, rather than non-finite values whose comparisons would all be meaningless
        std::ranges::fill(popvec, 0.);
        return false;
      }

      // if groundpop is below MINPOP then the NLTE solution fails
      const size_t index_ion_ground = get_nlte_vector_index(element, ion, 0, first_ion_used);
      if (index == index_ion_ground && population < MINPOP) {
        printlnlog(
            "  [warning] cell {} ts {}: NLTE solver gave ground pop less than MINPOP for index {} (Z={} ionstage {} "
            "level {}), pop = {:g}. Returning nltepop_matrix_solve failure",
            nltelog.modelgridindex, nltelog.timestep, index, get_atomicnumber(element), ionstage, level, population);
        return false;
      }

      if (population < 0.0) {
        printlnlog(
            "  [warning] cell {} ts {}: NLTE solver gave negative population for index {} (Z={} ionstage {} "
            "level {}), pop = {:g}",
            nltelog.modelgridindex, nltelog.timestep, index, get_atomicnumber(element), ionstage, level, population);
        if (population < -MINPOP) {
          printlnlog(
              "  [warning] negative pop = {:g} less than -MINPOP (-{:g}) unlikely a rounding error to zero so "
              "returning nltepop_matrix_solve failure",
              population, MINPOP);
          return false;
        }
        printlnlog(
            "  [warning] negative pop = {:g} greater than -MINPOP (-{:g}) likely a rounding error to zero so continue "
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
            const bool solve_failed = (inversion_factor > STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL);
            printlnlog(
                "  [warning] cell {} ts {}: pop inversion greater than factor {:g}: (g_pop {:g})/(e_pop {:g}) = {:g} "
                "is less than (g_sw {:g})/(e_sw {:g}) = {:g} for index {} Z={} ionstage {} level {} (factor {:g} "
                "inversion) - {}",
                nltelog.modelgridindex, nltelog.timestep, STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING,
                ground_pop, population, ground_pop / population, statweight_ground, statweight,
                statweight_ground / statweight, index, get_atomicnumber(element), ionstage, level, inversion_factor,
                solve_failed ? "large pop inversion so matrix solve failed"
                             : "relatively small pop inversion so continue with NLTE solution");
            if (solve_failed) {
              return false;
            }
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
            const bool solve_failed = (inversion_factor > STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL);
            printlnlog(
                "  [warning] cell {} ts {}: superlevel pop inversion greater than factor {:g}: (g_pop "
                "{:g})/(SL_first_level_pop {:g}) = {:g} is less than (g_sw {:g})/(SL_first_level_sw {:g}) = {:g} for "
                "index {} Z={} ionstage {} level {} (factor {:g} inversion) - {}",
                nltelog.modelgridindex, nltelog.timestep, STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING,
                ground_pop, sublevel_pop, ground_pop / sublevel_pop, statweight_ground, statweight,
                statweight_ground / statweight, index, get_atomicnumber(element), ionstage, level, inversion_factor,
                solve_failed ? "large pop inversion so matrix solve failed"
                             : "relatively small pop inversion so continue with NLTE solution");
            if (solve_failed) {
              return false;
            }
          }
        }
      }
    }
  }
  return true;
}

// row-major so that the rate matrix assembled in a flat std::vector can be mapped without a copy
using RowMajorMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

// Check the diagonal of an LU decomposition for zero or non-finite elements, which indicate levels that are
// disconnected from the rest of the NLTE rate matrix. Reports and returns true if the matrix is singular.
[[nodiscard]] auto lumatrix_is_singular(const RowMajorMatrixXd& lumatrix, const int element, const int first_ion_used,
                                        const int nions_used) -> bool {
  bool is_singular = false;
  for (Eigen::Index i = 0; i < lumatrix.rows(); i++) {
    // diagonal elements of LU matrix must be non-zero and finite
    const double matrixelement = lumatrix(i, i);
    if (matrixelement == 0. || !std::isfinite(matrixelement)) {
      if (!is_singular) {
        printlog("  [error] cell {} ts {}: NLTE matrix is singular for element Z={} (disconnected:",
                 nltelog.modelgridindex, nltelog.timestep, get_atomicnumber(element));
        is_singular = true;
      }
      printlog_nlte_vector_index(i, element, first_ion_used, nions_used);
    }
  }
  if (is_singular) {
    printlnlog(")");
  }
  return is_singular;
}

// The non-finite detection in this solver relies on std::isfinite and NaN comparisons keeping their IEEE semantics,
// which the FASTMATH build preserves only because -fno-finite-math-only is appended after -ffast-math (Makefile).
#if defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__
#error \
    "the NLTE solver requires std::isfinite to detect NaN/inf: compile without -ffinite-math-only (see FASTMATH in the Makefile)"
#endif

// A non-finite entry anywhere makes a solve meaningless. Report and return false, which lets the caller drop an
// ion stage and retry, instead of aborting the run. balance_vector may be empty (the GTH path has none).
[[nodiscard]] auto nlte_system_is_finite(const std::span<const double> rate_matrix,
                                         const std::span<const double> balance_vector, const int element) -> bool {
  const auto is_finite = [](const double x) { return std::isfinite(x); };
  if (!std::ranges::all_of(rate_matrix, is_finite) || !std::ranges::all_of(balance_vector, is_finite)) {
    printlnlog(
        "  [error] cell {} ts {}: NLTE rate matrix or balance vector contains a non-finite value for element Z={}",
        nltelog.modelgridindex, nltelog.timestep, get_atomicnumber(element));
    return false;
  }
  return true;
}

// warn above this componentwise relative backward error (from get_max_relative_residual) on either solver path
constexpr double RELATIVE_RESIDUAL_WARN_TOLERANCE = 1e-8;

// The largest componentwise relative backward error of a solution: max_i |b_i - (A*x)_i| / (|b_i| + sum_j |A_ij *
// x_j|), i.e. each row's residual measured against the size of the terms that were differenced to form it. The
// residual on its own is not a convergence measure: nltepop_matrix_normalise() rescales each row by an arbitrary
// accumulated balancing factor, so the raw residual elements carry no fixed physical magnitude, while in the
// per-row ratio those factors cancel. The measure must be per-row rather than a single normwise scale because the
// equilibration only balances each row's norm against its column's, never rows against each other, and with
// FORCE_SAHA_ION_BALANCE the constraint rows carry right hand sides that differ by many orders of magnitude - a
// normwise ratio would let a badly violated small-scale row hide behind the largest row. A row whose scale is
// exactly zero has an exactly zero residual and is skipped. An empty balance_vector is treated as all-zero (the
// homogeneous system solved on the GTH path).
[[nodiscard]] auto get_max_relative_residual(const std::span<const double> rate_matrix,
                                             const std::span<const double> balance_vector,
                                             const std::span<const double> vec_x) -> double {
  const auto nlte_dimension = vec_x.size();
  double max_relative_residual = 0.;
  for (auto i = 0ZU; i < nlte_dimension; i++) {
    const double b_i = balance_vector.empty() ? 0. : balance_vector[i];
    double rowdot = 0.;
    double rowscale = std::abs(b_i);
    for (auto j = 0ZU; j < nlte_dimension; j++) {
      const double term = rate_matrix[(i * nlte_dimension) + j] * vec_x[j];
      rowdot += term;
      rowscale += std::abs(term);
    }
    if (rowscale > 0.) {
      max_relative_residual = std::max(max_relative_residual, std::abs(b_i - rowdot) / rowscale);
    }
  }
  return max_relative_residual;
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

  if (!nlte_system_is_finite(rate_matrix, balance_vector, element)) {
    return false;
  }

  // solution vector for the matrix equation
  THREADLOCALONHOST std::vector<double> vec_x;
  vec_x.reserve(max_nlte_dimension);
  vec_x.resize(nlte_dimension);

  auto eigen_vec_x = Eigen::Map<Eigen::VectorX<double>>(vec_x.data(), nlte_dimension);
  const auto eigen_balance_vector = Eigen::Map<const Eigen::VectorX<double>>(balance_vector.data(), nlte_dimension);

  const auto eigen_rate_matrix = Eigen::Map<const RowMajorMatrixXd>(rate_matrix.data(), nlte_dimension, nlte_dimension);

  THREADLOCALONHOST Eigen::PartialPivLU<RowMajorMatrixXd> eigen_rate_matrix_lu;

  eigen_rate_matrix_lu.compute(eigen_rate_matrix);
  if (lumatrix_is_singular(eigen_rate_matrix_lu.matrixLU(), element, first_ion_used, nions_used)) {
    return false;
  }

  // a non-finite solution is not checked for here: every column of the rate matrix has a non-zero element (an
  // all-zero column would make the matrix singular, which is rejected above), so a non-finite element of the solution
  // vector makes at least one element of the residual below non-finite, wherever it lands, and the residual reduction
  // reports that as a non-finite error. It then leaves error_best negative, which is reported and returned as a
  // solver failure.
  eigen_vec_x = eigen_rate_matrix_lu.solve(eigen_balance_vector);

  // population solution vector with lowest error
  THREADLOCALONHOST std::vector<double> vec_x_best;
  vec_x_best.reserve(max_nlte_dimension);
  vec_x_best.resize(nlte_dimension);

  THREADLOCALONHOST std::vector<double> vec_residual;
  vec_residual.reserve(max_nlte_dimension);
  vec_residual.resize(nlte_dimension);

  // correction applied to the solution vector by the iterative refinement below. Solving into a named buffer and
  // adding it on avoids the evaluator temporary that Eigen heap-allocates for `vec_x += lu.solve(...)` on every
  // iteration (the permuted-rhs temporary inside solve itself remains either way).
  THREADLOCALONHOST std::vector<double> vec_correction;
  vec_correction.reserve(max_nlte_dimension);
  vec_correction.resize(nlte_dimension);

  auto eigen_vec_residual = Eigen::Map<Eigen::VectorX<double>>(vec_residual.data(), nlte_dimension);
  auto eigen_vec_correction = Eigen::Map<Eigen::VectorX<double>>(vec_correction.data(), nlte_dimension);

  constexpr int maxiterations = 10;
  double error_best = -1.;
  int iteration_best = -1;  // 1-based pass that produced vec_x_best (pass 1 is the unrefined LU solution)
  for (int iteration = 0; iteration < maxiterations; iteration++) {
    // one step of iterative refinement, x <- x + A^-1 (b - A x), using the residual left by the previous iteration
    if (iteration > 0) {
      eigen_vec_correction = eigen_rate_matrix_lu.solve(eigen_vec_residual);
      eigen_vec_x += eigen_vec_correction;
    }

    // residual = b - A x
    eigen_vec_residual = eigen_balance_vector - eigen_rate_matrix * eigen_vec_x;

    // PropagateNaN is required for the non-finite detection below: Eigen documents the default (PropagateFast) as
    // undefined in the presence of a NaN, and its scalar and SIMD paths genuinely differ there.
    const double error = eigen_vec_residual.cwiseAbs().maxCoeff<Eigen::PropagateNaN>();

    // only ever store a finite error, so that a non-finite iteration cannot latch error_best and block a later
    // iteration from being accepted. If every iteration is non-finite then error_best stays negative and the solve is
    // reported as a failure below.
    if (std::isfinite(error) && (error_best < 0. || error < error_best)) {
      std::ranges::copy(vec_x, vec_x_best.begin());
      error_best = error;
      iteration_best = iteration + 1;
    }

    // deliberately not converted into a usable relative tolerance, which would truncate the best-of-ten selection
    // above and change the solution. As an absolute threshold on residual rows that carry arbitrary equilibration
    // factors, it is unreachable in any realistic solve, so every solve runs all of the refinement passes.
    if (error < 1e-40) {
      break;
    }
  }

  if (!(error_best >= 0.)) {
    printlnlog("  [error] cell {} ts {}: NLTE solver failed to find a solution with finite error for element Z={}",
               nltelog.modelgridindex, nltelog.timestep, get_atomicnumber(element));
    return false;
  }
  std::ranges::copy(vec_x_best, vec_x.begin());

  // error_best carries arbitrary per-row equilibration factors, so convergence is judged by the componentwise
  // relative backward error rather than by comparison against a fixed absolute constant, which either fired on
  // every healthy solve or never fired at all depending only on the element number density.
  const double max_relative_residual = get_max_relative_residual(rate_matrix, balance_vector, vec_x);
  if (max_relative_residual > RELATIVE_RESIDUAL_WARN_TOLERANCE) {
    printlnlog(
        "  [warning] cell {} ts {}: NLTE solver iterative refinement: the best solution vector was from pass {} of "
        "up to {} (pass 1 is the unrefined LU solution), with a max residual of {:g} in the equilibrated system and "
        "a max componentwise relative backward error of {:g}",
        nltelog.modelgridindex, nltelog.timestep, iteration_best, maxiterations, error_best, max_relative_residual);
  }

  // get the unnormalised populations from the x solution vector and the normalisation factors
  for (auto i = 0ZU; i < nlte_dimension; i++) {
    popvec[i] = vec_x[i] * pop_normfactors[i];
  }

  return solution_pops_are_valid(nonemptymgi, element, popvec, first_ion_used, nions_used);
}

// GTH counterpart of nltepop_matrix_solve: solve for the stationary populations of the assembled rate matrix
// directly, with no normalisation row, balance vector, equilibration, or iterative refinement. popvec must be
// zero-filled on entry and holds a solution only on success.
[[nodiscard]] auto nltepop_matrix_solve_gth(const int element, const int nonemptymgi, RateMatrices& rate_matrices,
                                            std::span<double> popvec, const double nnelement, const int first_ion_used,
                                            const int nions_used) -> bool {
  const auto nlte_dimension = std::ssize(popvec);

  const auto rate_matrix = rate_matrices.get_summed_rate_matrix();
  assert_always(std::ssize(rate_matrix) == (nlte_dimension * nlte_dimension));

  if (!nlte_system_is_finite(rate_matrix, {}, element)) {
    return false;
  }

  // run the elimination on a scratch copy so that the summed matrix stays pristine for the residual check below
  THREADLOCALONHOST std::vector<double> elimination_matrix;
  elimination_matrix.assign(rate_matrix.begin(), rate_matrix.end());

  if (const auto disconnected_index = gth_stationary_distribution(elimination_matrix, popvec)) {
    // popvec is still zero-filled here, which is what can_remove_ion's triage expects
    printlog(
        "  [error] cell {} ts {}: NLTE matrix is singular for element Z={} (GTH elimination found no departure rate "
        "from:",
        nltelog.modelgridindex, nltelog.timestep, get_atomicnumber(element));
    printlog_nlte_vector_index(*disconnected_index, element, first_ion_used, nions_used);
    printlnlog(")");
    return false;
  }

  // scale the unit-sum stationary distribution to the element number density
  for (double& pop : popvec) {
    pop *= nnelement;
  }

  // componentwise relative backward error of the statistical equilibrium equations rate_matrix * popvec = 0 (with
  // no normalisation row on this path, every row is a rate equation)
  const double max_relative_residual = get_max_relative_residual(rate_matrix, {}, popvec);
  if (max_relative_residual > RELATIVE_RESIDUAL_WARN_TOLERANCE) {
    printlnlog(
        "  [warning] cell {} ts {}: NLTE GTH solution for element Z={} has a max componentwise relative backward "
        "error of {:g}",
        nltelog.modelgridindex, nltelog.timestep, get_atomicnumber(element), max_relative_residual);
  }

  return solution_pops_are_valid(nonemptymgi, element, popvec, first_ion_used, nions_used);
}

auto can_remove_ion(const int element, const int ion, const int first_ion_used, const int nions_used,
                    const double nnelement, std::span<const double> popvec) -> bool {
  if (nions_used <= 1) {
    // can't remove an ion if there is only one ion stage used in the NLTE matrix
    printlnlog(
        "  [warning] cell {} ts {}: can't remove ion stage {} from NLTE matrix for element Z={} as only one ion stage "
        "used",
        nltelog.modelgridindex, nltelog.timestep, get_ionstage(element, ion), get_atomicnumber(element));
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
        "  [warning] cell {} ts {}: {} ion (Z={} ionstage {}) ground state population too large to remove ion "
        "(ground_pop / nnelement) = ({}/{}) = {} > {}",
        nltelog.modelgridindex, nltelog.timestep, ionname, get_atomicnumber(element), get_ionstage(element, ion),
        ground_pop, nnelement, ground_pop / nnelement, NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION);

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
          "  [warning] cell {} ts {}: {} ion (Z={} ionstage {}) excited state (level {}) population too large to "
          "remove ion (nlte_excited_pop_{}_ion / nnelement) = ({}/{}) > ({})",
          nltelog.modelgridindex, nltelog.timestep, ionname, get_atomicnumber(element), get_ionstage(element, ion),
          level, ionname, levelpop, nnelement, NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION);
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
          "  [warning] cell {} ts {}: {} ion (Z={} ionstage {}) superlevel population too large to remove ion "
          "(superlevel_pop_{}_ion / nnelement) = ({:g}/{:g}) > ({:g})",
          nltelog.modelgridindex, nltelog.timestep, ionname, get_atomicnumber(element), get_ionstage(element, ion),
          ionname, superlevel_pop, nnelement, NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION);
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
// Return whether a solution was found and the range of ions included in it.
// Both solver paths normalise the solution so that the sum over every matrix slot equals the element number
// density: the LU path through its replaced zeroth matrix row, the GTH path by scaling its unit-sum stationary
// distribution.
auto nltepop_solve_matrix_with_ion_reduction(const int element, const int nonemptymgi, const double t_mid,
                                             const double nnelement,
                                             const std::vector<std::vector<double>>& s_renorm_allions,
                                             RateMatrices& rate_matrices, std::vector<double>& popvec)
    -> std::tuple<bool, int, int> {
  const int atomic_number = get_atomicnumber(element);
  const int nions = get_nions(element);
  const auto max_nlte_dimension = get_max_nlte_dimension();
  // the solver choice is fixed for the element: FORCE_SAHA_ION_BALANCE constraint rows are incompatible with the
  // Markov-generator structure that the GTH elimination requires, so those elements always use the LU path
  const bool use_gth_solver = NLTE_USE_GTH_SOLVER && !FORCE_SAHA_ION_BALANCE(atomic_number);
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

    if (use_gth_solver) {
      matrix_solve_success =
          nltepop_matrix_solve_gth(element, nonemptymgi, rate_matrices, popvec, nnelement, first_ion_used, nions_used);
    } else {
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

      // calculate the normalisation factors and apply them to the matrix columns and balance vector elements
      THREADLOCALONHOST std::vector<double> pop_normfactors;
      pop_normfactors.reserve(max_nlte_dimension);
      pop_normfactors.resize(nlte_dimension);
      std::ranges::fill(pop_normfactors, 1.0);
      nltepop_matrix_normalise(rate_matrix, balance_vector, pop_normfactors);

      matrix_solve_success = nltepop_matrix_solve(element, nonemptymgi, rate_matrix, balance_vector, popvec,
                                                  pop_normfactors, max_nlte_dimension, first_ion_used, nions_used);
    }

    matrix_solve_required = false;  // will be set to true if we need to retry with a different ion range
    if (matrix_solve_success) {
      if (nions_used < nions) {
        printlnlog(
            "successfully solved NLTE matrix when reducing ions used for element to Z={} ionstage={} to ionstage={}",
            atomic_number, get_ionstage(element, first_ion_used), get_ionstage(element, max_ion_used));
      }
    } else if (NLTE_LIMIT_ION_STAGES_AFTER_FAILURE) {
      printlnlog("  [warning] cell {} ts {}: NLTE matrix solution failed for element Z={} using ionstage {} to {}",
                 nltelog.modelgridindex, nltelog.timestep, atomic_number, get_ionstage(element, first_ion_used),
                 get_ionstage(element, max_ion_used));

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
            "  [warning] cell {} ts {}: can't remove top or bottom ion stage from NLTE equations, therefore unable to "
            "find an NLTE solution for element Z={}",
            nltelog.modelgridindex, nltelog.timestep, atomic_number);
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
      printlnlog(
          "  [warning] cell {} ts {}: Z={} ionstage {} removed from NLTE rate matrix. Setting all levelpops for ion to "
          "zero",
          modelgridindex, timestep, get_atomicnumber(element), get_ionstage(element, ion));

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
        "  [warning] timestep {} nlteiter {} Z={} element population is: {:g} (from abundance) and {:g} (from matrix "
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

// Grassmann-Taksar-Heyman (GTH) state-elimination solve for the stationary distribution of the continuous-time
// Markov chain whose generator is stored transposed in the NLTE rate matrix layout, i.e. rate_matrix[(to * n) +
// from] with n = vec_x.size(). The off-diagonal rates must be finite (the caller checks this) and non-negative.
// Only sums, products, and divisions of the rates are used
// and the diagonal is never read, so there is no subtractive cancellation and assembly rounding in the
// column sums cannot affect the result. Small negative autoionisation off-diagonals (warned about during
// assembly) only weaken that guarantee locally; any resulting invalid populations are policed by the caller.
// rate_matrix is overwritten. On success, returns std::nullopt and fills vec_x with the stationary distribution
// normalised to a sum of one. On failure, returns the index of a state with no departure rate into the remaining
// chain (a reducible/disconnected matrix) and leaves vec_x untouched. An off-diagonal that grows to infinity
// during the elimination is not a failure here: it is passed through as a non-finite distribution for the
// caller's population validity check to handle, rather than being rescued (see the back-substitution below).
// Defined with external linkage (declared in nltepop.h) so that unittests.cc can exercise it.
auto gth_stationary_distribution(std::span<double> rate_matrix, std::span<double> vec_x)
    -> std::optional<std::ptrdiff_t> {
  const auto n = std::ssize(vec_x);
  assert_always(n >= 1);
  assert_always(std::ssize(rate_matrix) == (n * n));

  // eliminate the highest-numbered remaining state by folding it into the transition rates among the states below
  for (auto k = n - 1; k >= 1; k--) {
    double departure_sum = 0.;
    for (auto i = 0Z; i < k; i++) {
      departure_sum += rate_matrix[(i * n) + k];
    }
    if (!std::isfinite(departure_sum) || !(departure_sum > 0.)) {
      return k;
    }
    // stash the departure sum in the diagonal slot, which GTH never otherwise reads
    rate_matrix[(k * n) + k] = departure_sum;
    for (auto i = 0Z; i < k; i++) {
      const double rate_k_to_i = rate_matrix[(i * n) + k];
      if (rate_k_to_i == 0.) {
        continue;
      }
      const double factor = rate_k_to_i / departure_sum;
      for (auto j = 0Z; j < k; j++) {
        rate_matrix[(i * n) + j] += factor * rate_matrix[(k * n) + j];
      }
    }
  }

  // back-substitution: stationary weights, starting relative to vec_x[0] = 1
  vec_x[0] = 1.;
  for (auto k = 1Z; k < n; k++) {
    const double departure_sum = rate_matrix[(k * n) + k];
    const auto get_inflow = [&] {
      double inflow = 0.;
      for (auto j = 0Z; j < k; j++) {
        inflow += vec_x[j] * rate_matrix[(k * n) + j];
      }
      return inflow;
    };
    double inflow = get_inflow();
    // subtraction-free overflow rescue, applied before the division so that even a single weight ratio beyond the
    // double range (e.g. rates of 1e200 up and 1e-200 back between two states) cannot overflow to infinity: while
    // the next weight would exceed 1e250, rescale the earlier weights and recompute the inflow. The ratios are
    // preserved, weights that underflow to zero are at least ~1e58 times smaller than the largest weight and
    // physically negligible, and capping every stored weight at 1e250 keeps the normalising sum below from
    // overflowing.
    //
    // The number of passes is capped because rescaling cannot rescue every non-finite inflow. Once a weight has
    // underflowed to zero, an infinite rate out of that state gives inflow = 0 * inf = NaN, which survives any
    // number of further passes, so an uncapped loop would spin here forever. Three passes already span the whole
    // dynamic range of a double, so the cap never cuts a rescue that would have succeeded. Giving up leaves a
    // non-finite population for solution_pops_are_valid() to reject or replace, as for any other unusable solve.
    constexpr int max_rescale_passes = 8;
    for (int pass = 0; pass < max_rescale_passes && (!std::isfinite(inflow) || inflow > departure_sum * 1e250);
         pass++) {
      for (auto j = 0Z; j < k; j++) {
        vec_x[j] *= 1e-200;
      }
      inflow = get_inflow();
    }
    vec_x[k] = inflow / departure_sum;
  }

  const double x_sum = std::accumulate(vec_x.begin(), vec_x.end(), 0.);
  for (double& x : vec_x) {
    x /= x_sum;
  }

  return std::nullopt;
}

// Solve the statistical balance equations to find NLTE level populations for all ions of an element
// (ionisation balance follows from this too). Failure modes can lead to quasi-LTE populations being set
// instead.
void solve_nlte_pops_element(const int element, const int nonemptymgi, const int timestep, const int nlte_iter) {
  const int atomic_number = get_atomicnumber(element);
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  nltelog = {.modelgridindex = modelgridindex, .timestep = timestep, .nlte_iter = nlte_iter};

  if (grid::get_elem_massfrac(nonemptymgi, element) <= 0.) {
    // abundance of this element is zero, so do not store any NLTE populations
    nltepop_reset_element(nonemptymgi, element);
    return;
  }

  const double cell_Te = grid::Te_allcells[nonemptymgi];

  if (cell_Te == MINTEMP) {
    printlnlog(
        "Not solving for NLTE populations in cell {} at timestep {} NLTE iteration {} for element Z={} due to low "
        "temperature Te=MINTEMP={:g} [K]",
        modelgridindex, timestep, nlte_iter, atomic_number, cell_Te);
    set_element_pops_lte(nonemptymgi, element);
    return;
  }

  const auto sys_time_start_nltesolver = std::chrono::steady_clock::now();

  const double t_mid = globals::timesteps[timestep].mid;
  const int nions = get_nions(element);
  const double nnelement = grid::get_elem_numberdens(nonemptymgi, element);

  if (nlte_iter == 0) {
    // the per-iteration progress is summarised by the NLTE convergence messages in update_grid.cc,
    // so this announcement is only made for the first iteration
    printlnlog(
        "Solving for NLTE populations in cell {} at timestep {} for element Z={} (mass fraction "
        "{:.2e}, nnelement {:.2e} [cm^-3])",
        modelgridindex, timestep, atomic_number, grid::get_elem_massfrac(nonemptymgi, element), nnelement);
  }

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
        "[warning] Can't solve for NLTE populations in cell {} at timestep {} for element Z={} due to singular matrix, "
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
  const double T_exc = LTEPOP_EXCITATION_USE_TJ ? grid::TJ_allcells[nonemptymgi] : grid::Te_allcells[nonemptymgi];
  const double E_level = epsilon(element, ion, level);
  const double E_superlevel = epsilon(element, ion, level_superlevel_start);

  return stat_weight(element, ion, level) / stat_weight(element, ion, level_superlevel_start) *
         exp(-(E_level - E_superlevel) / KB / T_exc);
}

void nltepop_open_file() {
  nlte_file = open_rank_outfile("nlte");
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
    printlnlog("[error] Expected {} NLTE levels but found {} in restart file", globals::total_nlte_levels,
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

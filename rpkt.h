#ifndef RPKT_H
#define RPKT_H

#include <cstddef>
#include <ctime>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "ltepop.h"
#include "packet.h"
#include "sn3d.h"

struct Phixslist {
  std::vector<double> groundcont_gamma_contr;  // for either USE_LUT_PHOTOION = true or USE_LUT_BFHEATING = true
  std::vector<double> chi_bf_sum;
  std::vector<double> gamma_contr;  // needed for DETAILED_BF_ESTIMATORS_ON
  int allcontend{-1};
  int allcontbegin{0};
  int bfestimend{-1};
  int bfestimbegin{0};
};

struct Rpkt_continuum_absorptioncoeffs {
  double nu{-1.};  // frequency at which opacity was calculated
  double total{0.};
  double ffescat{0.};
  double ffheat{0.};
  double bf{0.};
  int nonemptymgi{-1};
  int timestep{-1};
  Phixslist phixslist{};
};

void do_rpkt(Packet &pkt, double t2);
void emit_rpkt(Packet &pkt);
void calculate_chi_rpkt_cont(double nu_cmf, Rpkt_continuum_absorptioncoeffs &chi_rpkt_cont, int nonemptymgi);
[[nodiscard]] auto sample_planck_times_expansion_opacity(int nonemptymgi) -> double;
void allocate_expansionopacities();
void calculate_expansion_opacities(int nonemptymgi);
void MPI_Bcast_binned_opacities(ptrdiff_t nonemptymgi, int root_node_id);

[[nodiscard]] constexpr auto get_linedistance(const double prop_time, const double nu_cmf, const double nu_trans,
                                              const double d_nu_on_d_l) -> double {
  // distance from packet position to redshifting into line at frequency nu_trans

  if (nu_cmf <= nu_trans) {
    return 0;  // photon was propagated too far, make sure that we don't miss a line
  }

  if constexpr (USE_RELATIVISTIC_DOPPLER_SHIFT) {
    // With special relativity, the Doppler shift formula has an extra factor of 1/gamma in it,
    // which changes the distance reach a line resonance and creates a dependence
    // on packet position and direction

    // use linear interpolation of frequency along the path
    return (nu_trans - nu_cmf) / d_nu_on_d_l;
  }

  return CLIGHT * prop_time * (nu_cmf / nu_trans - 1);
}

constexpr auto closest_transition(const double nu_cmf, const int next_trans, const int nlines,
                                  const auto *const linelist) -> int
// for the propagation through non empty cells
// find the next transition lineindex redder than nu_cmf
// return -1 if no transition can be reached
{
  if (next_trans > (nlines - 1)) {
    // packet is tagged as having no more line interactions
    return -1;
  }
  // if nu_cmf is smaller than the lowest frequency in the linelist,
  // no line interaction is possible: return negative value as a flag
  if (nu_cmf < linelist[nlines - 1].nu) {
    return -1;
  }

  if (next_trans > 0) [[likely]] {
    // if next_trans > 0 we know the next line we should interact with, independent of the packets
    // current nu_cmf which might be smaller than globals::linelist[left].nu due to propagation errors
    return next_trans;
  }
  if (nu_cmf >= linelist[0].nu) {
    // if nu_cmf is larger than the highest frequency in the the linelist,
    // interaction with the first line occurs - no search
    return 0;
  }
  // otherwise go through the list until nu_cmf is located between two
  // entries in the line list and get the index of the closest line
  // to lower frequencies

  // will find the highest frequency (lowest index) line with nu_line <= nu_cmf
  // lower_bound matches the first element where the comparison function is false
  const int matchindex = static_cast<int>(
      std::lower_bound(linelist, linelist + nlines, nu_cmf,
                       [](const auto &line, const double find_nu_cmf) -> bool { return line.nu > find_nu_cmf; }) -
      linelist);

  if (matchindex >= nlines) [[unlikely]] {
    return -1;
  }

  return matchindex;
}

[[nodiscard]] inline auto get_ionestimindex_nonemptymgi(const int nonemptymgi, const int element, const int ion)
    -> int {
  assert_testmodeonly(ion >= 0);
  assert_testmodeonly(ion < get_nions(element) - 1);
  const int groundcontindex = globals::elements[element].ions[ion].groundcontindex;
  assert_always(groundcontindex >= 0);
  return (nonemptymgi * globals::nbfcontinua_ground) + groundcontindex;
}

inline auto keep_this_cont(int element, const int ion, const int level, const int nonemptymgi, const float nnetot)
    -> bool {
  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    return grid::get_elem_abundance(nonemptymgi, element) > 0;
  }
  return ((get_nnion(nonemptymgi, element, ion) / nnetot > 1.e-6) || (level == 0));
}

#endif  // RPKT_H

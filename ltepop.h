#ifndef LTEPOP_H
#define LTEPOP_H

#include <cmath>
#include <vector>

#include "constants.h"

[[gnu::pure]] [[nodiscard]] auto get_groundlevelpop(int nonemptymgi, int element, int ion) -> double;

[[gnu::pure]] [[nodiscard]] auto calculate_levelpop(int nonemptymgi, int element, int ion, int level) -> double;

[[gnu::pure]] [[nodiscard]] auto calculate_levelpop_boltzmann(int nonemptymgi, int element, int ion, int level)
    -> double;

[[gnu::pure]] [[nodiscard]] auto get_levelpop(int nonemptymgi, int uniquelevelindex) -> double;
[[gnu::pure]] [[nodiscard]] auto get_levelpop(int nonemptymgi, int element, int ion, int level) -> double;
[[gnu::pure]] [[nodiscard]] auto get_nnion(int nonemptymgi, int element, int ion) -> double;
[[nodiscard]] auto find_uppermost_ion(int nonemptymgi, int element, double nne_hi, bool force_saha) -> int;
void calculate_ion_balance_nne(int nonemptymgi);
void calculate_cellpartfuncts(int nonemptymgi, int element);
[[nodiscard]] auto calculate_ionfractions(int element, int nonemptymgi, double nne, bool use_phi_saha)
    -> std::vector<double>;
void set_groundlevelpops(int nonemptymgi, int element, float nne, bool force_saha);

// calculates saha factor in LTE: Phi_level,ion,element = nn_level,ion,element/(nne*nn_upper,ion+1,element)
[[gnu::const]] [[nodiscard]] constexpr auto calculate_sahafact(const double g_lower, const double g_upper,
                                                               const double T, const double E_threshold) -> double {
  const double sf = SAHACONST * g_lower / g_upper * std::pow(T, -1.5) * std::exp(E_threshold / KB / T);

  return sf;
}

#endif  // LTEPOP_H

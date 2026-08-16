// Declarations for photoionisation/recombination rate coefficients (ratecoeff.cc).

#ifndef RATECOEFF_H
#define RATECOEFF_H

#include <cstdint>

#include "constants.h"
#include "random.h"

void ratecoefficients_init();

void setup_photoion_luts();

[[nodiscard]] DEVICE_FUNC auto select_continuum_nu(int element, int lowerion, int lower, int phixstargetindex,
                                                   float T_e, rngstate_type& rngstate) -> double;
[[nodiscard]] DEVICE_FUNC auto get_spontrecombcoeff(int uniquelevelindex, int phixstargetindex, float T_e) -> double;

[[nodiscard]] DEVICE_FUNC auto get_bfcoolingcoeff(int element, int lowerion, int lowerionlevel, int phixstargetindex,
                                                  float T_e) -> double;

[[nodiscard]] DEVICE_FUNC auto get_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex,
                                                     int nonemptymgi, bool use_cellcache) -> double;
[[nodiscard]] auto get_corrphotoioncoeff_ana(int element, int ion, int level, int phixstargetindex, float T_R)
    -> double;

[[nodiscard]] auto iongamma_is_zero(int nonemptymgi, int element, int ion) -> bool;

[[nodiscard]] auto calculate_iongamma_per_gspop(int nonemptymgi, int element, int ion) -> double;
[[nodiscard]] auto calculate_iongamma_per_ionpop(int nonemptymgi, int element, int lowerion,
                                                 bool collisional_not_radiative, bool force_bfintegral) -> double;

// Population that the ion recombination coefficient returned by calculate_ionrecombcoeff() is normalised to
// (i.e. what it must be multiplied by, along with nne, to give a rate)
enum class IonRecombNorm : std::uint8_t {
  IONPOP,  // the total upper ion population (all upper levels)
  GROUNDMULTIPLETPOP,  // the population of the upper ion's ground multiplet
  TARGETLEVELPOP,  // for each lower level, the summed population of that level's photoionisation target
                   // levels, i.e. the normalisation of get_ion_spontrecombcoeff() / phi_rate_balance()
};

struct IonRecombCoeffOptions {
  bool assume_lte = false;
  bool collisional_not_radiative = false;
  bool lower_superlevel_only = false;
  IonRecombNorm norm = IonRecombNorm::IONPOP;
};

// nonemptymgi = -1 (no cell) with options.assume_lte = true is used by the startup recombination-rate
// calibration and ion alpha_sp table, which take LTE populations at the given temperature
[[nodiscard]] auto calculate_ionrecombcoeff(int nonemptymgi, float T_e, int element, int upperion,
                                            IonRecombCoeffOptions options) -> double;

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_ion_spontrecombcoeff(int uniqueionindex, float T_e) -> double;

#endif  // RATECOEFF_H

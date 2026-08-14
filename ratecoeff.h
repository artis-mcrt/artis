// Declarations for photoionisation/recombination rate coefficients (ratecoeff.cc).

#ifndef RATECOEFF_H
#define RATECOEFF_H

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

[[nodiscard]] auto calculate_ionrecombcoeff(int nonemptymgi, float T_e, int element, int upperion, bool assume_lte,
                                            bool collisional_not_radiative, bool lower_superlevel_only,
                                            bool per_groundmultipletpop) -> double;

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_ion_spontrecombcoeff(int uniqueionindex, float T_e) -> double;

#endif  // RATECOEFF_H

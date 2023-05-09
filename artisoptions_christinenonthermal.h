#ifndef ARTISOPTIONS_H  // NOLINT(llvm-header-guard)
#define ARTISOPTIONS_H
// NOLINTBEGIN(modernize*,misc-unused-parameters)

#include <cstdlib>

#include "constants.h"

constexpr int MPKTS = 10000;

constexpr int GRID_TYPE = GRID_UNIFORM;
constexpr int CUBOID_NCOORDGRID_X = 100;
constexpr int CUBOID_NCOORDGRID_Y = 100;
constexpr int CUBOID_NCOORDGRID_Z = 100;

constexpr bool NLTE_POPS_ON = true;

constexpr bool NLTE_POPS_ALL_IONS_SIMULTANEOUS = true;

constexpr int NLTEITER = 30;

constexpr bool LEVEL_IS_NLTE(int element_z, int ionstage, int level) {
  if (element_z < 22) {
    return (level <= 200);
  }
  return (level <= 300);
}

#define LTEPOP_EXCITATIONTEMPERATURE grid::get_Te(modelgridindex)

constexpr bool single_level_top_ion = false;

constexpr bool single_ground_level = false;

constexpr int NLEVELS_REQUIRETRANSITIONS(int Z, int ionstage) {
  return ((Z == 26 && ionstage == 5) && ionstage >= 1) ? 80 : 0;
}

constexpr bool UNIFORM_PELLET_ENERGIES = true;

constexpr bool DIRECT_COL_HEAT = true;
constexpr bool NO_INITIAL_PACKETS = false;
constexpr bool RECORD_LINESTAT = false;

constexpr bool USE_MODEL_INITIAL_ENERGY = true;

constexpr int TABLESIZE = 100;
constexpr double MINTEMP = 3000.;
constexpr double MAXTEMP = 140000.;

constexpr double RECOMBCALIBRATION_T_ELEC = 15000.;

constexpr bool DIPOLE = false;
constexpr bool POL_ON = false;

// #define VPKT_ON
constexpr size_t GSLWSIZE = 16384;

constexpr bool TRACK_ION_STATS = false;
constexpr bool TRACK_ION_MASTATS = false;

constexpr double MINPOP = 1e-40;

constexpr double NU_MIN_R = 1e13;
constexpr double NU_MAX_R = 5e16;

constexpr bool MULTIBIN_RADFIELD_MODEL_ON = true;

constexpr int RADFIELDBINCOUNT = 64;

constexpr int FIRST_NLTE_RADFIELD_TIMESTEP = 12;

constexpr double nu_lower_first_initial = (CLIGHT / (40000e-8));
constexpr double nu_upper_last_initial = (CLIGHT / (500e-8));
constexpr double nu_upper_superbin = (CLIGHT / (50e-8));

constexpr double T_R_min = 500;
constexpr double T_R_max = 250000;

constexpr bool DETAILED_LINE_ESTIMATORS_ON = false;

constexpr bool DETAILED_BF_ESTIMATORS_ON = true;

constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP = 13;

#define DETAILED_BF_ESTIMATORS_BYTYPE false

constexpr bool NO_LUT_PHOTOION = true;

constexpr bool NO_LUT_BFHEATING = true;

#define SEPARATE_STIMRECOMB false

constexpr bool NT_ON = true;

constexpr bool NT_SOLVE_SPENCERFANO = true;

constexpr int SFPTS = 4096;

constexpr double SF_EMAX = 16000;

constexpr double SF_EMIN = 0.1;

#define SF_USE_LOG_E_INCREMENT false

constexpr int SF_MAX_TIMESTEPS_BETWEEN_SOLUTIONS = 0;

constexpr double NT_MAX_FRACDIFF_NNEPERION_BETWEEN_SOLUTIONS = 0.05;

constexpr int NTEXCITATION_MAXNLEVELS_LOWER = 5;
constexpr int NTEXCITATION_MAXNLEVELS_UPPER = 250;

constexpr int MAX_NT_EXCITATIONS_STORED = 25000;

#define NT_EXCITATION_ON false

constexpr int NT_BLOCKSIZEEXCITATION = 5192;

constexpr bool NT_USE_VALENCE_IONPOTENTIAL = false;

constexpr int NT_MAX_AUGER_ELECTRONS = 2;

constexpr bool SF_AUGER_CONTRIBUTION_ON = true;

constexpr bool SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN = false;

constexpr double TEMPERATURE_SOLVER_ACCURACY = 1e-3;

constexpr double CONTINUUM_NU_INTEGRAL_ACCURACY = 1e-3;

constexpr double RATECOEFF_INTEGRAL_ACCURACY = 1e-3;

constexpr double IONGAMMA_POPFRAC_LEVELS_INCLUDED = 0.999;

constexpr bool USE_RELATIVISTIC_DOPPLER_SHIFT = false;

constexpr bool USE_CALCULATED_MEANATOMICWEIGHT = false;

constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC = false;

constexpr bool INSTANT_PARTICLE_DEPOSITION = true;

enum timestepsizemethods {
  TIMESTEP_SIZES_LOGARITHMIC = 0,
  TIMESTEP_SIZES_CONSTANT = 1,
  TIMESTEP_SIZES_LOGARITHMIC_THEN_CONSTANT = 2,
  TIMESTEP_SIZES_CONSTANT_THEN_LOGARITHMIC = 3,
};

constexpr enum timestepsizemethods TIMESTEP_SIZE_METHOD = TIMESTEP_SIZES_LOGARITHMIC;

constexpr double FIXED_TIMESTEP_WIDTH = -1.;

constexpr double TIMESTEP_TRANSITION_TIME = -1.;

constexpr bool USE_GSL_RANDOM = true;

constexpr bool KEEP_ALL_RESTART_FILES = false;

constexpr bool BFCOOLING_USELEVELPOPNOTIONPOP = false;

// NOLINTEND(modernize*,misc-unused-parameters)
#endif  // ARTISOPTIONS_H
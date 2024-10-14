#ifndef ARTISOPTIONS_H  // NOLINT(llvm-header-guard)
#define ARTISOPTIONS_H
// NOLINTBEGIN(modernize*,misc-unused-parameters)

#include <cstdlib>

#include "constants.h"

constexpr int MPKTS = 100000;

constexpr auto GRID_TYPE = GridType::CARTESIAN3D;
constexpr int CUBOID_NCOORDGRID_X = 100;
constexpr int CUBOID_NCOORDGRID_Y = 100;
constexpr int CUBOID_NCOORDGRID_Z = 100;
constexpr bool FORCE_SPHERICAL_ESCAPE_SURFACE = false;

constexpr int NLTEITER = 30;

constexpr bool LEVEL_IS_NLTE(int element_z, int ionstage, int level) { return false; }

constexpr bool LTEPOP_EXCITATION_USE_TJ = true;

constexpr bool FORCE_SAHA_ION_BALANCE(int element_z) { return false; }

constexpr bool single_level_top_ion = true;

constexpr bool single_ground_level = true;

constexpr int NLEVELS_REQUIRETRANSITIONS(int Z, int ionstage) { return 0; }

constexpr bool UNIFORM_PELLET_ENERGIES = true;

constexpr bool DIRECT_COL_HEAT = false;
constexpr bool INITIAL_PACKETS_ON = true;
constexpr bool RECORD_LINESTAT = false;

constexpr bool USE_MODEL_INITIAL_ENERGY = true;

constexpr int TABLESIZE = 100;
constexpr double MINTEMP = 3500.;
constexpr double MAXTEMP = 140000.;

constexpr double RECOMBCALIBRATION_T_ELEC = 6000.;

constexpr bool DIPOLE = true;
constexpr bool POL_ON = true;

constexpr bool VPKT_ON = false;
constexpr bool VPKT_WRITE_CONTRIBS = false;

constexpr bool TRACK_ION_STATS = false;

constexpr double MINPOP = 1e-30;

constexpr double NU_MIN_R = 1e14;
constexpr double NU_MAX_R = 5e15;

constexpr bool PHIXS_CLASSIC_NO_INTERPOLATION = true;

constexpr bool MULTIBIN_RADFIELD_MODEL_ON = false;

constexpr int RADFIELDBINCOUNT = 256;

constexpr int FIRST_NLTE_RADFIELD_TIMESTEP = 12;

constexpr double nu_lower_first_initial = (CLIGHT / (40000e-8));

constexpr double nu_upper_last_initial = (CLIGHT / (1085e-8));

constexpr double nu_upper_superbin = (CLIGHT / (10e-8));

constexpr double T_R_min = 500;
constexpr double T_R_max = 250000;

constexpr bool DETAILED_LINE_ESTIMATORS_ON = false;

constexpr bool DETAILED_BF_ESTIMATORS_ON = false;

constexpr bool LEVEL_HAS_BFEST(int element_z, int ionstage, int level) { return false; }

constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP = 13;

constexpr bool USE_LUT_PHOTOION = true;

constexpr bool USE_LUT_BFHEATING = true;

#define SEPARATE_STIMRECOMB false

constexpr bool NT_ON = false;

constexpr bool NT_SOLVE_SPENCERFANO = false;

constexpr int SFPTS = 4096;

constexpr double SF_EMAX = 16000;

constexpr double SF_EMIN = 0.1;

constexpr int SF_MAX_TIMESTEPS_BETWEEN_SOLUTIONS = 0;

constexpr double NT_MAX_FRACDIFF_NNEPERION_BETWEEN_SOLUTIONS = 0.05;

constexpr int NTEXCITATION_MAXNLEVELS_LOWER = 5;
constexpr int NTEXCITATION_MAXNLEVELS_UPPER = 250;

constexpr int MAX_NT_EXCITATIONS_STORED = 25000;

constexpr bool NT_EXCITATION_ON = false;

constexpr bool NT_USE_VALENCE_IONPOTENTIAL = false;

constexpr int NT_MAX_AUGER_ELECTRONS = 2;

constexpr bool SF_AUGER_CONTRIBUTION_ON = true;

constexpr bool SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN = false;

constexpr bool NT_WORKFUNCTION_USE_SHELL_OCCUPANCY_FILE = false;

constexpr double TEMPERATURE_SOLVER_ACCURACY = 1e-2;

constexpr double CONTINUUM_NU_INTEGRAL_ACCURACY = 1e-2;

constexpr double RATECOEFF_INTEGRAL_ACCURACY = 1e-2;

constexpr double IONGAMMA_POPFRAC_LEVELS_INCLUDED = 1.;

constexpr bool USE_RELATIVISTIC_DOPPLER_SHIFT = false;

constexpr bool USE_CALCULATED_MEANATOMICWEIGHT = false;

constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC = false;

constexpr auto TIMESTEP_SIZE_METHOD = TimeStepSizeMethod::LOGARITHMIC;

constexpr double FIXED_TIMESTEP_WIDTH = -1.;

constexpr double TIMESTEP_TRANSITION_TIME = -1.;

constexpr bool KEEP_ALL_RESTART_FILES = false;

constexpr bool BFCOOLING_USELEVELPOPNOTIONPOP = false;

constexpr bool EXPANSIONOPACITIES_ON = false;

constexpr float RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY = -1.;

constexpr bool USE_XCOM_GAMMAPHOTOION = false;

constexpr auto PARTICLE_THERMALISATION_SCHEME = ThermalisationScheme::INSTANT;

constexpr auto GAMMA_THERMALISATION_SCHEME = ThermalisationScheme::DETAILED;

// NOLINTEND(modernize*,misc-unused-parameters)
#endif  // ARTISOPTIONS_H

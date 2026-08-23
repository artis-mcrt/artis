// A compile-time options preset: symlink or copy one of the artisoptions_*.h files to
// artisoptions.h to select the run configuration. Every option is described in artisoptions_doc.md.

#ifndef ARTISOPTIONS_H
#define ARTISOPTIONS_H
// NOLINTBEGIN(modernize-use-trailing-return-type,misc-unused-parameters)

#include <cstdlib>
#include <optional>

#include "constants.h"

constexpr int MPKTS = 10000;

constexpr std::optional<GridType> GRID_TYPE_OVERRIDE;
constexpr int CUBOID_NCOORDGRID_X = 100;
constexpr int CUBOID_NCOORDGRID_Y = 100;
constexpr int CUBOID_NCOORDGRID_Z = 100;
constexpr bool FORCE_SPHERICAL_ESCAPE_SURFACE = false;

constexpr int NLTEITER = 30;

constexpr bool NLTE_TE_NNE_USE_ANDERSON_ACCEL = false;

constexpr double NLTE_OUTER_RELTOL = 0.04;

constexpr int ION_NLEVELS_EXCITED_NLTE(int element_z, int ionstage) {
  if (element_z == 26 && ionstage == 2) {
    return 197;
  }
  return 80;
}

constexpr bool LTEPOP_EXCITATION_USE_TJ = false;

constexpr bool FORCE_SAHA_ION_BALANCE(int element_z) { return false; }

constexpr bool SINGLE_LEVEL_TOP_ION = true;

constexpr int NLEVELS_REQUIRETRANSITIONS(int Z, int ionstage) { return 0; }

constexpr bool UNIFORM_PELLET_ENERGIES = true;

constexpr bool DIRECT_COL_HEAT = true;
constexpr bool INITIAL_PACKETS_ON = true;

constexpr bool USE_MODEL_INITIAL_ENERGY = true;

constexpr int TABLESIZE = 200;
constexpr double MINTEMP = 4000.;
constexpr double MAXTEMP = 140000.;

constexpr double RECOMBCALIBRATION_T_ELEC = 15000.;

constexpr bool DIPOLE = true;
constexpr bool POL_ON = true;

constexpr bool VPKT_ON = false;
constexpr bool VPKT_WRITE_CONTRIBS = false;

constexpr double MINPOP = 1e-40;

constexpr double NU_MIN_R = 1e14;
constexpr double NU_MAX_R = 5e16;

constexpr bool PHIXS_CLASSIC_NO_INTERPOLATION = false;

constexpr bool MULTIBIN_RADFIELD_MODEL_ON = true;

constexpr int RADFIELDBINCOUNT = 512;

constexpr int FIRST_NLTE_RADFIELD_TIMESTEP = 12;

constexpr double RADFIELDBINS_NU_MIN = (CLIGHT / 40000e-8);
constexpr double RADFIELDBINS_NU_MAX = (CLIGHT / 100e-8);
constexpr double RADFIELDBINS_T_E_SUPERBIN_NU_MAX = (CLIGHT / 10e-8);

constexpr bool DETAILED_LINE_ESTIMATORS_ON = false;

constexpr bool DETAILED_BF_ESTIMATORS_ON = true;

constexpr bool LEVEL_HAS_BFEST(int element_z, int ionstage, int level) { return true; }

constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP = 13;

constexpr bool USE_LUT_PHOTOION = false;

constexpr bool USE_ION_BFHEATING_ESTIMATORS = false;

constexpr bool STRICT_POPULATION_CHECKING = false;

constexpr bool NLTE_LIMIT_ION_STAGES_AFTER_FAILURE = false;

constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL = 1000.;

constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING = 10.;

constexpr double NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION = 1e-9;

constexpr bool NLTE_USE_GTH_SOLVER = false;

constexpr std::optional<int> NLTE_TIME_DEPENDENT_FIRST_TIMESTEP = std::nullopt;

constexpr bool NT_ON = true;

constexpr bool NT_SOLVE_SPENCERFANO = true;

constexpr int SF_MAX_TIMESTEPS_BETWEEN_SOLUTIONS = 0;

constexpr double NT_MAX_FRACDIFF_NNEPERION_BETWEEN_SOLUTIONS = 0.05;

constexpr int NTEXCITATION_MAXNLEVELS_LOWER = 5;
constexpr int NTEXCITATION_MAXNLEVELS_UPPER = 250;

constexpr int MAX_NT_EXCITATIONS_STORED = 25000;

constexpr bool NT_EXCITATION_ON = false;

constexpr bool NT_USE_VALENCE_IONPOTENTIAL = false;

constexpr int NT_MAX_AUGER_ELECTRONS = 2;

constexpr bool SF_AUGER_CONTRIBUTION_ON = true;

constexpr bool USE_RELATIVISTIC_DOPPLER_SHIFT = false;

constexpr bool USE_CALCULATED_MEANATOMICWEIGHT = false;

constexpr bool WRITE_EMISSIONABSORPTION_SPEC_AT_END = false;

constexpr bool KEEP_ESCAPED_GAMMAS = true;

constexpr TimeStepSizeMethod TIMESTEP_SIZE_METHOD = TimeStepSizeMethod::LOGARITHMIC;

constexpr double FIXED_TIMESTEP_WIDTH = 0.1;

constexpr double TIMESTEP_TRANSITION_TIME = 5;

constexpr bool KEEP_ALL_RESTART_FILES = false;

constexpr bool BFCOOLING_USELEVELPOPNOTIONPOP = true;

constexpr bool RPKT_USE_EXPANSION_OPACITIES = false;

constexpr bool VPKT_USE_EXPANSION_OPACITIES = false;

constexpr std::optional<float> RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY;

constexpr RpktGreyType RPKT_GREY_TYPE = RpktGreyType::FEGROUP_APPROX;

constexpr bool USE_XCOM_GAMMAPHOTOION = false;

constexpr auto PARTICLE_THERMALISATION_SCHEME = ParticleThermalisationScheme::INSTANTFULLDEPOSITION;

constexpr auto GAMMA_THERMALISATION_SCHEME = GammaThermalisationScheme::FREQUENCYDEPENDENT;

constexpr std::optional<double> GAMMA_USE_KAPPA_GREY;

constexpr bool ENABLE_CHARGE_TRANSFER_REACTIONS = false;

constexpr bool USE_MICROCLUMPING = false;

constexpr float clumping_factor([[maybe_unused]] double tmid, [[maybe_unused]] double rad_vel) { return 1.; }

// NOLINTEND(modernize-use-trailing-return-type,misc-unused-parameters)
#endif  // ARTISOPTIONS_H

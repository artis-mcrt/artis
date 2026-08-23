// A compile-time options preset: symlink or copy one of the artisoptions_*.h files to
// artisoptions.h to select the run configuration. Every option is described in artisoptions_doc.md.

#ifndef ARTISOPTIONS_H
#define ARTISOPTIONS_H
// NOLINTBEGIN(modernize-use-trailing-return-type,misc-unused-parameters)

#include <cstdlib>
#include <optional>

#include "constants.h"

constexpr int MPKTS = 1000000;

constexpr std::optional<GridType> GRID_TYPE_OVERRIDE;
constexpr int CUBOID_NCOORDGRID_X = 50;
constexpr int CUBOID_NCOORDGRID_Y = 50;
constexpr int CUBOID_NCOORDGRID_Z = 50;
constexpr bool FORCE_SPHERICAL_ESCAPE_SURFACE = false;

constexpr int NLTE_TE_NNE_MAXITER = 30;

constexpr bool NLTE_TE_NNE_USE_ANDERSON_ACCEL = true;

constexpr double NLTE_TE_NNE_RELTOL = 0.01;

constexpr int ION_NLEVELS_EXCITED_NLTE(int element_z, int ionstage) { return 0; }

constexpr bool LTEPOP_EXCITATION_USE_TJ = true;

constexpr bool FORCE_SAHA_ION_BALANCE(int element_z) { return true; }

constexpr bool SINGLE_LEVEL_TOP_ION = false;

constexpr int NLEVELS_REQUIRETRANSITIONS(int element_z, int ionstage) {
  return ((element_z == 26 || element_z == 28) && ionstage >= 1) ? 80 : 0;
}

constexpr bool UNIFORM_PELLET_ENERGIES = true;

constexpr bool COL_HEAT_FROM_LEVELPOPS = true;
constexpr bool INITIAL_PACKETS_ON = true;

constexpr int RATECOEFF_TABLESIZE = 200;
constexpr double MINTEMP = 500.;
constexpr double MAXTEMP = 100000.;

constexpr double RECOMBCALIBRATION_T_ELEC = 6000.;

constexpr bool DIPOLE = false;
constexpr bool POL_ON = false;

constexpr bool VPKT_ON = false;
constexpr bool VPKT_WRITE_CONTRIBS = false;

constexpr double MIN_LEVELPOP = 1e-40;

constexpr double NU_MIN_R = 1e13;
constexpr double NU_MAX_R = 5e16;

constexpr bool PHIXS_CLASSIC_NO_INTERPOLATION = false;

constexpr bool MULTIBIN_RADFIELD_MODEL_ON = false;

constexpr int RADFIELDBINCOUNT = 256;

constexpr int FIRST_NLTE_RADFIELD_TIMESTEP = 6;

constexpr double RADFIELDBINS_NU_MIN = (CLIGHT / 40000e-8);
constexpr double RADFIELDBINS_NU_MAX = (CLIGHT / 1085e-8);
constexpr double RADFIELDBINS_T_E_SUPERBIN_NU_MAX = (CLIGHT / 10e-8);

constexpr bool DETAILED_LINE_ESTIMATORS_ON = false;

constexpr bool DETAILED_BF_ESTIMATORS_ON = false;

constexpr bool LEVEL_HAS_BFEST(int element_z, int ionstage, int level) { return true; }

constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP = 6;

constexpr bool USE_LUT_PHOTOION = true;

constexpr bool USE_ION_BFHEATING_ESTIMATORS = true;

constexpr bool STRICT_POPULATION_CHECKING = true;

constexpr bool NLTE_LIMIT_ION_STAGES_AFTER_FAILURE = true;

constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL = 50.;

constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING = 2.;

constexpr double NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION = 1e-9;

constexpr bool NLTE_USE_GTH_SOLVER = false;

constexpr std::optional<int> NLTE_TIME_DEPENDENT_FIRST_TIMESTEP = std::nullopt;

constexpr NonThermalScheme NT_SCHEME = NonThermalScheme::NT_OFF;

constexpr int SF_MAX_TIMESTEPS_BETWEEN_SOLUTIONS = 0;

constexpr double NT_MAX_FRACDIFF_NNEPERION_BETWEEN_SOLUTIONS = 0.05;

constexpr int NTEXCITATION_MAXNLEVELS_LOWER = 5;
constexpr int NTEXCITATION_MAXNLEVELS_UPPER = 250;

constexpr int MAX_NT_EXCITATIONS_STORED = 25000;

constexpr bool NT_USE_VALENCE_IONPOTENTIAL = false;

constexpr int NT_MAX_AUGER_ELECTRONS = 2;

constexpr bool SF_AUGER_CONTRIBUTION_ON = true;

constexpr bool USE_RELATIVISTIC_DOPPLER_SHIFT = true;

constexpr bool USE_CALCULATED_MEANATOMICWEIGHT = true;

constexpr bool KEEP_ESCAPED_GAMMAS = true;

constexpr TimeStepSizeMethod TIMESTEP_SIZE_METHOD = TimeStepSizeMethod::LOGARITHMIC;

constexpr double FIXED_TIMESTEP_WIDTH = -1.;

constexpr double TIMESTEP_TRANSITION_TIME = -1.;

constexpr bool BFCOOLING_USELEVELPOPNOTIONPOP = false;

constexpr bool RPKT_USE_EXPANSION_OPACITIES = false;

constexpr bool VPKT_USE_EXPANSION_OPACITIES = false;

constexpr std::optional<float> RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY;

constexpr RpktGreyType RPKT_GREY_TYPE = RpktGreyType::TANAKA2020_ELECTRONFRAC;

constexpr bool USE_XCOM_GAMMAPHOTOION = false;

constexpr auto PARTICLE_THERMALISATION_SCHEME = ParticleThermalisationScheme::TIMEDEPENDENT;

constexpr auto GAMMA_THERMALISATION_SCHEME = GammaThermalisationScheme::FREQUENCYDEPENDENT;

constexpr std::optional<double> GAMMA_USE_KAPPA_GREY;

constexpr bool ENABLE_CHARGE_TRANSFER_REACTIONS = true;

constexpr bool USE_MICROCLUMPING = false;

constexpr float clumping_factor([[maybe_unused]] double tmid, [[maybe_unused]] double rad_vel) { return 1.; }

// NOLINTEND(modernize-use-trailing-return-type,misc-unused-parameters)
#endif  // ARTISOPTIONS_H

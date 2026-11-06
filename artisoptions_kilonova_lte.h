// The kilonova LTE preset: kilonovae in LTE, as Shingles et al. (2023), ApJL, 954, L41, doi:10.3847/2041-8213/acf29a
// describes the decay and the thermalisation. The preset has the time-dependent thermalisation of the decay
// particles, a grey opacity from a fit to the electron fraction, the relativistic Doppler shift, and a low MINTEMP
// for the late ejecta. The tests kilonova_1d, kilonova_2d, kilonova_2d_expansionopac,
// kilonova_2d_barnesthermalisation, and kilonova_2d_xcomgammaphotoion use this preset.
//
// A compile-time options preset: symlink one of the artisoptions_*.h files to artisoptions.h to select
// the run configuration. A copy through an existing symlink replaces the content of the tracked preset.
// Every option is described in artisoptions_doc.md.

#ifndef ARTISOPTIONS_H
#define ARTISOPTIONS_H
// NOLINTBEGIN(modernize-use-trailing-return-type,misc-unused-parameters)

#include <cstdint>
#include <cstdlib>
#include <optional>

#include "constants.h"

// --- Packets and timesteps ---

constexpr std::int64_t NUM_PACKETS = 2e8;

constexpr bool INITIAL_PACKETS_ON = true;

constexpr bool UNIFORM_PELLET_ENERGIES = true;

constexpr TimeStepSizeMethod TIMESTEP_SIZE_METHOD = TimeStepSizeMethod::LOGARITHMIC;

constexpr std::optional<double> TIMESTEP_SIZE_CONSTANT_MAX_WIDTH_DAYS;

constexpr std::optional<double> TIMESTEP_SIZE_TRANSITION_TIME_DAYS;

// --- Grid and composition ---

constexpr std::optional<GridType> GRID_TYPE_OVERRIDE;

constexpr int CUBOID_NCOORDGRID_X = 50;
constexpr int CUBOID_NCOORDGRID_Y = 50;
constexpr int CUBOID_NCOORDGRID_Z = 50;

constexpr bool FORCE_SPHERICAL_ESCAPE_SURFACE = false;

constexpr bool USE_CALCULATED_MEANATOMICWEIGHT = true;

constexpr bool USE_MICROCLUMPING = false;

constexpr float MICROCLUMPING_FACTOR([[maybe_unused]] double tmid, [[maybe_unused]] double rad_vel) { return 1.; }

// --- Atomic data and limits ---

constexpr bool SINGLE_LEVEL_TOP_ION = false;

constexpr bool PHIXS_CLASSIC_NO_INTERPOLATION = false;

constexpr int RATECOEFF_TABLESIZE = 200;
constexpr double MINTEMP = 500.;
constexpr double MAXTEMP = 150000.;

constexpr double RECOMBCALIBRATION_T_ELEC = 6000.;

constexpr double MINPOP = 1e-40;

// --- Ionisation and excitation ---

constexpr int ION_NLEVELS_EXCITED_NLTE(int element_z, int ionstage) { return 0; }

constexpr int NLEVELS_ADD_MISSING_COLLISIONAL_TRANSITIONS(int element_z, int ionstage) {
  return ((element_z == 26 || element_z == 28) && ionstage >= 1) ? 80 : 0;
}

constexpr bool LTEPOP_EXCITATION_USE_TJ = false;

constexpr bool FORCE_SAHA_ION_BALANCE(int element_z) { return element_z < 38 || element_z > 40; }

constexpr int NLTE_TE_NNE_MAXITER = 30;

constexpr bool NLTE_TE_NNE_USE_ANDERSON_ACCEL = true;

constexpr double NLTE_TE_NNE_RELTOL = 0.01;

constexpr bool NLTE_USE_GTH_SOLVER = true;

constexpr std::optional<int> NLTE_TIME_DEPENDENT_FIRST_TIMESTEP = std::nullopt;

constexpr bool STRICT_POPULATION_CHECKING = true;

constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL = 50.;

constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING = 2.;

constexpr std::optional<double> NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION = 1e-9;

constexpr bool ENABLE_CHARGE_TRANSFER_REACTIONS = false;

// --- Heating and cooling ---

constexpr bool COLLISIONAL_HEATING_FROM_LEVELPOPS = true;

constexpr bool BFCOOLING_USELEVELPOPNOTIONPOP = false;

constexpr bool USE_ION_BFHEATING_ESTIMATORS = false;

constexpr bool WRITE_ION_HEATING_COOLING_RATES = false;

// --- Radiation field and estimators ---

constexpr bool USE_LUT_PHOTOION = false;

constexpr std::optional<int> MULTIBIN_RADFIELD_FIRST_TIMESTEP;

constexpr int RADFIELDBINCOUNT = 256;

constexpr double RADFIELDBINS_NU_MIN = (CLIGHT / 40000e-8);
constexpr double RADFIELDBINS_NU_MAX = (CLIGHT / 1085e-8);

constexpr double RADFIELDBINS_T_E_SUPERBIN_NU_MAX = (CLIGHT / 10e-8);

constexpr bool DETAILED_LINE_ESTIMATORS_ON = false;

constexpr std::optional<int> DETAILED_BF_ESTIMATORS_FIRST_TIMESTEP;

constexpr bool DETAILED_BF_ESTIMATORS_LEVEL_SELECTED(int element_z, int ionstage, int level) { return false; }

// --- Non-thermal particles ---

constexpr NonThermalScheme NT_SCHEME = NonThermalScheme::NT_OFF;

constexpr double SF_RE_SOLVE_WITHIN_TIMESTEP_MIN_NNEPERION_FRACDIFF = 0.05;

constexpr int NTEXCITATION_MAXNLEVELS_LOWER = 5;
constexpr int NTEXCITATION_MAXNLEVELS_UPPER = 250;

constexpr int MAX_NT_EXCITATIONS_STORED = 25000;

constexpr bool NT_USE_VALENCE_IONPOTENTIAL = false;

constexpr int NT_MAX_AUGER_ELECTRONS = 2;

constexpr bool SF_AUGER_CONTRIBUTION_ON = true;

constexpr ParticleThermalisationScheme PARTICLE_THERMALISATION_SCHEME =
    ParticleThermalisationScheme::TIMEDEPENDENT_WITH_ADIABATIC_LOSS_ULTRARELATIVISTICLIMIT;

// --- Gamma rays ---

constexpr GammaThermalisationScheme GAMMA_THERMALISATION_SCHEME = GammaThermalisationScheme::FREQUENCYDEPENDENT;

constexpr bool GAMMA_USE_XCOM_PHOTOION = false;

constexpr std::optional<double> GAMMA_USE_KAPPA_GREY;

constexpr bool KEEP_ESCAPED_GAMMAS = false;

// --- Real packets ---

constexpr double NU_MIN_R = 1e13;
constexpr double NU_MAX_R = 5e16;

constexpr bool DIPOLE = false;

constexpr bool POL_ON = false;

constexpr bool USE_RELATIVISTIC_DOPPLER_SHIFT = true;

constexpr std::optional<float> RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY;

constexpr RpktGreyType RPKT_GREY_TYPE = RpktGreyType::TANAKA2020_ELECTRONFRAC;

constexpr bool SAMPLE_RPKT_EMISSION = false;

// --- Virtual packets ---

constexpr bool VPKT_ON = false;

constexpr bool VPKT_WRITE_CONTRIBS = false;

// --- Binned opacities ---

constexpr bool RPKT_USE_BINNED_OPACITIES = false;

constexpr bool VPKT_USE_BINNED_OPACITIES = false;

constexpr BinnedOpacityLineWeight BINNED_OPACITY_LINE_WEIGHT = BinnedOpacityLineWeight::EXPANSION;

constexpr double LINE_BINNED_CAPPED_TAU_MAX = 1.;

// --- Output files ---

constexpr bool WRITE_COMBINED_ALLRANK_OUT_FILES = true;

// NOLINTEND(modernize-use-trailing-return-type,misc-unused-parameters)
#endif  // ARTISOPTIONS_H

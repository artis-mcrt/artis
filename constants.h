#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <array>
#include <cstddef>
#include <numbers>
#include <string_view>
#include <vector>

// fundamental constants

constexpr double CLIGHT = 2.99792458e+10;  // Speed of light [cm/s]
constexpr double CLIGHT_PROP = CLIGHT;  // Speed of light for ray travel. Physically = CLIGHT but
                                        // can be changed for testing.
constexpr double H = 6.6260755e-27;  // Planck constant [erg s]
constexpr double MSUN = 1.98855e+33;  // Solar mass [g]
constexpr double LSUN = 3.826e+33;  // Solar luminosity [erg/s]
constexpr double MH = 1.67352e-24;  // Mass of hydrogen atom [g]
constexpr double ME = 9.1093897e-28;  // Mass of free electron [g]
constexpr double QE = 4.80325E-10;  // elementary charge in cgs units [statcoulomb]
constexpr double PI = std::numbers::pi;
constexpr double EV = 1.6021772e-12;  // eV to ergs [eV/erg]
constexpr double MEV = 1.6021772e-6;  // MeV to ergs [MeV/erg]
constexpr double DAY = 86400.;  // day to seconds [s/day]
constexpr double SIGMA_T = 6.6524e-25;  // Thomson cross-section [cm2]
constexpr double THOMSON_LIMIT = 1e-2;  // Limit below which e-scattering is Thomson
constexpr double PARSEC = 3.0857e+18;  // pc to cm [cm/pc]
constexpr double KB = 1.38064852e-16;  // Boltzmann constant [erg/K]
constexpr double STEBO = 5.670400e-5;  // Stefan-Boltzmann constant [erg cm^-2 s^-1 K^-4.]
                                       // (data from NIST http://physics.nist.gov/cgi-bin/cuu/Value?eqsigma)
constexpr double SAHACONST = 2.0706659e-16;  // Saha constant
constexpr double BARN = 1e-24;  // barn in cgs units [cm2]

constexpr double EULERGAMMA = std::numbers::egamma;
constexpr double PAIR_PROD_FREQUENCY =
    2.46636e+20;  // frequency of wave having the energy corresponding to the rest energy of an electron-positron pair

// numerical constants

constexpr double CLIGHTSQUARED = 8.9875518e+20;  // Speed of light squared [cm^2/s^2]
constexpr double TWOOVERCLIGHTSQUARED = 2.2253001e-21;
constexpr double TWOHOVERCLIGHTSQUARED = 1.4745007e-47;
constexpr double CLIGHTSQUAREDOVERTWOH = 6.7819570e+46;

constexpr double ONEOVERH = 1.509188961e+26;
constexpr double HOVERKB = 4.799243681748932e-11;
constexpr double FOURPI = 1.256637061600000e+01;
constexpr double ONEOVER4PI = 7.957747153555701e-02;
constexpr double HCLIGHTOVERFOURPI = 1.580764662876770e-17;
constexpr double OSCSTRENGTHCONVERSION = 1.3473837e+21;

constexpr double H_ionpot = 13.5979996 * EV;

enum class GridType {
  SPHERICAL1D,  // 1D radial shells (non-uniform dr)
  CYLINDRICAL2D,  // 2D cylindrical grid with uniform dz, drcyl
  CARTESIAN3D  // 3D Cartesian cubic grid with uniform dx=dy=dz
};

// constant for van-Regemorter approximation.
constexpr double C_0 = 5.465e-11;

constexpr size_t GSLWSIZE = 16384;  // GSL integration workspace size

enum class TimeStepSizeMethod { LOGARITHMIC, CONSTANT, LOGARITHMIC_THEN_CONSTANT, CONSTANT_THEN_LOGARITHMIC };

enum class ThermalisationScheme {
  INSTANT,
  DETAILED,
  DETAILEDWITHGAMMAPRODUCTS,
  BARNES,
  WOLLAEGER,
  GUTTMAN,
  TOT_THERM_FIT_BARNES_EQ34,
  TOT_THERM_FIT_BARNES_EQ32p33,
  TOT_THERM_FIT_MRW
};

// neutrino fractions for the e2e model, ONLY applicable for 50 timesteps from
// 0.1 to 20 days. TODO: improve
const std::vector<double> nu_frac_e2e = {
    0.37343442, 0.37479627, 0.37677179, 0.37912375, 0.38102108, 0.38384897, 0.38783418, 0.39096066, 0.39535378,
    0.39773747, 0.4006056,  0.40353313, 0.40539669, 0.40861926, 0.40818127, 0.40906144, 0.40954452, 0.40881866,
    0.40999089, 0.40937853, 0.41115956, 0.41162955, 0.41207718, 0.41291792, 0.41246466, 0.41429677, 0.41627379,
    0.41932937, 0.41619559, 0.41700607, 0.41727353, 0.41419981, 0.41302095, 0.40955967, 0.40755307, 0.40660224,
    0.4047616,  0.40244544, 0.4011143,  0.4013594,  0.40162403, 0.4025424,  0.40364484, 0.40789019, 0.4127288,
    0.41533031, 0.41929123, 0.42236189, 0.42632201, 0.43183402};

using Vec3d = std::array<double, 3>;
constexpr Vec3d syn_dir{0., 0., 1.};  // vector defining the theta=0 direction

constexpr std::string_view outdir_resfiles{"speclc_angle_res/"};

#endif

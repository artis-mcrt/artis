#ifndef CONSTANTS_H
#define CONSTANTS_H

/// fundamental constants
constexpr double CLIGHT = 2.99792458e+10;  /// Speed of light [cm/s]
constexpr double CLIGHT_PROP = CLIGHT;     // Speed of light for ray travel. Physically = CLIGHT but
                                           // can be changed for testing.
constexpr double H = 6.6260755e-27;        /// Planck constant [erg s]
constexpr double MSUN = 1.98855e+33;       /// Solar mass [g]
constexpr double LSUN = 3.826e+33;         /// Solar luminosity [erg/s]
constexpr double MH = 1.67352e-24;         /// Mass of hydrogen atom [g]
constexpr double ME = 9.1093897e-28;       /// Mass of free electron [g]
constexpr double QE = 4.80325E-10;         /// elementary charge in cgs units [statcoulomb]
constexpr double PI = 3.1415926535987;
constexpr double EV = 1.6021772e-12;         /// eV to ergs [eV/erg]
constexpr double MEV = 1.6021772e-6;         /// MeV to ergs [MeV/erg]
constexpr double DAY = 86400.0;              /// day to seconds [s/day]
constexpr double SIGMA_T = 6.6524e-25;       /// Thomson cross-section
constexpr double THOMSON_LIMIT = 1e-2;       /// Limit below which e-scattering is Thomson
constexpr double PARSEC = 3.0857e+18;        /// pc to cm [cm/pc]
constexpr double KB = 1.38064852e-16;        /// Boltzmann constant [erg/K]
constexpr double STEBO = 5.670400e-5;        /// Stefan-Boltzmann constant [erg cm^−2 s^−1 K^−4.]
                                             /// (data from NIST http://physics.nist.gov/cgi-bin/cuu/Value?eqsigma)
constexpr double SAHACONST = 2.0706659e-16;  /// Saha constant

/// numerical constants
constexpr double CLIGHTSQUARED = 8.9875518e+20;  /// Speed of light squared [cm^2/s^2]
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

constexpr int GRID_UNIFORM = 1;      // Simple cuboidal cells.
constexpr int GRID_SPHERICAL1D = 2;  // radial shells

// constant for van-Regemorter approximation.
constexpr double C_0 = 5.465e-11;

constexpr int MAXFILENAMELENGTH = 128;
constexpr size_t GSLWSIZE = 16384;  // GSL integration workspace size

enum timestepsizemethods {
  TIMESTEP_SIZES_LOGARITHMIC = 0,
  TIMESTEP_SIZES_CONSTANT = 1,
  TIMESTEP_SIZES_LOGARITHMIC_THEN_CONSTANT = 2,
  TIMESTEP_SIZES_CONSTANT_THEN_LOGARITHMIC = 3,
};

#endif
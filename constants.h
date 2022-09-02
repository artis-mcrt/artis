#ifndef CONSTANTS_H
#define CONSTANTS_H

/// fundamental constants
constexpr double CLIGHT = 2.99792458e+10;  /// Speed of light [cm/s]
#define H 6.6260755e-27                    /// Planck constant [erg s]
#define MSUN 1.98855e+33                   /// Solar mass [g]
#define LSUN 3.826e+33                     /// Solar luminosity [erg/s]
#define MH 1.67352e-24                     /// Mass of hydrogen atom [g]
#define ME 9.1093897e-28                   /// Mass of free electron [g]
#define QE 4.80325E-10                     /// elementary charge in cgs units [statcoulomb]
#define PI 3.1415926535987
#define EV 1.6021772e-12    /// eV to ergs [eV/erg]
#define MEV 1.6021772e-6    /// MeV to ergs [MeV/erg]
#define DAY 86400.0         /// day to seconds [s/day]
#define SIGMA_T 6.6524e-25  /// Thomson cross-section
#define THOMSON_LIMIT 1e-2  /// Limit below which e-scattering is Thomson
#define PARSEC 3.0857e+18   /// pc to cm [cm/pc]
#define KB 1.38064852e-16   /// Boltzmann constant [erg/K]
#define STEBO \
  5.670400e-5                    /// Stefan-Boltzmann constant [erg cm^−2 s^−1 K^−4.]
                                 /// (data from NIST http://physics.nist.gov/cgi-bin/cuu/Value?eqsigma)
#define SAHACONST 2.0706659e-16  /// Saha constant

/// numerical constants
#define CLIGHTSQUARED 8.9875518e+20  /// Speed of light squared [cm^2/s^2]
#define TWOOVERCLIGHTSQUARED 2.2253001e-21
#define TWOHOVERCLIGHTSQUARED 1.4745007e-47
#define CLIGHTSQUAREDOVERTWOH 6.7819570e+46

#define ONEOVERH 1.509188961e+26
#define HOVERKB 4.799243681748932e-11
#define FOURPI 1.256637061600000e+01
#define ONEOVER4PI 7.957747153555701e-02
#define HCLIGHTOVERFOURPI 1.580764662876770e-17
#define OSCSTRENGTHCONVERSION 1.3473837e+21

#define H_ionpot (13.5979996 * EV)

#define GRID_UNIFORM 1      // Simple cuboidal cells.
#define GRID_SPHERICAL1D 2  // radial shells

#endif
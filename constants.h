#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cstddef>
#include <numbers>

// fundamental constants
inline constexpr double CLIGHT = 2.99792458e+10;  // Speed of light [cm/s]
inline constexpr double CLIGHT_PROP = CLIGHT;     // Speed of light for ray travel. Physically = CLIGHT but
                                                  // can be changed for testing.
inline constexpr double H = 6.6260755e-27;        // Planck constant [erg s]
inline constexpr double MSUN = 1.98855e+33;       // Solar mass [g]
inline constexpr double LSUN = 3.826e+33;         // Solar luminosity [erg/s]
inline constexpr double MH = 1.67352e-24;         // Mass of hydrogen atom [g]
inline constexpr double ME = 9.1093897e-28;       // Mass of free electron [g]
inline constexpr double QE = 4.80325E-10;         // elementary charge in cgs units [statcoulomb]
inline constexpr double PI = std::numbers::pi;
inline constexpr double EV = 1.6021772e-12;    // eV to ergs [eV/erg]
inline constexpr double MEV = 1.6021772e-6;    // MeV to ergs [MeV/erg]
inline constexpr double DAY = 86400.;          // day to seconds [s/day]
inline constexpr double SIGMA_T = 6.6524e-25;  // Thomson cross-section [cm2]
inline constexpr double THOMSON_LIMIT = 1e-2;  // Limit below which e-scattering is Thomson
inline constexpr double PARSEC = 3.0857e+18;   // pc to cm [cm/pc]
inline constexpr double KB = 1.38064852e-16;   // Boltzmann constant [erg/K]
inline constexpr double STEBO = 5.670400e-5;   // Stefan-Boltzmann constant [erg cm^−2 s^−1 K^−4.]
                                               // (data from NIST http://physics.nist.gov/cgi-bin/cuu/Value?eqsigma)
inline constexpr double SAHACONST = 2.0706659e-16;  // Saha constant

// numerical constants
inline constexpr double CLIGHTSQUARED = 8.9875518e+20;  // Speed of light squared [cm^2/s^2]
inline constexpr double TWOOVERCLIGHTSQUARED = 2.2253001e-21;
inline constexpr double TWOHOVERCLIGHTSQUARED = 1.4745007e-47;
inline constexpr double CLIGHTSQUAREDOVERTWOH = 6.7819570e+46;

inline constexpr double ONEOVERH = 1.509188961e+26;
inline constexpr double HOVERKB = 4.799243681748932e-11;
inline constexpr double FOURPI = 1.256637061600000e+01;
inline constexpr double ONEOVER4PI = 7.957747153555701e-02;
inline constexpr double HCLIGHTOVERFOURPI = 1.580764662876770e-17;
inline constexpr double OSCSTRENGTHCONVERSION = 1.3473837e+21;

inline constexpr double H_ionpot = 13.5979996 * EV;

enum gridtypes {
  GRID_SPHERICAL1D = 1,    // 1D radial shells (non-uniform dr)
  GRID_CYLINDRICAL2D = 2,  // 2D cylindrical grid with uniform dz, drcyl
  GRID_CARTESIAN3D = 3     // 3D Cartesian cubic grid with uniform dx=dy=dz
};

inline constexpr int GRID_UNIFORM = GRID_CARTESIAN3D;  // deprecated alias for GRID_CARTESIAN3D

// constant for van-Regemorter approximation.
inline constexpr double C_0 = 5.465e-11;

inline constexpr int MAXFILENAMELENGTH = 128;
inline constexpr size_t GSLWSIZE = 16384;  // GSL integration workspace size

enum timestepsizemethods {
  TIMESTEP_SIZES_LOGARITHMIC = 0,
  TIMESTEP_SIZES_CONSTANT = 1,
  TIMESTEP_SIZES_LOGARITHMIC_THEN_CONSTANT = 2,
  TIMESTEP_SIZES_CONSTANT_THEN_LOGARITHMIC = 3,
};

#endif
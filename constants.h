#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <array>
#include <cstddef>
#include <numbers>
#include <string_view>

[[gnu::const]] constexpr auto pow2(auto x) { return x * x; }
[[gnu::const]] constexpr auto pow3(auto x) { return x * x * x; }
[[gnu::const]] constexpr auto pow4(auto x) { return pow2(x) * pow2(x); }
[[gnu::const]] constexpr auto pow5(auto x) { return pow4(x) * x; }

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

constexpr double EULERGAMMA = std::numbers::egamma;

// numerical constants

constexpr double CLIGHTSQUARED = pow2(CLIGHT);  // Speed of light squared [cm^2/s^2]
constexpr double TWOOVERCLIGHTSQUARED = 2 / CLIGHTSQUARED;
constexpr double TWOHOVERCLIGHTSQUARED = 2 * H / CLIGHTSQUARED;
constexpr double CLIGHTSQUAREDOVERTWOH = CLIGHTSQUARED / (2 * H);

constexpr double ONEOVERH = 1.0 / H;
constexpr double HOVERKB = H / KB;
constexpr double FOURPI = 4 * PI;
constexpr double ONEOVER4PI = 1 / (4 * PI);
constexpr double HCLIGHTOVERFOURPI = H * CLIGHT / FOURPI;

constexpr double H_ionpot = 13.5979996 * EV;

enum class GridType {
  SPHERICAL1D,  // 1D radial shells (non-uniform dr)
  CYLINDRICAL2D,  // 2D cylindrical grid with uniform dz, drcyl
  CARTESIAN3D  // 3D Cartesian cubic grid with uniform dx=dy=dz
};

// constant for van-Regemorter approximation.
constexpr double C_0 = 5.465e-11;

enum class TimeStepSizeMethod { LOGARITHMIC, CONSTANT, LOGARITHMIC_THEN_CONSTANT, CONSTANT_THEN_LOGARITHMIC };

enum class GammaThermalisationScheme { FREQUENCYDEPENDENT, BARNES, WOLLAEGER, GUTTMAN };
enum class ParticleThermalisationScheme {
  INSTANTFULLDEPOSITION,
  TIMEDEPENDENT,
  TIMEDEPENDENTWITHGAMMAPRODUCTS,
  BARNES,
  WOLLAEGER,
};

enum class RpktGreyType { FEGROUP_APPROX, TANAKA2020_ELECTRONFRAC, JUST2022_TEMP_LANTHANIDEFRAC };

using Vec3d = std::array<double, 3>;
constexpr Vec3d syn_dir{0., 0., 1.};  // vector defining the theta=0 direction

constexpr std::string_view outdir_resfiles{"speclc_angle_res/"};

constexpr std::array datafolders{"./", "data/", "artis/data/"};

#ifndef TESTMODE
#define TESTMODE false
#endif

#if defined(__NVCOMPILER_CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
#define THREADLOCALONHOST
#else
#define THREADLOCALONHOST thread_local static
#endif

#ifdef GPU_ON
#define DEVICE_FUNC __host__ __device__
#else
#define DEVICE_FUNC
#endif

#if defined REPRODUCIBLE && REPRODUCIBLE
#define SORT_OR_STABLE_SORT stable_sort
#else
#define SORT_OR_STABLE_SORT sort
#endif

#ifdef STDPAR_ON
#include <execution>
#define EXEC_PAR_UNSEQ std::execution::par_unseq,
#define EXEC_PAR std::execution::par,
#else
#define EXEC_PAR_UNSEQ
#define EXEC_PAR
#endif

#ifdef _OPENMP
#define atomicadd(var, val)                  \
  {                                          \
    _Pragma("omp atomic update") var += val; \
  }

#elifdef STDPAR_ON
#include <atomic>

template <typename T, typename U>
constexpr void atomicadd(T& var, U&& val) {
  std::atomic_ref<T>(var).fetch_add(std::forward<U>(val), std::memory_order_relaxed);
}

#else

#define atomicadd(var, val) var += (val);

#endif

#endif

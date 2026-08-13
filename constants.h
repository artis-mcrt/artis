// Physical constants (CGS units), unit conversions, and low-level helpers (atomicadd,
// cache-line alignment) used throughout the code.

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <new>
#include <numbers>
#include <string_view>
#include <version>

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
constexpr double EV = 1.6021772e-12;  // eV to ergs [erg/eV]
constexpr double MEV = 1.6021772e-6;  // MeV to ergs [erg/MeV]
constexpr double DAY = 86400.;  // day to seconds [s/day]
constexpr double HOUR = 3600.;  // hour to seconds [s/hour]
constexpr double SIGMA_T = 6.6524e-25;  // Thomson cross-section [cm2]
constexpr double THOMSON_LIMIT = 1e-2;  // Limit below which e-scattering is Thomson
constexpr double PARSEC = 3.0857e+18;  // pc to cm [cm/pc]
constexpr double KB = 1.38064852e-16;  // Boltzmann constant [erg/K]
constexpr double STEBO = 5.670400e-5;  // Stefan-Boltzmann constant [erg cm^-2 s^-1 K^-4.]
                                       // (data from NIST http://physics.nist.gov/cgi-bin/cuu/Value?eqsigma)
constexpr double SAHACONST = 2.0706659e-16;  // Saha constant

constexpr double EULERGAMMA = std::numbers::egamma;

// exp(EULERGAMMA) = 1.781072... This is the gamma of the Coulomb logarithm
// ln[m v^3 / (gamma e^2 omega_p)] in Schunk & Hays (1971) eq. 4 and Kozma & Fransson (1992) eq. 2.
// Schunk & Hays (p. 114) define it by "ln gamma is Euler's constant", i.e. gamma = exp(EULERGAMMA),
// not EULERGAMMA itself.
// Not calculated with std::exp() because that is not constexpr in libc++.
constexpr double EXP_EULERGAMMA = 1.7810724179901979852;

// numerical constants

constexpr double CLIGHTSQUARED = pow2(CLIGHT);  // Speed of light squared [cm^2/s^2]
constexpr double CLIGHTSQUAREDOVERTWOH = pow2(CLIGHT) / (2 * H);

constexpr double HOVERKB = H / KB;
constexpr double FOURPI = 4 * PI;
constexpr double HCLIGHTOVERFOURPI = H * CLIGHT / (4 * PI);

constexpr double H_ionpot = 13.5979996 * EV;

enum class GridType {
  SPHERICAL1D,  // 1D radial shells (non-uniform dr)
  CYLINDRICAL2D,  // 2D cylindrical grid with uniform dz, drcyl
  CARTESIAN3D,  // 3D Cartesian cubic grid with uniform dx=dy=dz
};

// constant for van-Regemorter approximation.
constexpr double C_0 = 5.465e-11;

enum class TimeStepSizeMethod { LOGARITHMIC, CONSTANT, LOGARITHMIC_THEN_CONSTANT, CONSTANT_THEN_LOGARITHMIC };

enum class GammaThermalisationScheme { FREQUENCYDEPENDENT, BARNES, WOLLAEGER, GUTTMAN };
enum class ParticleThermalisationScheme {
  INSTANTFULLDEPOSITION,
  TIMEDEPENDENT,
  TIMEDEPENDENT_WITH_ADIABATIC_LOSS,
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

#ifdef GPU_ON
constexpr bool cellcache_singleslot = false;
#else
constexpr bool cellcache_singleslot = true;
#endif

#ifdef __HIPCC__
#define THREADLOCALONHOST
#else
#define THREADLOCALONHOST thread_local
#endif

#ifdef __NVCOMPILER
// nvc++ : single-pass. Each macro becomes an independent "if target".
#include <nv/target>
#define MY_IF_DEVICE(...) NV_IF_TARGET(NV_IS_DEVICE, (__VA_ARGS__))
#define MY_IF_HOST(...) NV_IF_TARGET(NV_IS_HOST, (__VA_ARGS__))
#elif (defined(__HIP_DEVICE_COMPILE__) && __HIP_DEVICE_COMPILE__) || defined(__CUDA_ARCH__)
// Device pass of a multi-pass compiler (HIP-Clang / nvcc).
#define MY_IF_DEVICE(...) __VA_ARGS__
#define MY_IF_HOST(...)
#else
// Host pass of HIP/nvcc, OR a plain CPU-only compiler (gcc/clang/icc).
#define MY_IF_DEVICE(...)
#define MY_IF_HOST(...) __VA_ARGS__
#endif

#if defined(__device__) && defined(GPU_ON)
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

#ifdef __HIPCC__
#define EXEC_PAR_UNSEQ std::execution::par,
#else
#define EXEC_PAR_UNSEQ std::execution::par_unseq,
#endif
#define EXEC_PAR std::execution::par,
#else
#define EXEC_PAR_UNSEQ
#define EXEC_PAR
#endif

// Index of the lowest set bit. Precondition: bits != 0.
// std::countr_zero() cannot be translated into an nvc++ GPU accelerator region: libstdc++ implements it
// with __builtin_ctzll, and nvc++ 26.1 with the gcc-14 headers rejects that with
//   NVC++-F-0155-Compiler failed to translate accelerator region: Unsupported operation (bit: 418)
// Only this function is affected -- random.h calls std::rotl from the same header in device code -- so the
// fallback is limited to the builds that create accelerator regions. It is a binary search rather than a
// linear scan to stay cheap on that path. Note the two differ on zero (std:: gives 64, this gives 63),
// hence the precondition.
[[nodiscard]] DEVICE_FUNC constexpr auto lowest_set_bit(const std::uint64_t bits) -> int {
#if defined(__NVCOMPILER) && defined(GPU_ON)
  auto remaining = bits;
  int index = 0;
  for (auto width = 32U; width > 0U; width /= 2U) {
    if ((remaining & ((UINT64_C(1) << width) - 1U)) == 0U) {
      remaining >>= width;
      index += static_cast<int>(width);
    }
  }
  return index;
#else
  return std::countr_zero(bits);
#endif
}

// The nvc++ GPU job is the only build that instantiates the fallback, and it only compiles -- no test ever
// runs it -- so these are its whole coverage. They step through each halving width of the binary search.
static_assert(lowest_set_bit(1) == 0);
static_assert(lowest_set_bit(2) == 1);
static_assert(lowest_set_bit(0b1100) == 2);
static_assert(lowest_set_bit(~UINT64_C(0)) == 0);
static_assert(lowest_set_bit(UINT64_C(1) << 1) == 1);
static_assert(lowest_set_bit(UINT64_C(1) << 2) == 2);
static_assert(lowest_set_bit(UINT64_C(1) << 4) == 4);
static_assert(lowest_set_bit(UINT64_C(1) << 8) == 8);
static_assert(lowest_set_bit(UINT64_C(1) << 16) == 16);
static_assert(lowest_set_bit(UINT64_C(1) << 31) == 31);
static_assert(lowest_set_bit(UINT64_C(1) << 32) == 32);
static_assert(lowest_set_bit(UINT64_C(1) << 33) == 33);
static_assert(lowest_set_bit(UINT64_C(1) << 62) == 62);
static_assert(lowest_set_bit(UINT64_C(1) << 63) == 63);
static_assert(lowest_set_bit(UINT64_C(3) << 40) == 40);
static_assert(lowest_set_bit(~UINT64_C(0) << 17) == 17);

#if defined(__cpp_lib_hardware_interference_size) && __cpp_lib_hardware_interference_size >= 201603
constexpr std::size_t DESTRUCTIVE_INTERFERENCE_SIZE = std::hardware_destructive_interference_size;
#else
constexpr std::size_t DESTRUCTIVE_INTERFERENCE_SIZE = 64;
#endif

// Give a hot thread-shared accumulator variable its own cache line to prevent false sharing in the CPU multithreaded
// modes (OpenMP or C++ standard parallelism on multicore). Expands to nothing in single-threaded (MPI-only) and GPU
// modes, where the cache-line contention between host threads doesn't exist and padding would only waste memory.
#if (defined(_OPENMP) || defined(STDPAR_ON)) && !defined(GPU_ON)
constexpr bool align_to_avoid_false_sharing = true;
#define ALIGNAS_AVOID_FALSE_SHARING alignas(DESTRUCTIVE_INTERFERENCE_SIZE)
#else
constexpr bool align_to_avoid_false_sharing = false;
#define ALIGNAS_AVOID_FALSE_SHARING
#endif

#ifdef _OPENMP
#define atomicadd(var, val)                  \
  {                                          \
    _Pragma("omp atomic update") var += val; \
  }

// relaxed atomic read/write for values that concurrent workers may lazily fill with identical
// values (e.g. the continuum caches in the cell cache): the atomicity makes the unsynchronised
// same-value stores well-defined without imposing any ordering
template <typename T>
[[nodiscard]] auto atomicload(T& var) -> T {
  T val{};
#pragma omp atomic read
  val = var;
  return val;
}

template <typename T, typename U>
void atomicstore(T& var, const U val) {
#pragma omp atomic write
  var = val;
}

#elifdef STDPAR_ON
#include <atomic>

template <typename T, typename U>
constexpr void atomicadd(T& var, U&& val) {
  std::atomic_ref<T>(var).fetch_add(std::forward<U>(val), std::memory_order_relaxed);
}

// relaxed atomic read/write for values that concurrent workers may lazily fill with identical
// values (e.g. the continuum caches in the cell cache): the atomicity makes the unsynchronised
// same-value stores well-defined without imposing any ordering
template <typename T>
[[nodiscard]] constexpr auto atomicload(T& var) -> T {
  return std::atomic_ref<T>(var).load(std::memory_order_relaxed);
}

template <typename T, typename U>
constexpr void atomicstore(T& var, const U val) {
  std::atomic_ref<T>(var).store(val, std::memory_order_relaxed);
}

#else

#define atomicadd(var, val) var += (val);

template <typename T>
[[nodiscard]] constexpr auto atomicload(T& var) -> T {
  return var;
}

template <typename T, typename U>
constexpr void atomicstore(T& var, const U val) {
  var = val;
}

#endif

#endif

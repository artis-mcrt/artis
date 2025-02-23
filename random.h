// Cut-down version of MIT-licensed UTL random module Copyright (c) 2023 Dmitri Bogdanov
// https://github.com/DmitriBogdanov/UTL/blob/master/docs/module_random.md

#ifndef RANDOM_H
#define RANDOM_H

#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace utlrandom {

// Helper method to crush large uints to uint32_t,
// inspired by Melissa E. O'Neil's randutils https://gist.github.com/imneme/540829265469e673d045
template <class T>
  requires(std::is_integral_v<T> && sizeof(T) <= 8)
[[nodiscard]] constexpr auto _crush_to_uint32(T value) -> std::uint32_t {
  if constexpr (sizeof(value) <= 4) {
    return std::uint32_t(value);
  } else {
    std::uint64_t res = value;
    res *= 0xbc2ad017d719504d;
    return static_cast<std::uint32_t>(res ^ (res >> 32));
  }
}

template <class ResultType>
[[nodiscard]] constexpr auto _mix_seed(ResultType seed) -> ResultType {
  std::uint64_t state = (static_cast<std::uint64_t>(seed) + 0x9E3779B97f4A7C15);
  state = (state ^ (state >> 30)) * 0xBF58476D1CE4E5B9;
  state = (state ^ (state >> 27)) * 0x94D049BB133111EB;
  return static_cast<ResultType>(state ^ (state >> 31));
  // some of the 16/32-bit PRNGs have bad correlation on the successive seeds, this usually
  // can be alleviated by using a signle iteration of a "good" PRNG to pre-mix the seed
}

template <class T>
constexpr T _default_seed = (std::numeric_limits<T>::max() / 2) + 1;
// an "overall decent" default seed - doesn't gave too many zeroes,
// unlikely to accidentaly match with a user-defined seed

namespace generators {

// Implementation of 32-bit splitmix adopted from MurmurHash3, based on paper by Guy L. Steele,
// Doug Lea, and Christine H. Flood. 2014. "Fast splittable pseudorandom number generators"
// see http://marc-b-reynolds.github.io/shf/2017/09/27/LPRNS.html
//     https://gee.cs.oswego.edu/dl/papers/oopsla14.pdf
//     https://github.com/umireon/my-random-stuff/blob/e7b17f992955f4dbb02d4016682113b48b2f6ec1/xorshift/splitmix32.c
//
// Performance: Excellent
// Quality:     3/5
// State:       4 bytes
//
// One of the fastest 32-bit generators that requires only a single 'std::uint32_t' of state,
// making it the smallest state available. Some other PRNGs recommend using it for seeding their state.
// 32-bit version is somewhat lacking in terms of quality estimate data (relative to the widely used
// 64-bit version), however it still seems to be quite decent.
//
class SplitMix32 {
 public:
  using result_type = std::uint64_t;

 private:
  result_type s{};

 public:
  constexpr explicit SplitMix32(result_type seed = _default_seed<result_type>) noexcept { this->seed(seed); }

  [[nodiscard]] static constexpr auto min() noexcept -> result_type { return 0; }
  [[nodiscard]] static constexpr auto max() noexcept -> result_type { return std::numeric_limits<result_type>::max(); }

  constexpr void seed(result_type seed) noexcept {
    this->s = _mix_seed(seed);
    // naively seeded SplitMix32 has a horrible correlation between successive seeds, we can mostly alleviate
    // the issue by pre-mixing the seed with a single iteration of a "better" 64-bit algorithm
  }

  constexpr auto operator()() noexcept -> result_type {
    result_type result = (this->s += 0x9e3779b9);
    result = (result ^ (result >> 16)) * 0x21f0aaad;
    result = (result ^ (result >> 15)) * 0x735a2d97;
    return result ^ (result >> 15);
  }
};

// Implementation of Xoshiro128++ suggested by David Blackman and Sebastiano Vigna,
// see https://prng.di.unimi.it/
//     https://prng.di.unimi.it/xoshiro256plusplus.c
//
// Performance: Good
// Quality:     4/5
// State:       16 bytes
//
// Excellent choice as a general purpose 32-bit PRNG.
// Battle-tested and provides a good statistical quality at an excellent speed.
//
class Xoshiro128PP {
 public:
  using result_type = std::uint32_t;

 private:
  std::array<result_type, 4> s{};

 public:
  constexpr explicit Xoshiro128PP(result_type seed = _default_seed<result_type>) noexcept { this->seed(seed); }

  [[nodiscard]] static constexpr auto min() noexcept -> result_type { return 0; }
  // while zero-state is considered invalid, PRNG can still produce 0 as a result
  [[nodiscard]] static constexpr auto max() noexcept -> result_type { return std::numeric_limits<result_type>::max(); }

  constexpr void seed(result_type seed) noexcept {
    SplitMix32 splitmix{seed};
    this->s[0] = splitmix();  // Xoshiro family recommends using
    this->s[1] = splitmix();  // splitmix to initialize its state
    this->s[2] = splitmix();
    this->s[3] = splitmix();
  }

  constexpr auto operator()() noexcept -> result_type {
    const result_type result = std::rotl(this->s[0] + this->s[3], 7) + this->s[0];
    const result_type t = s[1] << 9;
    this->s[2] ^= this->s[0];
    this->s[3] ^= this->s[1];
    this->s[1] ^= this->s[2];
    this->s[0] ^= this->s[3];
    this->s[2] ^= t;
    this->s[3] = std::rotl(this->s[3], 11);
    return result;
  }
};

}  // namespace generators

template <class Gen>
constexpr auto generate_canonical_float(Gen& gen) noexcept(noexcept(gen())) -> float {
  using generated_type = typename Gen::result_type;

  constexpr generated_type prng_min = Gen::min();
  constexpr generated_type prng_max = Gen::max();

  static_assert(prng_min < prng_max);

  constexpr generated_type prng_range = prng_max - prng_min;
  constexpr generated_type type_range = std::numeric_limits<generated_type>::max();
  constexpr bool prng_is_bit_uniform = (prng_range == type_range);

  constexpr int exponent_bits_32 = 8;

  constexpr float mantissa_hex_32 = 0x1.0p-24F;  // == 2^-24, corresponds to 24 significant bits of float

  // Note 1: Note hexadecimal float literals, 'p' separates hex-base from the exponent
  // Note 2: Floats have 'mantissa_size + 1' significant bits due to having a sign bit

  static_assert(std::numeric_limits<float>::digits == 24, "Platform not supported, 'float' is expected to be 32-bit.");
  static_assert(prng_is_bit_uniform && sizeof(float) == 4 && sizeof(generated_type) == 4);

  // 32-bit float, 32-bit uniform PRNG
  // => multiplication algorithm tweaked for 32-bit
  return (gen() >> exponent_bits_32) * mantissa_hex_32;
}

}  // namespace utlrandom

#ifndef GPU_ON
#include <limits>
#include <random>
inline thread_local auto rng{std::mt19937{std::random_device{}()}};
#else
inline thread_local auto rng{utlrandom::generators::Xoshiro128PP{}};
#endif

auto rng_seed(auto&& seedval) { rng.seed(seedval); }

[[nodiscard]] inline auto get_rng_random_seed() -> std::int64_t {
#ifndef GPU_ON
  return std::random_device{}();
#else
  return utlrandom::_crush_to_uint32(std::chrono::high_resolution_clock::now().time_since_epoch().count());
#endif
}

inline auto rng_uniform() -> float {
  while (true) {
#ifndef GPU_ON
    const auto zrand = std::generate_canonical<float, std::numeric_limits<float>::digits>(rng);
#else
    const auto zrand = utlrandom::generate_canonical_float(rng);
#endif
    if (zrand != 1.) {
      return zrand;
    }
  }
}

inline auto rng_uniform_pos() -> float {
  while (true) {
    const auto zrand = rng_uniform();
    if (zrand > 0) {
      return zrand;
    }
  }
}

#endif  // RANDOM_H

// Standalone generator for the golden values in unittests.cc test_nonthermal_solve_upper_triangular().
// Compile in isolation against the vendored Eigen (from the repository root):
//   clang++ -std=c++20 -O2 -ffp-contract=off -isystem third_party scripts/generate_triangular_solver_golden.cc -o
//   goldengen && ./goldengen
// EIGEN_DONT_VECTORIZE and -ffp-contract=off make the solve run Eigen's scalar triangular_solve_vector code
// path with correctly-rounded IEEE arithmetic only, so the printed values are platform-independent and match
// what nonthermal::solve_upper_triangular() must produce in REPRODUCIBLE builds (see the contract at its
// definition). The matrix/rhs generation code must be kept identical to the unit test.

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <vector>

#ifndef EIGEN_DONT_VECTORIZE
#define EIGEN_DONT_VECTORIZE
#endif
#include <Eigen/Dense>

namespace {
constexpr int SFPTS = 4096;

uint64_t lcgstate = 20260808;

// uniform in [0, 1) from a 64-bit LCG (Knuth MMIX constants); fully deterministic across platforms
auto next_uniform() -> double {
  lcgstate = (lcgstate * 6364136223846793005ULL) + 1442695040888963407ULL;
  return std::ldexp(static_cast<double>(lcgstate >> 11), -53);
}

constexpr auto uppertriangular(const int i, const int j) -> int { return (SFPTS * i) - (i * (i + 1) / 2) + j; }
}  // namespace

auto main() -> int {
  // diagonally dominant upper triangular system: diagonal in [1, 2), off-diagonal magnitude < 2^-13
  std::vector<double> uppertri(static_cast<size_t>(SFPTS) * (SFPTS + 1) / 2);
  for (int i = 0; i < SFPTS; i++) {
    uppertri[uppertriangular(i, i)] = 1. + next_uniform();
    for (int j = i + 1; j < SFPTS; j++) {
      uppertri[uppertriangular(i, j)] = std::ldexp(next_uniform() - 0.5, -12);
    }
  }
  std::array<double, SFPTS> bvec{};
  for (auto& b : bvec) {
    b = 2. * (next_uniform() - 0.5);
  }
  // exact zeros in the last rows (solved first) exercise the skip-division-on-exact-zero path
  for (int i = SFPTS - 8; i < SFPTS; i++) {
    bvec[i] = 0.;
  }

  // expand to a full row-major matrix and solve with Eigen exactly as nonthermal.cc did before the rewrite
  // (a Map over flat storage, since a fixed-size 4096^2 Eigen matrix exceeds the stack-allocation limit)
  std::vector<double> fullmatrix(static_cast<size_t>(SFPTS) * SFPTS, 0.);
  for (int i = 0; i < SFPTS; i++) {
    for (int j = i; j < SFPTS; j++) {
      fullmatrix[(static_cast<size_t>(i) * SFPTS) + j] = uppertri[uppertriangular(i, j)];
    }
  }
  const Eigen::Map<const Eigen::Matrix<double, SFPTS, SFPTS, Eigen::RowMajor>> eigen_matrix{fullmatrix.data()};
  const Eigen::Map<const Eigen::Vector<double, SFPTS>> eigen_b{bvec.data()};
  std::array<double, SFPTS> xvec{};
  Eigen::Map<Eigen::Vector<double, SFPTS>> eigen_x{xvec.data()};
  eigen_x = eigen_matrix.triangularView<Eigen::Upper>().solve(eigen_b);

  // FNV-1a 64-bit hash over the solution, hashing each double's bits low byte first (endian-independent)
  uint64_t hash = 14695981039346656037ULL;
  for (const double x : xvec) {
    uint64_t bits = 0;
    std::memcpy(&bits, &x, sizeof(bits));
    for (int byteindex = 0; byteindex < 8; byteindex++) {
      hash ^= (bits >> (8 * byteindex)) & 0xffU;
      hash *= 1099511628211ULL;
    }
  }

  std::printf("hash 0x%016llx\n", static_cast<unsigned long long>(hash));
  std::printf("x[0]    %la\n", xvec[0]);
  std::printf("x[2048] %la\n", xvec[2048]);
  std::printf("x[4095] %la\n", xvec[4095]);
  return 0;
}

// Anderson acceleration of a fixed-point iteration x = G(x) for a small state vector (see
// NLTE_OUTER_ANDERSON_DEPTH). The accelerator keeps the last iterates and their map outputs, and
// combines them so that the residual of the combination is minimal in the least-squares sense.

#ifndef ANDERSON_H
#define ANDERSON_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

// N: dimension of the state. MAXDEPTH: the maximum number of previous residual differences. At most N
// residual differences are independent, so a larger MAXDEPTH only wastes solves.
template <std::size_t N, std::size_t MAXDEPTH = N>
class AndersonAccelerator {
  static_assert(MAXDEPTH >= 1);
  static_assert(MAXDEPTH <= N);

 public:
  using State = std::array<double, N>;
  static constexpr std::size_t max_depth = MAXDEPTH;

  explicit constexpr AndersonAccelerator(const std::size_t depth) : depth_(std::min(depth, MAXDEPTH)) {}

  // Forget the history. Call this when the map changes, e.g. when the set of solved ions changes.
  constexpr void reset() { nstored_ = 0; }

  // Give the iterate x and its map output g = G(x). Return the next iterate. With no history, or with a
  // depth of zero, the result is g (plain successive substitution).
  [[nodiscard]] constexpr auto next(const State& x, const State& g) -> State {
    State residual{};
    for (std::size_t i = 0; i < N; i++) {
      residual[i] = g[i] - x[i];
    }

    State result = g;
    const std::size_t m_max = std::min(nstored_, depth_);

    // consecutive differences: entry j is (pair k-j) - (pair k-j-1), with pair k = (g, residual). Each
    // difference of the residuals gets a unit norm, so the test of the solve below sees the angle between
    // the differences and not their size
    std::array<State, MAXDEPTH> dres_chain{};
    std::array<State, MAXDEPTH> dg_chain{};
    std::array<double, MAXDEPTH> dres_norm{};
    {
      State r_newer = residual;
      State g_newer = g;
      for (std::size_t j = 0; j < m_max; j++) {
        const auto& [gprev, rprev] = history_[(nstored_ - 1 - j) % MAXDEPTH];
        double sumsq = 0.;
        for (std::size_t i = 0; i < N; i++) {
          dres_chain[j][i] = r_newer[i] - rprev[i];
          dg_chain[j][i] = g_newer[i] - gprev[i];
          sumsq += dres_chain[j][i] * dres_chain[j][i];
        }
        dres_norm[j] = std::sqrt(sumsq);
        if (dres_norm[j] > 0.) {
          for (std::size_t i = 0; i < N; i++) {
            dres_chain[j][i] /= dres_norm[j];
          }
        }
        r_newer = rprev;
        g_newer = gprev;
      }
    }

    // a singular system (e.g. a 2-cycle, where two residual differences are exact negatives) makes
    // the solve fail. The loop then tries each smaller depth, down to the plain map output
    for (std::size_t m = m_max; m > 0; m--) {
      // normal equations (m x m) for the coefficients gamma: minimise |residual - sum_j gamma_j dres_chain[j]|
      std::array<std::array<double, MAXDEPTH>, MAXDEPTH> ata{};
      std::array<double, MAXDEPTH> atb{};
      bool zero_difference = false;
      for (std::size_t a = 0; a < m; a++) {
        zero_difference = zero_difference || !(dres_norm[a] > 0.);
        for (std::size_t b = 0; b < m; b++) {
          double sum = 0.;
          for (std::size_t i = 0; i < N; i++) {
            sum += dres_chain[a][i] * dres_chain[b][i];
          }
          ata[a][b] = sum;
        }
        double sum = 0.;
        for (std::size_t i = 0; i < N; i++) {
          sum += dres_chain[a][i] * residual[i];
        }
        atb[a] = sum;
      }
      if (zero_difference) {
        continue;
      }

      std::array<double, MAXDEPTH> gamma{};
      if (solve_small_system(ata, atb, m, gamma)) {
        for (std::size_t j = 0; j < m; j++) {
          for (std::size_t i = 0; i < N; i++) {
            result[i] -= gamma[j] / dres_norm[j] * dg_chain[j][i];
          }
        }
        break;
      }
    }

    history_[nstored_ % MAXDEPTH] = {g, residual};
    nstored_++;
    return result;
  }

 private:
  struct Entry {
    State g{};
    State residual{};
  };

  // Gaussian elimination with partial pivoting for the m x m normal equations of unit vectors. Return
  // false when the system is singular, so that the caller tries a smaller depth.
  static constexpr auto solve_small_system(std::array<std::array<double, MAXDEPTH>, MAXDEPTH> a,
                                           std::array<double, MAXDEPTH> b, const std::size_t m,
                                           std::array<double, MAXDEPTH>& solution) -> bool {
    for (std::size_t col = 0; col < m; col++) {
      std::size_t pivot = col;
      for (std::size_t row = col + 1; row < m; row++) {
        if (std::abs(a[row][col]) > std::abs(a[pivot][col])) {
          pivot = row;
        }
      }
      // the matrix holds the cosines of unit vectors, so a pivot at the round-off level means two
      // parallel residual differences (e.g. a 2-cycle) or a rank deficient system
      if (!(std::abs(a[pivot][col]) > 1e-8)) {
        return false;
      }
      std::swap(a[col], a[pivot]);
      std::swap(b[col], b[pivot]);
      for (std::size_t row = col + 1; row < m; row++) {
        const double factor = a[row][col] / a[col][col];
        for (std::size_t k = col; k < m; k++) {
          a[row][k] -= factor * a[col][k];
        }
        b[row] -= factor * b[col];
      }
    }
    for (std::size_t col = m; col-- > 0;) {
      double sum = b[col];
      for (std::size_t k = col + 1; k < m; k++) {
        sum -= a[col][k] * solution[k];
      }
      solution[col] = sum / a[col][col];
      if (!std::isfinite(solution[col])) {
        return false;
      }
    }
    return true;
  }

  std::size_t depth_;
  std::size_t nstored_{0};
  std::array<Entry, MAXDEPTH> history_{};
};

#endif  // ANDERSON_H

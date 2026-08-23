// Anderson acceleration of a fixed-point iteration x = G(x) for a small state vector (see
// NLTE_OUTER_ANDERSON_DEPTH). The accelerator keeps the last iterates and their map outputs, and
// combines them so that the residual of the combination is minimal in the least-squares sense.

#ifndef ANDERSON_H
#define ANDERSON_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

// N: dimension of the state. MAXDEPTH: the maximum number of previous residual differences.
template <std::size_t N, std::size_t MAXDEPTH>
class AndersonAccelerator {
 public:
  using State = std::array<double, N>;

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
    // a singular system (e.g. a 2-cycle, where two residual differences are exact negatives) makes
    // the solve fail. A smaller depth then gets a try, down to the plain map output
    for (std::size_t m = std::min(nstored_, depth_); m > 0; m--) {
      // consecutive differences: entry j is (pair k-j) - (pair k-j-1), with pair k = (x, g, residual)
      std::array<State, MAXDEPTH> dres_chain{};
      std::array<State, MAXDEPTH> dg_chain{};
      {
        State r_newer = residual;
        State g_newer = g;
        for (std::size_t j = 0; j < m; j++) {
          const auto& [xprev, gprev, rprev] = history_[(nstored_ - 1 - j) % MAXDEPTH];
          for (std::size_t i = 0; i < N; i++) {
            dres_chain[j][i] = r_newer[i] - rprev[i];
            dg_chain[j][i] = g_newer[i] - gprev[i];
          }
          r_newer = rprev;
          g_newer = gprev;
        }
      }

      // normal equations (m x m) for the coefficients gamma: minimise |residual - sum_j gamma_j dres_chain[j]|
      std::array<std::array<double, MAXDEPTH>, MAXDEPTH> ata{};
      std::array<double, MAXDEPTH> atb{};
      for (std::size_t a = 0; a < m; a++) {
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

      std::array<double, MAXDEPTH> gamma{};
      if (solve_small_system(ata, atb, m, gamma)) {
        for (std::size_t j = 0; j < m; j++) {
          for (std::size_t i = 0; i < N; i++) {
            result[i] -= gamma[j] * dg_chain[j][i];
          }
        }
        break;
      }
    }

    history_[nstored_ % MAXDEPTH] = {x, g, residual};
    nstored_++;
    return result;
  }

 private:
  struct Entry {
    State x{};
    State g{};
    State residual{};
  };

  // Gaussian elimination with partial pivoting for the m x m normal equations. Return false when the
  // system is singular, so that the caller falls back to the plain map output.
  static constexpr auto solve_small_system(std::array<std::array<double, MAXDEPTH>, MAXDEPTH> a,
                                           std::array<double, MAXDEPTH> b, const std::size_t m,
                                           std::array<double, MAXDEPTH>& solution) -> bool {
    double scale = 0.;
    for (std::size_t row = 0; row < m; row++) {
      for (std::size_t col = 0; col < m; col++) {
        scale = std::max(scale, std::abs(a[row][col]));
      }
    }
    if (!(scale > 0.)) {
      return false;
    }
    for (std::size_t col = 0; col < m; col++) {
      std::size_t pivot = col;
      for (std::size_t row = col + 1; row < m; row++) {
        if (std::abs(a[row][col]) > std::abs(a[pivot][col])) {
          pivot = row;
        }
      }
      // a pivot at the round-off level means a rank deficient system, e.g. more residual differences
      // than state components
      if (!(std::abs(a[pivot][col]) > 1e-12 * scale)) {
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

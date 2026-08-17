// Declarations for input file reading (input.cc).

#ifndef INPUT_H
#define INPUT_H

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <concepts>
#include <cstdlib>
#include <istream>
#include <limits>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include "constants.h"
#include "globals.h"
#include "packet.h"

void read_atomicdata();
void read_parameterfile(std::span<Packet> packets);
void update_parameterfile(int nts);
void setup_timesteps();

// calculate the timestep grid: ntimesteps entries with start, mid, and width set, plus a zero-width
// entry holding the end time. tmin and tmax are in seconds; the fixed width and transition time of the
// hybrid methods are in days
[[nodiscard]] auto calculate_timesteps(TimeStepSizeMethod method, double tmin, double tmax, int ntimesteps,
                                       double fixed_timestep_width_days, double timestep_transition_time_days)
    -> std::vector<globals::TimeStep>;

// Count the levels in the ground term of an ion. The ground term is the multiplet with the lowest energy.
// The caller gives the level energies and the statistical weights in the sequence of increasing energy.
//
// In an LS-coupled term, J changes by one from each level to the next level. Thus the statistical weight
// 2J+1 changes by two, in the same direction for the full term. A term with an integer J has an odd number
// of levels. The Landé interval rule gives a splitting E(J) - E(J-1) that is proportional to J. Thus the
// ratio of two splittings in a term is equal to the ratio of the higher J of each pair of levels. The term
// ends where a splitting is more than a tolerance factor larger than this rule permits, or where the
// sequence of statistical weights changes.
//
// The function uses two tolerances, because the reference splitting is not always in the same term:
//  - The splitting below the level is in the same term. The interval rule is accurate to approximately a
//    factor of two here. The largest error in the ARTIS data is 4.4, for Ni I 3F4 with 3D3,2 between its
//    levels.
//  - The splitting above the level is the reference only for the second level of the ion. If the ground
//    term has one level (for example 4S3/2, 6S5/2 or 7S3), this reference is the fine structure of the next
//    term. The fine structure is much smaller than the distance between the two terms, thus these
//    boundaries need a factor of 25 or more. The widest correct term needs 3.4 (the C-like 3P of Fe XXI,
//    where jj coupling makes the interval rule inaccurate).
//
// The function has these known limitations:
//  - The statistical weight must change by two from each listed level to the next level. Thus the function
//    stops too soon if the data does not list the levels in the sequence of J. This occurs for the 3P term
//    of the Po-like ions (Po I, At II, Rn III and Fr IV in the kilonova data). There jj coupling puts 3P0
//    below 3P1, the weights are 5, 1, 3, and the function counts only the ground level. To find these
//    terms, the function must examine the set of the weights and not their sequence. No rule that was
//    tested can separate this case from Mn II 7S3 / 5S2 / 5D4, where the weights are 7, 5, 9 and the
//    energies are similar.
//  - Some data does not resolve a fine-structure splitting. If the data gives a very small splitting in
//    place of zero, the function reads it as a term boundary and stops too soon. If the data gives two
//    levels with the same energy directly above the ground level, the function reads them as an unresolved
//    next term. This is correct for N I 4S3/2 with 2D5/2,3/2 above it, but incorrect for a 3P2,1,0 term
//    whose two upper levels have the same energy.
[[nodiscard]] constexpr auto count_groundterm_levels(const std::span<const double> energies,
                                                     const std::span<const float> statweights) -> int {
  assert_testmodeonly(energies.size() == statweights.size());
  const auto nlevels = static_cast<int>(energies.size());
  if (nlevels <= 1) {
    return nlevels;
  }

  // the largest factor by which a splitting can be more than the interval rule permits and stay in the term
  constexpr double MAX_DEPARTURE_WITHIN_MULTIPLET = 3.;
  constexpr double MAX_DEPARTURE_ACROSS_TERMS = 8.;

  // Two times the higher J of the pair of levels (a, a + 1), from g = 2J + 1. For a normal term this is the
  // upper level, and for an inverted term the lower level. The interval rule is proportional to this J.
  const auto twoj_upper = [&statweights](const int a) -> double {
    return static_cast<double>(std::max(statweights[a], statweights[a + 1])) - 1.;
  };

  // Return true if the interval rule permits the splitting of the pair (a, a + 1) and the splitting of the
  // reference pair (b, b + 1) in the same term. The test compares two products and not a ratio. Thus a
  // reference splitting of zero gives false and does not cause a division by zero.
  const auto splittings_consistent = [&](const int a, const int b, const double maxdeparture) -> bool {
    const double splitting_a = energies[a + 1] - energies[a];
    const double splitting_b = energies[b + 1] - energies[b];
    if (splitting_a < 0. || splitting_b < 0. || twoj_upper(b) <= 0.) {
      return false;  // energies not in sequence, or a reference pair with no change of J: no prediction
    }
    return (splitting_a * twoj_upper(b)) <= (maxdeparture * splitting_b * twoj_upper(a));
  };

  const bool increasing_j = statweights[1] > statweights[0];
  int nlevels_groundterm = 1;
  for (int level = 1; level < nlevels; level++) {
    const double delta_g = statweights[level] - statweights[level - 1];
    const double abs_delta_g = (delta_g < 0.) ? -delta_g : delta_g;  // std::abs is not constexpr in libc++
    // J must change by one, in the same direction for the full term. The limits accept a statistical weight
    // that the data does not give as an exact integer. They also reject a repeated statistical weight.
    if (abs_delta_g < 1.6 || abs_delta_g > 2.4 || ((delta_g > 0.) != increasing_j)) {
      break;
    }

    // Use the nearest splitting below the level that the data gives as more than zero. Some data rounds a
    // fine-structure splitting to zero, and a zero gives no scale for the comparison. If there is no such
    // splitting, use the splitting above the level, with the larger tolerance. A zero splitting above the
    // level shows that the next term is unresolved, which ends this term. Keep a level that has no reference
    // (the second level of the ion), because there is nothing to compare it with.
    int refpair = -1;
    for (int b = level - 2; b >= 0; b--) {
      if (energies[b + 1] > energies[b]) {
        refpair = b;
        break;
      }
    }
    if (refpair >= 0) {
      if (!splittings_consistent(level - 1, refpair, MAX_DEPARTURE_WITHIN_MULTIPLET)) {
        break;
      }
    } else if (level + 1 < nlevels) {
      if (!splittings_consistent(level - 1, level, MAX_DEPARTURE_ACROSS_TERMS)) {
        break;
      }
    } else if (level > 1 && !splittings_consistent(level - 1, level - 2, MAX_DEPARTURE_WITHIN_MULTIPLET)) {
      break;
    }
    nlevels_groundterm = level + 1;
  }

  // A term with an integer J has odd statistical weights and an odd number of levels. Thus an even count
  // shows that the last level is in the next term (for example 7S3 with 5S2 above it). Apply this rule only
  // if the loop found a level outside the term. If the list of levels ended first, keep all of them, because
  // the level limit in compositiondata.txt can cut an ion in the middle of its ground term.
  // Round the weight and do not truncate it: a weight of 6.999999 must also give an odd result.
  int g_ground = static_cast<int>(statweights[0]);
  if ((statweights[0] - static_cast<float>(g_ground)) > 0.5F) {
    g_ground++;
  }
  if (((g_ground % 2) == 1) && (nlevels_groundterm % 2) == 0 && nlevels_groundterm < nlevels) {
    nlevels_groundterm--;
  }

  return nlevels_groundterm;
}

// The energies are in eV and come from real atomic data. The unit tests give more examples.
static_assert(count_groundterm_levels(std::array{0.}, std::array{5.F}) == 1);  // single level
static_assert(count_groundterm_levels(std::array{0., 19.8196, 20.6158}, std::array{1.F, 3.F, 1.F}) ==
              1);  // He I 1S0 then 3S1 then 1S0
static_assert(count_groundterm_levels(std::array{0., 0.006}, std::array{1.F, 3.F}) ==
              2);  // 3P0,1 of an ion truncated to two levels: no boundary seen, so both are kept
static_assert(count_groundterm_levels(std::array{0., 0., 0.1, 1.9}, std::array{1.F, 3.F, 5.F, 5.F}) ==
              3);  // 3P0,1,2 with the first splitting rounded to zero in the data
static_assert(count_groundterm_levels(std::array{0., 0., 5.0}, std::array{1.F, 3.F, 5.F}) ==
              1);  // a 5 eV jump after a splitting rounded to zero is still a term boundary
static_assert(count_groundterm_levels(std::array{0., 2.38, 2.38, 3.58}, std::array{4.F, 6.F, 4.F, 2.F}) ==
              1);  // 4S3/2 then a 2D pair with zero splitting in the data
static_assert(count_groundterm_levels(std::array{0., 0.006, 0.0162, 1.899, 4.05},
                                      std::array{1.F, 3.F, 5.F, 5.F, 1.F}) == 3);  // N II 3P0,1,2 (normal)
static_assert(count_groundterm_levels(std::array{0., 0.0196, 0.0281, 1.9674, 4.1897},
                                      std::array{5.F, 3.F, 1.F, 5.F, 1.F}) ==
              3);  // O I 3P2,1,0 (inverted, splitting ratio 2.3)
static_assert(count_groundterm_levels(std::array{0., 9.1524, 14.5438, 30.3089, 46.0904},
                                      std::array{1.F, 3.F, 5.F, 5.F, 1.F}) ==
              3);  // Fe XXI 3P0,1,2: jj coupling shrinks the second splitting to 0.29 of the interval rule
static_assert(count_groundterm_levels(std::array{0., 2.3835, 2.3846, 3.5756, 3.5756},
                                      std::array{4.F, 6.F, 4.F, 2.F, 4.F}) == 1);  // N I 4S3/2 then 2D5/2,3/2
static_assert(count_groundterm_levels(std::array{0., 1.1745, 1.7762, 1.8094, 1.8326, 1.8475, 1.8548},
                                      std::array{7.F, 5.F, 9.F, 7.F, 5.F, 3.F, 1.F}) ==
              1);  // Mn II 7S3 then 5S2 then 5D4..0 (odd level count for integer J)
static_assert(count_groundterm_levels(std::array{0., 0.0477, 0.0828, 0.1070, 0.1211, 0.2322, 0.3013},
                                      std::array{10.F, 8.F, 6.F, 4.F, 2.F, 10.F, 8.F}) ==
              5);  // Fe II 6D9/2..1/2 then a4F
static_assert(count_groundterm_levels(std::array{0., 0.0176, 0.0517, 0.0996, 0.1590, 2.9825, 3.0912},
                                      std::array{1.F, 3.F, 5.F, 7.F, 9.F, 1.F, 9.F}) == 5);  // Fe V 5D0..4 (normal)
static_assert(count_groundterm_levels(std::array{0., 0.94, 0.96, 0.97}, std::array{7.F, 5.F, 1.F, 3.F}) ==
              1);  // Cr I 7S3 then 5S2 then 5D0,1
static_assert(count_groundterm_levels(std::array{0., 0., 0.0209, 1.4283}, std::array{4.F, 4.F, 6.F, 4.F}) ==
              1);  // repeated statistical weight ends the multiplet
static_assert(count_groundterm_levels(std::array{0., 0.0254, 0.1091, 0.1652, 0.2124},
                                      std::array{9.F, 7.F, 5.F, 7.F, 3.F}) ==
              1);  // Ni I 3F4 with 3D3,2 interleaved: the second splitting is far larger than the inverted 3F allows

// return true for whitespace-only lines, and lines that are exclusively whitespace up to a '#' character
[[nodiscard]] constexpr auto lineiscommentonly(const std::string_view line) -> bool {
  for (char const i : line) {
    if (i == '#') {  // anything to the right of a # character doesn't count
      return true;
    }
    if (i != ' ' && i != '\t' && i != '\r' && i != '\n') {
      return false;
    }
  }
  return true;
}

// read the next line, skipping any comment lines beginning with '#'
inline auto get_noncommentline(std::istream& input, std::string& line) -> bool {
  while (true) {
    const bool linefound = !(!std::getline(input, line));
    if (!linefound) {
      return false;
    }
    if (!lineiscommentonly(line)) {
      return true;
    }
  }
}

// parse the next whitespace-delimited token of the line as a number and advance past it.
// Return false if there is no token left or the token is not fully numeric.
// Accepts the same number spellings as stream extraction: leading plus signs are allowed, magnitudes below the
// range of T (or of double) read as zero, and out-of-range magnitudes and nan/inf spellings are rejected
template <typename T>
[[nodiscard]] inline auto parse_next_token(std::string_view& remainder, T& value) -> bool {
  constexpr std::string_view whitespace = " \t\r";
  const auto tokenstart = remainder.find_first_not_of(whitespace);
  if (tokenstart == std::string_view::npos) {
    return false;
  }
  const auto tokenend = std::min(remainder.find_first_of(whitespace, tokenstart), remainder.size());
  auto token = remainder.substr(tokenstart, tokenend - tokenstart);
  // stream extraction accepted a leading plus sign, but from_chars does not. Only remove the plus if it is not
  // followed by another sign, so that malformed tokens like "+-0" stay rejected
  if (token.size() >= 2 && token.front() == '+' && token[1] != '-' && token[1] != '+') {
    token.remove_prefix(1);
  }
  const auto* const tokenfirst = token.data();
  const auto* const tokenlast = std::to_address(token.cend());
  const auto [ptr, ec] = std::from_chars(tokenfirst, tokenlast, value);
  if (ptr != tokenlast) {
    return false;
  }
  if constexpr (std::floating_point<T>) {
    if (ec == std::errc{} && !std::isfinite(value)) {
      // reject nan and inf spellings, which from_chars accepts but stream extraction did not
      return false;
    }
    if (ec == std::errc::result_out_of_range) {
      // a syntactically valid number outside the range of T. Stream extraction stored zero for underflow (even
      // below double range, such as 1e-400) but rejected overflow, so use strtod, which returns zero for
      // underflow, to distinguish the two (the token views a null-terminated getline buffer, and strtod stops
      // at the whitespace or end of line that terminates the token)
      char* parseend{};
      const double dvalue = std::strtod(tokenfirst, &parseend);
      if (parseend == tokenlast && std::fabs(dvalue) <= std::numeric_limits<T>::max()) {
        value = static_cast<T>(dvalue);
        remainder.remove_prefix(tokenend);
        return true;
      }
    }
  }
  if (ec != std::errc{}) {
    return false;
  }
  remainder.remove_prefix(tokenend);
  return true;
}

#endif  // INPUT_H

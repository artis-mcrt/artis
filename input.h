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

// Count the levels of the ground term (the lowest fine-structure multiplet) of an ion from its level
// energies and statistical weights, both listed in order of increasing energy.
//
// Within an LS-coupled term, J changes by one from level to level, so the statistical weight 2J+1 changes by
// exactly two, in the same direction throughout (normal or inverted multiplets), and a term with integer J has
// an odd number of levels (2 min(L,S) + 1). By the Landé interval rule the fine-structure splittings
// E(J) - E(J-1) are proportional to J, so the ratio of two splittings within a multiplet is the ratio of the
// higher J of each pair of levels. A splitting that exceeds that prediction by more than a tolerance factor
// marks the boundary with the next term, as does a change in the statistical weight pattern.
//
// Two tolerances are needed because the reference splitting is not always a fine-structure splitting. Comparing
// with the splitting below stays within the multiplet, where the interval rule holds to about a factor of two
// (the largest departure in the ARTIS datasets is 4.4, for Ni I 3F4 with 3D3,2 interleaved). Comparing with the
// splitting above happens at the lowest level, and for a single-level ground term (4S3/2, 6S5/2, 7S3) that
// reference is the fine structure of the *next* term, which is far smaller than the term separation: those
// boundaries need a factor of 25 or more, while the widest real multiplet needs 3.4 (the C-like 3P of Fe XXI,
// which is deep enough into jj coupling that the interval rule barely applies).
//
// Known limitations:
//  - The statistical weight has to step by two from one listed level to the next, so a term whose levels are not
//    in J order is cut short. This happens to the 3P of the Po-like sequence (Po I, At II, Rn III, Fr IV in the
//    kilonova datasets), where jj coupling pushes 3P0 below 3P1 and the weights run 5, 1, 3: only the ground
//    level is counted. Detecting it needs the set of weights rather than their order, and no rule tried so far
//    separates it from Mn II 7S3 / 5S2 / 5D4, whose weights run 7, 5, 9 with a similar energy pattern.
//  - Datasets that do not resolve a fine-structure splitting: one written as a very small non-zero number
//    rather than rounded to zero looks like a term boundary and ends the multiplet early, and a pair of
//    degenerate levels directly above the ground level reads as an unresolved next term (which is what it is
//    for N I 4S3/2 followed by 2D5/2,3/2, but not for a 3P2,1,0 whose upper two levels are degenerate).
[[nodiscard]] constexpr auto count_groundterm_levels(const std::span<const double> energies,
                                                     const std::span<const float> statweights) -> int {
  assert_testmodeonly(energies.size() == statweights.size());
  const auto nlevels = static_cast<int>(energies.size());
  if (nlevels <= 1) {
    return nlevels;
  }

  // largest factor by which a splitting may exceed the interval-rule prediction and still count as the same term
  constexpr double MAX_DEPARTURE_WITHIN_MULTIPLET = 3.;
  constexpr double MAX_DEPARTURE_ACROSS_TERMS = 8.;

  // twice the higher J of the pair of levels (a, a + 1), from g = 2J + 1. Equal to 2J of the upper level of a
  // normal multiplet and of the lower level of an inverted one, which is the J the interval rule scales with.
  const auto twoj_upper = [&statweights](const int a) -> double {
    return static_cast<double>(std::max(statweights[a], statweights[a + 1])) - 1.;
  };

  // true if the splitting between levels (a, a + 1) is small enough compared to the splitting between levels
  // (b, b + 1) for the two pairs to belong to the same multiplet under the interval rule. Written as a pair of
  // products rather than a ratio so that a zero reference splitting fails rather than dividing by zero.
  const auto splittings_consistent = [&](const int a, const int b, const double maxdeparture) -> bool {
    const double splitting_a = energies[a + 1] - energies[a];
    const double splitting_b = energies[b + 1] - energies[b];
    if (splitting_a < 0. || splitting_b < 0. || twoj_upper(b) <= 0.) {
      return false;  // levels out of energy order, or a reference pair with no J step: no usable prediction
    }
    return (splitting_a * twoj_upper(b)) <= (maxdeparture * splitting_b * twoj_upper(a));
  };

  const bool increasing_j = statweights[1] > statweights[0];
  int nlevels_groundterm = 1;
  for (int level = 1; level < nlevels; level++) {
    const double delta_g = statweights[level] - statweights[level - 1];
    const double abs_delta_g = (delta_g < 0.) ? -delta_g : delta_g;  // std::abs is not constexpr in libc++
    // J must change by exactly one, in the same direction throughout. The window tolerates a statistical weight
    // that the dataset did not write as an exact integer; it also rejects a repeated statistical weight.
    if (abs_delta_g < 1.6 || abs_delta_g > 2.4 || ((delta_g > 0.) != increasing_j)) {
      break;
    }

    // Prefer the nearest splitting below that the data resolves as non-zero: datasets sometimes round a
    // fine-structure splitting to zero, and a zero carries no scale to compare against. Failing that, compare
    // with the splitting immediately above, which for a single-level ground term is the fine structure of the
    // next term rather than another splitting of this one, hence the looser tolerance (a zero there means the
    // next term is unresolved, which ends this one). A level with no reference at all, the second level of the
    // ion, is kept: there is nothing to compare it with.
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

  // A term with integer J (odd statistical weights) has an odd number of levels (2 min(L,S) + 1), so an even
  // count means the last level belongs to the next term (e.g. 7S3 followed by 5S2). This only applies when the
  // level after the counted ones was seen to be outside the term; if the level list itself ended (an ion
  // truncated by the level limit in compositiondata.txt) the retained levels are kept.
  // Round rather than truncate: a statistical weight written as 6.999999 must still test as odd.
  int g_ground = static_cast<int>(statweights[0]);
  if ((statweights[0] - static_cast<float>(g_ground)) > 0.5F) {
    g_ground++;
  }
  if (((g_ground % 2) == 1) && (nlevels_groundterm % 2) == 0 && nlevels_groundterm < nlevels) {
    nlevels_groundterm--;
  }

  return nlevels_groundterm;
}

// energies in eV from real atomic datasets (see the runtime unit tests for more)
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

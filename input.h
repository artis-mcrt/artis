// Declarations for input file reading (input.cc).

#ifndef INPUT_H
#define INPUT_H

#include <algorithm>
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

#include "packet.h"

void read_atomicdata();
void read_parameterfile(std::span<Packet> packets);
void update_parameterfile(int nts);
void setup_timesteps();

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

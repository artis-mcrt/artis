#ifndef INPUT_H
#define INPUT_H

#include <istream>
#include <string>
#include <string_view>

void input(int rank);
void read_parameterfile(int rank);
void update_parameterfile(int nts);
void setup_timesteps();
void write_timestep_file();
auto get_noncommentline(std::istream &input, std::string &line) -> bool;

// return true for whitespace-only lines, and lines that are exclusively whitespace up to a '#' character
[[nodiscard]] constexpr auto lineiscommentonly(const std::string_view line) -> bool {
  for (char const i : line) {
    if (i == '#') {  // anything to the right of a # character doesn't count
      return true;
    }
    if (i != ' ') {
      return false;
    }
  }
  return true;
}

#endif  // INPUT_H

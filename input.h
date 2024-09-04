#ifndef INPUT_H
#define INPUT_H

#include <string>

void input(int rank);
void read_parameterfile(int rank);
void update_parameterfile(int nts);
void time_init();
void write_timestep_file();
auto get_noncommentline(std::fstream &input, std::string &line) -> bool;

[[nodiscard]] constexpr auto lineiscommentonly(const std::string_view line) -> bool
// return true for whitespace-only lines, and lines that are exclusively whitespace up to a '#' character
{
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

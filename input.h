#ifndef INPUT_H
#define INPUT_H

#include <istream>
#include <string>
#include <string_view>

void read_atomicdata();
void read_parameterfile();
void update_parameterfile(int nts);
void setup_timesteps();

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

#endif  // INPUT_H

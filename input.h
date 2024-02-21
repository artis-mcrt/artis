#ifndef INPUT_H
#define INPUT_H

#include <string>

#include "exspec.h"

void input(int rank);
void read_parameterfile(int rank);
void update_parameterfile(int nts);
void time_init();
void write_timestep_file();
auto get_noncommentline(std::fstream &input, std::string &line) -> bool;
auto is_bfest(int element, int ion, int level) -> bool;

static inline auto lineiscommentonly(const std::string &line) -> bool
// return true for whitepace-only lines, and lines that are exclusively whitepace up to a '#' character
{
  for (size_t i = 0; i < line.length(); i++) {
    if (line[i] == '#') {  // anything to the right of a # character doesn't count
      return true;
    }
    if (line[i] != ' ') {
      return false;
    }
  }
  return true;
}

#endif  // INPUT_H

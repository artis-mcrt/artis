#ifndef INPUT_H
#define INPUT_H

#include <string>

#include "exspec.h"

void input(int rank);
void read_parameterfile(int rank);
void update_parameterfile(int nts);
void time_init();
void write_timestep_file();
bool get_noncommentline(std::fstream &input, std::string &line);

static inline bool lineiscommentonly(const std::string &line)
// return true for whitepace-only lines, and lines that are exclusively whitepace up to a '#' character
{
  int searchlength = line.find('#');  // ignore anything to the right of a # character
  if (searchlength < 0) {
    searchlength = line.length();
  }

  for (int i = 0; i < searchlength; i++) {
    if (line[i] != ' ') {
      return false;
    }
  }
  return true;
}

#endif  // INPUT_H

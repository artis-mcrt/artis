#ifndef INPUT_H
#define INPUT_H

#include <iostream>
#include <string>

#include "exspec.h"

void input(int rank);
void read_parameterfile(int rank);
void update_parameterfile(int nts);
void time_init(void);
void write_timestep_file(void);
bool get_noncommentline(std::istream &input, std::string &line);

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

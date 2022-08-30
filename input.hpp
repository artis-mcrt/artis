#ifndef INPUT_H
#define INPUT_H

#include <iostream>
#include <string>

#include "exspec.hpp"

void input(int rank);
void read_parameterfile(int rank);
void update_parameterfile(int nts);
void time_init(void);
void write_timestep_file(void);
bool get_noncommentline(std::istream &input, std::string &line);
bool lineiscommentonly(std::string &line);

#endif  // INPUT_H

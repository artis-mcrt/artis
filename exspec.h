// Constants for the exspec post-processing tool (exspec.cc).

#ifndef EXSPEC_H
#define EXSPEC_H

#include <cstddef>

constexpr ptrdiff_t MNUBINS = 1000;

constexpr int NPHIBINS = 10;
constexpr int NCOSTHETABINS = 10;
constexpr int MABINS = NPHIBINS * NCOSTHETABINS;

#endif  // EXSPEC_H

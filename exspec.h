// Constants for the exspec post-processing tool (exspec.cc).

#ifndef EXSPEC_H
#define EXSPEC_H

#include <cstddef>

#include "constants.h"

constexpr ptrdiff_t MNUBINS = 1000;

// frequency range of the spectrum of the escaped gamma packets
constexpr double NU_MIN_GAMMA = 0.05 * MEV / H;
constexpr double NU_MAX_GAMMA = 4. * MEV / H;

constexpr int NPHIBINS = 10;
constexpr int NCOSTHETABINS = 10;
constexpr int MABINS = NPHIBINS * NCOSTHETABINS;

#endif  // EXSPEC_H

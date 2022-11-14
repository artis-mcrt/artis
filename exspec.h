#ifndef EXSPEC_H
#define EXSPEC_H

/// Spectrum data structure
constexpr int MNUBINS = 1000;

constexpr int NPHIBINS = 10;
constexpr int NCOSTHETABINS = 10;
constexpr int MABINS = NPHIBINS * NCOSTHETABINS;

extern const bool do_exspec;

#endif  // EXSPEC_H

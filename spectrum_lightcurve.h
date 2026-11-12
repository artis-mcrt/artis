// Declarations for spectrum and light curve construction (spectrum_lightcurve.cc).

#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <span>

#include "packet.h"

constexpr int NPHIBINS = 10;
constexpr int NCOSTHETABINS = 10;
constexpr int MABINS = NPHIBINS * NCOSTHETABINS;

// packets_by_rank holds the packets of each rank. sn3d gives the packets of its own rank only.
void write_light_curves_and_spectra(int nts, std::span<const std::span<const Packet>> packets_by_rank);

// Remove the record before a simulation or a spectrum write changes its input or output files.
void remove_spectra_completion_record();

// Call after all ranks close the final packet files and the spectrum files.
void write_spectra_completion_record();

#endif  // SPECTRUM_H

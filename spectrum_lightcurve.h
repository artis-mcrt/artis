// Declarations for spectrum and light curve construction (spectrum_lightcurve.cc).

#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <span>

#include "packet.h"

// packets_by_rank holds the packets of each rank. sn3d gives the packets of its own rank only.
void write_light_curves_and_spectra(int nts, std::span<const std::span<const Packet>> packets_by_rank);

#endif  // SPECTRUM_H

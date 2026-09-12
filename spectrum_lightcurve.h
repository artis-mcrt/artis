// Declarations for spectrum and light curve construction (spectrum_lightcurve.cc).

#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <span>

#include "packet.h"

void write_light_curves_and_spectra(int nts, std::span<const Packet> packets);

#endif  // SPECTRUM_H

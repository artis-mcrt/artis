#pragma once
#ifndef LIGHT_CURVE_H
#define LIGHT_CURVE_H

#include <string>
#include <vector>

#include "exspec.h"
#include "packet.h"

void add_to_lc_res(const Packet &pkt_ptr, int current_abin, std::vector<double> &light_curve_lum,
                   std::vector<double> &light_curve_lumcmf);

void write_light_curve(const std::string &lc_filename, int current_abin, const std::vector<double> &light_curve_lum,
                       const std::vector<double> &light_curve_lumcmf, int numtimesteps);

#endif  // LIGHT_CURVE_H

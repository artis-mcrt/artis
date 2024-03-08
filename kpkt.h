#pragma once
#ifndef KPKT_H
#define KPKT_H

#include "globals.h"
#include "packet.h"
#include "thermalbalance.h"

constexpr double COOLING_UNDEFINED = -99;

namespace kpkt {

auto get_ncoolingterms() -> int;
void setup_coolinglist();
void set_kpktdiffusion(float kpktdiffusion_timescale_in, int n_kpktdiffusion_timesteps_in);
void calculate_cooling_rates(int modelgridindex, struct heatingcoolingrates *heatingcoolingrates);
void do_kpkt_blackbody(struct Packet &pkt_ptr);
void do_kpkt(struct Packet &pkt_ptr, double t2, int nts);

[[nodiscard]] inline auto get_coolinglistoffset(int element, int ion) -> int {
  return globals::elements[element].ions[ion].coolingoffset;
}

[[nodiscard]] inline auto get_ncoolingterms_ion(int element, int ion) -> int {
  return globals::elements[element].ions[ion].ncoolingterms;
}

}  // namespace kpkt

#endif  // KPKT_H

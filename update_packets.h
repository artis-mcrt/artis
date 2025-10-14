#ifndef UPDATE_PACKETS_H
#define UPDATE_PACKETS_H

#include <span>

#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "packet.h"

void update_packets(int nts, std::span<Packet> packets);

inline double neutrino_release_fraction = 0.;

inline auto calc_neutrino_release_fraction() -> double {
  double power_nonu = 0.;
  double power_total = 0.;
  const double tmid = globals::timesteps[globals::timestep].mid;
  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    const double cellmass = grid::get_rho_tmin(mgi) * grid::get_modelcell_assocvolume_tmin(mgi);

    power_nonu += cellmass * decay::get_particle_injection_rate(nonemptymgi, tmid, decay::DECAYTYPE_BETAPLUS);
    power_nonu += cellmass * decay::get_particle_injection_rate(nonemptymgi, tmid, decay::DECAYTYPE_BETAMINUS);
    power_nonu += cellmass * decay::get_particle_injection_rate(nonemptymgi, tmid, decay::DECAYTYPE_ALPHA);
    power_nonu += cellmass * (decay::get_gamma_emission_rate(nonemptymgi, tmid));

    for (const auto decaytype : decay::all_decaytypes) {
      power_total += cellmass * decay::get_qdot_modelcell(nonemptymgi, tmid, decaytype) * cellmass;
    }
  }
  return (1 - (power_nonu / power_total));
}

#endif  // UPDATE_PACKETS_H

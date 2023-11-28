#ifndef DECAY_H
#define DECAY_H

// #include <cstdio>

#include <array>
#include <string>
#include <vector>

#include "constants.h"
#include "packet.h"

namespace decay {
enum decaytypes {
  DECAYTYPE_ALPHA = 0,
  DECAYTYPE_ELECTRONCAPTURE = 1,
  DECAYTYPE_BETAPLUS = 2,
  DECAYTYPE_BETAMINUS = 3,
  DECAYTYPE_NONE = 4,
  DECAYTYPE_COUNT = 5,
};

constexpr std::array<enum decaytypes, 5> all_decaytypes = {
    decaytypes::DECAYTYPE_ALPHA, decaytypes::DECAYTYPE_ELECTRONCAPTURE, decaytypes::DECAYTYPE_BETAPLUS,
    decaytypes::DECAYTYPE_BETAMINUS, decaytypes::DECAYTYPE_NONE};

void init_nuclides(const std::vector<int> &zlist, const std::vector<int> &alist);
auto get_nucstring_z(const std::string &strnuc) -> int;
auto get_nucstring_a(const std::string &strnuc) -> int;
auto get_num_nuclides() -> int;
auto get_elname(int z) -> const char *;
auto get_nuc_z(int nucindex) -> int;
auto get_nuc_a(int nucindex) -> int;
auto get_nucindex(int z, int a) -> int;
auto nuc_exists(int z, int a) -> bool;
auto nucdecayenergygamma(int nucindex) -> double;
auto nucdecayenergygamma(int z, int a) -> double;
void set_nucdecayenergygamma(int nucindex, double value);
void update_abundances(int modelgridindex, int timestep, double t_current);
auto get_endecay_per_ejectamass_t0_to_time_withexpansion(int modelgridindex, double tstart) -> double;
auto get_modelcell_simtime_endecay_per_mass(int mgi) -> double;
void setup_decaypath_energy_per_mass();
void free_decaypath_energy_per_mass();
auto get_qdot_modelcell(int modelgridindex, double t, int decaytype) -> double;
auto get_particle_injection_rate(int modelgridindex, double t, int decaytype) -> double;
auto get_global_etot_t0_tinf() -> double;
void fprint_nuc_abundances(FILE *estimators_file, int modelgridindex, double t_current, int element);
void setup_radioactive_pellet(double e0, int mgi, struct packet *pkt_ptr);
void cleanup();

auto constexpr nucmass(int z, int a) -> double { return a * MH; }
}  // namespace decay

#endif  // DECAY_H

#ifndef DECAY_H
#define DECAY_H

#include <array>
#include <cstdio>
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

void init_nuclides(const std::vector<int> &custom_zlist, const std::vector<int> &custom_alist);
[[nodiscard]] auto get_nucstring_z(const std::string &strnuc) -> int;
[[nodiscard]] auto get_nucstring_a(const std::string &strnuc) -> int;
[[nodiscard]] auto get_num_nuclides() -> ptrdiff_t;
[[nodiscard]] auto get_elname(int z) -> std::string;
[[nodiscard]] auto get_nuc_z(int nucindex) -> int;
[[nodiscard]] auto get_nuc_a(int nucindex) -> int;
[[nodiscard]] auto get_nucindex(int z, int a) -> int;
[[nodiscard]] auto nuc_exists(int z, int a) -> bool;
[[nodiscard]] auto nucdecayenergygamma(int nucindex) -> double;
[[nodiscard]] auto nucdecayenergygamma(int z, int a) -> double;
void set_nucdecayenergygamma(int nucindex, double value);
void update_abundances(int nonemptymgi, int timestep, double t_current);
[[nodiscard]] auto get_endecay_per_ejectamass_t0_to_time_withexpansion(int nonemptymgi, double tstart) -> double;
[[nodiscard]] auto get_modelcell_simtime_endecay_per_mass(int nonemptymgi) -> double;
void setup_decaypath_energy_per_mass();
void free_decaypath_energy_per_mass();
[[nodiscard]] auto get_qdot_modelcell(int nonemptymgi, double t, int decaytype) -> double;
[[nodiscard]] auto get_particle_injection_rate(int nonemptymgi, double t, int decaytype) -> double;
[[nodiscard]] auto get_gamma_emission_rate(int nonemptymgi, double t) -> double;
[[nodiscard]] auto get_global_etot_t0_tinf() -> double;
void fprint_nuc_abundances(FILE *estimators_file, int nonemptmgi, double t_current, int element);
void setup_radioactive_pellet(double e0, int mgi, Packet &pkt);
void cleanup();

[[nodiscard]] auto constexpr nucmass(int /*z*/, int a) -> double { return a * MH; }
}  // namespace decay

#endif  // DECAY_H

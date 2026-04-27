#ifndef DECAY_H
#define DECAY_H

#include <array>
#include <cstddef>
#include <ostream>
#include <span>
#include <string>

#include "mpi_logging.h"
#include "packet.h"

namespace decay {

enum DecayType : int {
  DECAYTYPE_ALPHA = 0,
  DECAYTYPE_ELECTRONCAPTURE = 1,
  DECAYTYPE_BETAPLUS = 2,
  DECAYTYPE_BETAMINUS = 3,
  DECAYTYPE_NONE = 4,
  DECAYTYPE_SPONTFISSION = 5,
  DECAYTYPE_COUNT = 6,
};

constexpr std::array all_decaytypes{DecayType::DECAYTYPE_ALPHA, DecayType::DECAYTYPE_ELECTRONCAPTURE,
                                    DecayType::DECAYTYPE_BETAPLUS, DecayType::DECAYTYPE_BETAMINUS,
                                    DecayType::DECAYTYPE_SPONTFISSION};

void init_nuclides(std::span<const int> custom_zlist, std::span<const int> custom_alist);
[[nodiscard]] auto decaytype_is_used(DecayType decaytype) -> bool;
[[nodiscard]] auto get_nucstring_z(const std::string& strnuc) -> int;
[[nodiscard]] auto get_nucstring_a(const std::string& strnuc) -> int;
[[gnu::pure]] [[nodiscard]] auto get_num_nuclides() -> ptrdiff_t;
[[nodiscard]] auto get_elname(int z) -> std::string;
[[nodiscard]] auto get_nuc_z(int nucindex) -> int;
[[nodiscard]] auto get_nuc_a(int nucindex) -> int;
[[nodiscard]] auto get_nucindex(int z, int a) -> int;
[[nodiscard]] auto nuc_exists(int z, int a) -> bool;
[[nodiscard]] auto nucdecayenergygamma(int nucindex) -> double;
[[nodiscard]] auto nucdecayenergygamma(int z, int a) -> double;
[[nodiscard]] auto get_decay_neutrino_frac(int nucindex, DecayType decaytype) -> double;
void set_nucdecayenergygamma(int nucindex, double value);
void update_abundances(int nonemptymgi, double t_current);
[[nodiscard]] auto get_endecay_per_ejectamass_tmodel_to_time_withexpansion(int nonemptymgi, double tstart) -> double;
[[nodiscard]] auto get_modelcell_simtime_endecay_per_mass(int nonemptymgi,
                                                          std::span<const double> energy_per_mass_nonemptymgi_decaypath)
    -> double;
auto get_energy_per_mass_nonemptymgi_decaypath() -> MPI_shared_array<const double>;
[[nodiscard]] auto get_qdot_modelcell(int nonemptymgi, double t, DecayType decaytype) -> double;
[[nodiscard]] auto get_particle_injection_rate(int nonemptymgi, double t, DecayType decaytype) -> double;
[[nodiscard]] auto get_gamma_emission_rate(int nonemptymgi, double t) -> double;
[[nodiscard]] auto get_global_etot_tmodel_tinf() -> double;
void output_nuc_abundances(std::ostream& estimators_file, int nonemptymgi, double t_current, int element);
void setup_radioactive_pellet(double e_cmf_per_packet, int nonemptymgi, Packet& pkt,
                              std::span<const double> energy_per_mass_nonemptymgi_decaypath);

}  // namespace decay

#endif  // DECAY_H

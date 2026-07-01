#ifndef PACKET_H
#define PACKET_H

#include <cmath>
#include <span>
#include <string>
#include <vector>

#include "constants.h"

enum packet_type : int {
  TYPE_NONE = 0,
  TYPE_ESCAPE = 32,
  TYPE_RADIOACTIVE_PELLET = 100,
  TYPE_GAMMA = 10,
  TYPE_RPKT = 11,
  TYPE_KPKT = 12,
  TYPE_MA = 13,
  TYPE_NTLEPTON_DEPOSITED = 20,
  TYPE_NONTHERMAL_PREDEPOSIT_BETAMINUS = 21,
  TYPE_NONTHERMAL_PREDEPOSIT_BETAPLUS = 22,
  TYPE_NONTHERMAL_PREDEPOSIT_ALPHA = 23,
  TYPE_NTALPHA_FISPROD_DEPOSITED = 24,
  TYPE_PRE_KPKT = 120,
};

constexpr int EMTYPE_NOTSET{-9999000};
constexpr int EMTYPE_FREEFREE{-9999999};

struct MacroAtomState {
  int element;  // macro atom of type element (this is an element index)
  int ion;  // in ionstage ion (this is an ion index)
  int level;  // and level=level (this is a level index)
  int activatingline;  // Linelistindex of the activating line for bb activated MAs, -99 else.
};

#ifdef GPU_ON
#include "random.h"
#endif

struct Packet {
  double prop_time{-1.};  // internal clock to track how far in time the packet has been propagated
  Vec3d pos{};  // Position of the packet (x,y,z).
  Vec3d dir{};  // Direction of propagation. (x,y,z). Always a unit vector.
  double nu_cmf{0.};  // The frequency in the co-moving frame.
  double e_cmf{0.};  // The energy the packet carries in the co-moving frame.
  double nu_rf{0.};  // The frequency in the rest frame.
  double e_rf{0.};  // The energy the packet carries in the rest frame.
  int next_trans{-1};  // This keeps track of the next possible line interaction of a rpkt by storing
                       // its linelist index (to overcome numerical problems in propagating the rpkts).
  int nscatterings{0};  // records number of electron scatterings a r-pkt undergone since it was emitted
  int emissiontype{EMTYPE_NOTSET};  // records how the packet was emitted if it is a r-pkt
  Vec3d em_pos{NAN};  // Position of the last emission (x,y,z).
  float em_time{-1.};
  int absorptiontype{0};  // records linelistindex of the last absorption
                          // negative values give ff-abs (-1), bf-abs (-2), compton scattering of gammas (-3),
                          // photoelectric effect of gammas (-4), pair production of gammas (-5)
                          // decaying pellets of the 52Fe chain (-6) and pellets which decayed before the
                          // onset of the simulation (-7)
                          // decay of a positron pellet (-10)
  double absorptionfreq{};  // records nu_rf of packet at last absorption
  double stokes_q{0.};  // normalised Stokes q = Q/I
  double stokes_u{0.};  // normalised Stokes u = U/I
  int trueemissiontype = EMTYPE_NOTSET;  // emission type coming from a kpkt to rpkt (last thermal emission)
  Vec3d trueem_pos{NAN, NAN, NAN};
  float trueem_time{-1.};  // last thermal emission time [s]
  enum packet_type type {};  // type of packet (k-, r-, etc.)
  int cellindex{-1};  // The propagation grid cell that the packet is in.
  enum packet_type escape_type {};  // In which form when escaped from the grid.
  float escape_time{-1};  // time at which is passes out of the grid [s]
  double tdecay{-1.};  // Time at which pellet decays
  int number{-1};  // A unique number to identify the packet
  bool originated_from_particlenotgamma{false};  // first packet type after pellet decay
  int pellet_decaytype{-1};  // index into decay::decaytypes
  int pellet_nucindex{-1};  // nuclide index of the decaying species

#ifdef GPU_ON
  // per-packet RNG state so that GPU threads (which can't use thread_local)
  // don't share and race on a single global generator
  utlrandom::Xoshiro128PP rngstate{};
#endif

  auto operator<=>(const Packet& rhs) const = default;
};

#ifdef GPU_ON
constexpr DEVICE_FUNC auto get_rngstate([[maybe_unused]] Packet& packet) -> utlrandom::Xoshiro128PP& {
  return packet.rngstate;
}
#else
#include <random>

constexpr auto get_rngstate() -> std::mt19937& {
  thread_local std::mt19937 rng{std::random_device{}()};
  return rng;
}
constexpr auto get_rngstate([[maybe_unused]] const Packet& packet) -> std::mt19937& { return get_rngstate(); }
#endif

void packet_init(std::span<Packet> packets);
auto read_text_packets(const std::string& filename) -> std::vector<Packet>;
void write_text_packets(const std::string& filename, std::span<const Packet> packets);
void read_temp_packetsfile(int timestep, int my_rank, std::vector<Packet>& packets);
void write_temp_packetsfile(int timestep, int my_rank, std::span<const Packet> packets);

#endif  // PACKET_H

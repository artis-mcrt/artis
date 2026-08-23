// The Packet struct holding the full state of a Monte Carlo energy packet, and the enums for
// the packet types and interaction/emission types.

#ifndef PACKET_H
#define PACKET_H

#include <cmath>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "constants.h"

// Packet state in the indivisible energy packet scheme of Lucy (2002), doi:10.1051/0004-6361:20011756.
// do_packet() dispatches on this. Every packet starts as TYPE_RADIOACTIVE_PELLET; those that reach the grid
// surface end as TYPE_ESCAPE, while packets still in flight when the run ends keep whatever type they held,
// which is why exspec filters on TYPE_ESCAPE:
//
//   RADIOACTIVE_PELLET --(decay to gamma rays)--> GAMMA
//                      --(decay to a lepton/alpha)--> NONTHERMAL_PREDEPOSIT_{BETAMINUS,BETAPLUS,ALPHA}
//                      --(spontaneous fission)--> NTALPHA_FISPROD_DEPOSITED
//                      --(decay with no gamma spectrum at all, e.g. the 52Fe chain)--> KPKT
//                      --(decayed before tmin, or carrying the model's initial energy)--> PRE_KPKT
//   GAMMA --(Compton/photoelectric/pair production)--> NTLEPTON_DEPOSITED or a PREDEPOSIT type
//         --(leaves the grid)--> ESCAPE
//   NONTHERMAL_PREDEPOSIT_* --(thermalises)--> NTLEPTON_DEPOSITED / NTALPHA_FISPROD_DEPOSITED
//                           --(Barnes/Wollaeger schemes, no deposition)--> ESCAPE
//   NTLEPTON_DEPOSITED --(heating, or non-thermal ionisation/excitation)--> KPKT or RPKT
//   NTALPHA_FISPROD_DEPOSITED --(all to heating)--> KPKT
//   KPKT --(free-free, free-bound, or collisional emission)--> RPKT
//   RPKT --(free-free or bound-free absorption)--> KPKT
//        --(bound-bound absorption activates a macro-atom, which deactivates either
//          radiatively, the usual line scattering/fluorescence, or collisionally)--> RPKT or KPKT
//        --(leaves the grid)--> ESCAPE
//
// The values are written to packets*.out as type_id and escape_type_id and parsed by artistools, so they
// must not be renumbered.
enum packet_type : int {
  TYPE_NONE = 0,  // zero-init default; also the escape_type written for a packet that never escaped
  TYPE_GAMMA = 10,
  TYPE_RPKT = 11,

  // Normally destroyed by sampling a cooling channel, but do_kpkt() advances the packet by a small fraction of
  // the timestep first, so a k-packet can also survive to the next timestep. In optically thick cells, and
  // whenever RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY is set, that selection is skipped and
  // do_kpkt_blackbody() re-emits immediately: from a Planck function weighted by the expansion opacity when the
  // option is set and the cell is not thick, and from a plain Planck function otherwise.
  TYPE_KPKT = 12,

  // Never stored in Packet::type: do_macroatom() runs to deactivation within one call. Used only as a
  // provenance tag to vpkt::trace_vpkts(), marking an emission as a macro-atom deactivation.
  TYPE_MA = 13,

  TYPE_NTLEPTON_DEPOSITED = 20,  // awaiting partition into heating/ionisation/excitation by nonthermal.cc

  // Fast particles that have not yet given up their energy. PARTICLE_THERMALISATION_SCHEME decides whether
  // and when they do, so these states can persist across timesteps (time-dependent schemes) or escape
  // without depositing at all (Barnes/Wollaeger).
  TYPE_NONTHERMAL_PREDEPOSIT_BETAMINUS = 21,
  TYPE_NONTHERMAL_PREDEPOSIT_BETAPLUS = 22,
  TYPE_NONTHERMAL_PREDEPOSIT_ALPHA = 23,

  // All of this energy goes to heating, so it converts straight to a k-packet with no Spencer-Fano solve.
  TYPE_NTALPHA_FISPROD_DEPOSITED = 24,

  TYPE_ESCAPE = 32,  // inactive; binned into spectra and light curves by escape_type and escape_time
  TYPE_RADIOACTIVE_PELLET = 100,  // decays at Packet::tdecay, moving with the homologous flow until then

  // Energy released at or before tmin, re-emitted at tmin as a blackbody r-packet: pellets that decayed before
  // the simulation started (update_pellet() scales their energy for the work done on the ejecta in between),
  // and, under INITIAL_PACKETS_ON, the model's own initial thermal energy, which
  // enters as pellets with tdecay == tmin and so is not rescaled.
  TYPE_PRE_KPKT = 120,
};

// Sentinels for Packet::emissiontype: other values are a linelist index if non-negative, else a bound-free
// continuum encoded by get_emtype_continuum() (see get_bflistindex_from_emtype_continuum() in atomic.h).
constexpr int EMTYPE_NOTSET{-9999000};
constexpr int EMTYPE_FREEFREE{-9999999};

// negative absorptiontype values; non-negative values are linelist indices of bound-bound absorption
enum absorption_type : int {
  ABSTYPE_FREEFREE = -1,
  ABSTYPE_BOUNDFREE = -2,
  ABSTYPE_GAMMA_COMPTON = -3,
  ABSTYPE_GAMMA_PHOTOELECTRIC = -4,
  ABSTYPE_GAMMA_PAIRPRODUCTION = -5,
  ABSTYPE_PELLET_NOGAMMASPEC = -6,  // pellet decay with no known gamma spectrum (e.g. 52Fe chain)
  ABSTYPE_PELLET_BEFORESIMSTART = -7,  // pellet decayed before the onset of the simulation
  ABSTYPE_PELLET_PARTICLEDECAY = -10,  // pellet decay to non-thermal particle (beta+/-, alpha, fission fragment)
};

// The state a macro-atom is activated in. Local to a do_macroatom() call, not part of the packet's
// persistent state, since a macro-atom always runs to deactivation before the packet is handed back.
struct MacroAtomState {
  int element{-1};  // macro atom of type element (this is an element index)
  int ion{-1};  // in ionstage ion (this is an ion index)
  int level{-1};  // and level=level (this is a level index)
  // Linelist index of the activating line for bb activated MAs, -99 else. Does not affect which transition the
  // macro-atom samples, but it is not diagnostic-only: rpkt.cc copies it into Packet::absorptiontype, which is
  // the key for the bound-bound decomposition in absorption.out. It also feeds the resonance and
  // up/down-scattering counters and the macroatom.out activline column.
  int activatingline{-99};
};

#include "random.h"

struct Packet {
#ifdef GPU_ON
  // per-packet RNG state so that GPU threads (which can't use thread_local)
  // don't share and race on a single global generator
  rngstate_type rngstate{};
#endif
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

  // The process of the MOST RECENT emission, one of the two keys exspec decomposes the spectra by (see
  // trueemissiontype below). Overwritten by each emission rather than cleared when the packet re-enters the
  // thermal pool, except under RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY, where the thermal frequency
  // redistribution resets it to EMTYPE_NOTSET.
  int emissiontype{EMTYPE_NOTSET};
  Vec3d em_pos{NAN, NAN, NAN};  // Position of the last emission (x,y,z).
  float em_time{-1.};  // [s]
  int absorptiontype{0};  // records linelistindex of the last absorption
                          // or a negative absorption_type enum value
  double absorptionfreq{};  // records nu_rf of packet at last absorption
  double stokes_q{0.};  // normalised Stokes q = Q/I
  double stokes_u{0.};  // normalised Stokes u = U/I
  // The last emission out of the THERMAL POOL. Set when a k-packet emits and then carried unchanged through
  // subsequent scatterings and macro-atom deactivations, so it attributes escaped energy to where it was
  // thermalised rather than to the last surface it scattered off. Reset to EMTYPE_NOTSET at every site that
  // returns the packet to the thermal pool, so the next radiative emission starts a fresh record.
  int trueemissiontype = EMTYPE_NOTSET;
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

  auto operator<=>(const Packet& rhs) const = default;
};

#ifdef GPU_ON
constexpr DEVICE_FUNC auto get_rngstate([[maybe_unused]] Packet& packet) -> rngstate_type& { return packet.rngstate; }
#else
inline auto get_rngstate() -> rngstate_type& {
  // Every thread lazily seeds its own generator from a random source, so that OpenMP/stdpar worker
  // threads (which never run the seeding code in read_parameterfile) do not all share the identical
  // default-seeded sequence. The main thread is re-seeded deterministically in read_parameterfile()
  // to keep single-threaded runs reproducible.
  thread_local rngstate_type rng{static_cast<std::uint32_t>(get_rng_random_seed())};
  return rng;
}

inline auto get_rngstate([[maybe_unused]] const Packet& packet) -> rngstate_type& { return get_rngstate(); }
#endif

void packet_init(std::span<Packet> packets);
auto read_text_packets(const std::string& filename) -> std::vector<Packet>;
void write_text_packets(const std::string& filename, std::span<const Packet> packets);
void read_temp_packetsfile(int timestep, int my_rank, std::vector<Packet>& packets);
void write_temp_packetsfile(int timestep, int my_rank, std::span<const Packet> packets);

#endif  // PACKET_H

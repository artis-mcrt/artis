// The Packet struct holding the full state of a Monte Carlo energy packet, and the enums for
// the packet types and interaction/emission types.

#ifndef PACKET_H
#define PACKET_H

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include "artisoptions.h"
#include "constants.h"
#include "mpi_logging.h"

// Packet state in the indivisible energy packet scheme of Lucy (2002), A&A, 384, 725-735,
// doi:10.1051/0004-6361:20011756.
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
  // provenance tag to vpkt::trace_vpkts(), marking an emission as a macro-atom deactivation or, with
  // RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY, as a line scattering.
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
  // bound-bound absorption in a binned expansion opacity (RPKT_USE_EXPANSION_OPACITIES with
  // RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY), so no single line index is known
  ABSTYPE_BOUNDBOUND_EXPANSIONOPACITY = -11,
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

// The process of a sampled interaction (see Packet::sampled_interactions). An interaction is each event that sets
// Packet::em_time. The values are written to packets*.out as sampled<slot>_interactiontype and parsed by artistools,
// so they must not be renumbered.
enum interaction_type : int {
  INTERACTION_NONE = 0,  // an empty slot: the packet had fewer interactions than SAMPLED_INTERACTIONS_PER_PACKET
  INTERACTION_PELLET_PARTICLE_DECAY = 1,  // a pellet decays to a non-thermal particle
  INTERACTION_PELLET_GAMMA_DECAY = 2,  // a pellet decays to a gamma packet
  INTERACTION_PAIR_ANNIHILATION_GAMMA = 3,  // a pair production gives a 511 keV gamma packet
  INTERACTION_COMPTON_SCATTERING = 4,  // a gamma packet stays a gamma packet after a Compton scattering
  INTERACTION_KPKT_EMISSION = 5,  // a k-packet emits an r-packet (free-free, free-bound, or blackbody)
  INTERACTION_MACROATOM_BOUNDBOUND_EMISSION = 6,  // a radiative bound-bound deactivation of a macro-atom
  INTERACTION_MACROATOM_BOUNDFREE_EMISSION = 7,  // a radiative recombination of a macro-atom
  INTERACTION_ELECTRON_SCATTERING = 8,  // an electron scattering of an r-packet in a cell that is not thick
  INTERACTION_THICKCELL_GREY_SCATTERING = 9,  // a grey event of an r-packet in a thick cell
  // With RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY, a bound-bound event either scatters the r-packet at the same
  // comoving frequency or redistributes the frequency thermally.
  INTERACTION_BOUNDBOUND_SCATTERING = 10,
  INTERACTION_BOUNDBOUND_THERMALISATION = 11,
};

static_assert(SAMPLED_INTERACTIONS_PER_PACKET >= 0);

// The state of a packet directly after one interaction.
struct SampledInteraction {
  Vec3d pos{NAN, NAN, NAN};  // position of the interaction (x,y,z)
  double absorptionfreq{};  // Packet::absorptionfreq at the interaction
  float time{-1.};  // time of the interaction [s]
  enum interaction_type type { INTERACTION_NONE };
  int emissiontype{EMTYPE_NOTSET};  // Packet::emissiontype directly after the interaction
  int absorptiontype{0};  // Packet::absorptiontype at the interaction

  auto operator<=>(const SampledInteraction& rhs) const = default;
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
  // The number of electron scatterings of an r-packet since its last emission. A grey event in a thick cell also
  // adds one, because the code treats it as a coherent scattering. The grey opacity (see RPKT_GREY_TYPE) includes
  // the line opacity, so a grey event can also be a line interaction.
  int nscatterings{0};

  // The process of the MOST RECENT emission, one of the two keys exspec decomposes the spectra by (see
  // trueemissiontype below). Overwritten by each emission rather than cleared when the packet re-enters the
  // thermal pool, except under RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY, where the thermal frequency
  // redistribution resets it to EMTYPE_NOTSET.
  int emissiontype{EMTYPE_NOTSET};
  // Position of the last emission (x,y,z). A scattering also sets it: electron scattering of an r-packet, and
  // Compton scattering of a gamma packet.
  Vec3d em_pos{NAN, NAN, NAN};
  float em_time{-1.};  // [s]
  int absorptiontype{0};  // records linelistindex of the last absorption
                          // or a negative absorption_type enum value
  // nu_rf of the r-packet at its last absorption. A gamma-ray absorption and a pellet decay come before the first
  // r-packet absorption, so this value is 0 for them.
  double absorptionfreq{};
  double stokes_q{0.};  // normalised Stokes q = Q/I
  double stokes_u{0.};  // normalised Stokes u = U/I
  // The last emission out of the THERMAL POOL. A k-packet emission sets it. Scatterings and macro-atom
  // deactivations keep it, so it gives the escaped energy to the place of thermalisation, not to the last
  // scattering. Each site that hands the packet from the thermal pool to a macro-atom sets it to EMTYPE_NOTSET,
  // so the next radiative emission starts a fresh record.
  int trueemissiontype = EMTYPE_NOTSET;
  Vec3d trueem_pos{NAN, NAN, NAN};
  float trueem_time{-1.};  // last thermal emission time [s]
  enum packet_type type {};  // type of packet (k-, r-, etc.)
  int cellindex{-1};  // The propagation cell that the packet is in.
  enum packet_type escape_type {};  // In which form when escaped from the grid.
  float escape_time{-1};  // time at which is passes out of the grid [s]
  double tdecay{-1.};  // Time at which pellet decays
  int number{-1};  // A unique number to identify the packet
  bool originated_from_particlenotgamma{false};  // first packet type after pellet decay
  int pellet_decaytype{-1};  // decay::DecayType value of the pellet decay, or -1 for the initial-energy channel
  int pellet_nucindex{-1};  // nuclide index of the decaying species
  // The number of interactions of the packet since packet_init(). The count stays 0 if
  // SAMPLED_INTERACTIONS_PER_PACKET is 0.
  int ninteractions{0};
  // A uniform random sample without replacement of all the interactions of the packet. Each interaction of the
  // packet has the same probability min(1, SAMPLED_INTERACTIONS_PER_PACKET / ninteractions) to be in the sample.
  // A slot of the sample thus represents ninteractions / min(ninteractions, SAMPLED_INTERACTIONS_PER_PACKET)
  // interactions of the packet. The slots have no time order.
  std::array<SampledInteraction, SAMPLED_INTERACTIONS_PER_PACKET> sampled_interactions{};

  auto operator<=>(const Packet& rhs) const = default;
};

#ifdef GPU_ON
constexpr DEVICE_FUNC auto get_rngstate([[maybe_unused]] Packet& packet) -> rngstate_type& { return packet.rngstate; }
#else
inline auto get_rngstate() -> rngstate_type& {
  // Every thread lazily seeds its own generator from a random source, so that OpenMP/stdpar worker
  // threads (which never run the seeding code in read_parameterfile) do not all share the identical
  // default-seeded sequence. The main thread is re-seeded deterministically in read_parameterfile(), and a
  // resumed job restores its state in read_packet_restart_file(), to keep single-threaded runs reproducible.
  thread_local rngstate_type rng{static_cast<std::uint32_t>(get_rng_random_seed())};
  return rng;
}

inline auto get_rngstate([[maybe_unused]] const Packet& packet) -> rngstate_type& { return get_rngstate(); }
#endif

// Select the slot of a uniform sample of nslots items that item number nitems_seen (1-based) of a sequence goes into,
// or give -1 if the sample does not keep the item. After the last item, each item of the sequence has the same
// probability min(1, nslots / nitems_seen) to be in the sample. This is algorithm R of the reservoir sampling
// (Vitter, J. S. 1985, ACM Transactions on Mathematical Software, 11, 37-57, doi:10.1145/3147.3165).
//
// The uniform random integer comes from a hash of random_key and nitems_seen, and not from the random generator of
// the packets. The sample thus leaves the random sequence of the physics unchanged, and the sample of a packet is
// reproducible in a multithreaded run. The hash is the SplitMix64 finaliser (Steele, G. L., Lea, D., & Flood, C. H.
// 2014, Proceedings of the 2014 ACM International Conference on Object Oriented Programming Systems Languages &
// Applications, 453-472, doi:10.1145/2660193.2660195). The multiplication and the shift map the hash to the range
// [0, nitems_seen) (Lemire, D. 2019, ACM Transactions on Modeling and Computer Simulation, 29, 3:1-3:12,
// doi:10.1145/3230636). The bias of that map is less than nitems_seen / 2^32.
[[nodiscard]] constexpr DEVICE_FUNC auto get_reservoir_sample_slot(const int nitems_seen, const int nslots,
                                                                   const std::uint64_t random_key) -> int {
  if (nitems_seen <= nslots) {
    return nitems_seen - 1;
  }
  const auto splitmix64_finaliser = [](std::uint64_t state) {
    state += 0x9E3779B97F4A7C15ULL;
    state = (state ^ (state >> 30U)) * 0xBF58476D1CE4E5B9ULL;
    state = (state ^ (state >> 27U)) * 0x94D049BB133111EBULL;
    return state ^ (state >> 31U);
  };
  const std::uint64_t hash =
      splitmix64_finaliser(random_key ^ splitmix64_finaliser(static_cast<std::uint64_t>(nitems_seen)));
  const auto random_index = static_cast<int>(((hash >> 32U) * static_cast<std::uint64_t>(nitems_seen)) >> 32U);
  return (random_index < nslots) ? random_index : -1;
}

// Count the current interaction of the packet, and add it to Packet::sampled_interactions with the probability
// that keeps the sample uniform. Call this after the interaction has set emissiontype, absorptiontype, and
// absorptionfreq.
DEVICE_FUNC inline void sample_interaction(Packet& pkt, const enum interaction_type type) {
  if constexpr (SAMPLED_INTERACTIONS_PER_PACKET > 0) {
    assert_testmodeonly(pkt.ninteractions < std::numeric_limits<int>::max());
    pkt.ninteractions++;
    // Packet::number is unique only inside one rank, so the key also contains the rank
    const std::uint64_t random_key = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(globals::my_rank)) << 32U) |
                                     static_cast<std::uint32_t>(pkt.number);
    const int slot = get_reservoir_sample_slot(pkt.ninteractions, SAMPLED_INTERACTIONS_PER_PACKET, random_key);
    if (slot >= 0) {
      pkt.sampled_interactions[slot] = SampledInteraction{
          .pos = pkt.pos,
          .absorptionfreq = pkt.absorptionfreq,
          .time = static_cast<float>(pkt.prop_time),
          .type = type,
          .emissiontype = pkt.emissiontype,
          .absorptiontype = pkt.absorptiontype,
      };
    }
  }
}

void packet_init(std::span<Packet> packets);
auto read_text_packets(const std::string& filename) -> std::vector<Packet>;
void write_text_packets(const std::string& filename, std::span<const Packet> packets);
void read_packet_restart_file(int timestep, std::vector<Packet>& packets);
void write_packet_restart_file(int timestep, std::span<const Packet> packets);

#endif  // PACKET_H

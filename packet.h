#ifndef PACKET_H
#define PACKET_H

#include <compare>
#include <cstdio>

enum packet_type {
  TYPE_ESCAPE = 32,
  TYPE_RADIOACTIVE_PELLET = 100,
  TYPE_GAMMA = 10,
  TYPE_RPKT = 11,
  TYPE_KPKT = 12,
  TYPE_MA = 13,
  TYPE_NTLEPTON = 20,
  TYPE_NONTHERMAL_PREDEPOSIT = 21,
  TYPE_PRE_KPKT = 120,
};

constexpr int EMTYPE_NOTSET = -9999000;
constexpr int EMTYPE_FREEFREE = -9999999;

struct mastate {
  int element;         /// macro atom of type element (this is an element index)
  int ion;             /// in ionstage ion (this is an ion index)
  int level;           /// and level=level (this is a level index)
  int activatingline;  /// Linelistindex of the activating line for bb activated MAs, -99 else.
};

enum cell_boundary {
  COORD0_MIN = 101,
  COORD0_MAX = 102,
  COORD1_MIN = 103,
  COORD1_MAX = 104,
  COORD2_MIN = 105,
  COORD2_MAX = 106,
  BOUNDARY_NONE = 107,
};

struct packet {
  int where = -1;                                 // The propagation grid cell that the packet is in.
  enum packet_type type;                          // type of packet (k-, r-, etc.)
  enum cell_boundary last_cross = BOUNDARY_NONE;  // To avoid rounding errors on cell crossing.
  int interactions = 0;                           // number of interactions the packet undergone
  int nscatterings = 0;              // records number of electron scatterings a r-pkt undergone since it was emitted
  int last_event;                    // debug: stores information about the packets history
  double pos[3] = {0};               // Position of the packet (x,y,z).
  double dir[3] = {0};               // Direction of propagation. (x,y,z). Always a unit vector.
  double e_cmf;                      // The energy the packet carries in the co-moving frame.
  double e_rf;                       // The energy the packet carries in the rest frame.
  double nu_cmf;                     // The frequency in the co-moving frame.
  double nu_rf;                      // The frequency in the rest frame.
  int next_trans;                    // This keeps track of the next possible line interaction of a rpkt by storing
                                     // its linelist index (to overcome numerical problems in propagating the rpkts).
  int emissiontype = EMTYPE_NOTSET;  // records how the packet was emitted if it is a r-pkt
  double em_pos[3];                  // Position of the last emission (x,y,z).
  float em_time = -1.;
  double prop_time = -1.;  // internal clock to track how far in time the packet has been propagated
  int absorptiontype;      // records linelistindex of the last absorption
                           // negative values give ff-abs (-1), bf-abs (-2), compton scattering of gammas (-3),
                           // photoelectric effect of gammas (-4), pair production of gammas (-5)
                           // decaying pellets of the 52Fe chain (-6) and pellets which decayed before the
                           // onset of the simulation (-7)
                           // decay of a positron pellet (-10)
  int trueemissiontype = EMTYPE_NOTSET;  // emission type coming from a kpkt to rpkt (last thermal emission)
  float trueem_time = -1.;               // first thermal emission time [s]
  double absorptionfreq;                 // records nu_rf of packet at last absorption
  double absorptiondir[3] = {0.};  // Direction of propagation (x,y,z) when a packet was last absorbed in a line. Always
                                   // a unit vector.
  double stokes[3] = {0.};         // I, Q and U Stokes parameters
  double pol_dir[3] = {0.};      // unit vector which defines the coordinate system against which Q and U are measured;
                                 // should always be perpendicular to dir
  double tdecay = -1.;           // Time at which pellet decays
  enum packet_type escape_type;  // Flag to tell us in which form it escaped from the grid.
  float escape_time = -1;        // time at which is passes out of the grid [s]
  int number = -1;               // A unique number to identify the packet
  bool originated_from_particlenotgamma;  // first-non-pellet packet type was gamma
  int pellet_decaytype = -1;              // index into decay::decaytypes
  int pellet_nucindex = -1;               // nuclide index of the decaying species
  float trueemissionvelocity = -1;
  struct mastate mastate;

  inline bool operator==(const packet &rhs) {
    return (number == rhs.number && type == rhs.type &&
            (em_pos[0] == rhs.em_pos[0] && em_pos[1] == rhs.em_pos[1] && em_pos[2] == rhs.em_pos[2]) &&
            nu_cmf == rhs.nu_cmf && where == rhs.where && prop_time == rhs.prop_time &&
            mastate.activatingline == rhs.mastate.activatingline && tdecay == rhs.tdecay &&
            pellet_nucindex == rhs.pellet_nucindex);
  }
};

void packet_init(struct packet *pkt);
void write_packets(const char filename[], const struct packet *pkt);
void read_packets(const char filename[], struct packet *pkt);
void read_temp_packetsfile(int timestep, int my_rank, struct packet *pkt);
bool verify_temp_packetsfile(int timestep, int my_rank, const struct packet *pkt);

#endif  // PACKET_H

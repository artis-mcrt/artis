#ifndef PACKET_H
#define PACKET_H

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
  TYPE_GAMMA_KPKT = 121,
};

#include "boundary.h"
#include "macroatom.h"

struct packet {
  int where;                      /// The grid cell that the packet is in.
  enum packet_type type;          /// Identifies the type of packet (k-, r-, etc.)
  enum cell_boundary last_cross;  /// To avoid rounding errors on cell crossing.
  int interactions;               /// debug: number of interactions the packet undergone
  int nscatterings;               /// records number of electron scatterings a r-pkt undergone since it was emitted
  int last_event;                 /// debug: stores information about the packets history
  double pos[3];                  /// Position of the packet (x,y,z).
  double dir[3];                  /// Direction of propagation. (x,y,z). Always a unit vector.
  double e_cmf;                   /// The energy the packet carries in the co-moving frame.
  double e_rf;                    /// The energy the packet carries in the rest frame.
  double nu_cmf;                  /// The frequency in the co-moving frame.
  double nu_rf;                   /// The frequency in the rest frame.
  int next_trans;                 /// This keeps track of the next possible line interaction of a rpkt by storing
                                  /// its linelist index (to overcome numerical problems in propagating the rpkts).
  int emissiontype;               /// records how the packet was emitted if it is a r-pkt
  double em_pos[3];               /// Position of the packet (x,y,z).
  int em_time;
  double prop_time;      // internal clock to track how far in time the packet has been propagated
  int absorptiontype;    /// records linelistindex of the last absorption
                         /// negative values give ff-abs (-1), bf-abs (-2), compton scattering of gammas (-3),
                         /// photoelectric effect of gammas (-4), pair production of gammas (-5)
                         /// decaying pellets of the 52Fe chain (-6) and pellets which decayed before the
                         /// onset of the simulation (-7)
                         /// decay of a positron pellet (-10)
  int trueemissiontype;  // emission type coming from a kpkt to rpkt (last thermal emission)
  int trueem_time;
  double absorptionfreq;    /// records nu_cmf of packet at last absorption
  double absorptiondir[3];  /// Direction of propagation (x,y,z) when a packet was last absorbed in a line. Always a
                            /// unit vector.
  // short timestep;
  double stokes[3];   // I, Q and U Stokes parameters
  double pol_dir[3];  // unit vector which defines the coordinate system against which Q and U are measured; should
                      // always be perpendicular to dir
  double tdecay;                          /// Time at which pellet decays.
  enum packet_type escape_type;           /// Flag to tell us in which form it escaped from the grid.
  int escape_time;                        /// Time at which is passes out of the grid.
                                          /// Pos, dir, where, e_rf, nu_rf should all remain set at the exit point.
  int scat_count;                         /// WHAT'S THAT???
  int number;                             /// A unique number to identify which packet caused potential troubles.
  bool originated_from_particlenotgamma;  // first-non-pellet packet type was gamma
  int pellet_decaytype;                   // index into decay::decaytypes
  int pellet_nucindex;                    // nuclide index of the decaying species
  float trueemissionvelocity;
  struct mastate mastate;
};

void packet_init(int my_rank, struct packet *pkt);
void write_packets(char filename[], struct packet *pkt);
void read_packets(char filename[], struct packet *pkt);
void read_temp_packetsfile(const int timestep, const int my_rank, struct packet *const pkt);

#endif  // PACKET_H

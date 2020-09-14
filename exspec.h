#ifndef EXSPEC_H
#define EXSPEC_H

#include <stdbool.h>

/// Spectrum data structure
#define MNUBINS   1000
#define MABINS    100
#define MTBINS  200

typedef struct
{
  double *absorption;
  double *emission;
  double *trueemission;
} emstat_t;


///Specpol

struct specpol
{
    float lower_freq[MNUBINS];
    float delta_freq[MNUBINS];
    double flux[MNUBINS];
    float lower_time;
    float delta_t;
    emstat_t stat[MNUBINS];
};

typedef struct
{
  double dir[3];  /// Direction of propagation. (x,y,z). Always a unit vector.
  double e_rf;    /// The energy the packet carries in the rest frame.
  double e_cmf;    /// The energy the packet carries in the rest frame.
  double nu_rf;   /// The frequency in the rest frame.
  float arrive_time; /// Time at which is passes out of the grid.
  float arrive_time_cmf; /// Time at which is passes out of the grid.
  int emissiontype;   /// records how the packet was emitted if it is a r-pkt
  int trueemissiontype; /// records how the packet was emitted directly after coming from a kpkt
  int absorptiontype;   /// records how the packet was emitted if it is a r-pkt
  double absorptionfreq;   /// records how the packet was emitted if it is a r-pkt
  double em_pos[3]; /// Position of the packet (x,y,z) at last emission process
  int em_time;
  double stokes[3];
  float trueemissionvelocity;
} EPKT;

extern int nprocs_exspec;
extern bool do_emission_res;
extern int nepkts;

extern double dlognu;

extern double dlogtlc;
extern double dlogtlc_angle;

extern struct specpol stokes_i[MTBINS];
extern struct specpol stokes_q[MTBINS];
extern struct specpol stokes_u[MTBINS];

#endif //EXSPEC_H

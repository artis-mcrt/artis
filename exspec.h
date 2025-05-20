#pragma once
/// Spectrum data structure
#define MNUBINS   1000
#define MABINS    100
#define MTBINS  200

extern int nprocs_exspec;
extern int do_emission_res;
extern int nepkts;

typedef struct 
{
  double *absorption;
  double *emission;
} emstat_t;

struct spec
{
  float lower_freq[MNUBINS];
  float delta_freq[MNUBINS];
  double flux[MNUBINS];
  float lower_time;
  float delta_t;
  //int packets[MNUBINS];
  emstat_t stat[MNUBINS];
} extern spectra[MTBINS]; //, spectra_res[MTBINS][MABINS];

extern double dlognu;

/// Light curve data structure
#define MTLCBINS  200
#define MALCBINS  100
#define MANGLCBINS 100
extern double nu_min, nu_max; //limits on frequency range for gamma spectrum

struct lc
{
  float lower_time;
  float delta_t;
  double lum;
};
extern struct lc light_curve[MTLCBINS], light_curve_cmf[MTLCBINS], light_curve_angle[MTLCBINS][MANGLCBINS];

extern double dlogtlc;
extern double dlogtlc_angle;

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
extern struct specpol stokes_i[MTBINS], stokes_q[MTBINS], stokes_u[MTBINS];

extern struct spec spectra[MTBINS]; //, spectra_res[MTBINS][MABINS];
extern struct lc light_curve[MTLCBINS], light_curve_cmf[MTLCBINS], light_curve_angle[MTLCBINS][MANGLCBINS];
extern struct specpol stokes_i[MTBINS], stokes_q[MTBINS], stokes_u[MTBINS];

typedef struct 
{
  double dir[3];  /// Direction of propagation. (x,y,z). Always a unit vector.
  double e_rf;    /// The energy the packet carries in the rest frame.
  double e_cmf;    /// The energy the packet carries in the rest frame.
  double nu_rf;   /// The frequency in the rest frame.
  float arrive_time; /// Time at which is passes out of the grid.
  float arrive_time_cmf; /// Time at which is passes out of the grid.
  int emissiontype;   /// records how the packet was emitted if it is a r-pkt
  int absorptiontype;   /// records how the packet was emitted if it is a r-pkt
  double absorptionfreq;   /// records how the packet was emitted if it is a r-pkt
  double em_pos[3]; /// Position of the packet (x,y,z) at last emission process
  int em_time;
  double stokes[3];
} EPKT;
extern EPKT *epkts;

#include "exspec.h"

int nprocs_exspec;
int do_emission_res;
int nepkts;

struct spec spectra[MTBINS];

double dlognu;

/// Light curve data structure
#define MTLCBINS  200
#define MALCBINS  100
#define MANGLCBINS 100
double nu_min, nu_max; //limits on frequency range for gamma spectrum

struct lc light_curve[MTLCBINS], light_curve_cmf[MTLCBINS], light_curve_angle[MTLCBINS][MANGLCBINS];

double dlogtlc;
double dlogtlc_angle;


struct specpol stokes_i[MTBINS], stokes_q[MTBINS], stokes_u[MTBINS];

EPKT *epkts;

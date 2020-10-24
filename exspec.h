#ifndef EXSPEC_H
#define EXSPEC_H

#include <stdbool.h>

/// Spectrum data structure
#define MNUBINS   1000
#define MABINS    100
#define MTBINS    400

typedef struct
{
  double *absorption;
  double *emission;
  double *trueemission;
} emstat_t;


///Specpol

struct spec
{
  float lower_freq[MNUBINS];
  float delta_freq[MNUBINS];
  double flux[MNUBINS];
  emstat_t stat[MNUBINS];
};


extern int nprocs_exspec;
extern bool do_emission_res;

extern struct spec *stokes_i;
extern struct spec *stokes_q;
extern struct spec *stokes_u;


#endif //EXSPEC_H

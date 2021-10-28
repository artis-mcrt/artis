#ifndef EXSPEC_H
#define EXSPEC_H

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


struct spec
{
  float lower_freq[MNUBINS];
  float delta_freq[MNUBINS];
  double flux[MNUBINS];
  emstat_t stat[MNUBINS];
  double nu_min;
  double nu_max;
  bool do_emission_res;
};


extern int nprocs_exspec;
extern bool do_emission_res;
extern const bool do_exspec;


#endif //EXSPEC_H

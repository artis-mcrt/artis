#ifndef EXSPEC_H
#define EXSPEC_H

/// Spectrum data structure
#define MNUBINS   1000
#define MABINS    100


struct spec
{
  float *lower_freq;
  float *delta_freq;
  double *flux;
  double *absorption;
  double *emission;
  double *trueemission;
  double nu_min;
  double nu_max;
  bool do_emission_res;
};


extern const bool do_exspec;


#endif //EXSPEC_H

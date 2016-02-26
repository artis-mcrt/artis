#include "exspec.h"
#ifndef SPECTRUM_H
  #define SPECTRUM_H

  int write_spectrum(FILE *spec_file, FILE *emission_file, FILE *absorption_file);
  void init_spectrum(void);
  int gather_spectrum(int depth);
  int add_to_spec(EPKT *pkt_ptr);
  int gather_spectrum_res(int current_abin);
  int add_to_spec_res(EPKT *pkt_ptr, int current_abin);

#endif //SPECTRUM_H

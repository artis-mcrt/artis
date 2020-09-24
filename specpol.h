#ifndef SPECPOL_H
#define SPECPOL_H

#include <stdio.h>
#include "exspec.h"

int write_specpol(FILE *specpol_file, FILE *emissionpol_file, FILE *absorptionpol_file);
void init_specpol(void);
int gather_specpol(int depth);
int add_to_specpol(const EPKT *pkt_ptr);
int gather_specpol_res(int current_abin);
int add_to_specpol_res(const EPKT *pkt_ptr, int current_abin);

#endif //SPECPOL_H

#ifndef SPECPOL_H
#define SPECPOL_H

#include <stdio.h>
#include "exspec.h"

int write_specpol(FILE *specpol_file, FILE *emissionpol_file, FILE *absorptionpol_file);
void init_specpol(void);
int add_to_specpol(const PKT *pkt_ptr);
int add_to_specpol_res(const PKT *pkt_ptr, int current_abin);

#endif //SPECPOL_H

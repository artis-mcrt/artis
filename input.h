#ifndef INPUT_H
#define INPUT_H

#include <stdio.h>

int input(int rank);
void read_atomicdata(void);
void read_phixs_data(void);
void write_processed_modelatom(void);
void read_processed_modelatom(FILE *modelatom);
void read_parameterfile(int rank);
int read_1d_model(void);
int read_2d_model(void);
int read_3d_model(void);
int compare_linelistentry(const void *p1, const void *p2);
int search_groundphixslist(double nu_edge, int *index_in_groundlevelcontestimator,int element, int ion, int level);
void update_parameterfile(int nts);

#endif //INPUT_H

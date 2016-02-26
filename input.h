#ifndef INPUT_H
#define INPUT_H

int input(int rank);
void read_atomicdata();
void read_parameterfile();
int read_1d_model();
int read_2d_model();
int read_3d_model();
int compare_linelistentry(const void *p1, const void *p2);
int search_groundphixslist(double nu_edge, int *index_in_groundlevelcontestimator,int element, int ion, int level);
void update_parameterfile(int nts);

#endif //INPUT_H

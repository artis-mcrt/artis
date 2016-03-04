#ifndef INPUT_H
#define INPUT_H

int input(int rank);
static void read_atomicdata(void);
static void read_phixs_data(void);
static void write_processed_modelatom(void);
static void read_processed_modelatom(FILE *modelatom);
static void read_parameterfile(int rank);
static int read_1d_model(void);
static int read_2d_model(void);
static int read_3d_model(void);
static int compare_linelistentry(const void *p1, const void *p2);
static int search_groundphixslist(double nu_edge, int *index_in_groundlevelcontestimator,int element, int ion, int level);
void update_parameterfile(int nts);

#endif //INPUT_H

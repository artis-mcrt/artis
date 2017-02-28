#ifndef RADFIELD_H
#define RADFIELD_H

#include "sn3d.h"

void radfield_zero_estimators(int modelgridindex);
void radfield_init(int my_rank);
void radfield_write_to_file(int modelgridindex, int timestep);
void radfield_close_file(void);
void radfield_update_binned_estimators(int modelgridindex, double distance_e_cmf, double nu_cmf);
double radfield(double nu, int modelgridindex);
void radfield_fit_parameters(int modelgridindex, int timestep);
void radfield_set_J_normfactor(int modelgridindex, double normfactor);
void radfield_reduce_binned_estimators(void);
void radfield_MPI_Bcast(int root, int my_rank, int nstart, int ndo);
void radfield_write_restart_data(FILE *gridsave_file);
void radfield_read_restart_data(FILE *gridsave_file);


inline double radfield_dbb(double nu, float T, float W)
// returns J_nu for a diluted black body
{
  return W * TWOHOVERCLIGHTSQUARED * pow(nu,3) / expm1(HOVERKB * nu / T);
}


#endif //RADFIELD_H

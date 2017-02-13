#ifndef RADFIELD_H
#define RADFIELD_H

#include "sn3d.h"

void radfield_zero_estimators(int modelgridindex);
void radfield_init(int my_rank);
void radfield_write_to_file(int modelgridindex, int timestep);
void radfield_close_file(void);
void radfield_update_estimators(int modelgridindex, double distance_e_cmf, double nu_cmf);
double radfield(double nu, int modelgridindex);
void radfield_fit_parameters(int modelgridindex, int timestep);
void get_radfield_params_fullspec(double J, double nuJ, int modelgridindex, double *T_J, double *T_R, double *W);
void radfield_set_J_normfactor(int modelgridindex, double normfactor);
void radfield_reduce_estimators(int my_rank);
void radfield_MPI_Bcast(int root, int my_rank, int nstart, int ndo);
void radfield_write_restart_data(FILE *gridsave_file);
void radfield_read_restart_data(FILE *gridsave_file);
int radfield_select_bin(int modelgridindex, double nu);
double radfield_get_bin_nu_lower(int modelgridindex, int binindex);
double radfield_get_bin_nu_upper(int modelgridindex, int binindex);


inline double radfield2(double nu, double T, double W)
// returns J_nu for a diluted black body
{
  return W * TWOHOVERCLIGHTSQUARED * pow(nu,3) / (exp(HOVERKB * nu / T) - 1);
}


#endif //RADFIELD_H

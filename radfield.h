#ifndef RADFIELD_H
#define RADFIELD_H

void zero_radfield_estimators(int modelgridindex);
void init_radfield_file(void);
void write_to_radfield_file(int modelgridindex, int timestep);
void close_radfield_file(void);
void update_radfield(int modelgridindex, double distance, double e_cmf,
                     double nu_cmf);
double radfield(double nu, int modelgridindex);
void radfield_fit_parameters(int modelgridindex);
void radfield_set_J_normfactor(int modelgridindex, double normfactor);
double radfield_get_bin_J(int modelgridindex, int binindex);
int get_frequency_bin(int modelgridindex, double nu);
double get_bin_nu_lower(int modelgridindex, int binindex);

#endif //RADFIELD_H

#ifndef RADFIELD_H
#define RADFIELD_H

void zero_radfield_estimators(int modelgridindex);
void update_radfield(int modelgridindex, double distance, double e_cmf,
                     double nu_cmf);
double radfield(double nu, int modelgridindex);
void fit_radfield_parameters(int modelgridindex);
void radfield_set_J_normfactor(int modelgridindex, double normfactor);
int get_frequency_bin(int modelgridindex, double nu);
double get_bin_nu_lower(int modelgridindex, int binindex);

#endif //RADFIELD_H

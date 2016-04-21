#ifndef NLTEPOP_H
#define NLTEPOP_H

void nlte_pops_element(int element, int modelgridindex, int timestep);
double nlte_pops(int element, int ion, int modelgridindex, int timestep);
int get_nlte_vector_index(int element_in, int ion_in, int level_in);
void get_ion_level_of_nlte_vector_index(int index, int element, int *ion, int *level);
void filter_nlte_matrix(int element, int nlte_dimension, double *rate_matrix, double *balance_vector, const double *pop_norm_factors);
void eliminate_nlte_matrix_rowcol(int index, int gs_index, int nlte_dimension, double *rate_matrix, double *balance_vector);
double get_tot_nion(int modelgridindex);
double get_oneoverw(int element, int ion, int modelgridindex);
double get_mean_binding_energy(int element, int ion);
int read_binding_energies(void);
double nt_ionization_rate(int modelgridindex, int element, int ion);
double superlevel_boltzmann(int modelgridindex, int element, int ion, int level);

double nlte_matrix_rad_deexcitation(int modelgridindex, int lower, double epsilon_trans, int lineindex, double t_current);
double nlte_matrix_rad_excitation(int modelgridindex, int upper, double epsilon_trans, int lineindex, double t_current);

#endif //NLTEPOP_H

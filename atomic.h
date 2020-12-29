#ifndef ATOMIC_H
#define ATOMIC_H

#include <cassert>
#include "sn3d.h"

extern double last_phixs_nuovernuedge; // last photoion cross section point as a factor of nu_edge = last_phixs_nuovernuedge

int get_continuumindex(int element, int ion, int level, int upperionlevel);
int get_phixtargetindex(const int element, const int ion, const int level, const int upperionlevel);
double get_tau_sobolev(int modelgridindex, int lineindex, double t_current);
double get_nntot(int modelgridindex);
bool is_nlte(int element, int ion, int level);
bool level_isinsuperlevel(int element, int ion, int level);
double photoionization_crosssection_fromtable(float *photoion_xs, double nu_edge, double nu);
void set_nelements(const int nelements_in);
int get_nelements(void);
int get_element(int element);
int get_elementindex(int Z);
int get_nions(int element);
int get_ionstage(int element, int ion);
int get_nlevels(int element, int ion);
int get_nlevels_nlte(int element, int ion);
int get_nlevels_groundterm(int element, int ion);
int get_ionisinglevels(int element, int ion);
int get_uniqueionindex(int element, int ion);
void get_ionfromuniqueionindex(int allionsindex, int *element, int *ion);
double epsilon(int element, int ion, int level);
double stat_weight(int element, int ion, int level);
int get_maxrecombininglevel(int element, int ion);
bool ion_has_superlevel(int element, int ion);
int get_ndowntrans(int element, int ion, int level);
int get_nuptrans(int element, int ion, int level);
int get_nphixstargets(int element, int ion, int level);
int get_phixsupperlevel(int element, int ion, int level, int phixstargetindex);
double get_phixsprobability(int element, int ion, int level, int phixstargetindex);
void set_ndowntrans( int element,  int ion,  int level,  int ndowntrans);
void set_nuptrans( int element, int ion,  int level,  int nuptrans);
double einstein_spontaneous_emission(int lineindex);
double osc_strength(int lineindex);
double get_coll_str(int lineindex);
double statw_upper(int lineindex);
double statw_lower(int lineindex);
double photoionization_crosssection(int element, int ion, int level, double nu_edge, double nu);
double get_phixs_threshold(int element, int ion, int level, int phixstargetindex);


#endif //ATOMIC_H

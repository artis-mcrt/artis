#ifndef ATOMIC_H
#define ATOMIC_H

#include <array>

extern double
    last_phixs_nuovernuedge;  // last photoion cross section point as a factor of nu_edge = last_phixs_nuovernuedge
constexpr std::array<const char *, 3> phixsdata_filenames = {"version0ignore", "phixsdata.txt", "phixsdata_v2.txt"};
extern std::array<bool, 3> phixs_file_version_exists;  // first value in this array is not used but exists so the
                                                       // indexes match those of the phixsdata_filenames array

int get_continuumindex_phixstargetindex(int element, int ion, int level, int phixstargetindex);
int get_continuumindex(int element, int ion, int level, int upperionlevel);
int get_phixtargetindex(int element, int ion, int level, int upperionlevel);
double get_tau_sobolev(int modelgridindex, int lineindex, double t_current);
auto get_nnion_tot(int modelgridindex) -> double;
bool is_nlte(int element, int ion, int level);
bool level_isinsuperlevel(int element, int ion, int level);
double photoionization_crosssection_fromtable(const float *photoion_xs, double nu_edge, double nu);
void set_nelements(int nelements_in);
int get_nelements();
int get_atomicnumber(int element);
int get_elementindex(int Z);
int get_includedions();
void update_includedions_maxnions();
int get_max_nions();
int get_nions(int element);
int get_ionstage(int element, int ion);
int get_nlevels(int element, int ion);
int get_nlevels_nlte(int element, int ion);
int get_nlevels_groundterm(int element, int ion);
int get_ionisinglevels(int element, int ion);
int get_uniqueionindex(int element, int ion);
void get_ionfromuniqueionindex(int allionsindex, int *element, int *ion);
int get_uniquelevelindex(int element, int ion, int level);
void get_levelfromuniquelevelindex(int alllevelsindex, int *element, int *ion, int *level);
double epsilon(int element, int ion, int level);
double stat_weight(int element, int ion, int level);
int get_maxrecombininglevel(int element, int ion);
bool ion_has_superlevel(int element, int ion);
int get_ndowntrans(int element, int ion, int level);
int get_nuptrans(int element, int ion, int level);
int get_nphixstargets(int element, int ion, int level);
int get_phixsupperlevel(int element, int ion, int level, int phixstargetindex);
double get_phixsprobability(int element, int ion, int level, int phixstargetindex);
void set_ndowntrans(int element, int ion, int level, int ndowntrans);
void set_nuptrans(int element, int ion, int level, int nuptrans);
double einstein_spontaneous_emission(int lineindex);
double photoionization_crosssection(int element, int ion, int level, double nu_edge, double nu);
double get_phixs_threshold(int element, int ion, int level, int phixstargetindex);

#endif  // ATOMIC_H

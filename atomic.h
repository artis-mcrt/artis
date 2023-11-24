#ifndef ATOMIC_H
#define ATOMIC_H

#include <array>

extern double
    last_phixs_nuovernuedge;  // last photoion cross section point as a factor of nu_edge = last_phixs_nuovernuedge
constexpr std::array<const char *, 3> phixsdata_filenames = {"version0ignore", "phixsdata.txt", "phixsdata_v2.txt"};
extern std::array<bool, 3> phixs_file_version_exists;  // first value in this array is not used but exists so the
                                                       // indexes match those of the phixsdata_filenames array

auto get_continuumindex_phixstargetindex(int element, int ion, int level, int phixstargetindex) -> int;
auto get_continuumindex(int element, int ion, int level, int upperionlevel) -> int;
auto get_phixtargetindex(int element, int ion, int level, int upperionlevel) -> int;
auto get_tau_sobolev(int modelgridindex, int lineindex, double t_current) -> double;
auto get_nnion_tot(int modelgridindex) -> double;
auto is_nlte(int element, int ion, int level) -> bool;
auto level_isinsuperlevel(int element, int ion, int level) -> bool;
auto photoionization_crosssection_fromtable(const float *photoion_xs, double nu_edge, double nu) -> double;
void set_nelements(int nelements_in);
auto get_nelements() -> int;
auto get_atomicnumber(int element) -> int;
auto get_elementindex(int Z) -> int;
auto get_includedions() -> int;
void update_includedions_maxnions();
auto get_max_nions() -> int;
auto get_nions(int element) -> int;
auto get_ionstage(int element, int ion) -> int;
auto get_nlevels(int element, int ion) -> int;
auto get_nlevels_nlte(int element, int ion) -> int;
auto get_nlevels_groundterm(int element, int ion) -> int;
auto get_ionisinglevels(int element, int ion) -> int;
auto get_uniqueionindex(int element, int ion) -> int;
[[nodiscard]] auto get_ionfromuniqueionindex(int allionsindex) -> std::tuple<int, int>;
auto get_uniquelevelindex(int element, int ion, int level) -> int;
void get_levelfromuniquelevelindex(int alllevelsindex, int *element, int *ion, int *level);
auto epsilon(int element, int ion, int level) -> double;
auto stat_weight(int element, int ion, int level) -> double;
auto get_maxrecombininglevel(int element, int ion) -> int;
auto ion_has_superlevel(int element, int ion) -> bool;
auto get_ndowntrans(int element, int ion, int level) -> int;
auto get_nuptrans(int element, int ion, int level) -> int;
auto get_nphixstargets(int element, int ion, int level) -> int;
auto get_phixsupperlevel(int element, int ion, int level, int phixstargetindex) -> int;
auto get_phixsprobability(int element, int ion, int level, int phixstargetindex) -> double;
void set_ndowntrans(int element, int ion, int level, int ndowntrans);
void set_nuptrans(int element, int ion, int level, int nuptrans);
auto einstein_spontaneous_emission(int lineindex) -> double;
auto photoionization_crosssection(int element, int ion, int level, double nu_edge, double nu) -> double;
auto get_phixs_threshold(int element, int ion, int level, int phixstargetindex) -> double;
auto elem_has_nlte_levels(int element) -> bool;
auto elem_has_nlte_levels_search(int element) -> bool;

#endif  // ATOMIC_H

#ifndef ATOMIC_H
#define ATOMIC_H

#include "cuda.h"

extern __managed__ double
    last_phixs_nuovernuedge;  // last photoion cross section point as a factor of nu_edge = last_phixs_nuovernuedge
extern __managed__ int phixs_file_version;

__managed__ static const char *phixsdata_filenames[] = {"version0ignore", "phixsdata.txt", "phixsdata_v2.txt"};

__host__ __device__ int get_continuumindex(int element, int ion, int level, int upperionlevel);
__host__ __device__ int get_phixtargetindex(const int element, const int ion, const int level, const int upperionlevel);
__host__ __device__ double get_tau_sobolev(int modelgridindex, int lineindex, double t_current);
__host__ __device__ double get_nntot(int modelgridindex);
__host__ __device__ bool is_nlte(int element, int ion, int level);
__host__ __device__ bool level_isinsuperlevel(int element, int ion, int level);
__host__ __device__ double photoionization_crosssection_fromtable(float *photoion_xs, double nu_edge, double nu);
__host__ __device__ void set_nelements(const int nelements_in);
__host__ __device__ int get_nelements(void);
__host__ __device__ int get_element(int element);
__host__ __device__ int get_elementindex(int Z);
__host__ __device__ void increase_includedions(int nions);
__host__ __device__ int get_includedions(void);
__host__ __device__ void update_max_nions(const int nions);
__host__ __device__ int get_max_nions(void);
__host__ __device__ int get_nions(int element);
__host__ __device__ int get_ionstage(int element, int ion);
__host__ __device__ int get_nlevels(int element, int ion);
__host__ __device__ int get_nlevels_nlte(int element, int ion);
__host__ __device__ int get_nlevels_groundterm(int element, int ion);
__host__ __device__ int get_ionisinglevels(int element, int ion);
__host__ __device__ int get_uniqueionindex(int element, int ion);
__host__ __device__ void get_ionfromuniqueionindex(int allionsindex, int *element, int *ion);
__host__ __device__ int get_uniquelevelindex(const int element, const int ion, const int level);
__host__ __device__ void get_levelfromuniquelevelindex(const int alllevelsindex, int *element, int *ion, int *level);
__host__ __device__ double epsilon(int element, int ion, int level);
__host__ __device__ double stat_weight(int element, int ion, int level);
__host__ __device__ int get_maxrecombininglevel(int element, int ion);
__host__ __device__ bool ion_has_superlevel(int element, int ion);
__host__ __device__ int get_ndowntrans(int element, int ion, int level);
__host__ __device__ int get_nuptrans(int element, int ion, int level);
__host__ __device__ int get_nphixstargets(int element, int ion, int level);
__host__ __device__ int get_phixsupperlevel(int element, int ion, int level, int phixstargetindex);
__host__ __device__ double get_phixsprobability(int element, int ion, int level, int phixstargetindex);
__host__ __device__ void set_ndowntrans(int element, int ion, int level, int ndowntrans);
__host__ __device__ void set_nuptrans(int element, int ion, int level, int nuptrans);
__host__ __device__ double einstein_spontaneous_emission(int lineindex);
__host__ __device__ double osc_strength(int lineindex);
__host__ __device__ double get_coll_str(int lineindex);
__host__ __device__ double statw_upper(int lineindex);
__host__ __device__ double statw_lower(int lineindex);
__host__ __device__ double photoionization_crosssection(int element, int ion, int level, double nu_edge, double nu);
__host__ __device__ double get_phixs_threshold(int element, int ion, int level, int phixstargetindex);

#endif  // ATOMIC_H

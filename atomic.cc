#include "sn3d.h"
#include "atomic.h"
#include "grid_init.h"
#include "ltepop.h"
#include "update_grid.h"

__managed__ double last_phixs_nuovernuedge; // last photoion cross section point as a factor of nu_edge = last_phixs_nuovernuedge

extern __host__ __device__ inline int get_element(int element);
extern __host__ __device__ inline int get_elementindex(int Z);
extern __host__ __device__ inline int get_nions(int element);
extern __host__ __device__ inline int get_ionstage(int element, int ion);
extern __host__ __device__ inline int get_nlevels(int element, int ion);
extern __host__ __device__ inline int get_nlevels_nlte(int element, int ion);
extern __host__ __device__ inline int get_nlevels_groundterm(int element, int ion);
extern __host__ __device__ inline int get_ionisinglevels(int element, int ion);
extern __host__ __device__ inline int get_uniqueionindex(int element, int ion);
extern __host__ __device__ inline void get_ionfromuniqueionindex(int allionsindex, int *element, int *ion);
extern __host__ __device__ inline double epsilon(int element, int ion, int level);
extern __host__ __device__ inline bool level_isinsuperlevel(int element, int ion, int level);
extern __host__ __device__ inline double photoionization_crosssection_fromtable(float *photoion_xs, double nu_edge, double nu);
extern __host__ __device__ double stat_weight(int element, int ion, int level);

extern inline __host__ __device__ int get_maxrecombininglevel(int element, int ion);
extern inline __host__ __device__ bool ion_has_superlevel(int element, int ion);
extern inline __host__ __device__ int get_ndowntrans(int element, int ion, int level);
extern inline __host__ __device__ int get_nuptrans(int element, int ion, int level);
extern inline __host__ __device__ int get_nphixstargets(int element, int ion, int level);
extern inline __host__ __device__ int get_phixsupperlevel(int element, int ion, int level, int phixstargetindex);
extern inline __host__ __device__ double get_phixsprobability(int element, int ion, int level, int phixstargetindex);
extern inline __host__ __device__ void set_ndowntrans( int element,  int ion,  int level,  int ndowntrans);
extern inline __host__ __device__ void set_nuptrans( int element, int ion,  int level,  int nuptrans);
extern inline __host__ __device__ double einstein_spontaneous_emission(int lineindex);
extern inline __host__ __device__ double osc_strength(int lineindex);
extern inline __host__ __device__ double get_coll_str(int lineindex);
extern inline __host__ __device__ double statw_upper(int lineindex);
extern inline __host__ __device__ double statw_lower(int lineindex);
extern inline __host__ __device__ double photoionization_crosssection(int element, int ion, int level, double nu_edge, double nu);
extern inline __host__ __device__ double get_phixs_threshold(int element, int ion, int level, int phixstargetindex);


__host__ __device__
static int get_continuumindex_phixstargetindex(int element, int ion, int level, int phixstargetindex)
/// Returns the index of the continuum associated to the given level.
{
  return elements[element].ions[ion].levels[level].cont_index - phixstargetindex;
}


__host__ __device__
int get_phixtargetindex(const int element, const int ion, const int level, const int upperionlevel)
{
  for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++)
  {
    if (upperionlevel == get_phixsupperlevel(element, ion, level, phixstargetindex))
      return phixstargetindex;
  }
  assert(false);
  return -1;
}


__host__ __device__
int get_continuumindex(int element, int ion, int level, int upperionlevel)
/// Returns the index of the continuum associated to the given level.
{
  const int phixstargetindex = get_phixtargetindex(element, ion, level, upperionlevel);
  return get_continuumindex_phixstargetindex(element, ion, level, phixstargetindex);
}


__host__ __device__
double get_tau_sobolev(int modelgridindex, int lineindex, double t_current)
{
  const int element = linelist[lineindex].elementindex;
  const int ion = linelist[lineindex].ionindex;
  const int lower = linelist[lineindex].lowerlevelindex;
  const int upper = linelist[lineindex].upperlevelindex;

  const double n_l = calculate_exclevelpop(modelgridindex,element,ion,lower);
  const double n_u = calculate_exclevelpop(modelgridindex,element,ion,upper);

  const double nu_trans = (epsilon(element, ion, upper) - epsilon(element, ion, lower)) / H;
  const double A_ul = einstein_spontaneous_emission(lineindex);
  const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
  const double B_lu = stat_weight(element, ion, upper) / stat_weight(element, ion, lower) * B_ul;

  const double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current;
  return tau_sobolev;
}


__host__ __device__
double get_nntot(int modelgridindex)
// total ion (nuclei) density
{
  const double rho = get_rho(modelgridindex);
  int nntot = 0.;
  for (int element = 0; element < nelements; element++)
  {
    nntot += get_abundance(modelgridindex, element) / elements[element].mass * rho;
  }

  return nntot;
}

__host__ __device__
bool is_nlte(int element, int ion, int level)
// Returns true if (element,ion,level) is to be treated in nlte.
// (note this function returns true for the ground state,
//  although it is stored separately from the excited NLTE states)
{
  if (!NLTE_POPS_ON)
  {
    return false;
  }
  else if (get_element(element) == 26 && get_ionstage(element, ion) == 2)
    return (level <= 197);
  else
    return (level <= 80);
}



#ifndef DECAY_H
#define DECAY_H

// #include <cstdio>

enum decaypathways {
  DECAY_NI56 = 0,
  DECAY_NI56_CO56 = 1,
  DECAY_FE52 = 2,
  DECAY_FE52_MN52 = 3,
  DECAY_CR48 = 4,
  DECAY_CR48_V48 = 5,
  DECAY_CO56 = 6,
  DECAY_NI57 = 7,
  DECAY_NI57_CO57 = 8,
  DECAY_CO57 = 9,
  DECAYPATH_COUNT = 10,
};

#include "types.h"
#include "cuda.h"

namespace decay
{
  __host__ __device__ void init_nuclides(void);
  __host__ __device__ int get_num_nuclides(void);
  __host__ __device__ int get_nuc_z(int nucindex);
  __host__ __device__ int get_nuc_a(int nucindex);
  __host__ __device__ int get_nuc_index(int atomic_number, int mass_number);
  __host__ __device__ double nucdecayenergygamma(int z, int a);
  __host__ __device__ void set_nucdecayenergygamma(int z, int a, double value);
  __host__ __device__ double nucdecayenergypositrons(int z, int a);
  __host__ __device__ double nucdecayenergy(int z, int a);
  __host__ __device__ double meanlife(int z, int a);
  __host__ __device__ double nucmass(int z, int a);
  __host__ __device__ void update_abundances(const int modelgridindex, const int timestep, const double t_current);
  __host__ __device__ double get_simtime_endecay_per_ejectamass(const int mgi, enum decaypathways decaypath);
  __host__ __device__ double get_positroninjection_rate_density(const int modelgridindex, const double t);
  void setup_radioactive_pellet(const double e0, const int mgi, PKT *pkt_ptr);
}

#endif //DECAY_H

#ifndef DECAY_H
#define DECAY_H

// #include <cstdio>

enum radionuclides {
  NUCLIDE_NI57 = 0,
  NUCLIDE_NI56 = 1,
  NUCLIDE_CO56 = 2,
  FAKE_GAM_LINE_ID = 3,
  NUCLIDE_CR48 = 4,
  NUCLIDE_V48 = 5,
  NUCLIDE_CO57 = 6,
  NUCLIDE_FE52 = 7,
  NUCLIDE_MN52 = 8,
  RADIONUCLIDE_COUNT = 9,
};

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
  __host__ __device__ double nucdecayenergygamma(enum radionuclides nuclide_type);
  __host__ __device__ void set_nucdecayenergygamma(enum radionuclides nuclide_type, double value);
  __host__ __device__ double nucdecayenergypositrons(enum radionuclides nuclide_type);
  __host__ __device__ double nucdecayenergy(enum radionuclides nuclide_type);
  __host__ __device__ double meanlife(enum radionuclides nuclide_type);
  __host__ __device__ double nucmass(enum radionuclides nuclide_type);
  __host__ __device__ void update_abundances(const int modelgridindex, const int timestep, const double t_current);
  __host__ __device__ double get_simtime_endecay_per_ejectamass(const int mgi, enum decaypathways decaypath);
  __host__ __device__ double get_positroninjection_rate_density(const int modelgridindex, const double t);
  void setup_radioactive_pellet(const double e0, const int mgi, PKT *pkt_ptr);
}

#endif //DECAY_H

#ifndef DECAY_H
#define DECAY_H

// #include <cstdio>

const int FAKE_GAM_LINE_ID = 3;

#include "cuda.h"

#include <vector>


namespace decay
{
  enum decaytypes {
    DECAYTYPE_ALPHA = 0,
    DECAYTYPE_ELECTRONCAPTURE = 1,
    DECAYTYPE_BETAPLUS = 2,
    DECAYTYPE_BETAMINUS = 3,
    DECAYTYPE_NONE = 4,
    DECAYTYPE_COUNT = 5,
  };

  __host__ __device__ void init_nuclides(std::vector<int> zlist, std::vector<int> alist);
  __host__ __device__ int get_nucstring_z(const char *strnuc);
  __host__ __device__ int get_nucstring_a(const char *strnuc);
  __host__ __device__ int get_num_nuclides(void);
  const char *get_elname(const int z);
  __host__ __device__ int get_nuc_z(int nucindex);
  __host__ __device__ int get_nuc_a(int nucindex);
  __host__ __device__ int get_nuc_index(int z, int a);
  __host__ __device__ bool nuc_exists(int z, int a);
  __host__ __device__ double nucdecayenergygamma(int z, int a);
  __host__ __device__ double nucdecayenergyparticle(const int z_parent, const int a_parent, const int decaytype);
  __host__ __device__ void set_nucdecayenergygamma(int z, int a, double value);
  __host__ __device__ double nucdecayenergy(int z, int a);
  __host__ __device__ double get_meanlife(int z, int a);
  __host__ __device__ double nucmass(int z, int a);
  __host__ __device__ void update_abundances(const int modelgridindex, const int timestep, const double t_current);
  __host__ __device__ double get_endecay_per_ejectamass_t0_to_time_withexpansion(const int modelgridindex, const double tstart);
  __host__ __device__ double get_modelcell_simtime_endecay_per_mass(const int mgi);
  __host__ __device__ void setup_decaypath_energy_per_mass(void);
  __host__ __device__ void free_decaypath_energy_per_mass(void);
  __host__ __device__ double get_qdot_modelcell(const int modelgridindex, const double t, const int decaytype);
  __host__ __device__ double get_particle_injection_rate(const int modelgridindex, const double t, const int decaytype);
  __host__ __device__ double get_global_etot_t0_tinf(void);
  void fprint_nuc_abundances(FILE *estimators_file, const int modelgridindex, const double t_current, const int element);
  __host__ __device__ void setup_radioactive_pellet(const double e0, const int mgi, struct packet *pkt_ptr);
  void cleanup(void);
}

#endif //DECAY_H

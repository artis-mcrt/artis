#ifndef DECAY_H
#define DECAY_H

// #include <cstdio>

const int FAKE_GAM_LINE_ID = 3;

#include <string>
#include <vector>

#include "packet.h"

namespace decay {
enum decaytypes {
  DECAYTYPE_ALPHA = 0,
  DECAYTYPE_ELECTRONCAPTURE = 1,
  DECAYTYPE_BETAPLUS = 2,
  DECAYTYPE_BETAMINUS = 3,
  DECAYTYPE_NONE = 4,
  DECAYTYPE_COUNT = 5,
};

void init_nuclides(std::vector<int> zlist, std::vector<int> alist);
int get_nucstring_z(const std::string &strnuc);
int get_nucstring_a(const std::string &strnuc);
int get_num_nuclides(void);
const char *get_elname(int z);
int get_nuc_z(int nucindex);
int get_nuc_a(int nucindex);
int get_nucindex(int z, int a);
bool nuc_exists(int z, int a);
double nucdecayenergygamma(int nucindex);
double nucdecayenergygamma(int z, int a);
void set_nucdecayenergygamma(int nucindex, double value);
double nucmass(int z, int a);
void update_abundances(int modelgridindex, int timestep, double t_current);
double get_endecay_per_ejectamass_t0_to_time_withexpansion(int modelgridindex, double tstart);
double get_modelcell_simtime_endecay_per_mass(int mgi);
void setup_decaypath_energy_per_mass(void);
void free_decaypath_energy_per_mass(void);
double get_qdot_modelcell(int modelgridindex, double t, int decaytype);
double get_particle_injection_rate(int modelgridindex, double t, int decaytype);
double get_global_etot_t0_tinf(void);
void fprint_nuc_abundances(FILE *estimators_file, int modelgridindex, double t_current, int element);
void setup_radioactive_pellet(double e0, int mgi, struct packet *pkt_ptr);
void cleanup(void);
}  // namespace decay

#endif  // DECAY_H

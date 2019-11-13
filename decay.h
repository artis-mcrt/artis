#ifndef DECAY_H
#define DECAY_H

// #include <stdio.h>

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


double nucdecayenergygamma(enum radionuclides nuclide_type);
void set_nucdecayenergygamma(enum radionuclides nuclide_type, double value);
double nucdecayenergypositrons(enum radionuclides nuclide_type);
double nucdecayenergy(enum radionuclides nuclide_type);
double meanlife(enum radionuclides nuclide_type);
double nucmass(enum radionuclides nuclide_type);
double sample_decaytime(enum decaypathways decaypath, const double tdecaymin, const double tdecaymax);
enum packet_type get_decay_pellet_type(enum decaypathways decaypath, bool *originated_from_positron);
enum radionuclides decayparent(enum decaypathways decaypath);
enum radionuclides decaydaughter(enum decaypathways decaypath);
void update_abundances(const int modelgridindex, const int timestep, const double t_current);
double get_simtime_endecay_per_ejectamass(const int mgi, enum decaypathways decaypath);
double get_decay_power_density(enum decaypathways decaypath, const int modelgridindex, const double time);
double get_positroninjection_rate_density(const int modelgridindex, const double t);
double get_decayedenergy_per_ejectamass(const int n, const double tstart);


#endif //DECAY_H

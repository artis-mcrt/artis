#ifndef ENERGY_INPUT_H
#define ENERGY_INPUT_H

void energy_input_init(void);
double get_etot_fromenergyfile(void);
float get_modelcell_energydensity_init(int modelgridindex);
void set_modelcell_energydensity_init(int modelgridindex, float x);
static void setup_generic_pellet(const double e0, const int mgi, PKT *pkt_ptr);

extern __managed__ double *modelcell_energy; // energy in model grid cell read from energydistribution.txt

#endif //ENERGY_INPUT_H

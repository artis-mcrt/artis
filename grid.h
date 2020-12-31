#ifndef GRIDINIT_H
#define GRIDINIT_H

#include <cassert>
#include "sn3d.h"

int get_elements_uppermost_ion(const int modelgridindex, const int element);
void set_elements_uppermost_ion(const int modelgridindex, const int element, const int newvalue);
double wid_init(int cellindex);
double vol_init_modelcell(int modelgridindex);
double vol_init_gridcell(int cellindex);
double get_cellcoordmin(int cellindex, int axis);
int get_cellcoordpointnum(const int cellindex, const int axis);
int get_coordcellindexincrement(int axis);
int get_ngriddimensions(void);
float get_rhoinit(int modelgridindex);
float get_rho(int modelgridindex);
float get_nne(int modelgridindex);
float get_nnetot(int modelgridindex);
float get_ffegrp(int modelgridindex);
float get_elem_abundance(int modelgridindex, int element);
void set_elem_abundance(int modelgridindex, int element, float newabundance);
float get_kappagrey(int modelgridindex);
float get_Te(int modelgridindex);
float get_TR(int modelgridindex);
float get_TJ(int modelgridindex);
float get_W(int modelgridindex);
void set_rhoinit(int modelgridindex, float x);
void set_rho(int modelgridindex, float x);
void set_nne(int modelgridindex, float x);
void set_nnetot(int modelgridindex, float x);
void set_ffegrp(int modelgridindex, float x);
void set_kappagrey(int modelgridindex, float x);
void set_Te(int modelgridindex, float x);
void set_TR(int modelgridindex, float x);
void set_TJ(int modelgridindex, float x);
void set_W(int modelgridindex, float x);
void grid_init(int my_rank);
double get_cellradialpos(int cellindex);
float get_modelinitradioabund(int modelgridindex, enum radionuclides nuclide_type);
void set_modelinitradioabund(int modelgridindex, enum radionuclides nuclide_type, float abund);
float get_stable_abund(int mgi, int anumber);
int get_numassociatedcells(int modelgridindex);
enum model_types get_model_type(void);
void set_model_type(enum model_types model_type_value);


#endif //GRIDINIT_H

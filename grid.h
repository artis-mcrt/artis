#ifndef GRIDINIT_H
#define GRIDINIT_H

#include <cassert>
#include "sn3d.h"

__host__ __device__ int get_elements_uppermost_ion(const int modelgridindex, const int element);
__host__ __device__ void set_elements_uppermost_ion(const int modelgridindex, const int element, const int newvalue);
__host__ __device__ double wid_init(int cellindex);
__host__ __device__ double vol_init_modelcell(int modelgridindex);
__host__ __device__ double vol_init_gridcell(int cellindex);
__host__ __device__ double get_cellcoordmin(int cellindex, int axis);
__host__ __device__ int get_cellcoordpointnum(const int cellindex, const int axis);
__host__ __device__ int get_coordcellindexincrement(int axis);
__host__ __device__ int get_ngriddimensions(void);
__host__ __device__ float get_rhoinit(int modelgridindex);
__host__ __device__ float get_rho(int modelgridindex);
__host__ __device__ float get_nne(int modelgridindex);
__host__ __device__ float get_nnetot(int modelgridindex);
__host__ __device__ float get_ffegrp(int modelgridindex);
__host__ __device__ float get_elem_abundance(int modelgridindex, int element);
__host__ __device__ void set_elem_abundance(int modelgridindex, int element, float newabundance);
__host__ __device__ float get_kappagrey(int modelgridindex);
__host__ __device__ float get_Te(int modelgridindex);
__host__ __device__ float get_TR(int modelgridindex);
__host__ __device__ float get_TJ(int modelgridindex);
__host__ __device__ float get_W(int modelgridindex);
__host__ __device__ void set_rhoinit(int modelgridindex, float x);
__host__ __device__ void set_rho(int modelgridindex, float x);
__host__ __device__ void set_nne(int modelgridindex, float x);
__host__ __device__ void set_nnetot(int modelgridindex, float x);
__host__ __device__ void set_ffegrp(int modelgridindex, float x);
__host__ __device__ void set_kappagrey(int modelgridindex, float x);
__host__ __device__ void set_Te(int modelgridindex, float x);
__host__ __device__ void set_TR(int modelgridindex, float x);
__host__ __device__ void set_TJ(int modelgridindex, float x);
__host__ __device__ void set_W(int modelgridindex, float x);
void grid_init(int my_rank);
__host__ __device__ double get_cellradialpos(int cellindex);
__host__ __device__ float get_modelinitradioabund(int modelgridindex, enum radionuclides nuclide_type);
__host__ __device__ void set_modelinitradioabund(int modelgridindex, enum radionuclides nuclide_type, float abund);
__host__ __device__ float get_stable_abund(int mgi, int anumber);
__host__ __device__ int get_numassociatedcells(int modelgridindex);
__host__ __device__ enum model_types get_model_type(void);
__host__ __device__ void set_model_type(enum model_types model_type_value);


#endif //GRIDINIT_H

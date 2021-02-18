#ifndef GRIDINIT_H
#define GRIDINIT_H

#include <cassert>
#include "decay.h"
#include "types.h"

enum model_types {
  RHO_UNIFORM = 1,  // Constant density. NOT IN USE
  RHO_1D_READ = 2,  // Read model 1D
  RHO_2D_READ = 4,  // Read model 2D
  RHO_3D_READ = 3,  // Read model 3D
};


typedef struct modelgrid_t
{
  float Te;
  float TR;
  float TJ;
  float W;
  float nne;
  float initial_radial_pos;
  float rhoinit;
  float rho;
  //modelgrid nn_tot
  float nnetot;           // total electron density (free + bound).
  float initradioabund[RADIONUCLIDE_COUNT];
  float ffegrp;
  float fnistable;
  float fcostable;
  float ffestable;
  float fmnstable;
  float fcrstable;
  float fvstable;
  float ftistable;
  float kappagrey;
  float grey_depth;                      /// Grey optical depth to surface of the modelgridcell
                                         /// This is only stored to print it outside the OpenMP loop in update_grid to the estimatorsfile
                                         /// so there is no need to communicate it via MPI so far!
  int *elements_uppermost_ion; /// Highest ionisation stage which has a decent population for a particular element
                                                    /// in a given cell.
  compositionlist_entry *composition;    /// Pointer to an array which contains the time dependent abundances
                                        /// of all included elements and all the groundlevel
                                         /// populations and partition functions for their ions
  double *nlte_pops;                     /// Pointer to an array that contains the nlte-level
                                         /// populations for this cell

  double totalcooling;
  mgicooling_t *cooling;
  short thick;
} modelgrid_t;

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
__host__ __device__ int get_npts_model(void);
__host__ __device__ void set_npts_model(int npts_model);
__host__ __device__ int get_cell_modelgridindex(int cellindex);
void read_ejecta_model(enum model_types model_type);
void write_grid_restart_data(const int timestep);
void show_totmassradionuclides(void);
double get_totmassradionuclide(enum radionuclides nuc);

#endif //GRIDINIT_H

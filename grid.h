#ifndef GRIDINIT_H
#define GRIDINIT_H

#include "artisoptions.h"
#include "constants.h"
#include "cuda.h"

namespace grid {

struct compositionlist_entry {
  float abundance;        /// Abundance of the element (by mass!).
  float *groundlevelpop;  /// Pointer to an array of floats which contains the groundlevel populations
                          /// of all included ionisation stages for the element.
  float *partfunct;       /// Pointer to an array of floats which contains the partition functions
                          /// of all included ionisation stages for the element.
  // float *ltepartfunct;     /// Pointer to an array of floats which contains the LTE partition functions
  //                          /// of all included ionisation stages for the element.
};

struct gridcell {
  double pos_min[3];  // Initial co-ordinates of inner most corner of cell.
  // int xyz[3];       // Integer position of cell in grid.
  int modelgridindex;
};

enum model_types {
  RHO_UNIFORM = 1,  // Constant density. NOT IN USE
  RHO_1D_READ = 2,  // Read model 1D
  RHO_2D_READ = 4,  // Read model 2D
  RHO_3D_READ = 3,  // Read model 3D
};

struct modelgrid_t {
  float Te = 0.;
  float TR = 0.;
  float TJ = 0.;
  float W = 0.;
  float nne = 0.;
  float initial_radial_pos_sum = 0.;
  float rhoinit = 0.;
  float rho = 0.;
  // modelgrid nn_tot
  float nnetot = 0.;  // total electron density (free + bound).
  float *initradioabund = nullptr;
  float *initmassfracstable = nullptr;
  float *elem_meanweight = nullptr;
  float initelectronfrac;  // Ye: electrons (or protons) per nucleon
  float initenergyq;       // q: energy in the model at t_model to use with USE_MODEL_INITIAL_ENERGY [erg/g]
  float ffegrp;
  float kappagrey;
  float grey_depth;  /// Grey optical depth to surface of the modelgridcell
                     /// This is only stored to print it outside the OpenMP loop in update_grid to the estimatorsfile
                     /// so there is no need to communicate it via MPI so far!
  int *elements_uppermost_ion;  /// Highest ionisation stage which has a decent population for a particular element
                                /// in a given cell.
  compositionlist_entry *composition;  /// Pointer to an array which contains the time dependent abundances
                                       /// of all included elements and all the groundlevel
                                       /// populations and partition functions for their ions
  double *nlte_pops = nullptr;         /// Pointer to an array that contains the nlte-level
                                       /// populations for this cell

  double totalcooling;
  double **cooling_contrib_ion;
  short thick;
};

__host__ __device__ constexpr int get_ngriddimensions(void) { return (GRID_TYPE == GRID_SPHERICAL1D) ? 1 : 3; }

extern __managed__ struct modelgrid_t *modelgrid;

extern __managed__ int ncoordgrid[3];
extern __managed__ int ngrid;
extern __managed__ char coordlabel[3];

__host__ __device__ int get_elements_uppermost_ion(int modelgridindex, int element);
__host__ __device__ void set_elements_uppermost_ion(int modelgridindex, int element, int newvalue);
__host__ __device__ double wid_init(int cellindex);
__host__ __device__ double vol_init_modelcell(int modelgridindex);
__host__ __device__ double vol_init_gridcell(int cellindex);
__host__ __device__ double get_cellcoordmax(int cellindex, int axis);
__host__ __device__ double get_cellcoordmin(int cellindex, int axis);
__host__ __device__ int get_cellcoordpointnum(int cellindex, int axis);
__host__ __device__ int get_coordcellindexincrement(int axis);
__host__ __device__ float get_rhoinit(int modelgridindex);
__host__ __device__ float get_rho(int modelgridindex);
__host__ __device__ float get_nne(int modelgridindex);
__host__ __device__ float get_nnetot(int modelgridindex);
__host__ __device__ float get_ffegrp(int modelgridindex);
__host__ __device__ void set_elem_abundance(int modelgridindex, int element, float newabundance);
__host__ __device__ double get_elem_numberdens(int modelgridindex, int element);
__host__ __device__ double get_initelectronfrac(int modelgridindex);
__host__ __device__ double get_initenergyq(int modelgridindex);
__host__ __device__ float get_kappagrey(int modelgridindex);
__host__ __device__ float get_Te(int modelgridindex);
__host__ __device__ float get_TR(int modelgridindex);
__host__ __device__ float get_TJ(int modelgridindex);
__host__ __device__ float get_W(int modelgridindex);
__host__ __device__ void set_nne(int modelgridindex, float x);
__host__ __device__ void set_nnetot(int modelgridindex, float x);
__host__ __device__ void set_kappagrey(int modelgridindex, float x);
__host__ __device__ void set_Te(int modelgridindex, float x);
__host__ __device__ void set_TR(int modelgridindex, float x);
__host__ __device__ void set_TJ(int modelgridindex, float x);
__host__ __device__ void set_W(int modelgridindex, float x);
void grid_init(int my_rank);
__host__ __device__ double get_cellradialpos(int cellindex);
__host__ __device__ float get_modelinitradioabund(int modelgridindex, int z, int a);
__host__ __device__ float get_stable_initabund(int mgi, int anumber);
__host__ __device__ float get_element_meanweight(int mgi, int element);
__host__ __device__ void set_element_meanweight(int mgi, int element, float meanweight);
__host__ __device__ double get_electronfrac(int mgi);
__host__ __device__ int get_numassociatedcells(int modelgridindex);
__host__ __device__ int get_modelcell_nonemptymgi(int mgi);
__host__ __device__ int get_mgi_of_nonemptymgi(int nonemptymgi);
__host__ __device__ enum model_types get_model_type(void);
__host__ __device__ void set_model_type(enum model_types model_type_value);
__host__ __device__ int get_npts_model(void);
__host__ __device__ int get_nonempty_npts_model(void);
__host__ __device__ int get_t_model(void);
__host__ __device__ int get_cell_modelgridindex(int cellindex);
void read_ejecta_model(void);
void write_grid_restart_data(int timestep);
int get_maxndo(void);
int get_nstart(int rank);
int get_ndo(int rank);
int get_ndo_nonempty(int rank);
double get_totmassradionuclide(int z, int a);

__host__ __device__ static inline float get_elem_abundance(int modelgridindex, int element)
// mass fraction of an element (all isotopes combined)
{
  return modelgrid[modelgridindex].composition[element].abundance;
}

}  // namespace grid

#endif  // GRIDINIT_H

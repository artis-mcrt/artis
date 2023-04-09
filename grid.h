#ifndef GRIDINIT_H
#define GRIDINIT_H

#include <cinttypes>

#include "artisoptions.h"
#include "constants.h"

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
  float Te = -1.;
  float TR = -1.;
  float TJ = -1.;
  float W = -1.;
  float nne = -1.;
  float initial_radial_pos_sum = 0.;
  float rhoinit = -1.;
  float rho = -1.;
  // modelgrid nn_tot
  float nnetot = -1.;  // total electron density (free + bound).
  float *initradioabund = nullptr;
  float *initmassfracstable = nullptr;
  float *elem_meanweight = nullptr;
  float initelectronfrac = -1;  // Ye: electrons (or protons) per nucleon
  float initenergyq = 0.;       // q: energy in the model at tmin to use with USE_MODEL_INITIAL_ENERGY [erg/g]
  float ffegrp = 0.;
  float kappagrey = 0.;
  float grey_depth = 0.;  /// Grey optical depth to surface of the modelgridcell
                          /// This is only stored to print it outside the OpenMP loop in update_grid to the
                          /// estimatorsfile so there is no need to communicate it via MPI so far!
  int *elements_uppermost_ion = nullptr;  /// Highest ionisation stage which has a decent population for a particular
                                          /// element in a given cell.
  compositionlist_entry *composition = nullptr;  /// Pointer to an array which contains the time dependent abundances
                                                 /// of all included elements and all the groundlevel
                                                 /// populations and partition functions for their ions
  double *nlte_pops = nullptr;                   /// Pointer to an array that contains the nlte-level
                                                 /// populations for this cell

  double totalcooling = -1;
  double **cooling_contrib_ion = nullptr;
  uint_fast8_t thick = 0;
};

constexpr int get_ngriddimensions(void) { return (GRID_TYPE == GRID_SPHERICAL1D) ? 1 : 3; }

extern struct modelgrid_t *modelgrid;

extern int ncoordgrid[3];
extern int ngrid;
extern char coordlabel[3];

int get_elements_uppermost_ion(int modelgridindex, int element);
void set_elements_uppermost_ion(int modelgridindex, int element, int newvalue);
double wid_init(int cellindex);
double get_modelcell_assocvolume_tmin(int modelgridindex);
double get_gridcell_volume_tmin(int cellindex);
double get_cellcoordmax(int cellindex, int axis);
double get_cellcoordmin(int cellindex, int axis);
int get_cellcoordpointnum(int cellindex, int axis);
int get_coordcellindexincrement(int axis);
float get_rho_tmin(int modelgridindex);
float get_rho(int modelgridindex);
float get_nne(int modelgridindex);
float get_nnetot(int modelgridindex);
float get_ffegrp(int modelgridindex);
void set_elem_abundance(int modelgridindex, int element, float newabundance);
double get_elem_numberdens(int modelgridindex, int element);
double get_initelectronfrac(int modelgridindex);
double get_initenergyq(int modelgridindex);
float get_kappagrey(int modelgridindex);
float get_Te(int modelgridindex);
float get_TR(int modelgridindex);
float get_TJ(int modelgridindex);
float get_W(int modelgridindex);
void set_nne(int modelgridindex, float nne);
void set_nnetot(int modelgridindex, float x);
void set_kappagrey(int modelgridindex, float kappagrey);
void set_Te(int modelgridindex, float Te);
void set_TR(int modelgridindex, float TR);
void set_TJ(int modelgridindex, float TJ);
void set_W(int modelgridindex, float W);
void grid_init(int my_rank);
double get_cellradialpos(int cellindex);
float get_modelinitradioabund(int modelgridindex, int nucindex);
float get_stable_initabund(int mgi, int element);
float get_element_meanweight(int mgi, int element);
void set_element_meanweight(int mgi, int element, float meanweight);
double get_electronfrac(int modelgridindex);
int get_numassociatedcells(int modelgridindex);
int get_modelcell_nonemptymgi(int mgi);
int get_mgi_of_nonemptymgi(int nonemptymgi);
enum model_types get_model_type();
void set_model_type(enum model_types model_type_value);
int get_npts_model();
int get_nonempty_npts_model();
int get_t_model();
int get_cell_modelgridindex(int cellindex);
void read_ejecta_model();
void write_grid_restart_data(int timestep);
int get_maxndo();
int get_nstart(int rank);
int get_ndo(int rank);
int get_ndo_nonempty(int rank);
double get_totmassradionuclide(int z, int a);

static inline float get_elem_abundance(int modelgridindex, int element)
// mass fraction of an element (all isotopes combined)
{
  return modelgrid[modelgridindex].composition[element].abundance;
}

}  // namespace grid

#endif  // GRIDINIT_H

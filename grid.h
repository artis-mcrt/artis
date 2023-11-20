#ifndef GRIDINIT_H
#define GRIDINIT_H

#include <cinttypes>
#include <span>

#include "artisoptions.h"
#include "constants.h"
#include "packet.h"
#include "vectors.h"

namespace grid {

struct compositionlist_entry {
  float abundance;        /// Abundance of the element (by mass!).
  float *groundlevelpop;  /// Pointer to an array of floats which contains the groundlevel populations
                          /// of all included ionisation stages for the element.
  float *partfunct;       /// Pointer to an array of floats which contains the partition functions
                          /// of all included ionisation stages for the element.
};

struct gridcell {
  double pos_min[3];  // Initial co-ordinates of inner most corner of cell.
  int modelgridindex;
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

constexpr auto get_ngriddimensions() -> int {
  switch (GRID_TYPE) {
    case GRID_SPHERICAL1D:
      return 1;
    case GRID_CYLINDRICAL2D:
      return 2;
    case GRID_CARTESIAN3D:
      return 3;
    default:
      assert_always(false);
  }
}

extern struct modelgrid_t *modelgrid;

extern int ncoordgrid[3];
extern int ngrid;
extern char coordlabel[3];

auto get_elements_uppermost_ion(int modelgridindex, int element) -> int;
void set_elements_uppermost_ion(int modelgridindex, int element, int newvalue);
auto wid_init(int cellindex, int axis) -> double;
auto get_modelcell_assocvolume_tmin(int modelgridindex) -> double;
auto get_gridcell_volume_tmin(int cellindex) -> double;
auto get_cellcoordmax(int cellindex, int axis) -> double;
auto get_cellcoordmin(int cellindex, int axis) -> double;
auto get_cellcoordpointnum(int cellindex, int axis) -> int;
auto get_coordcellindexincrement(int axis) -> int;
auto get_rho_tmin(int modelgridindex) -> float;
auto get_rho(int modelgridindex) -> float;
auto get_nne(int modelgridindex) -> float;
auto get_nnetot(int modelgridindex) -> float;
auto get_ffegrp(int modelgridindex) -> float;
void set_elem_abundance(int modelgridindex, int element, float newabundance);
auto get_elem_numberdens(int modelgridindex, int element) -> double;
auto get_initelectronfrac(int modelgridindex) -> double;
auto get_initenergyq(int modelgridindex) -> double;
auto get_kappagrey(int modelgridindex) -> float;
auto get_Te(int modelgridindex) -> float;
auto get_TR(int modelgridindex) -> float;
auto get_TJ(int modelgridindex) -> float;
auto get_W(int modelgridindex) -> float;
void set_nne(int modelgridindex, float nne);
void set_nnetot(int modelgridindex, float x);
void set_kappagrey(int modelgridindex, float kappagrey);
void set_rho(int modelgridindex, float x);
void set_Te(int modelgridindex, float Te);
void set_TR(int modelgridindex, float TR);
void set_TJ(int modelgridindex, float TJ);
void set_W(int modelgridindex, float W);
void grid_init(int my_rank);
auto get_modelinitradioabund(int modelgridindex, int nucindex) -> float;
auto get_stable_initabund(int mgi, int element) -> float;
auto get_element_meanweight(int mgi, int element) -> float;
void set_element_meanweight(int mgi, int element, float meanweight);
auto get_electronfrac(int modelgridindex) -> double;
auto get_numassociatedcells(int modelgridindex) -> int;
auto get_modelcell_nonemptymgi(int mgi) -> int;
auto get_mgi_of_nonemptymgi(int nonemptymgi) -> int;
auto get_model_type() -> enum gridtypes;
void set_model_type(enum gridtypes model_type_value);
auto get_npts_model() -> int;
auto get_nonempty_npts_model() -> int;
auto get_t_model() -> double;
auto get_cell_modelgridindex(int cellindex) -> int;
auto get_cellindex_from_pos(std::span<const double, 3> pos, double time) -> int;
void read_ejecta_model();
void write_grid_restart_data(int timestep);
auto get_maxndo() -> int;
auto get_nstart(int rank) -> int;
auto get_ndo(int rank) -> int;
auto get_ndo_nonempty(int rank) -> int;
auto get_totmassradionuclide(int z, int a) -> double;
auto boundary_distance(std::span<const double, 3> dir, std::span<const double, 3> pos, double tstart, int cellindex,
                       int *snext, enum cell_boundary *pkt_last_cross) -> double;
void change_cell(struct packet *pkt_ptr, int snext);

static inline auto get_elem_abundance(int modelgridindex, int element) -> float
// mass fraction of an element (all isotopes combined)
{
  return modelgrid[modelgridindex].composition[element].abundance;
}

}  // namespace grid

#endif  // GRIDINIT_H

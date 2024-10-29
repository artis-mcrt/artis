#include <cstddef>
#include <span>

#ifndef GRIDINIT_H
#define GRIDINIT_H

#include <tuple>

#include "artisoptions.h"
#include "constants.h"
#include "globals.h"
#include "packet.h"
#include "sn3d.h"
#include "stats.h"

namespace grid {

struct ModelGridCell {
  float Te = -1.;
  float TR = -1.;
  float TJ = -1.;
  float W = -1.;
  float nne = -1.;
  float initial_radial_pos_sum = 0.;
  float rhoinit = -1.;
  float rho = -1.;
  // modelgrid nn_tot
  float nnetot = -1.;           // total electron density (free + bound).
  float initelectronfrac = -1;  // Ye: electrons (or protons) per nucleon
  float initenergyq = 0.;       // q: energy in the model at tmin to use with USE_MODEL_INITIAL_ENERGY [erg/g]
  float ffegrp = 0.;
  float kappagrey = 0.;
  float grey_depth = 0.;         // Grey optical depth to surface of the modelgridcell
                                 // This is only stored to print it outside the OpenMP loop in update_grid to the
                                 // estimatorsfile so there is no need to communicate it via MPI so far!
  float *ion_groundlevelpops{};  // groundlevel populations of all included ions
  float *ion_partfuncts{};       // partition functions for all included ions
  std::span<double> nlte_pops;   // Pointer to an array that contains the nlte-level populations for this cell
  double totalcooling = -1;
  double *ion_cooling_contribs{};
  int thick = 0;
};

consteval auto get_ngriddimensions() -> int {
  switch (GRID_TYPE) {
    case GridType::SPHERICAL1D:
      return 1;
    case GridType::CYLINDRICAL2D:
      return 2;
    case GridType::CARTESIAN3D:
      return 3;
    default:
      assert_always(false);
      return -1;
  }
}

inline std::span<ModelGridCell> modelgrid{};

inline int ngrid{0};

inline double mtot_input{0.};

inline float *elem_meanweight_allcells{};
inline float *elem_massfracs_allcells;  // mass fractions of elements in each cell for the current timestep

[[nodiscard]] auto get_elements_uppermost_ion(int modelgridindex, int element) -> int;
void set_elements_uppermost_ion(int modelgridindex, int element, int newvalue);
[[nodiscard]] auto wid_init(int cellindex, int axis) -> double;
[[nodiscard]] auto get_modelcell_assocvolume_tmin(int modelgridindex) -> double;
[[nodiscard]] auto get_gridcell_volume_tmin(int cellindex) -> double;
[[nodiscard]] auto get_cellcoordmax(int cellindex, int axis) -> double;
[[nodiscard]] auto get_cellcoordmin(int cellindex, int axis) -> double;
[[nodiscard]] auto get_cellcoordpointnum(int cellindex, int axis) -> int;
[[nodiscard]] auto get_cellradialposmid(int cellindex) -> double;
[[nodiscard]] auto get_coordcellindexincrement(int axis) -> int;
[[nodiscard]] auto get_rho_tmin(int modelgridindex) -> float;
[[nodiscard]] auto get_rho(int modelgridindex) -> float;
[[nodiscard]] auto get_nne(int modelgridindex) -> float;
[[nodiscard]] auto get_nnetot(int modelgridindex) -> float;
[[nodiscard]] auto get_ffegrp(int modelgridindex) -> float;
void set_elem_abundance(int nonemptymgi, int element, float newabundance);
[[nodiscard]] auto get_elem_numberdens(int modelgridindex, int element) -> double;
[[nodiscard]] auto get_initelectronfrac(int modelgridindex) -> double;
[[nodiscard]] auto get_initenergyq(int modelgridindex) -> double;
[[nodiscard]] auto get_kappagrey(int modelgridindex) -> float;
[[nodiscard]] auto get_Te(int modelgridindex) -> float;
[[nodiscard]] auto get_TR(int modelgridindex) -> float;
[[nodiscard]] auto get_TJ(int modelgridindex) -> float;
[[nodiscard]] auto get_W(int modelgridindex) -> float;
void set_nne(int modelgridindex, float nne);
void set_nnetot(int modelgridindex, float nnetot);
void set_kappagrey(int modelgridindex, float kappagrey);
void set_rho(int modelgridindex, float rho);
void set_Te(int modelgridindex, float Te);
void set_TR(int modelgridindex, float TR);
void set_TJ(int modelgridindex, float TJ);
void set_W(int modelgridindex, float W);
void grid_init(int my_rank);
[[nodiscard]] auto get_modelinitnucmassfrac(int modelgridindex, int nucindex) -> float;
[[nodiscard]] auto get_stable_initabund(int nonemptymgi, int element) -> float;
[[nodiscard]] auto get_element_meanweight(int mgi, int element) -> float;
[[nodiscard]] auto get_elem_abundance(int modelgridindex, int element) -> float;
void set_element_meanweight(int nonemptymgi, int element, float meanweight);
[[nodiscard]] auto get_electronfrac(int modelgridindex) -> double;
[[nodiscard]] auto get_numassociatedcells(int modelgridindex) -> int;
[[nodiscard]] auto get_modelcell_nonemptymgi(int mgi) -> int;
[[nodiscard]] auto get_mgi_of_nonemptymgi(int nonemptymgi) -> int;
[[nodiscard]] auto get_model_type() -> GridType;
void set_model_type(GridType model_type_value);
[[nodiscard]] auto get_npts_model() -> int;
[[nodiscard]] auto get_nonempty_npts_model() -> int;
[[nodiscard]] auto get_t_model() -> double;
[[nodiscard]] auto get_cell_modelgridindex(int cellindex) -> int;
[[nodiscard]] auto get_cellindex_from_pos(const std::array<double, 3> &pos, double time) -> int;
void read_ejecta_model();
void write_grid_restart_data(int timestep);
[[nodiscard]] auto get_maxndo() -> int;
[[nodiscard]] auto get_nstart(int rank) -> int;
[[nodiscard]] auto get_nstart_nonempty(int rank) -> int;
[[nodiscard]] auto get_ndo(int rank) -> int;
[[nodiscard]] auto get_ndo_nonempty(int rank) -> int;
[[nodiscard]] auto get_totmassradionuclide(int z, int a) -> double;
[[nodiscard]] auto boundary_distance(const std::array<double, 3> &dir, const std::array<double, 3> &pos, double tstart,
                                     int cellindex, enum cell_boundary *pkt_last_cross) -> std::tuple<double, int>;

void calculate_kappagrey();

inline void change_cell(Packet &pkt, const int snext)
// Routine to take a packet across a boundary.
{
  if (snext >= 0) {
    // Just need to update "where".
    pkt.where = snext;
  } else {
    // Then the packet is exiting the grid. We need to record
    // where and at what time it leaves the grid.
    pkt.escape_type = pkt.type;
    pkt.escape_time = pkt.prop_time;
    pkt.type = TYPE_ESCAPE;
    atomicadd(globals::nesc, 1);

    stats::increment(stats::COUNTER_CELLCROSSINGS);
  }
}

inline auto get_ejecta_kinetic_energy() {
  double E_kin = 0.;
  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    const int assoc_cells = grid::get_numassociatedcells(mgi);
    double M_cell = grid::get_rho_tmin(mgi) * grid::get_modelcell_assocvolume_tmin(mgi);
    const double radial_pos = grid::modelgrid[mgi].initial_radial_pos_sum / assoc_cells;
    E_kin += 0.5 * M_cell * std::pow(radial_pos / globals::tmin, 2);
  }

  return E_kin;
}

}  // namespace grid

#endif  // GRIDINIT_H

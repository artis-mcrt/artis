#include <cmath>
#include <cstddef>
#include <span>

#ifndef GRIDINIT_H
#define GRIDINIT_H

#include <tuple>

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
  float rho = -1.;
  // modelgrid nn_tot
  float nnetot = -1.;  // total electron density (free + bound).
  float kappagrey = 0.;
  float grey_depth = 0.;  // Grey optical depth to surface of the modelgridcell
                          // This is only stored to print it outside the OpenMP loop in update_grid to the
                          // estimatorsfile so there is no need to communicate it via MPI so far!
  double totalcooling = -1;
  int thick = 0;
};

inline std::span<ModelGridCell> modelgrid{};

inline ptrdiff_t ngrid{0};

inline double mtot_input{0.};

inline std::span<float> elem_meanweight_allcells{};
inline std::span<float> elem_massfracs_allcells;  // mass fractions of elements in each cell for the current timestep

inline std::span<double> nltepops_allcells{};
inline std::span<float> ion_groundlevelpops_allcells{};
inline std::span<float> ion_partfuncts_allcells{};
inline std::span<double> ion_cooling_contribs_allcells{};

[[nodiscard]] auto get_elements_uppermost_ion(int nonemptymgi, int element) -> int;
void set_elements_uppermost_ion(int nonemptymgi, int element, int uppermost_ion);
[[gnu::pure]] [[nodiscard]] auto wid_init(int cellindex, int axis) -> double;
[[gnu::pure]] [[nodiscard]] auto get_modelcell_assocvolume_tmin(int modelgridindex) -> double;
[[gnu::pure]] [[nodiscard]] auto get_propcell_volume_tmin(int cellindex) -> double;
[[gnu::pure]] [[nodiscard]] auto get_cellcoordmax(int cellindex, int axis) -> double;
[[gnu::pure]] [[nodiscard]] auto get_cellcoordmin(int cellindex, int axis) -> double;
[[gnu::pure]] [[nodiscard]] auto get_cellcoordpointnum(int cellindex, int axis) -> int;
[[gnu::pure]] [[nodiscard]] auto get_cellradialposmid(int cellindex) -> double;
[[gnu::pure]] [[nodiscard]] auto get_coordcellindexincrement(int axis) -> int;
[[gnu::pure]] [[nodiscard]] auto get_rho_tmin(int modelgridindex) -> float;
[[gnu::pure]] [[nodiscard]] auto get_rho(std::ptrdiff_t nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] auto get_nne(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] auto get_nnetot(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] auto get_ffegrp(int modelgridindex) -> float;
[[gnu::pure]] [[nodiscard]] auto get_modelcell_mean_radial_pos(int modelgridindex, double tratmid) -> double;
void set_elem_abundance(std::ptrdiff_t nonemptymgi, int element, float newabundance);
[[gnu::pure]] [[nodiscard]] auto get_elem_numberdens(std::ptrdiff_t nonemptymgi, int element) -> double;
[[gnu::pure]] [[nodiscard]] auto get_initenergyq(int modelgridindex) -> double;
[[gnu::pure]] [[nodiscard]] auto get_kappagrey(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] auto get_Te(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] auto get_TR(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] auto get_TJ(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] auto get_W(int nonemptymgi) -> float;
void set_nne(int nonemptymgi, float nne);
void set_nnetot(int nonemptymgi, float nnetot);
void set_kappagrey(int nonemptymgi, float kappagrey);
void set_rho(int nonemptymgi, float rho);
void set_Te(int nonemptymgi, float Te);
void set_TR(int nonemptymgi, float TR);
void set_TJ(int nonemptymgi, float TJ);
void set_W(int nonemptymgi, float W);
void init_grid(int my_rank);
[[gnu::pure]] [[nodiscard]] auto get_modelinitnucmassfrac(int modelgridindex, int nucindex) -> float;
[[gnu::pure]] [[nodiscard]] auto get_stable_initabund(std::ptrdiff_t nonemptymgi, int element) -> float;
[[gnu::pure]] [[nodiscard]] auto get_element_meanweight(std::ptrdiff_t nonemptymgi, int element) -> float;
[[gnu::pure]] [[nodiscard]] auto get_elem_abundance(std::ptrdiff_t nonemptymgi, int element) -> float;
void set_element_meanweight(std::ptrdiff_t nonemptymgi, int element, float meanweight);
[[gnu::pure]] [[nodiscard]] auto get_electronfrac(int nonemptymgi) -> double;
[[gnu::pure]] [[nodiscard]] auto get_numpropcells(int modelgridindex) -> int;
[[gnu::pure]] [[nodiscard]] auto get_nonemptymgi_of_mgi(int mgi) -> int;
[[gnu::pure]] [[nodiscard]] auto get_mgi_of_nonemptymgi(std::ptrdiff_t nonemptymgi) -> int;
[[gnu::pure]] [[nodiscard]] auto get_model_type() -> GridType;
void set_model_type(GridType model_type_value);
[[gnu::pure]] [[nodiscard]] auto get_npts_model() -> int;
[[gnu::pure]] [[nodiscard]] auto get_nonempty_npts_model() -> int;
[[gnu::pure]] [[nodiscard]] auto get_t_model() -> double;
[[gnu::pure]] [[nodiscard]] auto get_propcell_modelgridindex(int cellindex) -> int;
[[gnu::pure]] [[nodiscard]] auto get_propcell_nonemptymgi(int cellindex) -> int;
[[gnu::pure]] [[nodiscard]] auto get_cellindex_from_pos(const Vec3d& pos, double time) -> int;
void read_ejecta_model();
void write_grid_restart_data(int timestep);
[[gnu::pure]] [[nodiscard]] auto get_nstart(int rank) -> int;
[[gnu::pure]] [[nodiscard]] auto get_nstart_nonempty(int rank) -> int;
[[gnu::pure]] [[nodiscard]] auto get_ndo(int rank) -> int;
[[gnu::pure]] [[nodiscard]] auto get_ndo_nonempty(int rank) -> int;
[[gnu::pure]] [[nodiscard]] auto get_totmassradionuclide_tmodel(int z, int a) -> double;
[[nodiscard]] auto boundary_distance(const Vec3d& dir, const Vec3d& pos, double tstart, int cellindex)
    -> std::tuple<double, int>;

void calculate_kappagrey();

inline void change_cell(Packet& pkt, const int snext)
// Routine to take a packet across a boundary.
{
  if (snext >= 0) {
    // Just need to update "where".
    pkt.where = snext;
  } else {
    // Then the packet is exiting the grid. We need to record
    // where and at what time it leaves the grid.
    pkt.escape_type = pkt.type;
    pkt.escape_time = static_cast<float>(pkt.prop_time);
    pkt.type = TYPE_ESCAPE;
    atomicadd(globals::nesc, 1);

    stats::increment(stats::COUNTER_CELLCROSSINGS);
  }
}

inline auto get_ejecta_kinetic_energy() {
  double E_kin = 0.;
  for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);
    double const M_cell = get_rho_tmin(mgi) * grid::get_modelcell_assocvolume_tmin(mgi);
    const double radial_pos = get_modelcell_mean_radial_pos(mgi, 1.);
    E_kin += 0.5 * M_cell * std::pow(radial_pos / globals::tmin, 2);
  }

  return E_kin;
}

}  // namespace grid

#endif  // GRIDINIT_H

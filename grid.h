#include <cmath>
#include <cstddef>
#include <span>

#ifndef GRIDINIT_H
#define GRIDINIT_H

#include <tuple>

#include "constants.h"
#include "globals.h"
#include "packet.h"
#include "stats.h"

namespace grid {

// these arrays are indexed by nonemptymgi
inline std::span<float> rho_allcells;
inline std::span<float> Te_allcells;
inline std::span<float> TJ_allcells;
inline std::span<float> TR_allcells;
inline std::span<float> W_allcells;
inline std::span<float> nne_allcells;
inline std::span<float> nnetot_allcells;  // total electron density (free + bound).
inline std::span<float> kappagrey_allcells;
inline std::span<float> grey_depth_allcells;  // Grey optical depth to surface of the modelgridcell
inline std::span<int>
    thick_allcells;  // whether the cell is optically thick (1) or not (0), or (2) thick for vpkts only

inline ptrdiff_t ngrid{0};

inline double mtot_input{0.};

inline std::span<float> elem_meanweight_allcells{};
inline std::span<float> elem_massfracs_allcells{};  // mass fractions of elements in each cell for the current timestep

inline std::span<double> nltepops_allcells{};
inline std::span<float> ion_groundlevelpops_allcells{};
inline std::span<float> ion_partfuncts_allcells{};

[[nodiscard]] auto get_elements_uppermost_ion(int nonemptymgi, int element) -> int;
void set_elements_uppermost_ion(int nonemptymgi, int element, int uppermost_ion);
[[gnu::pure]] [[nodiscard]] auto wid_init(int cellindex, int axis) -> double;
[[gnu::pure]] [[nodiscard]] auto get_modelcell_assocvolume_tmin(int modelgridindex) -> double;
[[gnu::pure]] [[nodiscard]] auto get_propcell_volume_tmin(int cellindex) -> double;
[[gnu::pure]] [[nodiscard]] auto get_cellcoordmax(int cellindex, int axis) -> double;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellcoordmin(int cellindex, int axis) -> double;
[[gnu::pure]] [[nodiscard]] auto get_cellcoordpointnum(int cellindex, int axis) -> int;
[[gnu::pure]] [[nodiscard]] auto get_cellradialposmid(int cellindex) -> double;
[[gnu::pure]] [[nodiscard]] auto get_coordcellindexincrement(int axis) -> int;
[[gnu::pure]] [[nodiscard]] auto get_rho_tmin(int modelgridindex) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_rho(std::ptrdiff_t nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nne(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nnetot(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_ffegrp(int modelgridindex) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_modelcell_mean_radial_pos(int modelgridindex, double tratmid)
    -> double;
[[gnu::pure]] [[nodiscard]] auto get_modelcell_mean_radial_vel(int modegridindex, double tratmid) -> double;
void set_elem_abundance(std::ptrdiff_t nonemptymgi, int element, float newabundance);
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_elem_numberdens(std::ptrdiff_t nonemptymgi, int element) -> double;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_initenergyq(int modelgridindex) -> double;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_kappagrey(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_Te(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_TR(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_TJ(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_W(int nonemptymgi) -> float;
void set_nne(int nonemptymgi, float nne);
void set_nnetot(int nonemptymgi);
void set_kappagrey(int nonemptymgi, float kappagrey);
void set_rho(int nonemptymgi, float rho);
void set_Te(int nonemptymgi, float Te);
void set_TR(int nonemptymgi, float TR);
void set_TJ(int nonemptymgi, float TJ);
void set_W(int nonemptymgi, float W);
void init_grid(int my_rank);
[[gnu::pure]] [[nodiscard]] auto get_modelinitnucmassfrac(int modelgridindex, int nucindex) -> float;
[[gnu::pure]] [[nodiscard]] auto get_otherstable_initabund(std::ptrdiff_t nonemptymgi, int element) -> float;
[[gnu::pure]] [[nodiscard]] auto get_element_meanweight(std::ptrdiff_t nonemptymgi, int element) -> float;
[[gnu::pure]] [[nodiscard]] auto get_elem_abundance(std::ptrdiff_t nonemptymgi, int element) -> float;
void set_element_meanweight(std::ptrdiff_t nonemptymgi, int element, float meanweight);
[[gnu::pure]] [[nodiscard]] auto get_electronfrac(int nonemptymgi) -> double;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_numpropcells(int modelgridindex) -> int;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nonemptymgi_of_mgi(int mgi) -> int;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_mgi_of_nonemptymgi(std::ptrdiff_t nonemptymgi) -> int;
[[gnu::pure]] [[nodiscard]] auto get_model_type() -> GridType;
void set_model_type(GridType model_type_value);
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_npts_model() -> int;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nonempty_npts_model() -> int;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_t_model() -> double;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_propcell_modelgridindex(int cellindex) -> int;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_propcell_nonemptymgi(int cellindex) -> int;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellindex_from_pos(const Vec3d& pos, double time) -> int;
void read_ejecta_model();
void write_grid_restart_data(int timestep);
[[gnu::pure]] [[nodiscard]] auto get_nstart(int rank) -> int;
[[gnu::pure]] [[nodiscard]] auto get_nstart_nonempty(int rank) -> int;
[[gnu::pure]] [[nodiscard]] auto get_ndo(int rank) -> int;
[[gnu::pure]] [[nodiscard]] auto get_ndo_nonempty(int rank) -> int;
[[gnu::pure]] [[nodiscard]] auto get_totmassnuclide_tmodel(int z, int a) -> double;
[[nodiscard]] DEVICE_FUNC auto boundary_distance(const Vec3d& dir, const Vec3d& pos, double tstart, int cellindex)
    -> std::tuple<double, int>;

void calculate_kappagrey();

// Routine to take a packet across a boundary.
inline void change_cell(Packet& pkt, const int snext) {
  if (snext >= 0) {
    // Just need to update "where".
    pkt.where = snext;
    stats::increment(stats::Counter::CELLCROSSINGS);
  } else {
    // Then the packet is exiting the grid. We need to record
    // where and at what time it leaves the grid.
    pkt.escape_type = pkt.type;
    pkt.escape_time = static_cast<float>(pkt.prop_time);
    pkt.type = TYPE_ESCAPE;
    stats::increment(stats::Counter::PKTESCAPES);
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

// TODO: remove this function. i think it's not needed any more
// Applies the clumping factor to `val` `pow` times if microclumping is being used
inline auto apply_clumping(const float val, const float oneoverfv, const int pow) -> float {
  if constexpr (USE_MICROCLUMPING) {
    return val * std::pow(oneoverfv, pow);
  }
  return val;
}

// Applies the clumping factor to `val` if microclumping is being used
inline auto apply_clumping(const float val, const float oneoverfv) -> float {
  if constexpr (USE_MICROCLUMPING) {
    return val * oneoverfv;
  }
  return val;
}

}  // namespace grid

#endif  // GRIDINIT_H

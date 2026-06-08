#ifndef GRIDINIT_H
#define GRIDINIT_H

#include <cstddef>
#include <tuple>

#include "constants.h"
#include "globals.h"
#include "mpi_logging.h"
#include "packet.h"

namespace grid {

// these arrays are indexed by nonemptymgi
inline MPI_shared_array<float> rho_allcells;
inline MPI_shared_array<float> Te_allcells;
inline MPI_shared_array<float> TJ_allcells;
inline MPI_shared_array<float> TR_allcells;
inline MPI_shared_array<float> W_allcells;
inline MPI_shared_array<float> nne_allcells;
inline MPI_shared_array<float> nnetot_allcells;  // total electron density (free + bound).
inline MPI_shared_array<float> kappagrey_allcells;
inline MPI_shared_array<float> grey_depth_allcells;  // Grey optical depth to surface of the modelgridcell
inline MPI_shared_array<int>
    thick_allcells;  // whether the cell is optically thick (1) or not (0), or (2) thick for vpkts only

inline ptrdiff_t ngrid{0};

inline double mtot_input{0.};

inline MPI_shared_array<float> elem_meanweight_allcells;

// mass fractions of elements in each cell for the current timestep
inline MPI_shared_array<float> elem_massfracs_allcells;

inline MPI_shared_array<float> ion_groundlevelpops_allcells;
inline MPI_shared_array<float> ion_partfuncts_allcells;

[[nodiscard]] auto get_elements_uppermost_ion(int nonemptymgi, int element) -> int;
void set_elements_uppermost_ion(int nonemptymgi, int element, int uppermost_ion);
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto propcell_width_tmin(int cellindex, int axis) -> double;
[[gnu::pure]] [[nodiscard]] auto get_modelcell_assocvolume_tmin(int modelgridindex) -> double;
[[gnu::pure]] [[nodiscard]] auto get_propcell_volume_tmin(int cellindex) -> double;
[[gnu::pure]] [[nodiscard]] auto get_rho_tmin(int modelgridindex) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_rho(std::ptrdiff_t nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nne(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nnetot(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_ffegrp(int modelgridindex) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_modelcell_mean_radial_pos(int modelgridindex, double tratmid)
    -> double;
void set_elem_massfrac(std::ptrdiff_t nonemptymgi, int element, float newmassfrac);
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
void init_grid();
[[gnu::pure]] [[nodiscard]] auto get_modelinitnucmassfrac(int modelgridindex, int nucindex) -> float;
[[gnu::pure]] [[nodiscard]] auto get_otherstable_initabund(std::ptrdiff_t nonemptymgi, int element) -> float;
[[gnu::pure]] [[nodiscard]] auto get_elem_massfrac(std::ptrdiff_t nonemptymgi, int element) -> float;
void set_element_meanweight(std::ptrdiff_t nonemptymgi, int element, float meanweight);
[[gnu::pure]] [[nodiscard]] auto get_electronfrac(int nonemptymgi) -> double;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_numpropcells(int modelgridindex) -> int;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nonemptymgi_of_mgi(int mgi) -> int;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_mgi_of_nonemptymgi(std::ptrdiff_t nonemptymgi) -> int;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_modelgridtype() -> GridType;
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
[[nodiscard]] auto get_propcell_random_position_tmin(int cellindex) -> Vec3d;
[[nodiscard]] DEVICE_FUNC auto boundary_distance(const Vec3d& dir, const Vec3d& pos, double tstart, int cellindex)
    -> std::tuple<double, int>;
DEVICE_FUNC void change_cell(Packet& pkt, int next_cellindex);

[[nodiscard]] auto calculate_cell_kappagrey(int nonemptymgi) -> float;

inline auto get_ejecta_kinetic_energy() {
  double E_kin = 0.;
  for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);
    double const M_cell = get_rho_tmin(mgi) * grid::get_modelcell_assocvolume_tmin(mgi);
    const double radial_pos = get_modelcell_mean_radial_pos(mgi, 1.);
    E_kin += 0.5 * M_cell * pow2(radial_pos / globals::tmin);
  }

  return E_kin;
}

}  // namespace grid

#endif  // GRIDINIT_H

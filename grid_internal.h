// The state and the helper functions that the grid source files share: grid.cc, grid_geometry.cc, and
// grid_modelinput.cc. The other modules include grid.h.

#ifndef GRID_INTERNAL_H
#define GRID_INTERNAL_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <tuple>
#include <vector>

#include "constants.h"
#include "grid.h"
#include "mpi_logging.h"
#include "vectors.h"

namespace grid {

inline std::array<int, 3> ncoordgrid{0, 0, 0};  // propagation grid dimensions

inline std::optional<GridType> model_type{};

// The optional key mirror_plane z=0 of a model.txt of format version 2 gives mirror symmetry in the plane z = 0. The
// file then holds the cells with z >= 0, and each model cell also fills its mirror cell. The propagation grid has the
// cells of both halves, and the two cells of a mirror pair have the same model cell.
inline bool model_mirror_z{false};

inline ptrdiff_t npts_model = 0;  // number of model grid cells

inline double t_model = -1.;  // time at which densities in input model are correct.
inline std::vector<double> vout_model{};
inline std::array<int, 3> ncoord_model{};  // the model.txt input grid dimensions

// The cell edges of a 2D spherical or a 3D spherical model, which model.txt gives in its header. Each list has one
// value more than the count of the cells on that axis. A 2D spherical model has one phi cell from 0 to 2 pi.
inline std::vector<double> model_vel_r_edges{};  // radial velocity [cm/s], from 0 to vmax
// cos(theta) of the polar angle theta from the +z axis, from -1 to 1 (from 0 to 1 with a mirror)
inline std::vector<double> model_costheta_edges{};
inline std::vector<double> model_phi_edges{};  // azimuthal angle phi = atan2(y, x) [rad], from 0 to 2 pi

inline double mtot_input{0.};  // Total mass of the input model
inline double mfegroup{0.};  // Total mass of Fe group elements in ejecta

inline constexpr int max_atomic_number = 150;  // abundances.txt can give elements up to this atomic number

// with FORCE_SPHERICAL_ESCAPE_SURFACE: 1 for a propagation cell whose inner corner lies outside the
// spherical escape surface. boundary_distance() runs for every packet step, so
// allocate_nonemptymodelcells() precomputes the per-cell radius comparison.
inline std::vector<uint8_t> propcell_outside_escape_surface;

inline std::vector<int> modelgrid_numpropcells;
inline std::vector<int> nonemptymgi_of_mgi;

inline MPI_shared_array<float> initnucmassfrac_allcells{};

struct ModelGridCellInput {
  float rhoinit = -1.;
  float ffegrp = 0.;
  // sum of the midpoint radius per propagation cell. A double, because the sum runs over up to millions of
  // propagation cells.
  double initial_radial_pos_sum = 0.;
  // sum of the volume averaged r^2 per propagation cell, for the kinetic energy of the BARNES
  // thermalisation scheme. A double, because the sum runs over up to millions of propagation cells
  // with radii of ~1e15 cm.
  double initial_radial_pos_squared_sum = 0.;
  float initelectronfrac = -1.;  // Ye: electrons (or protons) per nucleon. Negative until the input sets it
  // q: the trapped radiation energy per mass [erg/g], scaled from t_model to tmin (INITIAL_PACKETS_ON)
  float initenergyq = 0.;
};
inline MPI_shared_array<ModelGridCellInput> modelgrid_input{};

// Get number of dimensions
constexpr auto get_ndim(const GridType gridtype) -> int {
  switch (gridtype) {
    case GridType::SPHERICAL1D:
      return 1;
    case GridType::CYLINDRICAL2D:
    case GridType::SPHERICAL2D:
      return 2;
    case GridType::CARTESIAN3D:
    case GridType::SPHERICAL3D:
      return 3;
  }
  fatal_crash("Unknown grid type {}", static_cast<int>(gridtype));
}

constexpr auto has_spherical_angle_axes(const GridType gridtype) -> bool {
  return gridtype == GridType::SPHERICAL2D || gridtype == GridType::SPHERICAL3D;
}

constexpr auto get_grid_type_name(const GridType gridtype) -> std::string_view {
  switch (gridtype) {
    case GridType::SPHERICAL1D:
    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D:
      return "spherical";
    case GridType::CYLINDRICAL2D:
      return "cylindrical";
    case GridType::CARTESIAN3D:
      return "uniform cuboidal";
  }
  return "unknown";
}

[[nodiscard]] DEVICE_FUNC constexpr auto get_coordlabel(const GridType gridtype, const int axis) -> std::string_view {
  assert_always(axis >= 0 && axis < get_ndim(gridtype));
  switch (gridtype) {
    case GridType::CARTESIAN3D:
      return std::array<std::string_view, 3>{"x", "y", "z"}[axis];
    case GridType::CYLINDRICAL2D:
      return std::array<std::string_view, 2>{"r", "z"}[axis];
    case GridType::SPHERICAL1D:
      return "r";
    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D:
      return std::array<std::string_view, 3>{"r", "costheta", "phi"}[axis];
  }
  fatal_crash("Unknown grid type {}", static_cast<int>(gridtype));
}

// grid.cc
void set_rho_tmin(int modelgridindex, float x);
[[nodiscard]] auto get_tracked_isotope_massfrac_sum_by_atomic_number(int mgi)
    -> std::array<double, max_atomic_number + 1>;
void calc_modelinit_totmassnuclides();

// grid_geometry.cc
[[gnu::pure]] DEVICE_FUNC auto get_propgridtype() -> GridType;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellcoordindex(int cellindex, int axis) -> int;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellcoordmin(int cellindex, int axis) -> double;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellcoordmax(int cellindex, int axis) -> double;
auto get_cell_r_inner(int cellindex, GridType prop_gridtype) -> double;
auto get_cellradialposmid(int cellindex) -> double;
auto get_cellradialposmeansquared(int cellindex) -> double;
[[nodiscard]] auto get_spherical_model_cellindex_tmin(const Vec3d& pos_tmin) -> int;
void require_increasing_shell_radii_at_tmin(std::span<const double> vel_r_edges);
void setup_grid_cartesian_3d();
void setup_grid_spherical_1d();
void setup_grid_cylindrical_2d();
void setup_grid_spherical();

// grid_modelinput.cc
void set_modelinitnucmassfrac(int modelgridindex, int nucindex, float abund);
auto read_elem_abundances() -> std::tuple<std::vector<float>, int>;
void set_spherical_model_cell_counts(GridType spherical_model_type);

}  // namespace grid

#endif  // GRID_INTERNAL_H

// Declarations for the model/propagation grid (grid.cc) and inline cell-property accessors.

#ifndef GRID_H
#define GRID_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <tuple>

#include "artisoptions.h"
#include "constants.h"
#include "globals.h"
#include "mpi_logging.h"
#include "packet.h"
#include "random.h"
#include "stats.h"
#include "vectors.h"

namespace grid {

enum class BoundaryType : std::uint8_t { LOWER, UPPER };

// Find the first forward crossing of an expanding circle or sphere.
// Return -1 when the ray has no crossing in the requested radial direction.
template <BoundaryType boundarytype, size_t S1>
[[gnu::pure]] [[nodiscard]] constexpr auto expanding_shell_intersection(const std::array<double, S1>& pos,
                                                                        const std::array<double, S1>& dir,
                                                                        const double speed,
                                                                        const double shellradiuststart,
                                                                        const double tstart) -> double {
  static_assert(S1 == 2 || S1 == 3);
  assert_testmodeonly(shellradiuststart > 0);

  // quadratic equation for intersection of ray with sphere
  // a*d^2 + b*d + c = 0
  const double a = dot(dir, dir) - pow2(shellradiuststart / tstart / speed);
  const double b = 2 * (dot(dir, pos) - (pow2(shellradiuststart) / tstart / speed));
  const double c = dot(pos, pos) - pow2(shellradiuststart);

  const double discriminant = pow2(b) - (4 * a * c);

  if (discriminant < 0) {
    // no intersection
    assert_testmodeonly(boundarytype == BoundaryType::LOWER);
    assert_testmodeonly(shellradiuststart < vec_len(pos));
    return -1;
  }

  if (discriminant > 0) {
    // two intersections
    // Use addition for one root and the product of the roots for the other.
    // Subtraction loses precision when the packet speed is close to the shell speed.
    // Equal speeds give a = 0 and a linear equation, with the single root c / distance_numerator = -c / b.
    const double distance_numerator = (-b - std::copysign(std::sqrt(discriminant), b)) / 2.;
    double dist1 = (a != 0.) ? distance_numerator / a : -1.;
    double dist2 = c / distance_numerator;

    std::array<double, S1> posfinal1{};
    std::array<double, S1> posfinal2{};
    for (auto d = 0ZU; d < S1; d++) {
      posfinal1[d] = pos[d] + (dist1 * dir[d]);
      posfinal2[d] = pos[d] + (dist2 * dir[d]);
    }

    const double v_rad_shell = shellradiuststart / tstart;
    const double v_rad_final1 = dot(dir, posfinal1) * speed / vec_len(posfinal1);
    const double v_rad_final2 = dot(dir, posfinal2) * speed / vec_len(posfinal2);

    // Keep only a crossing in the requested radial direction.
    if constexpr (boundarytype == BoundaryType::LOWER) {
      // The packet must enter through the inner shell, so the shell must overtake it.
      if (v_rad_final1 > v_rad_shell) {
        dist1 = -1;
      }
      if (v_rad_final2 > v_rad_shell) {
        dist2 = -1;
      }
    } else {
      // The packet must leave through the outer shell, so it must overtake the shell.
      if (v_rad_final1 < v_rad_shell) {
        dist1 = -1;
      }
      if (v_rad_final2 < v_rad_shell) {
        dist2 = -1;
      }
    }

#if (TESTMODE)
    if (dist1 >= 0) {
      const double shellradiusfinal1 = shellradiuststart / tstart * (tstart + (dist1 / speed));
      assert_testmodeonly(fabs((vec_len(posfinal1) / shellradiusfinal1) - 1.) < 1e-3);
    }
    if (dist2 >= 0) {
      const double shellradiusfinal2 = shellradiuststart / tstart * (tstart + (dist2 / speed));
      assert_testmodeonly(fabs((vec_len(posfinal2) / shellradiusfinal2) - 1.) < 1e-3);
    }
#endif

    // A negative distance lies behind the packet. Keep the first forward crossing.
    if (dist1 < 0 && dist2 < 0) {
      return -1;
    }
    if (dist2 < 0) {
      return dist1;
    }
    if (dist1 < 0) {
      return dist2;
    }
    return std::min(dist1, dist2);
  }

  // A tangent touches the shell and keeps the packet in its current cell.
  assert_testmodeonly(shellradiuststart <= vec_len(pos));
  return -1.;
}

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
// optical-thickness treatment of a cell. Do not renumber: the values appear in the estimators and gridsave files
enum class CellThickness : int {
  THIN = 0,  // detailed opacity treatment
  THICK = 1,  // optically thick: grey opacity treatment
  THICK_VPKT_ONLY = 2,  // treated as optically thick for virtual packets only
};

inline MPI_shared_array<CellThickness> thick_allcells;
inline MPI_shared_array<float> clumpfactor_allcells;

inline ptrdiff_t ngrid{0};

inline MPI_shared_array<float> elem_meanweight_allcells;

// mass fractions of elements in each cell for the current timestep
inline MPI_shared_array<float> elem_massfracs_allcells;

inline MPI_shared_array<float> ion_groundlevelpops_allcells;
inline MPI_shared_array<float> ion_partfuncts_allcells;

[[nodiscard]] auto get_elements_uppermost_ion(std::ptrdiff_t nonemptymgi, int element) -> int;
void set_elements_uppermost_ion(std::ptrdiff_t nonemptymgi, int element, int uppermost_ion);
// The charge transfer reactions and the Anderson acceleration read the NLTE solution ranges, so the
// lowermost ion of the range then needs storage.
constexpr bool NLTE_TRACK_SOLUTION_RANGES = ENABLE_CHARGE_TRANSFER_REACTIONS || NLTE_TE_NNE_USE_ANDERSON_ACCEL;
[[nodiscard]] auto get_elements_lowermost_ion(std::ptrdiff_t nonemptymgi, int element) -> int;
void set_elements_lowermost_ion(std::ptrdiff_t nonemptymgi, int element, int lowermost_ion);
// Exchange the NLTE solution ranges between the node leaders.
void do_MPI_Bcast_nlte_solution_ranges(ptrdiff_t nstart_nonempty, ptrdiff_t ndo_nonempty, int root_node_id);
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto propcell_width_tmin(int cellindex, int axis) -> double;
[[gnu::pure]] [[nodiscard]] auto get_modelcell_assocvolume_tmin(int modelgridindex) -> double;
[[gnu::pure]] [[nodiscard]] auto get_propcell_volume_tmin(int cellindex) -> double;
[[gnu::pure]] [[nodiscard]] auto get_rho_tmin(int modelgridindex) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_rho(std::ptrdiff_t nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nne(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nnetot(int nonemptymgi) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_ffegrp(int modelgridindex) -> float;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_modelcell_mean_radial_pos_tmin(int modelgridindex) -> double;
[[gnu::pure]] [[nodiscard]] auto get_modelcell_mean_radial_pos_squared_tmin(int modelgridindex) -> double;
void set_elem_massfrac(std::ptrdiff_t nonemptymgi, int element, float newmassfrac);
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_elem_numberdens(std::ptrdiff_t nonemptymgi, int element) -> double;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_initenergyq(int modelgridindex) -> double;
void set_nne(int nonemptymgi, float nne);
void set_nnetot(int nonemptymgi);
void set_rho(int nonemptymgi, float rho);
void init_grid();
[[gnu::pure]] [[nodiscard]] auto get_modelinitnucmassfrac(int modelgridindex, int nucindex) -> float;
[[gnu::pure]] [[nodiscard]] auto get_elem_untrackedstable_initmassfrac(std::ptrdiff_t nonemptymgi, int element)
    -> float;
[[gnu::pure]] [[nodiscard]] auto get_elem_massfrac(std::ptrdiff_t nonemptymgi, int element) -> float;
void set_element_meanweight(std::ptrdiff_t nonemptymgi, int element, float meanweight);
[[gnu::pure]] [[nodiscard]] auto get_electronfrac(int nonemptymgi) -> double;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_numpropcells(int modelgridindex) -> int;
// must not be called for an empty model cell
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nonemptymgi_of_mgi(int mgi) -> int;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_mgi_of_nonemptymgi(std::ptrdiff_t nonemptymgi) -> int;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_modelgridtype() -> GridType;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_npts_model() -> int;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nonempty_npts_model() -> int;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_t_model() -> double;
// Three cell index spaces are used throughout ARTIS. They are all plain integers, so mixing them up compiles
// cleanly but silently reads the wrong cell:
//  - cellindex: propagation cell, [0, ngrid). The geometry packets move through; Packet::cellindex.
//  - modelgridindex (mgi): input model grid cell, [0, get_npts_model()). Where the ejecta model file defines
//    density and abundances. Several propagation cells can share one model cell (get_numpropcells()).
//  - nonemptymgi: model cells containing matter, [0, get_nonempty_npts_model()). The per-cell physics arrays
//    are allocated over only these, so this is the index the physics modules use.
// returns -1 for a propagation cell containing no matter, so check before indexing a per-cell array
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_propcell_modelgridindex(int cellindex) -> int;
// returns -1 for a propagation cell containing no matter, so check before indexing a per-cell array
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_propcell_nonemptymgi(int cellindex) -> int;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellindex_from_pos(const Vec3d& pos, double time) -> int;
void read_ejecta_model();
void write_grid_restart_data(int timestep);
// init_grid() sets up the rank ranges of the grid update and the output files
[[nodiscard]] auto get_nstart(int rank) -> int;
[[nodiscard]] auto get_nstart_nonempty(int rank) -> int;
[[nodiscard]] auto get_ndo(int rank) -> int;
[[nodiscard]] auto get_ndo_nonempty(int rank) -> int;
[[gnu::pure]] [[nodiscard]] auto get_totmassnuclide_tmodel(int z, int a) -> double;
[[nodiscard]] auto get_propcell_random_xyz_position_tmin(int cellindex, rngstate_type& rngstate) -> Vec3d;
[[nodiscard]] DEVICE_FUNC auto boundary_distance(const Vec3d& dir, const Vec3d& pos, double tstart, int cellindex)
    -> std::tuple<double, int>;
void set_clumpfactor(int nonemptymgi, float clumpfactor);
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_clumpfactor(int nonemptymgi) -> float;
DEVICE_FUNC void snap_pos_to_cell(Vec3d& pos, double time, int cellindex);

[[nodiscard]] auto calculate_cell_kappagrey(int nonemptymgi) -> float;

// pass tally_stats = false for virtual ray traces on packet copies, so that the global
// cell-crossing and escape counters only reflect real packet propagation
inline void change_cell_or_escape(Packet& pkt, const int next_cellindex, const bool tally_stats = true) {
  if (next_cellindex >= 0) {
    if (next_cellindex != pkt.cellindex) {
      // make the position exactly consistent with the new cell by snapping it onto the
      // crossed boundary, so that rounding errors cannot accumulate over many crossings.
      // (when boundary_distance() returns early due to max_path_step, next_cellindex equals
      // pkt.cellindex and the packet is mid-cell, so it must not be snapped)
      snap_pos_to_cell(pkt.pos, pkt.prop_time, next_cellindex);
    }
    pkt.cellindex = next_cellindex;
    if (tally_stats) {
      stats::increment(stats::Counter::CELLCROSSINGS);
    }
  } else {
    // Then the packet is exiting the grid. We need to record
    // where and at what time it leaves the grid.
    pkt.escape_type = pkt.type;
    pkt.escape_time = static_cast<float>(pkt.prop_time);
    pkt.type = TYPE_ESCAPE;
    if (tally_stats) {
      stats::increment(stats::Counter::PKTESCAPES);
    }
  }
}

inline auto get_ejecta_kinetic_energy() {
  // Fixed by the cell masses and velocities at tmin, so compute once and cache. The Barnes thermalisation of
  // the particles and of the gamma packets reads it once per packet, after the grid setup.
  static const double E_kin = [] {
    double e_kin = 0.;
    for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
      const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);
      const double M_cell = get_rho_tmin(mgi) * grid::get_modelcell_assocvolume_tmin(mgi);
      // the volume averaged mean of r^2 gives the exact kinetic energy of a uniform-density cell.
      // The square of the mean radius underestimates it.
      e_kin += 0.5 * M_cell * get_modelcell_mean_radial_pos_squared_tmin(mgi) / pow2(globals::tmin);
    }
    return e_kin;
  }();

  return E_kin;
}

// The grid-mapped mass [g] of the ejecta: the sum over the nonempty cells, consistent with
// get_ejecta_kinetic_energy(). The staircase mapping and the FORCE_SPHERICAL_ESCAPE_SURFACE
// blanking make this differ from mtot_input, the analytic mass of the input model.
inline auto get_ejecta_mass() {
  static const double M_ejecta = [] {
    double mtot = 0.;
    for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
      const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);
      mtot += get_rho_tmin(mgi) * grid::get_modelcell_assocvolume_tmin(mgi);
    }
    return mtot;
  }();

  return M_ejecta;
}

}  // namespace grid

#endif  // GRID_H

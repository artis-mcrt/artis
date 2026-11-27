// Declarations for the model/propagation grid (grid.cc, grid_geometry.cc, and grid_modelinput.cc) and inline
// cell-property accessors.

#ifndef GRID_H
#define GRID_H

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
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

// The inner (LOWER) and the outer (UPPER) boundary of a cell along one coordinate. The grid expands homologously, so
// each boundary moves with a speed that is proportional to its coordinate value.
enum class BoundaryType : std::uint8_t { LOWER, UPPER };

// Get the path length d from pos along dir to the point where a packet crosses an expanding shell. The shell is a
// sphere around the origin for a cell boundary of a 1D grid and for the spherical escape surface. For the cylindrical
// radius of a 2D grid, the shell is a circle in the xy plane. Then pos and dir are the projections into that plane,
// and packet_speed is the speed of the packet in that plane. In each case, dir must be a unit vector.
//
// The shell radius increases homologously with the time of flight d / packet_speed:
//   R(d) = R_0 * (t_0 + d / packet_speed) / t_0, with R_0 = shell_radius_at_tstart and t_0 = tstart.
// For an UPPER boundary, the result is the distance to the point where the packet leaves the shell from inside. For
// a LOWER boundary, the result is the distance to the point where the shell overtakes a packet outside the shell.
// The function returns -1 if no such crossing lies ahead of the packet.
template <BoundaryType boundarytype, size_t VECDIM>
[[gnu::pure]] [[nodiscard]] constexpr auto expanding_shell_intersection(const std::array<double, VECDIM>& pos,
                                                                        const std::array<double, VECDIM>& dir,
                                                                        const double packet_speed,
                                                                        const double shell_radius_at_tstart,
                                                                        const double tstart) -> double {
  static_assert(VECDIM == 2 || VECDIM == 3);
  assert_testmodeonly(shell_radius_at_tstart > 0);
  assert_testmodeonly(std::fabs(dot(dir, dir) - 1.) < 1e-6);

  // The packet is at pos + d * dir. The crossing condition |pos + d * dir|^2 = R(d)^2 gives the quadratic equation
  // a * d^2 + b * d + c = 0. With the shell speed v_shell = R_0 / t_0, the coefficients are:
  //   a = |dir|^2 - (v_shell / packet_speed)^2, which is positive when the packet is faster than the shell;
  //   b = 2 * (pos . dir - R_0 * v_shell / packet_speed);
  //   c = |pos|^2 - R_0^2, which is negative for a packet inside the shell and positive for a packet outside it.
  const double a = dot(dir, dir) - pow2(shell_radius_at_tstart / tstart / packet_speed);
  const double b = 2 * (dot(dir, pos) - (pow2(shell_radius_at_tstart) / tstart / packet_speed));
  const double c = dot(pos, pos) - pow2(shell_radius_at_tstart);

  const double discriminant = pow2(b) - (4 * a * c);

  if (discriminant < 0) {
    // The path of the packet does not meet the shell. This occurs only for a packet that passes outside an inner shell.
    assert_testmodeonly(boundarytype == BoundaryType::LOWER);
    assert_testmodeonly(shell_radius_at_tstart < vec_len(pos));
    return -1.;
  }

  if (discriminant == 0) {
    // The path of the packet touches the shell at one point and does not cross it, so the packet stays in its cell.
    assert_testmodeonly(shell_radius_at_tstart <= vec_len(pos));
    return -1.;
  }

  // Select the root for the direction of the crossing. The function f(d) = |pos + d * dir|^2 - R(d)^2 is negative
  // inside the shell and positive outside it. At a root, f(d) = 0 and the derivative f'(d) = 2 * a * d + b is
  // +sqrt(discriminant) or -sqrt(discriminant). At a crossing, |pos + d * dir| = R(d), and
  //   f'(d) = (2 * R(d) / packet_speed) * (v_rad - v_shell),
  // with the radial speed v_rad = packet_speed * dir . (pos + d * dir) / R(d) of the packet. So, at the root with
  // +sqrt(discriminant), the packet moves outward faster than the shell and leaves it (UPPER). At the root with
  // -sqrt(discriminant), the shell moves faster than the packet and overtakes it (LOWER).
  //
  // Each root has two equal forms, for example (-b + sqrt(discriminant)) / (2 * a) = 2 * c / (-b - sqrt(discriminant)).
  // The calculation uses the form that adds two numbers of the same sign, because the difference of two almost equal
  // numbers loses precision. For equal packet and shell speeds, a = 0. One form then gives the root -c / b of the
  // linear equation, and the other form gives an infinite distance.
  const double sqrt_discriminant = std::sqrt(discriminant);
  double distance_to_crossing = -1.;
  if constexpr (boundarytype == BoundaryType::UPPER) {
    distance_to_crossing = (b < 0.) ? (-b + sqrt_discriminant) / (2 * a) : 2 * c / (-b - sqrt_discriminant);
  } else {
    distance_to_crossing = (b > 0.) ? (-b - sqrt_discriminant) / (2 * a) : 2 * c / (-b + sqrt_discriminant);
  }

  // A negative distance is behind the packet. An infinite distance is a crossing that never occurs.
  if (!std::isfinite(distance_to_crossing) || distance_to_crossing < 0.) {
    return -1.;
  }

#if (TESTMODE)
  // At the crossing, the distance of the packet from the centre must equal the shell radius at that time.
  std::array<double, VECDIM> pos_at_crossing{};
  for (auto axis = 0ZU; axis < VECDIM; axis++) {
    pos_at_crossing[axis] = pos[axis] + (distance_to_crossing * dir[axis]);
  }
  const double shell_radius_at_crossing =
      shell_radius_at_tstart / tstart * (tstart + (distance_to_crossing / packet_speed));
  assert_testmodeonly(fabs((vec_len(pos_at_crossing) / shell_radius_at_crossing) - 1.) < 1e-3);
#endif

  return distance_to_crossing;
}

// Get the path length d from pos along dir to the point where a packet crosses a cell boundary of cos(theta) of a
// spherical grid, where theta is the polar angle from the +z axis. The boundary is the cone with the half angle
// theta_edge around the z axis and its apex at the centre. The homologous expansion does not move this cone, thus the
// distance does not depend on the time. dir must be a unit vector.
//
// For an UPPER boundary, the result is the distance to the point where cos(theta) of the packet increases through
// cos_theta_edge. For a LOWER boundary, cos(theta) decreases through cos_theta_edge there. The function returns -1 if
// no such crossing lies ahead of the packet. cos(theta) along a straight path has at most one extremum, thus a packet
// can cross each of the two boundaries of its cell, whatever its present direction is.
template <BoundaryType boundarytype>
[[gnu::pure]] [[nodiscard]] constexpr auto costheta_cone_intersection(const Vec3d& pos, const Vec3d& dir,
                                                                      const double cos_theta_edge) -> double {
  assert_testmodeonly(std::fabs(dot(dir, dir) - 1.) < 1e-6);

  double distance_to_crossing = -1.;
  if (cos_theta_edge == 0.) {
    // The cone is the plane z = 0. cos(theta) increases where z increases.
    const bool moves_towards_boundary = (boundarytype == BoundaryType::UPPER) ? (dir[2] > 0.) : (dir[2] < 0.);
    if (!moves_towards_boundary) {
      return -1.;
    }
    distance_to_crossing = -pos[2] / dir[2];
  } else {
    // The packet is at pos + d * dir. The function f(d) = z(d)^2 - cos(theta_edge)^2 * r(d)^2 is zero on the cone
    // and on its mirror image below the plane z = 0. It is positive nearer to the z axis than the cone, and negative
    // nearer to the plane z = 0. f(d) = a * d^2 + 2 * half_b * d + c.
    const double cos_squared = pow2(cos_theta_edge);
    const double a = pow2(dir[2]) - cos_squared;
    const double half_b = (pos[2] * dir[2]) - (cos_squared * dot(pos, dir));
    const double c = pow2(pos[2]) - (cos_squared * dot(pos, pos));

    // half_b^2 - a * c with the angular momentum L = pos x dir. This form has no difference of two almost equal
    // numbers for a cone near the plane z = 0.
    const auto angmom = cross_prod(pos, dir);
    const double quarter_discriminant =
        cos_squared * (pow2(angmom[0]) + pow2(angmom[1]) - (cos_squared * dot(angmom, angmom)));
    if (quarter_discriminant <= 0.) {
      // The path does not meet the cone or touches it at one point, e.g. a radial path.
      return -1.;
    }
    const double half_sqrt_discriminant = std::sqrt(quarter_discriminant);

    // At a root, the derivative f'(d) is +2 * half_sqrt_discriminant or -2 * half_sqrt_discriminant. cos(theta)
    // increases through the cone where the packet moves towards the +z axis. That is an increase of f for a cone
    // above the plane z = 0, and a decrease of f for a cone below it. Each root has two equal forms, and the
    // calculation uses the form that adds two numbers of the same sign (see expanding_shell_intersection).
    const bool costheta_increases = (boundarytype == BoundaryType::UPPER);
    const bool crossing_has_increasing_f = (cos_theta_edge > 0.) == costheta_increases;
    if (crossing_has_increasing_f) {
      distance_to_crossing =
          (half_b < 0.) ? (-half_b + half_sqrt_discriminant) / a : c / (-half_b - half_sqrt_discriminant);
    } else {
      distance_to_crossing =
          (half_b > 0.) ? (-half_b - half_sqrt_discriminant) / a : c / (-half_b + half_sqrt_discriminant);
    }
  }

  // A negative distance is behind the packet. An infinite distance is a crossing that never occurs.
  if (!std::isfinite(distance_to_crossing) || distance_to_crossing < 0.) {
    return -1.;
  }
  // a root on the other side of the plane z = 0 is a crossing of the mirror image of the cone
  if (cos_theta_edge != 0. && ((pos[2] + (distance_to_crossing * dir[2])) * cos_theta_edge) <= 0.) {
    return -1.;
  }
  return distance_to_crossing;
}

// Get the path length d from pos along dir to the point where a packet crosses a cell boundary of the azimuthal
// angle phi = atan2(y, x) of a spherical grid. The boundary is the half of a plane through the z axis with phi =
// phi_edge. The homologous expansion does not move this half-plane.
//
// For an UPPER boundary, the result is the distance to the point where phi of the packet increases through phi_edge.
// For a LOWER boundary, phi decreases through phi_edge there. The function returns -1 if no such crossing lies ahead
// of the packet. phi along a straight path changes in one sense only, which is the sign of the z component of the
// angular momentum pos x dir.
template <BoundaryType boundarytype>
[[gnu::pure]] [[nodiscard]] constexpr auto phi_halfplane_intersection(const Vec3d& pos, const Vec3d& dir,
                                                                      const double cos_phi_edge,
                                                                      const double sin_phi_edge) -> double {
  // The signed distance of the packet from the full plane is positive on the side of larger phi near the half-plane.
  // It changes at a constant rate along the path.
  const double distance_from_plane = (cos_phi_edge * pos[1]) - (sin_phi_edge * pos[0]);
  const double distance_from_plane_rate = (cos_phi_edge * dir[1]) - (sin_phi_edge * dir[0]);
  const double angmom_z = (pos[0] * dir[1]) - (pos[1] * dir[0]);

  // A path meets the full plane one time. That point is on the half with phi = phi_edge only if the sign of the rate
  // agrees with the sense of phi. If it does not, the point is on the half with phi = phi_edge + pi.
  const bool crosses_halfplane = (boundarytype == BoundaryType::UPPER)
                                     ? (angmom_z > 0. && distance_from_plane_rate > 0.)
                                     : (angmom_z < 0. && distance_from_plane_rate < 0.);
  if (!crosses_halfplane) {
    return -1.;
  }
  const double distance_to_crossing = -distance_from_plane / distance_from_plane_rate;
  if (!std::isfinite(distance_to_crossing) || distance_to_crossing < 0.) {
    return -1.;
  }
  return distance_to_crossing;
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
[[gnu::pure]] [[nodiscard]] auto get_thermal_input_energy(int modelgridindex) -> double;
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
// Set up a 2D spherical or a 3D spherical propagation grid from the cell edges of the radial velocity [cm/s], of
// cos(theta) of the polar angle from -1 to 1, and of the azimuthal angle [rad], with no model.txt. globals::tmin must
// have its value. A 2D spherical grid ignores phi_edges. The unit tests of the cell geometry use this function.
void setup_spherical_grid_from_edges(GridType spherical_grid_type, std::span<const double> vel_r_edges,
                                     std::span<const double> costheta_edges, std::span<const double> phi_edges);
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
